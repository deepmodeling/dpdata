#!/usr/bin/env python3
# %%
from __future__ import annotations

import re
import warnings
from collections import OrderedDict

import numpy as np

from dpdata.periodic_table import Element

# Possible keys for the energy field in the extxyz comment line,
# checked in order of priority.
_ENERGY_KEYS = ("energy", "Energy", "free_energy", "REF_energy", "energies")

# Accepted per-atom property names for forces.
_FORCE_KEYS = ("force", "forces")

# Accepted header keys for virial tensor.
_VIRIAL_KEYS = ("virial", "virials")

# Accepted header keys for stress tensor.
_STRESS_KEYS = ("stress", "stresses")


def _parse_stress_to_virials(stress_str, cell, stress_sign=-1):
    """Convert a stress field string to virial tensor.

    Parameters
    ----------
    stress_str : str
        Space-separated stress values.  Accepts either 9 values (3x3 matrix,
        row-major) or 6 values (Voigt notation: xx yy zz yz xz xy).
    cell : np.ndarray
        3x3 cell matrix (angstrom).
    stress_sign : int
        Sign convention for ``virial = stress_sign * volume * stress``.
        Default ``-1`` follows the ASE convention where
        ``virial = -V * stress`` (stress in eV/angstrom^3).

    Returns
    -------
    np.ndarray
        Virial tensor with shape ``(1, 3, 3)`` in eV.
    """
    vals = list(filter(bool, stress_str.split(" ")))
    vals = np.array(vals, dtype=np.float64)
    if len(vals) == 9:
        stress = vals.reshape(3, 3)
    elif len(vals) == 6:
        # Voigt order: xx yy zz yz xz xy
        xx, yy, zz, yz, xz, xy = vals
        stress = np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
    else:
        raise ValueError(
            f"stress field must have 6 (Voigt) or 9 (3x3) values, got {len(vals)}"
        )
    volume = abs(np.linalg.det(cell))
    virials = stress_sign * volume * stress
    return np.array([virials])


class QuipGapxyzSystems:
    """Parse an extended XYZ (QUIP/GAP) file frame by frame."""

    def __init__(self, file_name, **kwargs):
        self.file_object = open(file_name)
        self.kwargs = kwargs
        self.block_generator = self.get_block_generator()

    def __iter__(self):
        return self

    def __next__(self):
        return self.handle_single_xyz_frame(next(self.block_generator), **self.kwargs)

    def __del__(self):
        self.file_object.close()

    def get_block_generator(self):
        p3 = re.compile(r"^\s*(\d+)\s*")
        while True:
            line = self.file_object.readline()
            if not line:
                break
            if p3.match(line):
                atom_num = int(p3.match(line).group(1))
                lines = []
                lines.append(line)
                for ii in range(atom_num + 1):
                    lines.append(self.file_object.readline())
                if not lines[-1]:
                    raise RuntimeError(
                        f"this xyz file may lack of lines, should be {atom_num + 2};lines:{lines}"
                    )
                yield lines

    @staticmethod
    def handle_single_xyz_frame(lines, stress_sign=-1, **kwargs):
        """Parse a single extended XYZ frame.

        Parameters
        ----------
        lines : list[str]
            Raw lines for one frame (atom count + comment + atom lines).
        stress_sign : int, optional
            Sign convention for stress→virial conversion.
            ``-1`` (default) follows the ASE convention:
            ``virial = -V * stress``.
        **kwargs : dict
            Additional keyword arguments (reserved for future use).
        """
        atom_num = int(lines[0].strip("\n").strip())
        if len(lines) != atom_num + 2:
            raise RuntimeError(
                f"format error, atom_num=={atom_num}, {len(lines)}!=atom_num+2"
            )
        data_format_line = lines[1].strip("\n").strip() + " "
        field_value_pattern = re.compile(
            r"(?P<key>\S+)=(?P<quote>[\'\"]?)(?P<value>.*?)(?P=quote)\s+"
        )
        prop_pattern = re.compile(
            r"(?P<key>\w+?):(?P<datatype>[a-zA-Z]):(?P<value>\d+)"
        )

        data_format_list = [
            kv_dict.groupdict()
            for kv_dict in field_value_pattern.finditer(data_format_line)
        ]
        field_dict = {}
        for item in data_format_list:
            field_dict[item["key"]] = item["value"]

        Properties = field_dict["Properties"]
        prop_list = [
            kv_dict.groupdict() for kv_dict in prop_pattern.finditer(Properties)
        ]

        data_lines = []
        for line in lines[2:]:
            data_lines.append(list(filter(bool, line.strip().split())))
        data_array = np.array(data_lines)
        used_colomn = 0

        type_array = None
        coords_array = None
        Z_array = None
        force_array = None
        for kv_dict in prop_list:
            field_length = int(kv_dict["value"])
            key = kv_dict["key"]

            if key == "species":
                if kv_dict["datatype"] != "S":
                    raise RuntimeError(
                        f"datatype for species must be 'S' instead of {kv_dict['datatype']}"
                    )
                type_array = data_array[
                    :, used_colomn : used_colomn + field_length
                ].flatten()
                used_colomn += field_length
            elif key == "pos":
                if kv_dict["datatype"] != "R":
                    raise RuntimeError(
                        f"datatype for pos must be 'R' instead of {kv_dict['datatype']}"
                    )
                coords_array = data_array[:, used_colomn : used_colomn + field_length]
                used_colomn += field_length
            elif key == "Z":
                if kv_dict["datatype"] != "I":
                    raise RuntimeError(
                        f"datatype for Z must be 'I' instead of {kv_dict['datatype']}"
                    )
                Z_array = data_array[
                    :, used_colomn : used_colomn + field_length
                ].flatten()
                used_colomn += field_length
            elif key in _FORCE_KEYS:
                if kv_dict["datatype"] != "R":
                    raise RuntimeError(
                        f"datatype for {key} must be 'R' instead of {kv_dict['datatype']}"
                    )
                force_array = data_array[:, used_colomn : used_colomn + field_length]
                used_colomn += field_length
            else:
                # Skip unknown per-atom properties (e.g. magmom, charges,
                # tags, local_energy) instead of crashing.
                warnings.warn(
                    f"Skipping unknown per-atom property '{key}' "
                    f"(type={kv_dict['datatype']}, width={field_length})",
                    stacklevel=2,
                )
                used_colomn += field_length

        # --- atom type bookkeeping ---
        type_num_dict = OrderedDict()
        atom_type_list = []
        type_map = {}
        temp_atom_max_index = 0
        if type_array is None:
            raise RuntimeError("type_array can't be None type, check .xyz file")
        for ii in type_array:
            if ii not in type_map:
                type_map[ii] = temp_atom_max_index
                temp_atom_max_index += 1
                temp_atom_index = type_map[ii]
                atom_type_list.append(temp_atom_index)
                type_num_dict[ii] = 1
            else:
                temp_atom_index = type_map[ii]
                atom_type_list.append(temp_atom_index)
                type_num_dict[ii] += 1
        type_num_list = []
        for atom_type, atom_num in type_num_dict.items():
            type_num_list.append((atom_type, atom_num))
        type_num_array = np.array(type_num_list)

        # --- cells / Lattice (parsed early so volume is available for stress→virial) ---
        info_dict = {}
        if "Lattice" in field_dict and field_dict["Lattice"].strip():
            lattice_values = list(filter(bool, field_dict["Lattice"].split(" ")))
            cells = np.array(lattice_values, dtype=np.float64).reshape(3, 3)
            info_dict["cells"] = np.array([cells])
            info_dict["nopbc"] = False
        else:
            cells = np.diag([100.0, 100.0, 100.0])
            info_dict["cells"] = np.array([cells])
            info_dict["nopbc"] = True

        # Override nopbc if explicit pbc field is present
        if "pbc" in field_dict:
            pbc_flags = field_dict["pbc"].replace('"', "").replace("'", "").split()
            if all(f.upper() in ("F", "FALSE", "0") for f in pbc_flags):
                info_dict["nopbc"] = True
            elif all(f.upper() in ("T", "TRUE", "1") for f in pbc_flags):
                info_dict["nopbc"] = False

        # --- virial / stress ---
        virials = None
        virial_raw = None
        for vkey in _VIRIAL_KEYS:
            if field_dict.get(vkey):
                virial_raw = field_dict[vkey]
                break
        stress_raw = None
        for skey in _STRESS_KEYS:
            if field_dict.get(skey):
                stress_raw = field_dict[skey]
                break

        if virial_raw is not None:
            virials = np.array(
                [np.array(list(filter(bool, virial_raw.split(" ")))).reshape(3, 3)]
            ).astype(np.float64)
        elif stress_raw is not None:
            virials = _parse_stress_to_virials(
                stress_raw, cells, stress_sign=stress_sign
            )

        # --- energy (try several common keys) ---
        energy_value = None
        for ekey in _ENERGY_KEYS:
            if ekey in field_dict:
                energy_value = field_dict[ekey]
                break
        if energy_value is None:
            raise ValueError(
                f"No energy field found in extxyz comment line. "
                f"Tried: {_ENERGY_KEYS}. Available keys: {list(field_dict.keys())}"
            )

        # --- assemble output ---
        info_dict["atom_names"] = list(type_num_array[:, 0])
        info_dict["atom_numbs"] = list(type_num_array[:, 1].astype(int))
        info_dict["atom_types"] = np.array(atom_type_list).astype(int)
        info_dict["coords"] = np.array([coords_array]).astype(np.float64)
        info_dict["energies"] = np.array([energy_value]).astype(np.float64)
        info_dict["forces"] = np.array([force_array]).astype(np.float64)
        if virials is not None:
            info_dict["virials"] = virials
        info_dict["orig"] = np.zeros(3)
        return info_dict


def format_single_frame(data, frame_idx):
    """Format a single frame of system data into extended XYZ format lines.

    Parameters
    ----------
    data : dict
        system data
    frame_idx : int
        frame index

    Returns
    -------
    list[str]
        lines for the frame
    """
    # Number of atoms
    natoms = len(data["atom_types"])

    # Build header line with metadata
    header_parts = []

    # Energy
    energy = data["energies"][frame_idx]
    header_parts.append(f"energy={energy:.12e}")

    # Virial and stress (if present)
    if "virials" in data:
        virial = data["virials"][frame_idx]
        virial_str = "    ".join(f"{v:.12e}" for v in virial.flatten())
        header_parts.append(f'virial="{virial_str}"')
        # Also write stress for ASE compatibility: stress = -virial / volume
        cell = data["cells"][frame_idx]
        volume = abs(np.linalg.det(cell))
        if volume > 0:
            stress = -virial / volume
            stress_str = "    ".join(f"{s:.12e}" for s in stress.flatten())
            header_parts.append(f'stress="{stress_str}"')

    # Lattice
    cell = data["cells"][frame_idx]
    lattice_str = "   ".join(f"{c:.12e}" for c in cell.flatten())
    header_parts.append(f'Lattice="{lattice_str}"')

    # pbc
    if data.get("nopbc", False):
        header_parts.append('pbc="F F F"')
    else:
        header_parts.append('pbc="T T T"')

    # Properties — use "forces" for ASE compatibility (not "force")
    header_parts.append("Properties=species:S:1:pos:R:3:Z:I:1:forces:R:3")

    header_line = "    ".join(header_parts)

    # Format atom lines
    atom_lines = []
    coords = data["coords"][frame_idx]
    forces = data["forces"][frame_idx]
    atom_names = np.array(data["atom_names"])
    atom_types = data["atom_types"]

    for i in range(natoms):
        atom_type_idx = atom_types[i]
        species = atom_names[atom_type_idx]
        x, y, z = coords[i]
        fx, fy, fz = forces[i]
        atomic_number = Element(species).Z

        atom_line = f"{species}    {x:.11e}   {y:.11e}   {z:.11e}   {atomic_number}    {fx:.11e}  {fy:.11e}   {fz:.11e}"
        atom_lines.append(atom_line)

    # Combine all lines for this frame
    frame_lines = [str(natoms), header_line] + atom_lines
    return frame_lines
