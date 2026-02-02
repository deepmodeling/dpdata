#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

ptr_float_fmt = "%15.10f"
ptr_int_fmt = "%6d"
ptr_key_fmt = "%15s"

# Mapping of LAMMPS atom styles to their column layouts
# Format: (atom_id_col, atom_type_col, x_col, y_col, z_col, has_molecule_id, has_charge, charge_col)
ATOM_STYLE_COLUMNS = {
    "atomic": (0, 1, 2, 3, 4, False, False, None),
    "angle": (0, 2, 3, 4, 5, True, False, None),
    "bond": (0, 2, 3, 4, 5, True, False, None),
    "charge": (0, 1, 3, 4, 5, False, True, 2),
    "full": (0, 2, 4, 5, 6, True, True, 3),
    "molecular": (0, 2, 3, 4, 5, True, False, None),
    "dipole": (0, 1, 3, 4, 5, False, True, 2),
    "sphere": (0, 1, 4, 5, 6, False, False, None),
}


def detect_atom_style(lines: list[str]) -> str | None:
    """Detect LAMMPS atom style from data file content.

    Parameters
    ----------
    lines : list
        Lines from LAMMPS data file

    Returns
    -------
    str or None
        Detected atom style, or None if not detected
    """
    # Look for atom style in comments after "Atoms" section header
    atom_lines = get_atoms(lines)
    if not atom_lines:
        return None

    # Find the "Atoms" line
    for idx, line in enumerate(lines):
        if "Atoms" in line:
            # Check if there's a comment with atom style after "Atoms"
            if "#" in line:
                comment_part = line.split("#")[1].strip().lower()
                for style in ATOM_STYLE_COLUMNS:
                    if style in comment_part:
                        return style
            break

    # If no explicit style found, try to infer from first data line
    if atom_lines:
        first_line = atom_lines[0].split()
        num_cols = len(first_line)

        # Try to match based on number of columns and content patterns
        # This is a heuristic approach
        if num_cols == 5:
            # Could be atomic style: atom-ID atom-type x y z
            return "atomic"
        elif num_cols == 6:
            # Could be charge or bond/molecular style
            # Try to determine if column 2 (index 2) looks like a charge (float) or type (int)
            try:
                val = float(first_line[2])
                # If it's a small float, likely a charge
                if abs(val) < 10 and val != int(val):
                    return "charge"
                else:
                    # Likely molecule ID (integer), so bond/molecular style
                    return "bond"
            except ValueError:
                return "atomic"  # fallback
        elif num_cols == 7:
            # Could be full style: atom-ID molecule-ID atom-type charge x y z
            return "full"
        elif num_cols >= 8:
            # Could be dipole or sphere style
            # For now, default to dipole if we have enough columns
            return "dipole"

    return None  # Unable to detect


def _get_block(lines, keys):
    for idx in range(len(lines)):
        if keys in lines[idx]:
            break
    if idx == len(lines) - 1:
        return None
    idx_s = idx + 2
    idx = idx_s
    ret = []
    while True:
        if idx == len(lines) or len(lines[idx].split()) == 0:
            break
        else:
            ret.append(lines[idx])
        idx += 1
    return ret


def lmpbox2box(lohi, tilt):
    xy = tilt[0]
    xz = tilt[1]
    yz = tilt[2]
    orig = np.array([lohi[0][0], lohi[1][0], lohi[2][0]])
    lens = []
    for dd in range(3):
        lens.append(lohi[dd][1] - lohi[dd][0])
    xx = [lens[0], 0, 0]
    yy = [xy, lens[1], 0]
    zz = [xz, yz, lens[2]]
    return orig, np.array([xx, yy, zz])


def box2lmpbox(orig, box):
    lohi = np.zeros([3, 2])
    for dd in range(3):
        lohi[dd][0] = orig[dd]
    tilt = np.zeros(3)
    tilt[0] = box[1][0]
    tilt[1] = box[2][0]
    tilt[2] = box[2][1]
    lens = np.zeros(3)
    lens[0] = box[0][0]
    lens[1] = box[1][1]
    lens[2] = box[2][2]
    for dd in range(3):
        lohi[dd][1] = lohi[dd][0] + lens[dd]
    return lohi, tilt


def get_atoms(lines):
    return _get_block(lines, "Atoms")


def get_natoms(lines):
    for ii in lines:
        if "atoms" in ii:
            return int(ii.split()[0])
    return None


def get_natomtypes(lines):
    for ii in lines:
        if "atom types" in ii:
            return int(ii.split()[0])
    return None


def _atom_info_mol(line):
    vec = line.split()
    # idx, mole_type, atom_type, charge, x, y, z
    return (
        int(vec[0]),
        int(vec[1]),
        int(vec[2]),
        float(vec[3]),
        float(vec[4]),
        float(vec[5]),
        float(vec[6]),
    )


def _atom_info_atom(line):
    vec = line.split()
    # idx, atom_type, x, y, z
    return int(vec[0]), int(vec[1]), float(vec[2]), float(vec[3]), float(vec[4])


def _atom_info_style(line: str, atom_style: str = "atomic") -> dict[str, int | float]:
    """Parse atom information based on the specified atom style.

    Parameters
    ----------
    line : str
        The atom line from LAMMPS data file
    atom_style : str
        The LAMMPS atom style (atomic, full, charge, etc.)

    Returns
    -------
    dict
        Dictionary containing parsed atom information with keys:
        'atom_id', 'atom_type', 'x', 'y', 'z', 'molecule_id' (if present), 'charge' (if present)
    """
    if atom_style not in ATOM_STYLE_COLUMNS:
        raise ValueError(
            f"Unsupported atom style: {atom_style}. Supported styles: {list(ATOM_STYLE_COLUMNS.keys())}"
        )

    vec = line.split()
    columns = ATOM_STYLE_COLUMNS[atom_style]

    result = {
        "atom_id": int(vec[columns[0]]),
        "atom_type": int(vec[columns[1]]),
        "x": float(vec[columns[2]]),
        "y": float(vec[columns[3]]),
        "z": float(vec[columns[4]]),
    }

    # Add molecule ID if present
    if columns[5]:  # has_molecule_id
        result["molecule_id"] = int(
            vec[1]
        )  # molecule ID is always in column 1 when present

    # Add charge if present
    if columns[6]:  # has_charge
        result["charge"] = float(vec[columns[7]])  # charge_col

    return result


def get_natoms_vec(lines: list[str], atom_style: str = "atomic") -> list[int]:
    """Get number of atoms for each atom type.

    Parameters
    ----------
    lines : list
        Lines from LAMMPS data file
    atom_style : str
        The LAMMPS atom style

    Returns
    -------
    list
        Number of atoms for each atom type
    """
    atype = get_atype(lines, atom_style=atom_style)
    natoms_vec = []
    natomtypes = get_natomtypes(lines)
    for ii in range(natomtypes):
        natoms_vec.append(sum(atype == ii + 1))
    assert sum(natoms_vec) == get_natoms(lines)
    return natoms_vec


def get_atype(
    lines: list[str], type_idx_zero: bool = False, atom_style: str = "atomic"
) -> np.ndarray:
    """Get atom types from LAMMPS data file.

    Parameters
    ----------
    lines : list
        Lines from LAMMPS data file
    type_idx_zero : bool
        Whether to use zero-based indexing for atom types
    atom_style : str
        The LAMMPS atom style

    Returns
    -------
    np.ndarray
        Array of atom types
    """
    alines = get_atoms(lines)
    atype = []
    for ii in alines:
        atom_info = _atom_info_style(ii, atom_style)
        at = atom_info["atom_type"]
        if type_idx_zero:
            atype.append(at - 1)
        else:
            atype.append(at)
    return np.array(atype, dtype=int)


def get_posi(lines: list[str], atom_style: str = "atomic") -> np.ndarray:
    """Get atomic positions from LAMMPS data file.

    Parameters
    ----------
    lines : list
        Lines from LAMMPS data file
    atom_style : str
        The LAMMPS atom style

    Returns
    -------
    np.ndarray
        Array of atomic positions
    """
    atom_lines = get_atoms(lines)
    posis = []
    for ii in atom_lines:
        atom_info = _atom_info_style(ii, atom_style)
        posis.append([atom_info["x"], atom_info["y"], atom_info["z"]])
    return np.array(posis)


def get_charges(lines: list[str], atom_style: str = "atomic") -> np.ndarray | None:
    """Get atomic charges from LAMMPS data file if the atom style supports charges.

    Parameters
    ----------
    lines : list
        Lines from LAMMPS data file
    atom_style : str
        The LAMMPS atom style

    Returns
    -------
    np.ndarray or None
        Array of atomic charges if atom style has charges, None otherwise
    """
    if atom_style not in ATOM_STYLE_COLUMNS:
        raise ValueError(f"Unsupported atom style: {atom_style}")

    # Check if this atom style has charges
    if not ATOM_STYLE_COLUMNS[atom_style][6]:  # has_charge
        return None

    atom_lines = get_atoms(lines)
    charges = []
    for ii in atom_lines:
        atom_info = _atom_info_style(ii, atom_style)
        charges.append(atom_info["charge"])
    return np.array(charges)


def get_spins(lines: list[str], atom_style: str = "atomic") -> np.ndarray | None:
    atom_lines = get_atoms(lines)
    if len(atom_lines[0].split()) < 8:
        return None
    spins_ori = []
    spins_norm = []
    for ii in atom_lines:
        iis = ii.split()
        spins_ori.append([float(jj) for jj in iis[5:8]])
        spins_norm.append([float(iis[-1])])
    return np.array(spins_ori) * np.array(spins_norm)


def get_lmpbox(lines):
    box_info = []
    tilt = np.zeros(3)
    for ii in lines:
        if "xlo" in ii and "xhi" in ii:
            box_info.append([float(ii.split()[0]), float(ii.split()[1])])
            break
    for ii in lines:
        if "ylo" in ii and "yhi" in ii:
            box_info.append([float(ii.split()[0]), float(ii.split()[1])])
            break
    for ii in lines:
        if "zlo" in ii and "zhi" in ii:
            box_info.append([float(ii.split()[0]), float(ii.split()[1])])
            break
    for ii in lines:
        if "xy" in ii and "xz" in ii and "yz" in ii:
            tilt = np.array([float(jj) for jj in ii.split()[0:3]])
    return box_info, tilt


def system_data(
    lines: list[str],
    type_map: list[str] | None = None,
    type_idx_zero: bool = True,
    atom_style: str = "atomic",
) -> dict:
    """Parse LAMMPS data file to system data format.

    Parameters
    ----------
    lines : list
        Lines from LAMMPS data file
    type_map : list, optional
        Mapping from atom types to element names
    type_idx_zero : bool
        Whether to use zero-based indexing for atom types
    atom_style : str
        The LAMMPS atom style (atomic, full, charge, etc.)

    Returns
    -------
    dict
        System data dictionary
    """
    system = {}
    system["atom_numbs"] = get_natoms_vec(lines, atom_style=atom_style)
    system["atom_names"] = []
    if type_map is None:
        for ii in range(len(system["atom_numbs"])):
            system["atom_names"].append("Type_%d" % ii)  # noqa: UP031
    else:
        assert len(type_map) >= len(system["atom_numbs"])
        for ii in range(len(system["atom_numbs"])):
            system["atom_names"].append(type_map[ii])
    lohi, tilt = get_lmpbox(lines)
    orig, cell = lmpbox2box(lohi, tilt)
    system["orig"] = np.array(orig)
    system["cells"] = [np.array(cell)]
    natoms = sum(system["atom_numbs"])
    system["atom_types"] = get_atype(
        lines, type_idx_zero=type_idx_zero, atom_style=atom_style
    )
    system["coords"] = [get_posi(lines, atom_style=atom_style)]
    system["cells"] = np.array(system["cells"])
    system["coords"] = np.array(system["coords"])

    # Add charges if the atom style supports them
    charges = get_charges(lines, atom_style=atom_style)
    if charges is not None:
        system["charges"] = np.array([charges])

    spins = get_spins(lines, atom_style=atom_style)
    if spins is not None:
        system["spins"] = np.array([spins])

    return system


def to_system_data(
    lines: list[str],
    type_map: list[str] | None = None,
    type_idx_zero: bool = True,
    atom_style: str = "atomic",
) -> dict:
    """Parse LAMMPS data file to system data format.

    Parameters
    ----------
    lines : list
        Lines from LAMMPS data file
    type_map : list, optional
        Mapping from atom types to element names
    type_idx_zero : bool
        Whether to use zero-based indexing for atom types
    atom_style : str
        The LAMMPS atom style. If "auto", attempts to detect automatically
        from file. Default is "atomic".

    Returns
    -------
    dict
        System data dictionary
    """
    # Attempt automatic detection if requested
    if atom_style == "auto":
        detected_style = detect_atom_style(lines)
        if detected_style:
            atom_style = detected_style
        else:
            atom_style = "atomic"  # fallback to default

    return system_data(
        lines, type_map=type_map, type_idx_zero=type_idx_zero, atom_style=atom_style
    )


def rotate_to_lower_triangle(
    cell: np.ndarray, coord: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Rotate the cell to lower triangular and ensure the diagonal elements are non-negative.

    Args:
        cell (np.ndarray): The original cell matrix.
        coord (np.ndarray): The coordinates of the atoms.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]: The rotated cell and adjusted coordinates.
    """
    q, _ = np.linalg.qr(cell.T)
    cell = np.matmul(cell, q)
    coord = np.matmul(coord, q)

    # Ensure the diagonal elements of the cell are non-negative
    rot = np.eye(3)
    if cell[0][0] < 0:
        rot[0][0] = -1
    if cell[1][1] < 0:
        rot[1][1] = -1
    if cell[2][2] < 0:
        rot[2][2] = -1
    cell = np.matmul(cell, rot)
    coord = np.matmul(coord, rot)
    return cell, coord


def from_system_data(system, f_idx=0):
    ret = ""
    ret += "\n"
    natoms = sum(system["atom_numbs"])
    ntypes = len(system["atom_numbs"])
    cell, coord = rotate_to_lower_triangle(
        system["cells"][f_idx], system["coords"][f_idx]
    )
    ret += "%d atoms\n" % natoms  # noqa: UP031
    ret += "%d atom types\n" % ntypes  # noqa: UP031
    ret += (ptr_float_fmt + " " + ptr_float_fmt + " xlo xhi\n") % (
        0,
        cell[0][0],
    )  # noqa: UP031
    ret += (ptr_float_fmt + " " + ptr_float_fmt + " ylo yhi\n") % (
        0,
        cell[1][1],
    )  # noqa: UP031
    ret += (ptr_float_fmt + " " + ptr_float_fmt + " zlo zhi\n") % (
        0,
        cell[2][2],
    )  # noqa: UP031
    ret += (
        ptr_float_fmt + " " + ptr_float_fmt + " " + ptr_float_fmt + " xy xz yz\n"
    ) % (
        cell[1][0],
        cell[2][0],
        cell[2][1],
    )  # noqa: UP031
    ret += "\n"
    ret += "Atoms # atomic\n"
    ret += "\n"
    coord_fmt = (
        ptr_int_fmt
        + " "
        + ptr_int_fmt
        + " "
        + ptr_float_fmt
        + " "
        + ptr_float_fmt
        + " "
        + ptr_float_fmt
        + "\n"
    )  # noqa: UP031

    if "spins" in system:
        coord_fmt = (
            coord_fmt.strip("\n")
            + " "
            + ptr_float_fmt
            + " "
            + ptr_float_fmt
            + " "
            + ptr_float_fmt
            + " "
            + ptr_float_fmt
            + "\n"
        )  # noqa: UP031
        spins_norm = np.linalg.norm(system["spins"][f_idx], axis=1)
    for ii in range(natoms):
        if "spins" in system:
            if spins_norm[ii] != 0:
                ret += coord_fmt % (
                    ii + 1,
                    system["atom_types"][ii] + 1,
                    coord[ii][0] - system["orig"][0],
                    coord[ii][1] - system["orig"][1],
                    coord[ii][2] - system["orig"][2],
                    system["spins"][f_idx][ii][0] / spins_norm[ii],
                    system["spins"][f_idx][ii][1] / spins_norm[ii],
                    system["spins"][f_idx][ii][2] / spins_norm[ii],
                    spins_norm[ii],
                )  # noqa: UP031
            else:
                ret += coord_fmt % (
                    ii + 1,
                    system["atom_types"][ii] + 1,
                    coord[ii][0] - system["orig"][0],
                    coord[ii][1] - system["orig"][1],
                    coord[ii][2] - system["orig"][2],
                    system["spins"][f_idx][ii][0],
                    system["spins"][f_idx][ii][1],
                    system["spins"][f_idx][ii][2] + 1,
                    spins_norm[ii],
                )  # noqa: UP031
        else:
            ret += coord_fmt % (
                ii + 1,
                system["atom_types"][ii] + 1,
                coord[ii][0] - system["orig"][0],
                coord[ii][1] - system["orig"][1],
                coord[ii][2] - system["orig"][2],
            )  # noqa: UP031
    return ret


if __name__ == "__main__":
    fname = "water-SPCE.data"
    lines = open(fname).read().split("\n")
    bonds, tilt = get_lmpbox(lines)
    # print(bonds, tilt)
    orig, box = lmpbox2box(bonds, tilt)
    # print(orig, box)
    bonds1, tilt1 = box2lmpbox(orig, box)
    # print(bonds1, tilt1)
    print(bonds1 - bonds)
    print(tilt1 - tilt)
    print(box)
    print(get_atype(lines))
    print(get_posi(lines))
