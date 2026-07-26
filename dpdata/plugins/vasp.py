from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np

import dpdata.formats.vasp.outcar
import dpdata.formats.vasp.poscar
import dpdata.formats.vasp.xml
from dpdata.data_type import Axis, DataType
from dpdata.format import Format
from dpdata.utils import open_file, uniq_atom_names

if TYPE_CHECKING:
    from dpdata.utils import FileType


def register_move_data(data):
    if "move" in data:
        dt = DataType(
            "move",
            np.ndarray,
            (Axis.NFRAMES, Axis.NATOMS, 3),
            required=False,
            deepmd_name="move",
        )
        dpdata.System.register_data_type(dt)


@Format.register("poscar")
@Format.register("contcar")
@Format.register("vasp/poscar")
@Format.register("vasp/contcar")
class VASPPoscarFormat(Format):
    """VASP POSCAR or CONTCAR structure file.

    `VASP <https://www.vasp.at/>`_ (Vienna Ab initio Simulation Package) is
    a computer program for atomic scale materials modelling.

    POSCAR/CONTCAR stores one periodic configuration and optional selective
    dynamics flags. It does not contain energies or forces, so it maps to
    :class:`dpdata.System` rather than :class:`dpdata.LabeledSystem`.
    """

    @Format.post("rot_lower_triangular")
    def from_system(self, file_name: FileType, **kwargs):
        """Load a VASP POSCAR or CONTCAR file.

        Parameters
        ----------
        file_name : str or os.PathLike or file-like object
            POSCAR/CONTCAR input.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            System data, including ``move`` selective-dynamics flags when
            present.
        """
        with open_file(file_name) as fp:
            lines = [line.rstrip("\n") for line in fp]
        data = dpdata.formats.vasp.poscar.to_system_data(lines)
        data = uniq_atom_names(data)
        register_move_data(data)
        return data

    def to_system(self, data, file_name: FileType, frame_idx=0, **kwargs):
        """Dump the system in vasp POSCAR format.

        Parameters
        ----------
        data : dict
            The system data
        file_name : str
            The output file name
        frame_idx : int
            The index of the frame to dump
        **kwargs : dict
            other parameters
        """
        w_str = VASPStringFormat().to_system(data, frame_idx=frame_idx)
        with open_file(file_name, "w") as fp:
            fp.write(w_str)


@Format.register("vasp/string")
class VASPStringFormat(Format):
    """In-memory VASP POSCAR text representation.

    `VASP <https://www.vasp.at/>`_ (Vienna Ab initio Simulation Package) is
    a computer program for atomic scale materials modelling.

    Unlike ``vasp/poscar``, this write-only helper returns the POSCAR content
    as a string instead of writing it to a file.
    """

    def to_system(self, data, frame_idx=0, **kwargs):
        """Dump the system in vasp POSCAR format string.

        Parameters
        ----------
        data : dict
            The system data
        frame_idx : int
            The index of the frame to dump
        **kwargs : dict
            other parameters
        """
        assert frame_idx < len(data["coords"])
        return dpdata.formats.vasp.poscar.from_system_data(data, frame_idx)


# rotate the system to lammps convention
@Format.register("outcar")
@Format.register("vasp/outcar")
class VASPOutcarFormat(Format):
    """VASP ``OUTCAR`` labeled trajectory.

    `VASP <https://www.vasp.at/>`_ (Vienna Ab initio Simulation Package) is
    a computer program for atomic scale materials modelling.

    The reader extracts ionic-step cells, coordinates, energies, forces, and
    virials. It supports frame subsampling, convergence filtering, and
    recursive loading of conventionally named ``OUTCAR`` files into
    :class:`dpdata.MultiSystems`.
    """

    def from_multi_systems(self, directory, **kwargs):
        """Find conventionally named OUTCAR files below ``directory``.

        VASP calculations are commonly stored one calculation per directory,
        so the objects consumed by :meth:`from_labeled_system` are the OUTCAR
        files themselves rather than the calculation directories.  Searching
        recursively also supports grouping calculations below intermediate
        directories such as workflow stages or temperatures.

        Parameters
        ----------
        directory : str or os.PathLike
            Root directory containing VASP calculation directories.
        **kwargs : dict
            Additional format options. They are consumed later when each
            discovered OUTCAR is loaded.

        Returns
        -------
        list[str]
            Deterministically ordered paths to files named ``OUTCAR``.
        """
        outcar_files = []
        for root, _, files in os.walk(directory):
            # Match VASP's canonical output name so unrelated files in the
            # calculation tree are not accidentally parsed as OUTCAR data.
            if "OUTCAR" in files:
                outcar_files.append(os.path.join(root, "OUTCAR"))
        return sorted(outcar_files)

    @Format.post("rot_lower_triangular")
    def from_labeled_system(
        self, file_name, begin=0, step=1, convergence_check=True, **kwargs
    ):
        """Load labeled ionic steps from a VASP OUTCAR.

        Parameters
        ----------
        file_name : str or os.PathLike
            VASP ``OUTCAR`` file.
        begin : int, default=0
            Index of the first ionic step to load.
        step : int, default=1
            Load every ``step``-th ionic step.
        convergence_check : bool, default=True
            Exclude unconverged electronic or ionic steps when enabled.
        **kwargs : dict
            Additional options. ``ml=True`` reads labels from VASP's machine-
            learning force-field output blocks.

        Returns
        -------
        dict
            Labeled trajectory data. Forces or virials are omitted when the
            corresponding records are unavailable.
        """
        data = {}
        ml = kwargs.get("ml", False)
        (
            data["atom_names"],
            data["atom_numbs"],
            data["atom_types"],
            data["cells"],
            data["coords"],
            data["energies"],
            tmp_force,
            tmp_virial,
        ) = dpdata.formats.vasp.outcar.get_frames(
            file_name,
            begin=begin,
            step=step,
            ml=ml,
            convergence_check=convergence_check,
        )
        if tmp_force is not None:
            data["forces"] = tmp_force
        if tmp_virial is not None:
            data["virials"] = tmp_virial
        # scale virial to the unit of eV
        if "virials" in data:
            v_pref = 1 * 1e3 / 1.602176621e6
            for ii in range(data["cells"].shape[0]):
                vol = np.linalg.det(np.reshape(data["cells"][ii], [3, 3]))
                data["virials"][ii] *= v_pref * vol
        data = uniq_atom_names(data)
        register_move_data(data)
        return data


# rotate the system to lammps convention
@Format.register("xml")
@Format.register("vasp/xml")
class VASPXMLFormat(Format):
    """VASP ``vasprun.xml`` labeled trajectory.

    `VASP <https://www.vasp.at/>`_ (Vienna Ab initio Simulation Package) is
    a computer program for atomic scale materials modelling.

    XML output contains structured ionic-step cells, coordinates, energies,
    forces, and stresses and is useful when text ``OUTCAR`` parsing is not
    desired.
    """

    @Format.post("rot_lower_triangular")
    def from_labeled_system(
        self, file_name, begin=0, step=1, convergence_check=True, **kwargs
    ):
        """Load labeled ionic steps from ``vasprun.xml``.

        Parameters
        ----------
        file_name : str or os.PathLike
            VASP XML output file.
        begin : int, default=0
            Index of the first ionic step to load.
        step : int, default=1
            Load every ``step``-th ionic step.
        convergence_check : bool, default=True
            Exclude unconverged calculations when enabled.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled trajectory data.
        """
        data = {}
        (
            data["atom_names"],
            data["atom_types"],
            data["cells"],
            data["coords"],
            data["energies"],
            data["forces"],
            tmp_virial,
        ) = dpdata.formats.vasp.xml.analyze(
            file_name,
            type_idx_zero=True,
            begin=begin,
            step=step,
            convergence_check=convergence_check,
        )
        data["atom_numbs"] = []
        for ii in range(len(data["atom_names"])):
            data["atom_numbs"].append(sum(data["atom_types"] == ii))
        # the vasp xml assumes the direct coordinates
        # apply the transform to the cartesan coordinates
        for ii in range(data["cells"].shape[0]):
            data["coords"][ii] = np.matmul(data["coords"][ii], data["cells"][ii])
        # scale virial to the unit of eV
        if tmp_virial.size > 0:
            data["virials"] = tmp_virial
            v_pref = 1 * 1e3 / 1.602176621e6
            for ii in range(data["cells"].shape[0]):
                vol = np.linalg.det(np.reshape(data["cells"][ii], [3, 3]))
                data["virials"][ii] *= v_pref * vol
        data = uniq_atom_names(data)
        register_move_data(data)
        return data
