from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

import dpdata.formats.pwmat.atomconfig
import dpdata.formats.pwmat.movement
from dpdata.format import Format
from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType


@Format.register("movement")
@Format.register("mlmd")
@Format.register("pwmat/movement")
@Format.register("pwmat/mlmd")
@Format.register("pwmat/output")
class PwmatOutputFormat(Format):
    """PWmat ``MOVEMENT``/``OUT.MLMD`` labeled trajectory.

    `PWmat <https://www.pwmat.com/>`_ is a plane-wave based DFT
    electronic-structure calculation software using GPU acceleration.

    These output files contain a sequence of cells and coordinates together
    with energies and optional force or virial labels. The reader supports
    frame subsampling and convergence filtering.
    """

    @Format.post("rot_lower_triangular")
    def from_labeled_system(
        self, file_name, begin=0, step=1, convergence_check=True, **kwargs
    ):
        """Load a labeled PWmat trajectory.

        Parameters
        ----------
        file_name : str or os.PathLike
            PWmat ``MOVEMENT`` or ``OUT.MLMD`` file.
        begin : int, default=0
            Index of the first frame to load.
        step : int, default=1
            Load every ``step``-th frame.
        convergence_check : bool, default=True
            Exclude frames marked as unconverged when enabled.
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
            data["atom_numbs"],
            data["atom_types"],
            data["cells"],
            data["coords"],
            data["energies"],
            tmp_force,
            tmp_virial,
        ) = dpdata.formats.pwmat.movement.get_frames(
            file_name, begin=begin, step=step, convergence_check=convergence_check
        )
        if tmp_force is not None:
            data["forces"] = tmp_force
        if tmp_virial is not None:
            data["virials"] = tmp_virial
        # scale virial to the unit of eV
        if "virials" in data:
            v_pref = 1 * 1e3 / 1.602176621e6
            for ii in range(data["coords"].shape[0]):
                vol = np.linalg.det(np.reshape(data["cells"][ii], [3, 3]))
                data["virials"][ii] *= v_pref * vol
        return data


@Format.register("atom.config")
@Format.register("final.config")
@Format.register("pwmat/atom.config")
@Format.register("pwmat/final.config")
class PwmatAtomconfigFormat(Format):
    """PWmat ``atom.config`` or ``final.config`` structure file.

    `PWmat <https://www.pwmat.com/>`_ is a plane-wave DFT software using GPU
    acceleration.

    The format stores a cell and one atomic configuration. Reading normalizes
    the cell to dpdata's lower-triangular convention.
    """

    @Format.post("rot_lower_triangular")
    def from_system(self, file_name: FileType, **kwargs):
        """Load one PWmat configuration.

        Parameters
        ----------
        file_name : str or os.PathLike or file-like object
            ``atom.config`` or ``final.config`` input.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            System data for one configuration.
        """
        with open_file(file_name) as fp:
            lines = [line.rstrip("\n") for line in fp]
        return dpdata.formats.pwmat.atomconfig.to_system_data(lines)

    def to_system(self, data, file_name: FileType, frame_idx=0, *args, **kwargs):
        """Dump the system in pwmat atom.config format.

        Parameters
        ----------
        data : dict
            The system data
        file_name : str
            The output file name
        frame_idx : int
            The index of the frame to dump
        *args : list
            other parameters
        **kwargs : dict
            other parameters
        """
        assert frame_idx < len(data["coords"])
        w_str = dpdata.formats.pwmat.atomconfig.from_system_data(data, frame_idx)
        with open_file(file_name, "w") as fp:
            fp.write(w_str)
