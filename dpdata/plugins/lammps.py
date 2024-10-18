from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

import dpdata.lammps.dump
import dpdata.lammps.lmp
from dpdata.data_type import Axis, DataType
from dpdata.format import Format
from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType


def register_spin(data):
    if "spins" in data:
        dt = DataType(
            "spins",
            np.ndarray,
            (Axis.NFRAMES, Axis.NATOMS, 3),
            required=False,
            deepmd_name="spin",
        )
        dpdata.System.register_data_type(dt)


@Format.register("lmp")
@Format.register("lammps/lmp")
class LAMMPSLmpFormat(Format):
    @Format.post("shift_orig_zero")
    def from_system(self, file_name: FileType, type_map=None, **kwargs):
        with open_file(file_name) as fp:
            lines = [line.rstrip("\n") for line in fp]
        data = dpdata.lammps.lmp.to_system_data(lines, type_map)
        register_spin(data)
        return data

    def to_system(self, data, file_name: FileType, frame_idx=0, **kwargs):
        """Dump the system in lammps data format.

        Parameters
        ----------
        data : dict
            System data
        file_name : str
            The output file name
        frame_idx : int
            The index of the frame to dump
        **kwargs : dict
            other parameters
        """
        assert frame_idx < len(data["coords"])
        w_str = dpdata.lammps.lmp.from_system_data(data, frame_idx)
        with open_file(file_name, "w") as fp:
            fp.write(w_str)


@Format.register("dump")
@Format.register("lammps/dump")
class LAMMPSDumpFormat(Format):
    @Format.post("shift_orig_zero")
    def from_system(
        self,
        file_name: str,
        type_map: list[str] = None,
        begin: int = 0,
        step: int = 1,
        unwrap: bool = False,
        input_file: str = None,
        **kwargs,
    ):
        """Read the data from a lammps dump file.

        Parameters
        ----------
        file_name : str
            The dump file name
        type_map : List[str], optional
            The atom type list
        begin : int, optional
            The begin step
        step : int, optional
            The step
        unwrap : bool, optional
            Whether to unwrap the coordinates
        input_file : str, optional
            The input file name

        Returns
        -------
        dict
            The system data
        """
        lines = dpdata.lammps.dump.load_file(file_name, begin=begin, step=step)
        data = dpdata.lammps.dump.system_data(
            lines, type_map, unwrap=unwrap, input_file=input_file
        )
        register_spin(data)
        return data
