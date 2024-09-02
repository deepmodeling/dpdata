from __future__ import annotations

import dpdata.lammps.dump
import dpdata.lammps.lmp
from dpdata.format import Format
from dpdata.utils import FileType, open_file


@Format.register("lmp")
@Format.register("lammps/lmp")
class LAMMPSLmpFormat(Format):
    @Format.post("shift_orig_zero")
    def from_system(self, file_name: FileType, type_map=None, **kwargs):
        with open_file(file_name) as fp:
            lines = [line.rstrip("\n") for line in fp]
        return dpdata.lammps.lmp.to_system_data(lines, type_map)

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
        self, file_name, type_map=None, begin=0, step=1, unwrap=False, **kwargs
    ):
        lines = dpdata.lammps.dump.load_file(file_name, begin=begin, step=step)
        return dpdata.lammps.dump.system_data(lines, type_map, unwrap=unwrap)
