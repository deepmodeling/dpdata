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


def register_charge(data: dict) -> None:
    if "charges" in data:
        dt = DataType(
            "charges",
            np.ndarray,
            (Axis.NFRAMES, Axis.NATOMS),
            required=False,
            deepmd_name="charge",
        )
        dpdata.System.register_data_type(dt)


@Format.register("lmp")
@Format.register("lammps/lmp")
class LAMMPSLmpFormat(Format):
    @Format.post("shift_orig_zero")
    def from_system(
        self, file_name: FileType, type_map=None, atom_style="auto", **kwargs
    ):
        """Load LAMMPS data file to system data format.

        This method supports multiple LAMMPS atom styles with automatic charge extraction
        and maintains backward compatibility. The parser can automatically detect the atom
        style from the LAMMPS data file header when possible.

        Parameters
        ----------
        file_name : str or Path
            Path to LAMMPS data file
        type_map : list, optional
            Mapping from atom types to element names
        atom_style : str, optional
            The LAMMPS atom style. Default is "auto" which attempts to detect
            the style automatically from the file. Can also be explicitly set to:
            atomic, full, charge, bond, angle, molecular, dipole, sphere
        **kwargs : dict
            Other parameters

        Returns
        -------
        dict
            System data dictionary with additional data based on atom style:
            - charges: For styles with charge information (full, charge, dipole)
            - molecule_ids: For styles with molecule information (full, bond, angle, molecular)

        Examples
        --------
        Load LAMMPS data with automatic detection:

        >>> system = dpdata.System("data.lmp", type_map=["O", "H"])

        Load with specific atom styles:

        >>> # Full style with charges and molecule IDs
        >>> system = dpdata.System("data.lmp", type_map=["O", "H"], atom_style="full")
        >>> print(system["charges"])  # Access extracted charges

        >>> # Charge style with charges only
        >>> system = dpdata.System("data.lmp", type_map=["O", "H"], atom_style="charge")

        >>> # Bond/molecular styles with molecule IDs
        >>> system = dpdata.System("data.lmp", type_map=["O", "H"], atom_style="bond")

        Notes
        -----
        Atom Style Column Layouts:
        - atomic: atom-ID atom-type x y z (default)
        - full: atom-ID molecule-ID atom-type charge x y z
        - charge: atom-ID atom-type charge x y z
        - bond: atom-ID molecule-ID atom-type x y z
        - angle: atom-ID molecule-ID atom-type x y z
        - molecular: atom-ID molecule-ID atom-type x y z
        - dipole: atom-ID atom-type charge x y z mux muy muz
        - sphere: atom-ID atom-type diameter density x y z
        """
        with open_file(file_name) as fp:
            lines = [line.rstrip("\n") for line in fp]
        data = dpdata.lammps.lmp.to_system_data(lines, type_map, atom_style=atom_style)
        register_spin(data)
        register_charge(data)
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
