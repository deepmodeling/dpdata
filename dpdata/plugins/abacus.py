from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np

import dpdata.abacus.md
import dpdata.abacus.relax
import dpdata.abacus.scf
from dpdata.data_type import Axis, DataType
from dpdata.format import Format
from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType


@Format.register("abacus/stru")
@Format.register("stru")
class AbacusSTRUFormat(Format):
    def from_system(self, file_name, **kwargs):
        data = dpdata.abacus.scf.get_frame_from_stru(file_name)
        register_mag_data(data)
        return data

    def to_system(self, data, file_name: FileType, frame_idx=0, **kwargs):
        """Dump the system into ABACUS STRU format file.

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
        stru_string = dpdata.abacus.scf.make_unlabeled_stru(
            data=data,
            frame_idx=frame_idx,
            dest_dir=os.path.dirname(file_name),
            **kwargs,
        )
        with open_file(file_name, "w") as fp:
            fp.write(stru_string)


def register_mag_data(data):
    if "spins" in data:
        dt = DataType(
            "spins",
            np.ndarray,
            (Axis.NFRAMES, Axis.NATOMS, 3),
            required=False,
            deepmd_name="spin",
        )
        dpdata.System.register_data_type(dt)
        dpdata.LabeledSystem.register_data_type(dt)
    if "force_mags" in data:
        dt = DataType(
            "force_mags",
            np.ndarray,
            (Axis.NFRAMES, Axis.NATOMS, 3),
            required=False,
            deepmd_name="force_mag",
        )
        dpdata.System.register_data_type(dt)
        dpdata.LabeledSystem.register_data_type(dt)


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


@Format.register("abacus/scf")
@Format.register("abacus/pw/scf")
@Format.register("abacus/lcao/scf")
class AbacusSCFFormat(Format):
    # @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, **kwargs):
        data = dpdata.abacus.scf.get_frame(file_name)
        register_mag_data(data)
        register_move_data(data)
        return data


@Format.register("abacus/md")
@Format.register("abacus/pw/md")
@Format.register("abacus/lcao/md")
class AbacusMDFormat(Format):
    # @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, **kwargs):
        data = dpdata.abacus.md.get_frame(file_name)
        register_mag_data(data)
        register_move_data(data)
        return data


@Format.register("abacus/relax")
@Format.register("abacus/pw/relax")
@Format.register("abacus/lcao/relax")
class AbacusRelaxFormat(Format):
    # @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, **kwargs):
        data = dpdata.abacus.relax.get_frame(file_name)
        register_mag_data(data)
        register_move_data(data)
        return data
