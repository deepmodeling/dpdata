from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np

import dpdata.formats.abacus.md
import dpdata.formats.abacus.relax
import dpdata.formats.abacus.scf
from dpdata.data_type import Axis, DataType
from dpdata.format import Format
from dpdata.formats.abacus.stru import get_frame_from_stru, make_unlabeled_stru
from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType


@Format.register("abacus/stru")
@Format.register("stru")
class AbacusSTRUFormat(Format):
    """ABACUS structure file.

    `ABACUS <https://abacus.ustc.edu.cn/>`_ (Atomic-orbital Based Ab-initio
    Computation at UStc) is an open-source DFT package based on LCAO and
    plane-wave basis sets.

    ``STRU`` stores the cell, species, coordinates, pseudopotential/orbital
    references, and optional movement or magnetic-moment fields for one
    ABACUS configuration. This format reads and writes unlabeled
    :class:`dpdata.System` objects.
    """

    def from_system(self, file_name, **kwargs):
        """Load one ABACUS ``STRU`` file.

        Parameters
        ----------
        file_name : str or os.PathLike
            Input ``STRU`` file.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            System data, including optional ``move`` and magnetic fields when
            present in the file.
        """
        data = get_frame_from_stru(file_name)
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
            Additional STRU fields described below.

        Other Parameters
        ----------------
        pp_file : list[str] or dict[str, str], optional
            Pseudopotential file for each atom type.
        numerical_orbital : list[str] or dict[str, str], optional
            Numerical orbital file for each atom type.
        numerical_descriptor : str, optional
            Numerical descriptor file used by ABACUS.
        mass : list[float], optional
            Atomic mass for each atom type.
        move : array-like, optional
            Per-frame, per-atom Cartesian movement flags.
        velocity : array-like, optional
            Initial Cartesian velocity for each atom.
        mag : array-like, optional
            Scalar or vector magnetic moment for each atom.
        angle1, angle2 : array-like, optional
            Polar and azimuthal magnetic-moment angles for noncollinear spins.
        sc : array-like, optional
            Spin-constraint flags.
        lambda_ : array-like, optional
            Spin-constraint lambda values.
        link_file : bool, default=False
            Write basenames and create symbolic links for referenced files.
        """
        stru_string = make_unlabeled_stru(
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
    """ABACUS self-consistent-field calculation directory.

    `ABACUS <https://abacus.ustc.edu.cn/>`_ is an open-source DFT package.
    The reader combines the calculation's ``INPUT`` and ``STRU`` files with
    the corresponding ``OUT.<suffix>/running_scf.log`` output and returns the
    final labeled configuration.
    """

    # @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, **kwargs):
        """Load an ABACUS SCF calculation.

        Parameters
        ----------
        file_name : str or os.PathLike
            Calculation directory containing ``INPUT``, ``STRU``, and the
            ABACUS output directory.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled system data with energy and any available forces, virial,
            movement flags, or magnetic fields.
        """
        data = dpdata.formats.abacus.scf.get_frame(file_name)
        register_mag_data(data)
        register_move_data(data)
        return data


@Format.register("abacus/md")
@Format.register("abacus/pw/md")
@Format.register("abacus/lcao/md")
class AbacusMDFormat(Format):
    """ABACUS molecular-dynamics calculation directory.

    `ABACUS <https://abacus.ustc.edu.cn/>`_ is an open-source DFT package.
    This format reads the trajectory and labels emitted by an ABACUS MD run,
    including optional force, virial, movement, and magnetic data.
    """

    # @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, **kwargs):
        """Load frames from an ABACUS MD calculation.

        Parameters
        ----------
        file_name : str or os.PathLike
            Calculation directory containing ``INPUT``, ``STRU``, and
            ``OUT.<suffix>`` MD output files.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled trajectory data.
        """
        data = dpdata.formats.abacus.md.get_frame(file_name)
        register_mag_data(data)
        register_move_data(data)
        return data


@Format.register("abacus/relax")
@Format.register("abacus/pw/relax")
@Format.register("abacus/lcao/relax")
class AbacusRelaxFormat(Format):
    """ABACUS ionic- or cell-relaxation calculation directory.

    `ABACUS <https://abacus.ustc.edu.cn/>`_ is an open-source DFT package.
    The reader reconstructs relaxation frames from the ABACUS log and saved
    ``STRU_ION*_D`` structures and attaches the available energies, forces,
    virials, movement flags, and magnetic data.
    """

    # @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, **kwargs):
        """Load an ABACUS relaxation trajectory.

        Parameters
        ----------
        file_name : str or os.PathLike
            Calculation directory containing the ABACUS input and relaxation
            output files.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled relaxation frames.
        """
        data = dpdata.formats.abacus.relax.get_frame(file_name)
        register_mag_data(data)
        register_move_data(data)
        return data
