from __future__ import annotations

import io
from typing import TYPE_CHECKING

import numpy as np

from dpdata.format import Format
from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType
from dpdata.formats.xyz.quip_gap_xyz import QuipGapxyzSystems, format_single_frame
from dpdata.formats.xyz.xyz import coord_to_xyz, xyz_to_coord


@Format.register("xyz")
class XYZFormat(Format):
    """Plain XYZ molecular structure file.

    Plain XYZ stores element symbols and Cartesian coordinates but no cell or
    labels. dpdata therefore treats it as nonperiodic and assigns a placeholder
    cell. Use ``extxyz`` when energies, forces, virials, or multiple chemical
    formulas must be preserved.

    Examples
    --------
    >>> import dpdata
    >>> system = dpdata.System("POSCAR", fmt="vasp/poscar")
    >>> system.to("xyz", "a.xyz")
    """

    def to_system(self, data, file_name: FileType, **kwargs):
        """Write all frames as concatenated plain XYZ records.

        Parameters
        ----------
        data : dict
            System data. Cell and label fields are not written.
        file_name : str or os.PathLike or file-like object
            Destination XYZ file.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.
        """
        buff = []
        types = np.array(data["atom_names"])[data["atom_types"]]
        for cc in data["coords"]:
            buff.append(coord_to_xyz(cc, types))
        with open_file(file_name, "w") as fp:
            fp.write("\n".join(buff))

    def from_system(self, file_name: FileType, **kwargs):
        """Load the first structure from a plain XYZ file.

        Parameters
        ----------
        file_name : str or os.PathLike or file-like object
            Input XYZ file.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Nonperiodic System data with a placeholder cell.
        """
        with open_file(file_name) as fp:
            coords, types = xyz_to_coord(fp.read())
        atom_names, atom_types, atom_numbs = np.unique(
            types, return_inverse=True, return_counts=True
        )
        return {
            "atom_names": list(atom_names),
            "atom_numbs": list(atom_numbs),
            "atom_types": atom_types,
            "coords": coords.reshape((1, *coords.shape)),
            "cells": np.eye(3).reshape((1, 3, 3)) * 100,
            "nopbc": True,
            "orig": np.zeros(3),
        }


@Format.register("quip/gap/xyz")
@Format.register("quip/gap/xyz_file")
@Format.register("extxyz")
@Format.register("gpumd/xyz")
@Format.register("nequip/xyz")
@Format.register("mace/xyz")
class QuipGapXYZFormat(Format):
    """Extended XYZ used by QUIP/GAP and atomistic ML tools.

    `QUIP/GAP <https://libatoms.github.io/QUIP/>`_ provides a Gaussian
    Approximation Potential framework, while `MACE
    <https://github.com/ACEsuit/mace>`_, `NequIP
    <https://github.com/mir-group/nequip>`_, and `GPUMD
    <https://github.com/brucefan1983/GPUMD>`_ are modern machine-learning
    interatomic potential packages.

    The comment-line ``Lattice`` and ``Properties`` metadata can store cells,
    energies, forces, virials, and per-atom fields. A single file may contain
    multiple frames and formulas, so the format supports
    :class:`dpdata.MultiSystems`. The aliases ``extxyz``, ``mace/xyz``,
    ``nequip/xyz``, and ``gpumd/xyz`` share this implementation.
    """

    def from_labeled_system(self, data, **kwargs):
        """Load the first labeled frame from an extended XYZ source.

        Parameters
        ----------
        data : str, os.PathLike, or dict
            Input extended XYZ file, or an already parsed frame supplied by
            :meth:`from_multi_systems`.
        **kwargs : dict
            Extended-XYZ parsing options described below.

        Other Parameters
        ----------------
        stress_sign : int, default=-1
            Sign in ``virial = stress_sign * volume * stress``. The default
            follows ASE's ``virial = -V * stress`` convention.

        Returns
        -------
        dict
            Labeled data for the first frame.
        """
        # When called via from_multi_systems iteration, data is already
        # a parsed info_dict — return as-is.
        if isinstance(data, dict):
            return data
        # When called directly with a filename, read the first frame.
        file_name = data
        for frame in QuipGapxyzSystems(file_name, **kwargs):
            return frame
        raise RuntimeError(f"No frames found in {file_name}")

    def from_multi_systems(self, file_name, **kwargs):
        """Iterate over all frames and formulas in an extended XYZ file.

        Parameters
        ----------
        file_name : str or os.PathLike
            Input extended XYZ file.
        **kwargs : dict
            Extended-XYZ parsing options described below.

        Other Parameters
        ----------------
        stress_sign : int, default=-1
            Sign in ``virial = stress_sign * volume * stress``. The default
            follows ASE's ``virial = -V * stress`` convention.

        Returns
        -------
        collections.abc.Iterable[dict]
            Parsed labeled frame dictionaries.
        """
        # here directory is the file_name
        return QuipGapxyzSystems(file_name, **kwargs)

    def to_labeled_system(self, data, file_name: FileType, **kwargs):
        """Write LabeledSystem data to QUIP/GAP XYZ format file.

        Parameters
        ----------
        data : dict
            system data
        file_name : FileType
            output file name or file handler
        **kwargs : dict
            additional arguments
        """
        frames = []
        nframes = len(data["energies"])

        for frame_idx in range(nframes):
            frame_lines = format_single_frame(data, frame_idx)
            frames.append("\n".join(frame_lines))

        content = "\n".join(frames)

        if isinstance(file_name, io.IOBase):
            file_name.write(content)
            if not content.endswith("\n"):
                file_name.write("\n")
        else:
            with open_file(file_name, "w") as fp:
                fp.write(content)

    def to_multi_systems(self, formulas, directory, **kwargs):
        """Return single filename for all systems in QUIP/GAP XYZ format.

        For QUIP/GAP XYZ format, all systems are written to a single file.

        Parameters
        ----------
        formulas : list[str]
            list of system names/formulas
        directory : str
            output filename
        **kwargs : dict
            additional arguments

        Yields
        ------
        file handler
            file handler for all systems
        """
        with open_file(directory, "w") as f:
            # Just create/truncate the file, then yield file handlers
            for _ in formulas:
                yield f
