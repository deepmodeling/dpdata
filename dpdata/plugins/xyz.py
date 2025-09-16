from __future__ import annotations

import io
from typing import TYPE_CHECKING

import numpy as np

from dpdata.format import Format
from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType
from dpdata.xyz.quip_gap_xyz import QuipGapxyzSystems, format_single_frame
from dpdata.xyz.xyz import coord_to_xyz, xyz_to_coord


@Format.register("xyz")
class XYZFormat(Format):
    """XYZ foramt.

    Examples
    --------
    >>> s.to("xyz", "a.xyz")
    """

    def to_system(self, data, file_name: FileType, **kwargs):
        buff = []
        types = np.array(data["atom_names"])[data["atom_types"]]
        for cc in data["coords"]:
            buff.append(coord_to_xyz(cc, types))
        with open_file(file_name, "w") as fp:
            fp.write("\n".join(buff))

    def from_system(self, file_name: FileType, **kwargs):
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
    def from_labeled_system(self, data, **kwargs):
        return data

    def from_multi_systems(self, file_name, **kwargs):
        # here directory is the file_name
        return QuipGapxyzSystems(file_name)

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
