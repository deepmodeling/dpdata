import numpy as np

from dpdata.format import Format
from dpdata.xyz.quip_gap_xyz import QuipGapxyzSystems
from dpdata.xyz.xyz import coord_to_xyz, xyz_to_coord


@Format.register("xyz")
class XYZFormat(Format):
    """XYZ foramt.

    Examples
    --------
    >>> s.to("xyz", "a.xyz")
    """

    def to_system(self, data, file_name, **kwargs):
        buff = []
        types = np.array(data["atom_names"])[data["atom_types"]]
        for cc in data["coords"]:
            buff.append(coord_to_xyz(cc, types))
        with open(file_name, "w") as fp:
            fp.write("\n".join(buff))

    def from_system(self, file_name, **kwargs):
        with open(file_name) as fp:
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
class QuipGapXYZFormat(Format):
    def from_labeled_system(self, data, **kwargs):
        return data

    def from_multi_systems(self, file_name, **kwargs):
        # here directory is the file_name
        return QuipGapxyzSystems(file_name)
