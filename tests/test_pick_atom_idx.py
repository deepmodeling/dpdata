from __future__ import annotations

import unittest

import numpy as np
from comp_sys import CompSys, IsNoPBC
from context import dpdata

try:
    exist_module = True
except Exception:
    exist_module = False


class TestPickAtomIdx(unittest.TestCase, CompSys, IsNoPBC):
    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6
        self.system_1 = dpdata.LabeledSystem(
            "gaussian/methane_reordered.gaussianlog", fmt="gaussian/log"
        ).pick_atom_idx(slice(4))
        self.system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )


class TestPickAtomOwnership(unittest.TestCase):
    def test_picked_atoms_are_independent(self):
        for idx in (slice(0, 2), slice(None, None, 2), np.array([0, 2])):
            with self.subTest(idx=idx):
                source = dpdata.LabeledSystem(
                    data={
                        "atom_names": ["H", "O"],
                        "atom_numbs": [2, 2],
                        "atom_types": np.array([0, 1, 0, 1]),
                        "orig": np.zeros(3),
                        "cells": np.tile(np.eye(3), (2, 1, 1)),
                        "coords": np.arange(24, dtype=float).reshape(2, 4, 3),
                        "energies": np.array([1.0, 2.0]),
                        "forces": np.arange(24, dtype=float).reshape(2, 4, 3),
                    }
                )
                picked = source.pick_atom_idx(idx)
                picked.check_data()
                names = ("atom_types", "coords", "forces", "cells", "energies")
                original = {name: source[name].copy() for name in names}
                for name in names:
                    picked[name].flat[0] += 10
                    np.testing.assert_array_equal(source[name], original[name])
                picked["atom_names"][0] = "He"
                self.assertEqual(source["atom_names"], ["H", "O"])

                selected = {name: picked[name].copy() for name in names}
                for name in names:
                    source[name].flat[0] -= 10
                    np.testing.assert_array_equal(picked[name], selected[name])
                source["atom_names"][1] = "Ne"
                self.assertEqual(picked["atom_names"], ["He", "O"])


@unittest.skipIf(not exist_module, "skip")
class TestPickByAmberMask(unittest.TestCase, CompSys, IsNoPBC):
    def setUp(self):
        parmfile = "amber/corr/qmmm.parm7"
        ep = r"@%EP"
        target = ":1"
        cutoff = 6.0
        interactwith = f"({target})<:{cutoff:f}&!{ep}"
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6
        self.system_1 = dpdata.LabeledSystem(
            "amber/corr/dp_corr", fmt="deepmd/npy"
        ).pick_by_amber_mask(parmfile, interactwith, pass_coords=True, nopbc=True)[
            "C6EP0H11HW192O6OW96P1"
        ]
        self.system_2 = dpdata.LabeledSystem(
            "amber/corr/dp_amber_mask/C6EP0H11HW192O6OW96P1", fmt="deepmd/npy"
        )


if __name__ == "__main__":
    unittest.main()
