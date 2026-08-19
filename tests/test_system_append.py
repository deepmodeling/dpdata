from __future__ import annotations

import unittest

import numpy as np
from comp_sys import CompLabeledSys, IsNoPBC, IsPBC
from context import dpdata


class TestFailedAppend(unittest.TestCase):
    def test_failed_append(self):
        sys1 = dpdata.System("poscars/POSCAR.h2o.md", fmt="vasp/poscar")
        sys2 = dpdata.System("poscars/POSCAR.h4o3", fmt="vasp/poscar")
        with self.assertRaises(Exception) as c:
            sys1.append(sys2)
        self.assertTrue(
            "systems with inconsistent formula could not be append" in str(c.exception)
        )


class TestAppendOwnership(unittest.TestCase):
    """Regression tests for copy-on-slice and copy-on-first-append semantics."""

    def test_sub_system_does_not_alias_metadata(self):
        system = dpdata.System(
            data={
                "atom_names": ["H"],
                "atom_numbs": [1],
                "atom_types": np.array([0]),
                "orig": np.zeros(3),
                "cells": np.eye(3).reshape(1, 3, 3),
                "coords": np.zeros((1, 1, 3)),
            }
        )
        sub = system[0:1]
        sub.data["atom_names"][0] = "X"
        sub.data["atom_types"][0] = 1
        sub.data["coords"][0, 0, 0] = 123.0
        sub.data["cells"][0, 0, 0] = 456.0
        self.assertEqual(system.data["atom_names"], ["H"])
        np.testing.assert_array_equal(system.data["atom_types"], [0])
        self.assertEqual(system.data["coords"][0, 0, 0], 0.0)
        self.assertEqual(system.data["cells"][0, 0, 0], 1.0)

    def test_first_append_does_not_alias_source(self):
        source = dpdata.System("poscars/POSCAR.oh.d", fmt="vasp/poscar")
        target = dpdata.System()
        target.append(source)

        source.data["atom_names"][0] = "X"
        source.data["coords"][0, 0, 0] = 123.0
        self.assertEqual(target.data["atom_names"][0], "O")
        self.assertNotEqual(target.data["coords"][0, 0, 0], 123.0)


class TestVaspXmlAppend(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.places = 6
        # rotated vasp computation, subject to numerical error
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6
        begin = 2
        end = 10
        step = 3
        self.system_1 = dpdata.LabeledSystem("poscars/vasprun.h2o.md.10.xml")
        self.system_2 = dpdata.LabeledSystem("poscars/vasprun.h2o.md.10.xml")
        self.system_1.append(self.system_2)

        self.system_1 = self.system_1.sub_system([0, 12, 4, 16, 8])
        self.system_2 = dpdata.LabeledSystem(
            "poscars/vasprun.h2o.md.10.xml"
        ).sub_system(np.arange(0, 10, 2))


class TestDifferentOrderAppend(unittest.TestCase, CompLabeledSys, IsNoPBC):
    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        self.system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_reordered.gaussianlog", fmt="gaussian/log"
        )
        self.system_1.append(system_2)

        self.system_2 = self.system_1.sub_system([0, 0])


if __name__ == "__main__":
    unittest.main()
