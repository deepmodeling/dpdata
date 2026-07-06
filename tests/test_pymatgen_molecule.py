from __future__ import annotations

import os
import unittest

import numpy as np
from context import dpdata

try:
    import pymatgen  # noqa: F401
except ModuleNotFoundError:
    skip_pymatgen = True
else:
    skip_pymatgen = False


@unittest.skipIf(skip_pymatgen, "skip pymatgen related test. install pymatgen to fix")
class TestPOSCARCart(unittest.TestCase):
    def setUp(self):
        self.system = dpdata.System()
        self.system.from_pymatgen_molecule(os.path.join("pymatgen_data", "FA-001.xyz"))
        self.assertEqual(list(self.system["atom_types"]), [0, 1, 2, 1, 1, 2, 1, 1])

    def test_poscar_to_molecule(self):
        tmp_system = dpdata.System()
        tmp_system.from_vasp_poscar(os.path.join("pymatgen_data", "mol2.vasp"))
        natoms = len(tmp_system["coords"][0])
        tmpcoord = tmp_system["coords"][0]
        cog = np.average(tmpcoord, axis=0)
        dist = tmpcoord - np.tile(cog, [natoms, 1])
        max_dist_0 = np.max(np.linalg.norm(dist, axis=1))

        mols = tmp_system.to("pymatgen/molecule")
        cog = np.average(mols[-1].cart_coords, axis=0)
        dist = mols[-1].cart_coords - np.tile(cog, [natoms, 1])
        max_dist_1 = np.max(np.linalg.norm(dist, axis=1))
        self.assertAlmostEqual(max_dist_0, max_dist_1)

    def test_ungrouped_atom_types_species_match_coords(self):
        # Regression for gh-994: species must follow atom_types order, not the
        # grouped atom_names/atom_numbs order, otherwise sites are assigned the
        # wrong element when atoms are not grouped by type.
        system = dpdata.System(
            data={
                "atom_names": ["H", "O"],
                "atom_numbs": [2, 1],
                "atom_types": np.array([0, 1, 0]),
                "orig": np.zeros(3),
                "cells": np.eye(3).reshape(1, 3, 3) * 20,
                "coords": np.array(
                    [[[0.0, 0.0, 0.0], [9.0, 0.0, 0.0], [1.0, 0.0, 0.0]]]
                ),
            }
        )
        mol = system.to("pymatgen/molecule")[0]
        self.assertEqual([str(site.specie) for site in mol], ["H", "O", "H"])
        for site, coord in zip(mol, system["coords"][0]):
            np.testing.assert_allclose(site.coords, coord)


if __name__ == "__main__":
    unittest.main()
