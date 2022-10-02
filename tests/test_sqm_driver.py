import unittest
import shutil

import numpy as np
from context import dpdata
from comp_sys import CompSys, IsNoPBC


@unittest.skipIf(shutil.which("sqm") is None, "sqm is not installed")
class TestSQMdriver(unittest.TestCase, CompSys, IsNoPBC):
    """Test sqm with a hydrogen ion."""
    @classmethod
    def setUpClass(cls):
        cls.system_1 = dpdata.System(data={
            "atom_names": ["H"],
            "atom_numbs": [1],
            "atom_types": np.zeros((1,), dtype=int),
            "coords": np.zeros((1, 1, 3), dtype=np.float32),
            "cells": np.zeros((1, 3, 3), dtype=np.float32),
            "orig": np.zeros(3, dtype=np.float32),
            "nopbc": True,
        })
        cls.system_2 = cls.system_1.predict(theory="DFTB3", charge=1, driver="sqm")
        cls.places = 6
    
    def test_energy(self):
        self.assertAlmostEqual(self.system_2['energies'].ravel()[0], 6.549447)
    
    def test_forces(self):
        forces = self.system_2['forces']
        np.testing.assert_allclose(forces, np.zeros_like(forces))
