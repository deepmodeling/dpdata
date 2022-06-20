import unittest
import shutil
import importlib

import numpy as np
from context import dpdata
from comp_sys import CompSys, IsNoPBC


@unittest.skipIf(shutil.which("g16") is None, "g16 is not installed")
@unittest.skipIf(importlib.util.find_spec("openbabel") is None, "openbabel is not installed")
class TestGaussianDriver(unittest.TestCase, CompSys, IsNoPBC):
    """Test Gaussian with a hydrogen ion."""
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
        cls.system_2 = cls.system_1.predict(keywords="force B3LYP", charge=1, driver="gaussian")
        cls.places = 6
    
    def test_energy(self):
        self.assertAlmostEqual(self.system_2['energies'].ravel()[0], 0.)
    
    def test_forces(self):
        forces = self.system_2['forces']
        np.testing.assert_allclose(forces, np.zeros_like(forces))
