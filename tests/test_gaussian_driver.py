import unittest
import shutil
import importlib
import os

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


class TestMakeGaussian(unittest.TestCase):
    """This class will not check if the output is correct, but only see if there is any errors."""
    def setUp(self):
        self.system = dpdata.System(data={
            "atom_names": ["H"],
            "atom_numbs": [1],
            "atom_types": np.zeros((1,), dtype=int),
            "coords": np.zeros((1, 1, 3), dtype=np.float32),
            "cells": np.zeros((1, 3, 3), dtype=np.float32),
            "orig": np.zeros(3, dtype=np.float32),
            "nopbc": True,
        })

    @unittest.skipIf(importlib.util.find_spec("openbabel") is None, "requires openbabel")
    def test_make_fp_gaussian(self):
        self.system.to_gaussian_gjf("gaussian/tmp.gjf", keywords="wb97x/6-31g* force")

    def test_make_fp_gaussian_multiplicity_one(self):
        self.system.to_gaussian_gjf("gaussian/tmp.gjf", keywords="wb97x/6-31g* force", multiplicity=1)

    def test_detect_multiplicity(self):
        # oxygen O2 3
        self._check_multiplicity(['O', 'O'], 3)
        # methane CH4 1
        self._check_multiplicity(['C', 'H', 'H', 'H', 'H'], 1)
        # CH3 2
        self._check_multiplicity(['C', 'H', 'H', 'H'], 2)
        # CH2 1
        self._check_multiplicity(['C', 'H', 'H'], 1)
        # CH 2
        self._check_multiplicity(['C', 'H'], 2)

    def _check_multiplicity(self, symbols, multiplicity):
        self.assertEqual(dpdata.gaussian.gjf.detect_multiplicity(np.array(symbols)), multiplicity)

    def tearDown(self):
        if os.path.exists('gaussian/tmp.gjf'):
            os.remove('gaussian/tmp.gjf')


class TestDumpGaussianGjf(unittest.TestCase):
    def setUp(self):
        self.system = dpdata.LabeledSystem('gaussian/methane.gaussianlog',
                                           fmt='gaussian/log')

    def test_dump_to_gjf(self):
        self.system.to_gaussian_gjf("gaussian/tmp.gjf", keywords="force B3LYP/6-31G(d)", multiplicity=1)
        with open("gaussian/tmp.gjf") as f:
            f.readline()
            header = f.readline().strip()
            f.readline()
            title = f.readline().strip()
            f.readline()
            charge, mult = (int(x) for x in f.readline().strip().split())
            atoms = []
            coords = []
            for ii in range(5):
                line = f.readline().strip().split()
                atoms.append(line[0])
                coords.append([float(x) for x in line[1:]])

        self.assertEqual(header, "#force B3LYP/6-31G(d)")
        self.assertEqual(title, self.system.formula)
        self.assertEqual(charge, 0)
        self.assertEqual(mult, 1)
        self.assertEqual(atoms, ['C', 'H', 'H', 'H', 'H'])
        for i in range(self.system['coords'].shape[1]):
            for j in range(3):
                self.assertAlmostEqual(coords[i][j], self.system['coords'][0][i][j])

    def tearDown(self):
        if os.path.exists('gaussian/tmp.gjf'):
            os.remove('gaussian/tmp.gjf')
