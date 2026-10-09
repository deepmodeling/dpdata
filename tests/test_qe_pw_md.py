from __future__ import annotations

import os
import tempfile
import unittest

import numpy as np
from context import dpdata


class TestQEPWMD(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.prefix = os.path.join(self.tmpdir.name, "md")
        with open(self.prefix + ".in", "w") as fp:
            fp.write(
                """&CONTROL
 calculation = 'md'
/
&SYSTEM
 ibrav = 0, nat = 2, ntyp = 1
/
CELL_PARAMETERS (angstrom)
 4.0 0.0 0.0
 0.0 4.0 0.0
 0.0 0.0 4.0
ATOMIC_SPECIES
H 1.0 H.UPF
ATOMIC_POSITIONS (angstrom)
H 0.0 0.0 0.0
H 1.0 1.0 1.0
"""
            )
        with open(self.prefix + ".out", "w") as fp:
            fp.write(
                """     Program PWSCF
ATOMIC_POSITIONS (angstrom)
H 0.0 0.0 0.0
H 1.0 1.0 1.0

!    total energy              =   -1.00000000 Ry
     Forces acting on atoms (Ry/au):
     atom 1 type 1 force = 0.10000000 0.00000000 0.00000000
     atom 2 type 1 force = 0.00000000 0.10000000 0.00000000

ATOMIC_POSITIONS (bohr)
H 0.5 0.5 0.5
H 1.5 1.5 1.5

!    total energy              =   -2.00000000 Ry
     Forces acting on atoms (Ry/au):
     atom 1 type 1 force = 0.20000000 0.00000000 0.00000000
     atom 2 type 1 force = 0.00000000 0.20000000 0.00000000
"""
            )

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_pw_md_frames_and_sampling(self):
        system = dpdata.LabeledSystem(
            self.prefix + ".out", fmt="qe/pw/md", begin=1, step=1
        )
        self.assertEqual(system.get_nframes(), 1)
        self.assertEqual(system["coords"].shape, (1, 2, 3))
        self.assertAlmostEqual(
            system["energies"][0], -2.0 * dpdata.formats.qe.scf.ry2ev
        )
        np.testing.assert_allclose(
            system["coords"][0, 0], [0.5 * dpdata.formats.qe.scf.bohr2ang] * 3
        )
        self.assertNotIn("virials", system.data)
        npy_dir = os.path.join(self.tmpdir.name, "npy")
        system.to("deepmd/npy", npy_dir)
        self.assertTrue(os.path.isfile(os.path.join(npy_dir, "set.000", "coord.npy")))

    def test_variable_cell_units(self):
        output = self.prefix + ".vc.out"
        with open(self.prefix + ".vc.in", "w") as fp:
            fp.write(open(self.prefix + ".in").read())
        with open(output, "w") as fp:
            fp.write(
                """     Program PWSCF
CELL_PARAMETERS (bohr)
 5.0 0.0 0.0
 0.0 5.0 0.0
 0.0 0.0 5.0
ATOMIC_POSITIONS (crystal)
H 0.0 0.0 0.0
H 0.5 0.5 0.5

!    total energy = -1.00000000 Ry
     Forces acting on atoms (Ry/au):
     atom 1 type 1 force = 0.0 0.0 0.0
     atom 2 type 1 force = 0.0 0.0 0.0
"""
            )
        system = dpdata.LabeledSystem(output, fmt="qe/pw/md")
        cell = 5.0 * dpdata.formats.qe.scf.bohr2ang
        np.testing.assert_allclose(system["cells"][0], np.eye(3) * cell)
        np.testing.assert_allclose(system["coords"][0, 1], [cell / 2] * 3)
