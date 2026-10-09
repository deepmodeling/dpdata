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

ATOMIC_POSITIONS (bohr)
H 1.0 1.0 1.0
H 2.0 2.0 2.0
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

        # The original issue used the existing qe/pw/scf entry point for a
        # pw.x AIMD output. It should now retain every labeled frame too.
        scf_system = dpdata.LabeledSystem(self.prefix + ".out", fmt="qe/pw/scf")
        self.assertEqual(scf_system.get_nframes(), 2)
        np.testing.assert_allclose(scf_system["coords"][0, 0], [0.0, 0.0, 0.0])
        self.assertAlmostEqual(
            scf_system["energies"][0], -1.0 * dpdata.formats.qe.scf.ry2ev
        )

    def test_variable_cell_units(self):
        output = self.prefix + ".vc.out"
        with open(self.prefix + ".vc.in", "w") as fp:
            with open(self.prefix + ".in") as source:
                fp.write(source.read())
        with open(output, "w") as fp:
            fp.write(
                """     Program PWSCF
!    total energy = -1.00000000 Ry
     Forces acting on atoms (Ry/au):
     atom 1 type 1 force = 0.0 0.0 0.0
     atom 2 type 1 force = 0.0 0.0 0.0

CELL_PARAMETERS (bohr)
 5.0 0.0 0.0
 0.0 5.0 0.0
 0.0 0.0 5.0
ATOMIC_POSITIONS (crystal)
H 0.0 0.0 0.0
H 0.5 0.5 0.5

!    total energy = -2.00000000 Ry
     Forces acting on atoms (Ry/au):
     atom 1 type 1 force = 0.0 0.0 0.0
     atom 2 type 1 force = 0.0 0.0 0.0

CELL_PARAMETERS (bohr)
 6.0 0.0 0.0
 0.0 6.0 0.0
 0.0 0.0 6.0
ATOMIC_POSITIONS (crystal)
H 0.0 0.0 0.0
H 0.25 0.25 0.25
"""
            )
        system = dpdata.LabeledSystem(output, fmt="qe/pw/md", begin=1)
        cell = 5.0 * dpdata.formats.qe.scf.bohr2ang
        np.testing.assert_allclose(system["cells"][0], np.eye(3) * cell)
        np.testing.assert_allclose(system["coords"][0, 1], [cell / 2] * 3)

    def test_incomplete_output_positions_are_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "incomplete ATOMIC_POSITIONS"):
            dpdata.formats.qe.scf._output_position_frames(
                ["ATOMIC_POSITIONS (angstrom)", "H 0.0 0.0 0.0", ""],
                natoms=2,
                cell=np.eye(3),
            )
