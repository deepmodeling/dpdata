from __future__ import annotations

import os
import tempfile
import unittest

import dpdata


class TestQESingleLineEvp(unittest.TestCase):
    """Test for QE CP trajectory with single-line .evp file (issue #899)."""

    def setUp(self):
        # Create a temporary directory for test files
        self.test_dir = tempfile.mkdtemp()

        # Create a simple input file
        self.input_content = """&CONTROL
  calculation = 'cp'
  prefix = 'test'
/
&SYSTEM
  ibrav = 1
  celldm(1) = 10.0
  nat = 2
  ntyp = 1
/
&ELECTRONS
/
&IONS
/
ATOMIC_SPECIES
H 1.0 H.pbe-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS {alat}
H 0.0 0.0 0.0
H 0.5 0.5 0.5
"""

        # Create a single-line .evp file (this used to cause errors)
        self.evp_content = "# comments here\n    195  9.433649E-03  1.290500E-05  0.000000E+00  1.057182E-02      -1100.03076389      -1100.03076389      -1100.03075434      -1100.03073209  1.300606E+04      -3.31610"

        # Create a single-line .for file
        self.for_content = "    195\n0.0 0.0 0.0\n0.0 0.0 0.0"

        # Create a single-line .pos file
        self.pos_content = "    195\n0.0 0.0 0.0\n1.0 1.0 1.0"

        # Create a .cel file
        self.cel_content = "    195\n10.0 0.0 0.0\n0.0 10.0 0.0\n0.0 0.0 10.0"

        # Write files
        self.prefix = os.path.join(self.test_dir, "test")
        with open(self.prefix + ".in", "w") as f:
            f.write(self.input_content)
        with open(self.prefix + ".evp", "w") as f:
            f.write(self.evp_content)
        with open(self.prefix + ".for", "w") as f:
            f.write(self.for_content)
        with open(self.prefix + ".pos", "w") as f:
            f.write(self.pos_content)
        with open(self.prefix + ".cel", "w") as f:
            f.write(self.cel_content)

    def tearDown(self):
        # Clean up temporary files
        import shutil

        shutil.rmtree(self.test_dir)

    def test_single_line_evp_load(self):
        """Test that single-line .evp files can be loaded without errors."""
        # This should not raise an exception with the fix
        system = dpdata.LabeledSystem(self.prefix, fmt="qe/cp/traj")

        # Verify basic properties
        self.assertEqual(system.get_nframes(), 1)
        self.assertEqual(len(system["energies"]), 1)

        # The energy should be converted properly
        # -1100.03076389 Hartree * energy_convert factor
        # Let's check the actual value first and use a more reasonable tolerance
        self.assertAlmostEqual(system["energies"][0], -29933.361998680055, places=3)


if __name__ == "__main__":
    unittest.main()
