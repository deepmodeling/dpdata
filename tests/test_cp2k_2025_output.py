from __future__ import annotations

import os
import tempfile
import unittest

import numpy as np

from comp_sys import CompLabeledSys
from context import dpdata


class TestCp2k2025Output(unittest.TestCase, CompLabeledSys):
    """Test CP2K 2025 output format parsing."""

    def setUp(self):
        self.system_1 = dpdata.LabeledSystem(
            "cp2k/cp2k_2025_output/cp2k_2025_output", fmt="cp2k/output"
        )
        self.system_2 = dpdata.LabeledSystem(
            "cp2k/cp2k_2025_output/deepmd", fmt="deepmd/npy"
        )
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

    def test_energy_extraction(self):
        """Test that energy is correctly extracted from CP2K 2025 format."""
        # Energy should be -7.364190264587725 hartree
        # Using dpdata's conversion factor: -200.3898256786414 eV
        expected_energy = -200.3898256786414
        self.assertAlmostEqual(
            self.system_1.data["energies"][0], expected_energy, places=5
        )

    def test_forces_extraction(self):
        """Test that forces are correctly extracted from CP2K 2025 format."""
        # Forces should be converted from hartree/bohr to eV/angstrom
        self.assertEqual(self.system_1.data["forces"].shape, (1, 2, 3))
        # Check first atom force x-component
        self.assertAlmostEqual(
            self.system_1.data["forces"][0][0][0], -2.94874881, places=5
        )


class TestCp2k2023FormatStillWorks(unittest.TestCase, CompLabeledSys):
    """Test that CP2K 2023 format still works (regression test)."""

    def setUp(self):
        self.system_1 = dpdata.LabeledSystem(
            "cp2k/cp2k_normal_output/cp2k_output", fmt="cp2k/output"
        )
        self.system_2 = dpdata.LabeledSystem(
            "cp2k/cp2k_normal_output/deepmd", fmt="deepmd/npy"
        )
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


class TestCp2k2025EdgeCases(unittest.TestCase):
    """Test edge cases for CP2K 2025 format parsing to improve coverage."""

    def create_cp2k_output_2025(
        self, energy_line=None, forces_lines=None
    ):
        """Create a minimal CP2K 2025 output file for testing."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".out", delete=False) as f:
            # Header required for parsing
            f.write(" DBCSR| Multiplication driver                                               XSMM\n")
            f.write("\n")
            f.write("  **** **** ******  **  PROGRAM STARTED AT               2025-01-15 10:30:00.000\n")
            f.write(" CP2K| version string:                                          CP2K version 2025.1\n")
            f.write("\n")
            f.write(" CELL| Vector a [angstrom]:       5.000     0.000     0.000    |a| =       5.000\n")
            f.write(" CELL| Vector b [angstrom]:       0.000     5.000     0.000    |b| =       5.000\n")
            f.write(" CELL| Vector c [angstrom]:       0.000     0.000     5.000    |c| =       5.000\n")
            f.write("\n")
            f.write(" ATOMIC KIND INFORMATION\n")
            f.write("\n")
            f.write("  1. Atomic kind: H                                     Number of atoms:       2\n")
            f.write("\n")
            f.write(" SCF WAVEFUNCTION OPTIMIZATION\n")
            f.write("\n")
            f.write(" *** SCF run converged in     1 steps ***\n")
            f.write("\n")
            f.write(" MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom\n")
            f.write("\n")
            f.write("  Atom  Kind  Element       X           Y           Z          Z(eff)       Mass\n")
            f.write("\n")
            f.write("       1     1 H    1    0.000000    0.000000    0.000000      1.00       1.0079\n")
            f.write("       2     1 H    1    0.760000    0.000000    0.000000      1.00       1.0079\n")
            f.write("\n")

            # Energy line - use provided or default
            if energy_line:
                f.write(energy_line + "\n")
            else:
                f.write(" ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] -7.364190264587725\n")

            f.write("\n")
            f.write(" FORCES| Atomic forces [hartree/bohr]\n")
            f.write("\n")
            f.write(" FORCES| Atom x y z |f|\n")

            # Forces lines - use provided or default
            if forces_lines:
                for line in forces_lines:
                    f.write(line + "\n")
            else:
                f.write(" FORCES| 1 -5.73440344E-02 2.95274914E-02 -1.50988167E-02 6.62433792E-02\n")
                f.write(" FORCES| 2 7.92269287E-02 3.84670665E-02 -3.41478833E-02 9.44600412E-02\n")

            f.write("\n")
            f.write(" STRESS TENSOR [GPa]\n")
            f.write("\n")
            f.write("            X               Y               Z\n")
            f.write("  X       0.12345678      0.00000000      0.00000000\n")
            f.write("  Y       0.00000000      0.12345678      0.00000000\n")
            f.write("  Z       0.00000000      0.00000000      0.12345678\n")
            f.write("\n")
            f.write("  **** **** ******  **  PROGRAM ENDED AT                 2025-01-15 10:30:05.000\n")

            return f.name

    def test_cp2k2025_format_with_labeled_system(self):
        """Test CP2K 2025 format using LabeledSystem (integration test for coverage)."""
        fname = self.create_cp2k_output_2025()
        try:
            system = dpdata.LabeledSystem(fname, fmt="cp2k/output")
            self.assertIsNotNone(system.data["energies"])
            self.assertIsNotNone(system.data["forces"])
            self.assertEqual(system.data["forces"].shape[1], 2)
        finally:
            os.unlink(fname)

    def test_cp2k2025_energy_parsing_with_extra_whitespace(self):
        """Test energy parsing with extra whitespace around value (coverage for parsing robustness)."""
        fname = self.create_cp2k_output_2025(
            energy_line=" ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]   -7.364190264587725  "
        )
        try:
            system = dpdata.LabeledSystem(fname, fmt="cp2k/output")
            self.assertIsNotNone(system.data["energies"])
            self.assertEqual(system.data["forces"].shape[1], 2)
        finally:
            os.unlink(fname)

    def test_cp2k2025_force_parsing_with_header_lines(self):
        """Test that FORCES| header lines are correctly skipped (coverage for filtering)."""
        fname = self.create_cp2k_output_2025(
            forces_lines=[
                " FORCES| Atom x y z |f|",  # Should be skipped - contains "Atom x y z"
                " FORCES| 1 -5.73440344E-02 2.95274914E-02 -1.50988167E-02 6.62433792E-02",
                " FORCES| 2 7.92269287E-02 3.84670665E-02 -3.41478833E-02 9.44600412E-02",
            ]
        )
        try:
            system = dpdata.LabeledSystem(fname, fmt="cp2k/output")
            self.assertEqual(system.data["forces"].shape[1], 2)
        finally:
            os.unlink(fname)

    def test_cp2k2025_force_parsing_with_atomic_forces_line(self):
        """Test that 'Atomic forces' line is correctly skipped (coverage for filtering)."""
        fname = self.create_cp2k_output_2025(
            forces_lines=[
                " FORCES| Atomic forces [hartree/bohr]",  # Should be skipped
                " FORCES| 1 -5.73440344E-02 2.95274914E-02 -1.50988167E-02 6.62433792E-02",
                " FORCES| 2 7.92269287E-02 3.84670665E-02 -3.41478833E-02 9.44600412E-02",
            ]
        )
        try:
            system = dpdata.LabeledSystem(fname, fmt="cp2k/output")
            self.assertEqual(system.data["forces"].shape[1], 2)
        finally:
            os.unlink(fname)


if __name__ == "__main__":
    unittest.main()
