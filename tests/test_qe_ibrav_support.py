from __future__ import annotations

import os
import tempfile
import unittest

import numpy as np
from context import dpdata


class TestQEIbravSupport(unittest.TestCase):
    """Test support for various ibrav values in Quantum Espresso files."""

    def setUp(self):
        # Create temporary directory for test files
        self.temp_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        # Clean up temporary files
        import shutil
        shutil.rmtree(self.temp_dir)

    def _create_qe_files(self, ibrav, params_str, expected_cell):
        """Create QE input and output files for testing."""
        input_content = f"""&control
calculation='scf'
/
&system
ibrav={ibrav}
{params_str}
nat=1
ntyp=1
/
ATOMIC_SPECIES
H 1.0 H.pbe.UPF
ATOMIC_POSITIONS {{angstrom}}
H 0.0 0.0 0.0
"""
        
        output_content = """&control
calculation='scf'
/

!    total energy              =     -10.0 Ry

Forces acting on atoms (cartesian axes, Ry/au):

     atom    1 type  1   force =     0.00000000    0.00000000    0.00000000
"""
        
        input_file = os.path.join(self.temp_dir, "test.in")
        output_file = os.path.join(self.temp_dir, "test.out")
        
        with open(input_file, 'w') as f:
            f.write(input_content)
        with open(output_file, 'w') as f:
            f.write(output_content)
            
        return input_file, output_file

    def test_ibrav_1_cubic(self):
        """Test ibrav=1 (simple cubic) - existing functionality."""
        _, output_file = self._create_qe_files(1, "a=10", np.eye(3) * 10)
        
        sys = dpdata.LabeledSystem(output_file, fmt='qe/pw/scf')
        expected_cell = np.array([[10, 0, 0], [0, 10, 0], [0, 0, 10]])
        np.testing.assert_allclose(sys.data['cells'][0], expected_cell, rtol=1e-10)

    def test_ibrav_8_simple_orthorhombic(self):
        """Test ibrav=8 (simple orthorhombic) - main issue case."""
        _, output_file = self._create_qe_files(8, "a=10\nb=12\nc=8", None)
        
        sys = dpdata.LabeledSystem(output_file, fmt='qe/pw/scf')
        expected_cell = np.array([[10, 0, 0], [0, 12, 0], [0, 0, 8]])
        np.testing.assert_allclose(sys.data['cells'][0], expected_cell, rtol=1e-10)

    def test_ibrav_8_with_celldm(self):
        """Test ibrav=8 with celldm parameters."""
        # celldm in atomic units (bohr)
        bohr2ang = dpdata.qe.scf.bohr2ang
        a_bohr = 10 / bohr2ang
        _, output_file = self._create_qe_files(8, f"celldm(1)={a_bohr}\ncelldm(2)=1.2\ncelldm(3)=0.8", None)
        
        sys = dpdata.LabeledSystem(output_file, fmt='qe/pw/scf')
        expected_cell = np.array([[10, 0, 0], [0, 12, 0], [0, 0, 8]])
        np.testing.assert_allclose(sys.data['cells'][0], expected_cell, rtol=1e-6)

    def test_ibrav_2_fcc(self):
        """Test ibrav=2 (face-centered cubic)."""
        _, output_file = self._create_qe_files(2, "a=10", None)
        
        sys = dpdata.LabeledSystem(output_file, fmt='qe/pw/scf')
        # After rot_lower_triangular transformation
        expected_cell = np.array([[7.071068e+00, 0.000000e+00, 0.000000e+00],
                                 [3.535534e+00, 6.123724e+00, 0.000000e+00],
                                 [3.535534e+00, 2.041241e+00, 5.773503e+00]])
        np.testing.assert_allclose(sys.data['cells'][0], expected_cell, atol=1e-6)

    def test_ibrav_3_bcc(self):
        """Test ibrav=3 (body-centered cubic)."""
        _, output_file = self._create_qe_files(3, "a=10", None)
        
        sys = dpdata.LabeledSystem(output_file, fmt='qe/pw/scf')
        # After rot_lower_triangular transformation
        expected_cell = np.array([[8.660254e+00, 0.000000e+00, 0.000000e+00],
                                 [2.886751e+00, 8.164966e+00, 0.000000e+00],
                                 [-2.886751e+00, 4.082483e+00, 7.071068e+00]])
        np.testing.assert_allclose(sys.data['cells'][0], expected_cell, atol=1e-6)

    def test_ibrav_6_tetragonal(self):
        """Test ibrav=6 (simple tetragonal)."""
        _, output_file = self._create_qe_files(6, "a=10\nc=8", None)
        
        sys = dpdata.LabeledSystem(output_file, fmt='qe/pw/scf')
        expected_cell = np.array([[10, 0, 0], [0, 10, 0], [0, 0, 8]])
        np.testing.assert_allclose(sys.data['cells'][0], expected_cell, rtol=1e-10)

    def test_ibrav_missing_parameters(self):
        """Test error handling for missing required parameters."""
        _, output_file = self._create_qe_files(8, "a=10", None)  # Missing b and c
        
        with self.assertRaises(RuntimeError) as cm:
            dpdata.LabeledSystem(output_file, fmt='qe/pw/scf')
        self.assertIn("parameter 'b'", str(cm.exception))

    def test_unsupported_ibrav(self):
        """Test error handling for unsupported ibrav values."""
        _, output_file = self._create_qe_files(99, "a=10", None)
        
        with self.assertRaises(RuntimeError) as cm:
            dpdata.LabeledSystem(output_file, fmt='qe/pw/scf')
        self.assertIn("ibrav = 99 is not supported", str(cm.exception))


if __name__ == '__main__':
    unittest.main()