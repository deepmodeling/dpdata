from __future__ import annotations

import os
import shutil
import unittest
import tempfile
import numpy as np
from comp_sys import (
    CompSys,
    MultiSystems,
    IsNoPBC,
)
from context import dpdata


class TestMixedMultiSystemsUnlabeled(unittest.TestCase, MultiSystems, IsNoPBC):
    """Test MultiSystems with unlabeled systems using deepmd/npy/mixed format."""

    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        # Create unlabeled systems (no energies/forces)
        nframes = 2
        
        # H2O system without labels
        system1_data = {
            'atom_names': ['H', 'O'],
            'atom_numbs': [2, 1],
            'atom_types': np.array([0, 0, 1]),
            'orig': np.array([0., 0., 0.]),
            'cells': np.random.random((nframes, 3, 3)),
            'coords': np.random.random((nframes, 3, 3))
        }
        
        # CH4 system without labels
        system2_data = {
            'atom_names': ['C', 'H'],
            'atom_numbs': [1, 4],
            'atom_types': np.array([0, 1, 1, 1, 1]),
            'orig': np.array([0., 0., 0.]),
            'cells': np.random.random((nframes, 3, 3)),
            'coords': np.random.random((nframes, 5, 3))
        }

        # Create System objects (unlabeled)
        system_1 = dpdata.System(data=system1_data)
        system_2 = dpdata.System(data=system2_data)

        self.ms = dpdata.MultiSystems(system_1, system_2)
        
        # Create temporary directory for testing
        self.tmpdir = tempfile.mkdtemp()
        self.mixed_dir = os.path.join(self.tmpdir, "mixed_unlabeled")
        
        # Dump to mixed format
        self.ms.to_deepmd_npy_mixed(self.mixed_dir)
        
        # Load back using the fixed method
        self.systems = dpdata.MultiSystems.from_file(self.mixed_dir, fmt="deepmd/npy/mixed")
        
        self.ms_1 = self.ms
        self.ms_2 = self.systems
        
        # For compatibility with inherited tests
        self.system_1 = self.ms
        self.system_2 = self.systems
        
        self.system_names = ["H4O0C1", "H2O1C0"]
        self.system_sizes = {"H4O0C1": 1, "H2O1C0": 1}
        self.atom_names = ["H", "O", "C"]

    def tearDown(self):
        if os.path.exists(self.tmpdir):
            shutil.rmtree(self.tmpdir)

    def test_systems_are_unlabeled(self):
        """Test that loaded systems are System objects, not LabeledSystem."""
        for name, system in self.systems.systems.items():
            self.assertIsInstance(system, dpdata.System)
            self.assertNotIsInstance(system, dpdata.LabeledSystem)
            # Verify no energy/force data
            self.assertNotIn("energies", system.data)
            self.assertNotIn("forces", system.data)

    def test_len(self):
        self.assertEqual(len(self.ms), 2)
        self.assertEqual(len(self.systems), 2)

    def test_get_nframes(self):
        self.assertEqual(self.ms.get_nframes(), 4)
        self.assertEqual(self.systems.get_nframes(), 4)

    def test_str(self):
        self.assertEqual(str(self.ms), "MultiSystems (2 systems containing 4 frames)")
        self.assertEqual(str(self.systems), "MultiSystems (2 systems containing 4 frames)")


if __name__ == "__main__":
    unittest.main()