import unittest

import numpy as np
from context import dpdata


class TestSplitDataset(unittest.TestCase):
    def setUp(self):
        self.systems = dpdata.MultiSystems()
        sing_sys = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")
        for ii in range(10):
            self.systems.append(sing_sys.copy())

    def test_split_dataset(self):
        train, test, test_idx = self.systems.train_test_split(0.2)
        self.assertEqual(
            train.get_nframes(), int(np.floor(self.systems.get_nframes() * 0.8))
        )
        self.assertEqual(
            test.get_nframes(), int(np.floor(self.systems.get_nframes() * 0.2))
        )
        self.assertEqual(
            sum([np.count_nonzero(x) for x in test_idx.values()]),
            int(np.floor(self.systems.get_nframes() * 0.2)),
        )
