import unittest

import h5py
import numpy as np

import dpdata
from dpdata.system import Axis, DataType


class TestDeepmdLoadDumpComp(unittest.TestCase):
    def setUp(self):
        self.system = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")
        self.foo = np.ones((len(self.system), 2, 4))
        self.system.data["foo"] = self.foo
        self.system.check_data()

    def test_to_deepmd_raw(self):
        self.system.to_deepmd_raw("data_foo")
        foo = np.loadtxt("data_foo/foo.raw")
        np.testing.assert_allclose(foo.reshape(self.foo.shape), self.foo)

    def test_from_deepmd_raw(self):
        self.system.to_deepmd_raw("data_foo")
        x = dpdata.LabeledSystem("data_foo", fmt="deepmd/raw")
        np.testing.assert_allclose(x.data["foo"], self.foo)

    def test_to_deepmd_npy(self):
        self.system.to_deepmd_npy("data_foo")
        foo = np.load("data_foo/set.000/foo.npy")
        np.testing.assert_allclose(foo.reshape(self.foo.shape), self.foo)

    def test_from_deepmd_npy(self):
        self.system.to_deepmd_npy("data_foo")
        x = dpdata.LabeledSystem("data_foo", fmt="deepmd/npy")
        np.testing.assert_allclose(x.data["foo"], self.foo)

    def test_to_deepmd_hdf5(self):
        self.system.to_deepmd_hdf5("data_foo.h5")
        with h5py.File("data_foo.h5") as f:
            foo = f["set.000/foo.npy"][:]
        np.testing.assert_allclose(foo.reshape(self.foo.shape), self.foo)

    def test_from_deepmd_hdf5(self):
        self.system.to_deepmd_hdf5("data_foo.h5")
        x = dpdata.LabeledSystem("data_foo.h5", fmt="deepmd/hdf5")
        np.testing.assert_allclose(x.data["foo"], self.foo)
