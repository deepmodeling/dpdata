from __future__ import annotations

import unittest

import h5py  # noqa: TID253
import numpy as np

import dpdata
from dpdata.data_type import Axis, DataType


class DeepmdLoadDumpCompTest:
    def setUp(self):
        self.system = self.cls(
            data=dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar").data
        )
        self.foo = np.ones((len(self.system), *self.shape))
        self.system.data["foo"] = self.foo
        self.system.check_data()

    def test_to_deepmd_raw(self):
        self.system.to_deepmd_raw("data_foo")
        foo = np.loadtxt("data_foo/foo.raw")
        np.testing.assert_allclose(foo.reshape(self.foo.shape), self.foo)

    def test_from_deepmd_raw(self):
        self.system.to_deepmd_raw("data_foo")
        x = self.cls("data_foo", fmt="deepmd/raw")
        np.testing.assert_allclose(x.data["foo"], self.foo)

    def test_to_deepmd_npy(self):
        self.system.to_deepmd_npy("data_foo")
        foo = np.load("data_foo/set.000/foo.npy")
        np.testing.assert_allclose(foo.reshape(self.foo.shape), self.foo)

    def test_from_deepmd_npy(self):
        self.system.to_deepmd_npy("data_foo")
        x = self.cls("data_foo", fmt="deepmd/npy")
        np.testing.assert_allclose(x.data["foo"], self.foo)

    def test_to_deepmd_hdf5(self):
        self.system.to_deepmd_hdf5("data_foo.h5")
        with h5py.File("data_foo.h5") as f:
            foo = f["set.000/foo.npy"][:]
        np.testing.assert_allclose(foo.reshape(self.foo.shape), self.foo)

    def test_from_deepmd_hdf5(self):
        self.system.to_deepmd_hdf5("data_foo.h5")
        x = self.cls("data_foo.h5", fmt="deepmd/hdf5")
        np.testing.assert_allclose(x.data["foo"], self.foo)

    def test_duplicated_data_type(self):
        dt = DataType("foo", np.ndarray, (Axis.NFRAMES, *self.shape), required=False)
        n_dtypes_old = len(self.cls.DTYPES)
        with self.assertWarns(UserWarning):
            self.cls.register_data_type(dt)
        n_dtypes_new = len(self.cls.DTYPES)
        self.assertEqual(n_dtypes_old, n_dtypes_new)


class TestDeepmdLoadDumpCompUnlabeled(unittest.TestCase, DeepmdLoadDumpCompTest):
    cls = dpdata.System
    shape = (3, 3)

    def setUp(self):
        DeepmdLoadDumpCompTest.setUp(self)


class TestDeepmdLoadDumpCompLabeled(unittest.TestCase, DeepmdLoadDumpCompTest):
    cls = dpdata.LabeledSystem
    shape = (2, 4)

    def setUp(self):
        DeepmdLoadDumpCompTest.setUp(self)


class TestDeepmdLoadDumpCompAny(unittest.TestCase):
    def setUp(self):
        self.system = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")
        self.bar = np.ones((len(self.system), self.system.get_natoms(), 2))
        self.system.data["bar"] = self.bar
        self.system.check_data()

    def test_to_deepmd_raw(self):
        self.system.to_deepmd_raw("data_bar")
        bar = np.loadtxt("data_bar/bar.raw")
        np.testing.assert_allclose(bar.reshape(self.bar.shape), self.bar)

    def test_from_deepmd_raw(self):
        self.system.to_deepmd_raw("data_bar")
        x = dpdata.LabeledSystem("data_bar", fmt="deepmd/raw")
        np.testing.assert_allclose(x.data["bar"], self.bar)

    def test_to_deepmd_npy(self):
        self.system.to_deepmd_npy("data_bar")
        bar = np.load("data_bar/set.000/bar.npy")
        np.testing.assert_allclose(bar.reshape(self.bar.shape), self.bar)

    def test_from_deepmd_npy(self):
        self.system.to_deepmd_npy("data_bar")
        x = dpdata.LabeledSystem("data_bar", fmt="deepmd/npy")
        np.testing.assert_allclose(x.data["bar"], self.bar)

    def test_to_deepmd_hdf5(self):
        self.system.to_deepmd_hdf5("data_bar.h5")
        with h5py.File("data_bar.h5") as f:
            bar = f["set.000/bar.npy"][:]
        np.testing.assert_allclose(bar.reshape(self.bar.shape), self.bar)

    def test_from_deepmd_hdf5(self):
        self.system.to_deepmd_hdf5("data_bar.h5")
        x = dpdata.LabeledSystem("data_bar.h5", fmt="deepmd/hdf5")
        np.testing.assert_allclose(x.data["bar"], self.bar)
