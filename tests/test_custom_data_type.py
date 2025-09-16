from __future__ import annotations

import unittest
import warnings

import h5py  # noqa: TID253
import numpy as np

import dpdata
from dpdata.data_type import Axis, DataType


class TestDataType(unittest.TestCase):
    """Test DataType class methods."""

    def setUp(self):
        # Store original DTYPES to restore later
        self.original_dtypes = dpdata.System.DTYPES

    def tearDown(self):
        # Restore original DTYPES
        dpdata.System.DTYPES = self.original_dtypes

    def test_eq(self):
        """Test equality method."""
        dt1 = DataType("test", np.ndarray, shape=(Axis.NFRAMES, 3))
        dt2 = DataType("test", np.ndarray, shape=(Axis.NFRAMES, 3))
        dt3 = DataType("other", np.ndarray, shape=(Axis.NFRAMES, 3))

        self.assertTrue(dt1 == dt2)
        self.assertFalse(dt1 == dt3)
        self.assertFalse(dt1 == "not a DataType")

    def test_repr(self):
        """Test string representation."""
        dt = DataType("test", np.ndarray, shape=(Axis.NFRAMES, 3))
        expected = (
            "DataType(name='test', dtype=ndarray, "
            "shape=(<Axis.NFRAMES: 'nframes'>, 3), required=True, "
            "deepmd_name='test')"
        )
        self.assertEqual(repr(dt), expected)

    def test_register_same_data_type_no_warning(self):
        """Test registering identical DataType instances should not warn."""
        dt1 = DataType("test_same", np.ndarray, shape=(Axis.NFRAMES, 3))
        dt2 = DataType("test_same", np.ndarray, shape=(Axis.NFRAMES, 3))

        # Register first time
        dpdata.System.register_data_type(dt1)

        # Register same DataType again - should not warn
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            dpdata.System.register_data_type(dt2)
            # Check no warnings were issued
            self.assertEqual(len(w), 0)

    def test_register_different_data_type_with_warning(self):
        """Test registering different DataType instances with same name should warn."""
        dt1 = DataType("test_diff", np.ndarray, shape=(Axis.NFRAMES, 3))
        dt2 = DataType(
            "test_diff", list, shape=(Axis.NFRAMES, 4)
        )  # Different dtype and shape

        # Register first time
        dpdata.System.register_data_type(dt1)

        # Register different DataType with same name - should warn
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            dpdata.System.register_data_type(dt2)
            # Check warning was issued
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, UserWarning))
            self.assertIn(
                "registered twice with different definitions", str(w[-1].message)
            )


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

    def test_to_deepmd_npy_mixed(self):
        ms = dpdata.MultiSystems(self.system)
        ms.to_deepmd_npy_mixed("data_foo_mixed")
        x = dpdata.MultiSystems().load_systems_from_file(
            "data_foo_mixed",
            fmt="deepmd/npy/mixed",
            labeled=issubclass(self.cls, dpdata.LabeledSystem),
        )
        np.testing.assert_allclose(list(x.systems.values())[0].data["foo"], self.foo)


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
        print(self.system.data.keys())
        print(x.data.keys())
        np.testing.assert_allclose(x.data["bar"], self.bar)
