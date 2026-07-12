from __future__ import annotations

import os
import shutil
import unittest

import lmdb
import msgpack
import numpy as np

import dpdata
from dpdata.data_type import Axis, DataType
from dpdata.formats.lmdb.format import LMDBError, LMDBFrameError


class TestLMDBFrameData(unittest.TestCase):
    """Frame-dependent custom data must round-trip through LMDB."""

    def setUp(self):
        self.original_system_dtypes = dpdata.System.DTYPES
        self.original_labeled_system_dtypes = dpdata.LabeledSystem.DTYPES

        self.dt_frame = DataType(
            "frame_data", np.ndarray, shape=(Axis.NFRAMES, 2), required=False
        )
        dpdata.System.register_data_type(self.dt_frame)
        dpdata.LabeledSystem.register_data_type(self.dt_frame)

        self.lmdb_path = "tmp_custom_dtype.lmdb"

        try:
            self.system = dpdata.LabeledSystem(
                "poscars/OUTCAR.h2o.md", fmt="vasp/outcar"
            )
        except FileNotFoundError:
            self.system = dpdata.LabeledSystem(
                "tests/poscars/OUTCAR.h2o.md", fmt="vasp/outcar"
            )

        nframes = self.system.get_nframes()
        self.system.data["frame_data"] = np.random.rand(nframes, 2)
        self.system.check_data()

    def tearDown(self):
        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)

    def test_frame_data_preservation(self):
        self.system.to("lmdb", self.lmdb_path)
        loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        np.testing.assert_allclose(
            loaded.data["frame_data"], self.system.data["frame_data"]
        )

    def test_frame_data_auto_registration(self):
        self.system.to("lmdb", self.lmdb_path)

        # simulate a clean session
        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes
        self.assertNotIn("frame_data", [dt.name for dt in dpdata.LabeledSystem.DTYPES])

        loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        self.assertIn("frame_data", [dt.name for dt in dpdata.LabeledSystem.DTYPES])
        np.testing.assert_allclose(
            loaded.data["frame_data"], self.system.data["frame_data"]
        )


class TestLMDBFparamAparam(unittest.TestCase):
    def setUp(self):
        self.original_system_dtypes = dpdata.System.DTYPES
        self.original_labeled_system_dtypes = dpdata.LabeledSystem.DTYPES

        new_datatypes = [
            DataType("fparam", np.ndarray, shape=(Axis.NFRAMES, 2), required=False),
            DataType(
                "aparam",
                np.ndarray,
                shape=(Axis.NFRAMES, Axis.NATOMS, 3),
                required=False,
            ),
        ]
        for datatype in new_datatypes:
            dpdata.System.register_data_type(datatype)
            dpdata.LabeledSystem.register_data_type(datatype)

        self.lmdb_path = "tmp_fparam_aparam.lmdb"

        try:
            self.system = dpdata.LabeledSystem(
                "poscars/OUTCAR.h2o.md", fmt="vasp/outcar"
            )
        except FileNotFoundError:
            self.system = dpdata.LabeledSystem(
                "tests/poscars/OUTCAR.h2o.md", fmt="vasp/outcar"
            )

        nframes = self.system.get_nframes()
        natoms = self.system.get_natoms()
        self.system.data["fparam"] = np.random.rand(nframes, 2)
        self.system.data["aparam"] = np.random.rand(nframes, natoms, 3)
        self.system.check_data()

    def tearDown(self):
        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)

    def test_fparam_aparam_preservation(self):
        self.system.to("lmdb", self.lmdb_path)
        loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        np.testing.assert_allclose(loaded.data["fparam"], self.system.data["fparam"])
        np.testing.assert_allclose(loaded.data["aparam"], self.system.data["aparam"])

    def test_fparam_aparam_auto_registration(self):
        self.system.to("lmdb", self.lmdb_path)

        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes

        loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        self.assertIn("fparam", [dt.name for dt in dpdata.LabeledSystem.DTYPES])
        self.assertIn("aparam", [dt.name for dt in dpdata.LabeledSystem.DTYPES])
        np.testing.assert_allclose(loaded.data["fparam"], self.system.data["fparam"])
        np.testing.assert_allclose(loaded.data["aparam"], self.system.data["aparam"])

    def test_symbolic_axis_natoms_preservation(self):
        # natoms == 3 collides with the trailing coordinate dim; the stored
        # symbolic shape must still recover (NFRAMES, NATOMS, 3), not
        # (NFRAMES, NATOMS, NATOMS).
        self.system.to("lmdb", self.lmdb_path)

        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes

        dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        aparam_dt = next(
            dt for dt in dpdata.LabeledSystem.DTYPES if dt.name == "aparam"
        )
        self.assertEqual(aparam_dt.shape, (Axis.NFRAMES, Axis.NATOMS, 3))

        data_diff = {
            "atom_numbs": [5],
            "atom_names": ["H"],
            "atom_types": np.array([0, 0, 0, 0, 0]),
            "coords": np.random.rand(1, 5, 3),
            "cells": np.random.rand(1, 3, 3),
            "orig": np.array([0, 0, 0]),
            "energies": np.array([1.0]),
            "forces": np.random.rand(1, 5, 3),
            "aparam": np.random.rand(1, 5, 3),
        }
        try:
            dpdata.LabeledSystem(data=data_diff)
        except dpdata.data_type.DataError as e:
            self.fail(f"DataError raised despite symbolic NATOMS: {e}")


class TestLMDBFieldProtocol(unittest.TestCase):
    def setUp(self):
        self.original_system_dtypes = dpdata.System.DTYPES
        self.original_labeled_system_dtypes = dpdata.LabeledSystem.DTYPES
        self.lmdb_path = "tmp_field_protocol.lmdb"
        try:
            self.system = dpdata.LabeledSystem(
                "poscars/OUTCAR.h2o.md", fmt="vasp/outcar"
            )
        except FileNotFoundError:
            self.system = dpdata.LabeledSystem(
                "tests/poscars/OUTCAR.h2o.md", fmt="vasp/outcar"
            )

    def tearDown(self):
        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)

    def _register(self, *dtypes):
        dpdata.System.register_data_type(*dtypes)
        dpdata.LabeledSystem.register_data_type(*dtypes)

    def test_deepmd_name_used_on_disk_and_restored(self):
        spin = DataType(
            "spins",
            np.ndarray,
            (Axis.NFRAMES, Axis.NATOMS, 3),
            required=False,
            deepmd_name="spin",
        )
        self._register(spin)
        self.system.data["spins"] = np.arange(
            self.system.get_nframes() * self.system.get_natoms() * 3,
            dtype=float,
        ).reshape(self.system.get_nframes(), self.system.get_natoms(), 3)
        expected = self.system.data["spins"].copy()
        self.system.to("lmdb", self.lmdb_path)

        with lmdb.open(self.lmdb_path, readonly=True, lock=False) as env:
            with env.begin() as txn:
                metadata = msgpack.unpackb(
                    txn.get(b"__metadata__"), raw=False
                )
                frame = msgpack.unpackb(
                    txn.get(b"000000000000"), raw=False
                )
        self.assertIn("spin", frame)
        self.assertNotIn("spins", frame)
        self.assertEqual(metadata["dp_data_names"], {"spin": "spins"})

        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes
        loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        np.testing.assert_array_equal(loaded["spins"], expected)
        registered = next(
            dt for dt in dpdata.LabeledSystem.DTYPES if dt.name == "spins"
        )
        self.assertEqual(registered.deepmd_name, "spin")

    def test_duplicate_deepmd_name_rejected(self):
        first = DataType(
            "first_field",
            np.ndarray,
            (Axis.NFRAMES, 1),
            required=False,
            deepmd_name="shared",
        )
        second = DataType(
            "second_field",
            np.ndarray,
            (Axis.NFRAMES, 1),
            required=False,
            deepmd_name="shared",
        )
        self._register(first, second)
        nframes = self.system.get_nframes()
        self.system.data["first_field"] = np.zeros((nframes, 1))
        self.system.data["second_field"] = np.ones((nframes, 1))
        with self.assertRaisesRegex(
            LMDBError,
            "both map to",
        ):
            self.system.to("lmdb", self.lmdb_path)

    def test_deepmd_core_name_collision_rejected(self):
        malicious = DataType(
            "malicious_coords",
            np.ndarray,
            (Axis.NFRAMES, Axis.NATOMS, 3),
            required=False,
            deepmd_name="coord",
        )
        self._register(malicious)
        self.system.data["malicious_coords"] = np.full(
            (
                self.system.get_nframes(),
                self.system.get_natoms(),
                3,
            ),
            101.0,
        )
        with self.assertRaisesRegex(LMDBError, "reserved LMDB key"):
            self.system.to("lmdb", self.lmdb_path)

    def test_additional_protocol_alias_rejected(self):
        alias = DataType(
            "magmom_alias",
            np.ndarray,
            (Axis.NFRAMES, Axis.NATOMS, 3),
            required=False,
            deepmd_name="spin",
        )
        self._register(alias)
        self.system.data["magmom_alias"] = np.zeros(
            (
                self.system.get_nframes(),
                self.system.get_natoms(),
                3,
            )
        )
        with self.assertRaisesRegex(
            LMDBError, "belongs to 'spins'"
        ):
            self.system.to("lmdb", self.lmdb_path)

    def test_atom_types_key_collision_rejected(self):
        malicious = DataType(
            "malicious_types",
            np.ndarray,
            (Axis.NFRAMES, Axis.NATOMS),
            required=False,
            deepmd_name="atom_types",
        )
        self._register(malicious)
        self.system.data["malicious_types"] = np.zeros(
            (self.system.get_nframes(), self.system.get_natoms())
        )
        with self.assertRaisesRegex(LMDBError, "reserved LMDB key"):
            self.system.to("lmdb", self.lmdb_path)

    def test_static_field_roundtrip(self):
        static = DataType(
            "static_data", np.ndarray, shape=(2,), required=False
        )
        self._register(static)
        self.system.data["static_data"] = np.array([1.0, 2.0])
        self.system.to("lmdb", self.lmdb_path)
        loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        np.testing.assert_array_equal(loaded["static_data"], [1.0, 2.0])

    def test_inconsistent_static_field_raises(self):
        static = DataType(
            "static_data", np.ndarray, shape=(2,), required=False
        )
        self._register(static)
        self.system.data["static_data"] = np.array([1.0, 2.0])
        self.system.to("lmdb", self.lmdb_path)

        env = lmdb.open(self.lmdb_path, map_size=1 << 30)
        with env.begin(write=True) as txn:
            key = b"000000000001"
            frame = msgpack.unpackb(txn.get(key), raw=False)
            replacement = np.array([3.0, 4.0])
            frame["static_data"] = {
                "type": str(replacement.dtype),
                "shape": list(replacement.shape),
                "data": replacement.tobytes(),
            }
            txn.put(key, msgpack.packb(frame, use_bin_type=True))
        env.close()

        with self.assertRaisesRegex(
            LMDBFrameError,
            "Static field",
        ):
            dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")

    def test_nonleading_frame_axis_roundtrip(self):
        transposed = DataType(
            "transposed_frames",
            np.ndarray,
            shape=(2, Axis.NFRAMES),
            required=False,
        )
        self._register(transposed)
        values = np.arange(
            2 * self.system.get_nframes(), dtype=float
        ).reshape(2, self.system.get_nframes())
        self.system.data["transposed_frames"] = values
        self.system.to("lmdb", self.lmdb_path)
        loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        np.testing.assert_array_equal(loaded["transposed_frames"], values)

    def test_unused_shape_none_data_type_does_not_block_write(self):
        undefined = DataType(
            "unused_undefined_shape",
            np.ndarray,
            shape=None,
            required=False,
        )
        self._register(undefined)
        self.system.to("lmdb", self.lmdb_path)
        loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        self.assertEqual(loaded.get_nframes(), self.system.get_nframes())

    def test_used_shape_none_data_type_rejected(self):
        undefined = DataType(
            "used_undefined_shape",
            np.ndarray,
            shape=None,
            required=False,
        )
        self._register(undefined)
        self.system.data["used_undefined_shape"] = np.zeros(
            (self.system.get_nframes(), 2)
        )
        with self.assertRaisesRegex(LMDBError, "no declared shape"):
            self.system.to("lmdb", self.lmdb_path)

    def test_multiple_atom_axes_remove_virtual_atoms(self):
        hessian = DataType(
            "hessian",
            np.ndarray,
            (
                Axis.NFRAMES,
                Axis.NATOMS,
                3,
                Axis.NATOMS,
                3,
            ),
            required=False,
            deepmd_name="hessian",
        )
        self._register(hessian)
        data = {
            "atom_numbs": [3],
            "atom_names": ["MIXED_TOKEN"],
            "atom_types": np.zeros(3, dtype=int),
            "real_atom_names": ["H", "O"],
            "real_atom_types": np.array([[1, 0, -1]], dtype=int),
            "orig": np.zeros(3),
            "cells": np.eye(3)[None],
            "coords": np.zeros((1, 3, 3)),
            "energies": np.array([-1.0]),
            "hessian": np.zeros((1, 3, 3, 3, 3)),
        }
        from dpdata.formats.lmdb import dump_systems

        dump_systems([data], self.lmdb_path)
        with lmdb.open(self.lmdb_path, readonly=True, lock=False) as env:
            with env.begin() as txn:
                frame = msgpack.unpackb(
                    txn.get(b"000000000000"), raw=False
                )
        self.assertEqual(frame["hessian"]["shape"], [2, 3, 2, 3])


if __name__ == "__main__":
    unittest.main()
