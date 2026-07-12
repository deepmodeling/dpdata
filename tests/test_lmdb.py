from __future__ import annotations

import os
import shutil
import threading
import typing
import unittest
from pathlib import Path
from unittest import mock

import lmdb
import msgpack
import numpy as np
from comp_sys import (
    CompLabeledMultiSys,
    CompLabeledSys,
    CompSys,
    IsPBC,
    MSAllIsNoPBC,
)
from context import dpdata

from dpdata.data_type import Axis, DataError, DataType
from dpdata.formats.lmdb.format import (
    LMDBError,
    LMDBFrameError,
    LMDBMetadataError,
    _close_read_env,
    _open_read_env,
)


def _reference_encode_array(array):
    """Encode an array independently of the production implementation."""
    value = np.asarray(array)
    return {
        "type": str(value.dtype),
        "shape": list(value.shape),
        "data": value.tobytes(order="C"),
    }


def _reference_packb(value):
    """Pack a protocol value and assert the LMDB byte contract."""
    packed = msgpack.packb(value, use_bin_type=True)
    if not isinstance(packed, bytes):
        raise TypeError("msgpack.packb did not return bytes")
    return packed


class TestLMDBLabeledSystem(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")
        self.lmdb_path = "tmp_labeled.lmdb"
        self.system_1.to("lmdb", self.lmdb_path)
        self.system_2 = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

    def tearDown(self):
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)


class TestLMDBSystem(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.System("poscars/POSCAR.h2o.md", fmt="vasp/poscar")
        self.lmdb_path = "tmp_system.lmdb"
        self.system_1.to("lmdb", self.lmdb_path)
        self.system_2 = dpdata.System(self.lmdb_path, fmt="lmdb")
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

    def tearDown(self):
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)


class TestLMDBMultiSystems(unittest.TestCase, CompLabeledMultiSys, MSAllIsNoPBC):
    def setUp(self):
        self.lmdb_path = "tmp_multi.lmdb"
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_reordered.gaussianlog", fmt="gaussian/log"
        )
        system_3 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        self.ms_1 = dpdata.MultiSystems(system_1, system_2, system_3)
        for system in self.ms_1:
            system.sort_atom_types()
        self.ms_1.to("lmdb", self.lmdb_path)
        self.ms_2 = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")

        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

    def tearDown(self):
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)


class TestLMDBOnDiskFormat(unittest.TestCase):
    """The on-disk layout must match the DeePMD-kit / reference converters."""

    def setUp(self):
        self.lmdb_path = "tmp_format.lmdb"
        self.system = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")
        self.system.to("lmdb", self.lmdb_path)

    def tearDown(self):
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)

    def test_metadata_schema(self):
        with lmdb.open(self.lmdb_path, readonly=True, lock=False) as env:
            with env.begin() as txn:
                meta = msgpack.unpackb(txn.get(b"__metadata__"), raw=False)
        self.assertEqual(meta["nframes"], self.system.get_nframes())
        self.assertEqual(meta["frame_idx_fmt"], "012d")
        self.assertEqual(meta["type_map"], self.system.data["atom_names"])
        self.assertEqual(len(meta["frame_nlocs"]), self.system.get_nframes())
        self.assertTrue(all(n == self.system.get_natoms() for n in meta["frame_nlocs"]))
        self.assertEqual(len(meta["frame_system_ids"]), self.system.get_nframes())
        self.assertTrue(all(s == 0 for s in meta["frame_system_ids"]))
        # legacy field must be gone
        self.assertNotIn("system_info", meta)

    def test_frame_encoding(self):
        with lmdb.open(self.lmdb_path, readonly=True, lock=False) as env:
            with env.begin() as txn:
                raw = txn.get(b"000000000000")
                meta = msgpack.unpackb(txn.get(b"__metadata__"), raw=False)
        frame = msgpack.unpackb(raw, raw=False)
        # string keys, manual encoding
        for key in ("atom_types", "coords", "cells", "energies", "forces"):
            self.assertIn(key, frame)
            self.assertIsInstance(frame[key], dict)
            self.assertIn("type", frame[key])
            self.assertIn("shape", frame[key])
            self.assertIn("data", frame[key])
        # atom_types are int32 global indices
        self.assertEqual(frame["atom_types"]["type"], "int32")
        # energies stored as a 0-d array
        self.assertEqual(list(frame["energies"]["shape"]), [])
        # atom_numbs is a plain list over the full type_map
        self.assertIsInstance(frame["atom_numbs"], list)
        self.assertEqual(len(frame["atom_numbs"]), len(meta["type_map"]))

    def test_atom_types_are_global(self):
        with lmdb.open(self.lmdb_path, readonly=True, lock=False) as env:
            with env.begin() as txn:
                raw = txn.get(b"000000000000")
                meta = msgpack.unpackb(txn.get(b"__metadata__"), raw=False)
        frame = msgpack.unpackb(raw, raw=False)
        atype = np.frombuffer(frame["atom_types"]["data"], dtype=np.int32).reshape(
            frame["atom_types"]["shape"]
        )
        type_map = meta["type_map"]
        decoded_names = [type_map[i] for i in atype]
        expected = [
            self.system.data["atom_names"][t] for t in self.system.data["atom_types"]
        ]
        self.assertEqual(decoded_names, expected)


class TestLMDBMixedTypeRead(unittest.TestCase):
    """mixed_type=True keeps the full global type_map on every system."""

    def setUp(self):
        self.lmdb_path = "tmp_mixed.lmdb"
        s1 = dpdata.LabeledSystem("gaussian/methane.gaussianlog", fmt="gaussian/log")
        s2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )
        self.ms = dpdata.MultiSystems(s1, s2)
        self.type_map = ["H", "C", "N", "O"]
        self.ms.to("lmdb", self.lmdb_path, type_map=self.type_map)

    def tearDown(self):
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)

    def test_mixed_preserves_full_type_map(self):
        # via MultiSystems the element order is normalized (sorted), but the
        # full type_map set must be preserved on every system.
        ms = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb", mixed_type=True)
        for ss in ms:
            names = [str(n) for n in ss.data["atom_names"]]
            self.assertEqual(sorted(names), sorted(self.type_map))
            self.assertEqual(len(ss.data["atom_numbs"]), len(self.type_map))

    def test_mixed_single_system_preserves_order(self):
        # a direct single-system load of a single-composition file is not
        # reordered, so the global type_map order is preserved exactly.
        path = "tmp_mixed_single.lmdb"
        try:
            s = dpdata.LabeledSystem("gaussian/methane.gaussianlog", fmt="gaussian/log")
            s.to("lmdb", path, type_map=self.type_map)
            ls = dpdata.LabeledSystem(path, fmt="lmdb", mixed_type=True)
            self.assertEqual([str(n) for n in ls.data["atom_names"]], self.type_map)
        finally:
            if os.path.exists(path):
                shutil.rmtree(path)

    def test_standard_compresses_type_map(self):
        ms = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")
        # methane only contains C and H, not N/O
        for ss in ms:
            self.assertNotIn("N", ss.data["atom_names"])
            self.assertNotIn("O", ss.data["atom_names"])


class TestLMDBTypeMapOverride(unittest.TestCase):
    def setUp(self):
        self.lmdb_path = "tmp_type_map.lmdb"
        self.system = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")

    def tearDown(self):
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)

    def test_explicit_type_map(self):
        type_map = ["H", "He", "Li", "Be", "B", "C", "N", "O"]
        self.system.to("lmdb", self.lmdb_path, type_map=type_map)
        with lmdb.open(self.lmdb_path, readonly=True, lock=False) as env:
            with env.begin() as txn:
                meta = msgpack.unpackb(txn.get(b"__metadata__"), raw=False)
                raw = txn.get(b"000000000000")
        self.assertEqual(meta["type_map"], type_map)
        frame = msgpack.unpackb(raw, raw=False)
        atype = np.frombuffer(frame["atom_types"]["data"], dtype=np.int32).reshape(
            frame["atom_types"]["shape"]
        )
        decoded = [type_map[i] for i in atype]
        expected = [
            self.system.data["atom_names"][t] for t in self.system.data["atom_types"]
        ]
        self.assertEqual(decoded, expected)

    def test_missing_element_raises(self):
        from dpdata.formats.lmdb.format import LMDBError

        # water needs O and H; this type_map omits O.
        with self.assertRaises(LMDBError):
            self.system.to("lmdb", self.lmdb_path, type_map=["H", "He"])


class TestLMDBReferenceFormatInterop(unittest.TestCase):
    """dpdata must read an LMDB written exactly like the reference converters."""

    def setUp(self):
        self.lmdb_path = "tmp_reference.lmdb"
        self.type_map = ["H", "C", "N", "O"]
        # two frames: H2O (3 atoms) and CH4 (5 atoms)
        self.frames = [
            {
                "atom_types": np.array([3, 0, 0], dtype=np.int32),  # O H H
                "coords": np.random.rand(3, 3).astype(np.float32),
                "cells": (np.eye(3) * 10).astype(np.float32),
                "energies": np.array(-10.0, dtype=np.float32),
                "forces": np.random.rand(3, 3).astype(np.float32),
            },
            {
                "atom_types": np.array([1, 0, 0, 0, 0], dtype=np.int32),  # C H H H H
                "coords": np.random.rand(5, 3).astype(np.float32),
                "cells": (np.eye(3) * 12).astype(np.float32),
                "energies": np.array(-20.0, dtype=np.float32),
                "forces": np.random.rand(5, 3).astype(np.float32),
            },
        ]
        env = lmdb.open(self.lmdb_path, map_size=1 << 30)
        frame_nlocs = []
        with env.begin(write=True) as txn:
            for i, fr in enumerate(self.frames):
                out: dict[str, object] = {}
                for k, v in fr.items():
                    out[k] = _reference_encode_array(v)
                ntypes = len(self.type_map)
                counts = np.bincount(fr["atom_types"], minlength=ntypes)[:ntypes]
                out["atom_numbs"] = [int(c) for c in counts]
                txn.put(format(i, "012d").encode(), _reference_packb(out))
                frame_nlocs.append(len(fr["atom_types"]))
            meta: dict[str, object] = {
                "nframes": len(self.frames),
                "frame_idx_fmt": "012d",
                "type_map": self.type_map,
                "frame_nlocs": frame_nlocs,
                "frame_system_ids": list(range(len(self.frames))),
            }
            txn.put(b"__metadata__", _reference_packb(meta))
        env.close()

    def tearDown(self):
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)

    def test_read_reference(self):
        ms = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")
        self.assertEqual(ms.get_nframes(), 2)
        self.assertEqual(len(ms), 2)
        # composition-based grouping: water (3 atoms) + methane (5 atoms)
        self.assertEqual(sorted(s.get_natoms() for s in ms), [3, 5])

    def test_read_reference_values(self):
        ms = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")
        water = next(s for s in ms if s.get_natoms() == 3)
        # Frames are canonicalized by global type while coordinates and labels
        # follow the same stable permutation.
        names = [water.data["atom_names"][t] for t in water.data["atom_types"]]
        self.assertEqual(names, ["H", "H", "O"])
        permutation = np.argsort(self.frames[0]["atom_types"], kind="stable")
        np.testing.assert_allclose(
            water.data["coords"][0],
            self.frames[0]["coords"][permutation],
            rtol=1e-6,
            atol=1e-6,
        )
        np.testing.assert_allclose(
            water.data["energies"][0], self.frames[0]["energies"], rtol=1e-6
        )


def _make_labeled_system(
    atom_types: list[int], atom_names: list[str], nframes: int
) -> dpdata.LabeledSystem:
    atype = np.array(atom_types, dtype=int)
    natoms = len(atype)
    numbs = [int(np.count_nonzero(atype == i)) for i in range(len(atom_names))]
    data = {
        "atom_numbs": numbs,
        "atom_names": list(atom_names),
        "atom_types": atype,
        "orig": np.zeros(3),
        "cells": np.tile(np.eye(3) * 10.0, (nframes, 1, 1)),
        "coords": np.random.rand(nframes, natoms, 3),
        "energies": np.arange(nframes, dtype=float),
        "forces": np.random.rand(nframes, natoms, 3),
    }
    return dpdata.LabeledSystem(data=data)


def _write_raw_lmdb(path, frames, type_map, metadata_extra=None):
    """Write frames using the reference encoding (one record per frame)."""
    if os.path.exists(path):
        shutil.rmtree(path)
    env = lmdb.open(path, map_size=1 << 30)
    nlocs = []
    with env.begin(write=True) as txn:
        for i, fr in enumerate(frames):
            out: dict[str, object] = {
                k: _reference_encode_array(np.asarray(v)) for k, v in fr.items()
            }
            at = np.asarray(fr["atom_types"])
            counts = np.bincount(at, minlength=len(type_map))[: len(type_map)]
            out["atom_numbs"] = [int(c) for c in counts]
            txn.put(format(i, "012d").encode(), _reference_packb(out))
            nlocs.append(len(at))
        meta: dict[str, object] = {
            "nframes": len(frames),
            "frame_idx_fmt": "012d",
            "type_map": type_map,
            "frame_nlocs": nlocs,
            "frame_system_ids": list(range(len(frames))),
        }
        if metadata_extra:
            meta.update(metadata_extra)
        txn.put(b"__metadata__", _reference_packb(meta))
    env.close()


class TestLMDBReadRobustness(unittest.TestCase):
    """Regression tests for the read path."""

    def setUp(self):
        self.lmdb_path = "tmp_robust.lmdb"
        self.original_system_dtypes = dpdata.System.DTYPES
        self.original_labeled_system_dtypes = dpdata.LabeledSystem.DTYPES

    def tearDown(self):
        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)

    def test_type_map_remap_by_name(self):
        # file stores type_map [H, O] (so O has global index 1); reading with a
        # reordered type_map must remap by element name, not by index.
        frames = [
            {
                "atom_types": np.array([1, 0, 0], dtype=np.int32),  # O H H
                "coords": np.random.rand(3, 3).astype(np.float32),
                "cells": (np.eye(3) * 10).astype(np.float32),
                "energies": np.array(-1.0, dtype=np.float32),
            }
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["H", "O"])
        ms = dpdata.MultiSystems.from_file(
            self.lmdb_path, fmt="lmdb", type_map=["O", "N", "C", "H"]
        )
        s = ms[0]
        names = [str(s.data["atom_names"][t]) for t in s.data["atom_types"]]
        self.assertEqual(names, ["O", "H", "H"])

    def test_type_map_missing_element_raises(self):
        frames = [
            {
                "atom_types": np.array([0, 1, 1], dtype=np.int32),
                "coords": np.random.rand(3, 3).astype(np.float32),
                "cells": (np.eye(3) * 10).astype(np.float32),
                "energies": np.array(-1.0, dtype=np.float32),
            }
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["O", "H"])
        from dpdata.formats.lmdb.format import LMDBError

        with self.assertRaises(LMDBError):
            dpdata.MultiSystems.from_file(
                self.lmdb_path, fmt="lmdb", type_map=["C", "N"]
            )

    def test_inconsistent_keys_raise(self):
        # same composition, but only the first frame carries virials.
        frames = [
            {
                "atom_types": np.array([0, 1, 1], dtype=np.int32),
                "coords": np.random.rand(3, 3).astype(np.float32),
                "cells": (np.eye(3) * 10).astype(np.float32),
                "energies": np.array(-1.0, dtype=np.float32),
                "forces": np.random.rand(3, 3).astype(np.float32),
                "virials": np.random.rand(3, 3).astype(np.float32),
            },
            {
                "atom_types": np.array([0, 1, 1], dtype=np.int32),
                "coords": np.random.rand(3, 3).astype(np.float32),
                "cells": (np.eye(3) * 10).astype(np.float32),
                "energies": np.array(-2.0, dtype=np.float32),
                "forces": np.random.rand(3, 3).astype(np.float32),
            },
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["O", "H"])
        with self.assertRaisesRegex(LMDBFrameError, "must contain identical fields"):
            dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")

    def test_single_system_multi_composition_warns(self):
        frames = [
            {
                "atom_types": np.array([0, 1, 1], dtype=np.int32),
                "coords": np.random.rand(3, 3).astype(np.float32),
                "cells": (np.eye(3) * 10).astype(np.float32),
                "energies": np.array(-1.0, dtype=np.float32),
            },
            {
                "atom_types": np.array([0, 1, 1, 1, 1], dtype=np.int32),
                "coords": np.random.rand(5, 3).astype(np.float32),
                "cells": (np.eye(3) * 10).astype(np.float32),
                "energies": np.array(-2.0, dtype=np.float32),
            },
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["C", "H"])
        with self.assertWarns(UserWarning):
            ls = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        # only the first composition is returned
        self.assertEqual(ls.get_natoms(), 3)

    def test_same_composition_different_atom_order_is_canonicalized(self):
        first_types = np.array([1, 0, 1], dtype=np.int32)
        second_types = np.array([0, 1, 1], dtype=np.int32)
        first_values = np.array([10.0, 20.0, 30.0])
        second_values = np.array([40.0, 50.0, 60.0])
        frames = [
            {
                "atom_types": first_types,
                "coords": np.repeat(first_values[:, None], 3, axis=1),
                "cells": np.eye(3),
                "energies": np.array(-1.0),
                "forces": np.repeat(first_values[:, None], 3, axis=1),
            },
            {
                "atom_types": second_types,
                "coords": np.repeat(second_values[:, None], 3, axis=1),
                "cells": np.eye(3),
                "energies": np.array(-2.0),
                "forces": np.repeat(second_values[:, None], 3, axis=1),
            },
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["H", "O"])
        ms = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")
        self.assertEqual(len(ms), 1)
        system = ms[0]
        self.assertEqual(system.get_nframes(), 2)
        np.testing.assert_array_equal(system["atom_types"], [0, 1, 1])
        np.testing.assert_array_equal(
            system["coords"][:, :, 0],
            [[20.0, 10.0, 30.0], [40.0, 50.0, 60.0]],
        )
        np.testing.assert_array_equal(
            system["forces"][:, :, 0],
            [[20.0, 10.0, 30.0], [40.0, 50.0, 60.0]],
        )

    def test_same_composition_mixed_pbc_raises(self):
        frames = [
            {
                "atom_types": np.array([0, 1, 1], dtype=np.int32),
                "coords": np.zeros((3, 3)),
                "cells": np.eye(3),
                "energies": np.array(-1.0),
            },
            {
                "atom_types": np.array([0, 1, 1], dtype=np.int32),
                "coords": np.ones((3, 3)),
                "energies": np.array(-2.0),
            },
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["O", "H"])
        with self.assertRaisesRegex(LMDBFrameError, "must contain identical fields"):
            dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")

    def test_max_frames_guard(self):
        frames = [
            {
                "atom_types": np.array([0], dtype=np.int32),
                "coords": np.zeros((1, 3)),
            },
            {
                "atom_types": np.array([0], dtype=np.int32),
                "coords": np.ones((1, 3)),
            },
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["H"])
        with self.assertRaisesRegex(LMDBError, "exceeding max_frames"):
            dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb", max_frames=1)
        loaded = dpdata.MultiSystems.from_file(
            self.lmdb_path, fmt="lmdb", max_frames=None, labeled=False
        )
        self.assertEqual(loaded.get_nframes(), 2)

    def test_big_endian_dtype_preserved(self):
        frames = [
            {
                "atom_types": np.array([0], dtype=np.int32),
                "coords": np.array([[1.0, 2.0, 3.0]], dtype=">f8"),
                "cells": np.eye(3, dtype=">f8"),
                "energies": np.array(-1.0, dtype=">f8"),
            },
            {
                "atom_types": np.array([0], dtype=np.int32),
                "coords": np.array([[4.0, 5.0, 6.0]], dtype=">f8"),
                "cells": np.eye(3, dtype=">f8"),
                "energies": np.array(-2.0, dtype=">f8"),
            },
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["H"])
        system = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")[0]
        self.assertEqual(system["coords"].dtype.str, ">f8")
        self.assertEqual(system["cells"].dtype.str, ">f8")
        self.assertEqual(system["energies"].dtype.str, ">f8")

    def test_reference_fparam_shape_does_not_collide_with_natoms(self):
        frames = [
            {
                "atom_types": np.array([0, 0], dtype=np.int32),
                "coords": np.zeros((2, 3)),
                "fparam": np.array([1.0, 2.0]),
            },
            {
                "atom_types": np.array([0, 0], dtype=np.int32),
                "coords": np.ones((2, 3)),
                "fparam": np.array([3.0, 4.0]),
            },
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["H"])
        dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb", labeled=False)
        dtype = next(dt for dt in dpdata.System.DTYPES if dt.name == "fparam")
        shape = dtype.shape
        self.assertIsNotNone(shape)
        assert shape is not None
        self.assertEqual(shape[0], Axis.NFRAMES)
        self.assertNotIn(Axis.NATOMS, shape)

    def test_flattened_reference_aparam_is_normalized_and_reordered(self):
        fixed_aparam = DataType(
            "aparam",
            np.ndarray,
            (Axis.NFRAMES, Axis.NATOMS, 3),
            required=False,
        )
        dpdata.System.register_data_type(fixed_aparam)
        dpdata.LabeledSystem.register_data_type(fixed_aparam)
        for raw, expected in (
            (
                np.array([10.0, 20.0, 30.0]),
                np.array([[[20.0], [10.0], [30.0]]]),
            ),
            (
                np.array([10.0, 11.0, 20.0, 21.0, 30.0, 31.0]),
                np.array([[[20.0, 21.0], [10.0, 11.0], [30.0, 31.0]]]),
            ),
        ):
            with self.subTest(size=raw.size):
                frames = [
                    {
                        "atom_types": np.array([1, 0, 1], dtype=np.int32),
                        "coords": np.zeros((3, 3)),
                        "aparam": raw,
                    }
                ]
                _write_raw_lmdb(self.lmdb_path, frames, ["H", "O"])
                system = dpdata.MultiSystems.from_file(
                    self.lmdb_path, fmt="lmdb", labeled=False
                )[0]
                np.testing.assert_array_equal(system["aparam"], expected)
                registered = next(
                    dt for dt in dpdata.System.DTYPES if dt.name == "aparam"
                )
                self.assertEqual(
                    registered.shape,
                    (Axis.NFRAMES, Axis.NATOMS, -1),
                )

    def test_reference_spin_key_maps_to_dpdata_name(self):
        frames = [
            {
                "atom_types": np.array([0, 0], dtype=np.int32),
                "coords": np.zeros((2, 3)),
                "spin": np.ones((2, 3)),
            }
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["H"])
        system = dpdata.MultiSystems.from_file(
            self.lmdb_path, fmt="lmdb", labeled=False
        )[0]
        self.assertIn("spins", system.data)
        self.assertNotIn("spin", system.data)
        np.testing.assert_array_equal(system["spins"], np.ones((1, 2, 3)))

    def test_known_field_shape_hint_cannot_remove_atom_axis(self):
        frames = [
            {
                "atom_types": np.array([0, 0, 0], dtype=np.int32),
                "coords": np.zeros((3, 3)),
                "spin": np.ones((3, 3)),
            }
        ]
        _write_raw_lmdb(
            self.lmdb_path,
            frames,
            ["H"],
            metadata_extra={
                "dp_data_shapes": {
                    "spin": ["nframes", 3, 3],
                }
            },
        )
        with self.assertRaisesRegex(LMDBMetadataError, "changes its protocol shape"):
            dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb", labeled=False)

    def test_data_name_hint_must_match_registered_protocol(self):
        frames = [
            {
                "atom_types": np.array([0], dtype=np.int32),
                "coords": np.zeros((1, 3)),
                "foo": np.ones((1, 3)),
            }
        ]
        _write_raw_lmdb(
            self.lmdb_path,
            frames,
            ["H"],
            metadata_extra={"dp_data_names": {"foo": "spins"}},
        )
        with self.assertRaisesRegex(LMDBMetadataError, "protocol key is 'spin'"):
            dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb", labeled=False)

    def test_custom_deepmd_core_key_rejected_on_read(self):
        frames = [
            {
                "atom_types": np.array([0], dtype=np.int32),
                "coords": np.zeros((1, 3)),
                "coord": np.full((1, 3), 101.0),
            }
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["H"])
        with self.assertRaisesRegex(LMDBFrameError, "reserved LMDB key 'coord'"):
            dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb", labeled=False)

    def test_unknown_atomic_axis_inferred_across_atom_counts(self):
        frames = [
            {
                "atom_types": np.array([0, 0], dtype=np.int32),
                "coords": np.zeros((2, 3)),
                "mystery": np.zeros((2, 1)),
            },
            {
                "atom_types": np.array([0, 0, 0], dtype=np.int32),
                "coords": np.ones((3, 3)),
                "mystery": np.ones((3, 1)),
            },
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["H"])
        systems = dpdata.MultiSystems.from_file(
            self.lmdb_path, fmt="lmdb", labeled=False
        )
        self.assertEqual(systems.get_nframes(), 2)
        dtype = next(dt for dt in dpdata.System.DTYPES if dt.name == "mystery")
        self.assertEqual(
            dtype.shape,
            (
                Axis.NFRAMES,
                Axis.NATOMS,
                1,
            ),
        )

    def test_unknown_same_nloc_atom_axis_is_rejected_as_ambiguous(self):
        frames = [
            {
                "atom_types": np.array([1, 0, 1], dtype=np.int32),
                "coords": np.zeros((3, 3)),
                "mystery_atomic": np.array([10.0, 20.0, 30.0]),
            }
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["H", "O"])
        with self.assertRaisesRegex(LMDBFrameError, "ambiguous without dp_data_shapes"):
            dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb", labeled=False)

    def test_process_global_schema_conflict_raises(self):
        second_path = "tmp_robust_second.lmdb"
        try:
            first_frames = [
                {
                    "atom_types": np.array([0], dtype=np.int32),
                    "coords": np.zeros((1, 3)),
                    "custom": np.array([1.0, 2.0]),
                }
            ]
            second_frames = [
                {
                    "atom_types": np.array([0], dtype=np.int32),
                    "coords": np.zeros((1, 3)),
                    "custom": np.array([1.0, 2.0, 3.0]),
                }
            ]
            _write_raw_lmdb(self.lmdb_path, first_frames, ["H"])
            _write_raw_lmdb(second_path, second_frames, ["H"])
            dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb", labeled=False)
            with self.assertRaisesRegex(LMDBError, "process-global definition"):
                dpdata.MultiSystems.from_file(second_path, fmt="lmdb", labeled=False)
        finally:
            if os.path.exists(second_path):
                shutil.rmtree(second_path)

    def test_failed_read_does_not_partially_register_data_types(self):
        frames = [
            {
                "atom_types": np.array([0], dtype=np.int32),
                "coords": np.zeros((1, 3)),
                "new_scalar": np.array(1.0),
            },
            {
                "atom_types": np.array([0, 0], dtype=np.int32),
                "coords": np.zeros((2, 3)),
                "forces": np.zeros((2, 3)),
            },
            {
                "atom_types": np.array([0, 0], dtype=np.int32),
                "coords": np.ones((2, 3)),
            },
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["H"])
        with self.assertRaisesRegex(LMDBFrameError, "must contain identical fields"):
            dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb", labeled=False)
        self.assertNotIn("new_scalar", [dt.name for dt in dpdata.System.DTYPES])

    def test_caller_construction_failure_rolls_back_registration(self):
        frames = [
            {
                "atom_types": np.array([0], dtype=np.int32),
                "coords": np.zeros((1, 3)),
                "unlabeled_custom": np.array(1.0),
            }
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["H"])
        with self.assertRaises(DataError):
            dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")
        self.assertNotIn(
            "unlabeled_custom",
            [dt.name for dt in dpdata.System.DTYPES],
        )
        self.assertNotIn(
            "unlabeled_custom",
            [dt.name for dt in dpdata.LabeledSystem.DTYPES],
        )

    def test_direct_labeled_read_rejects_unlabeled_before_registration(self):
        frames = [
            {
                "atom_types": np.array([0], dtype=np.int32),
                "coords": np.zeros((1, 3)),
                "unlabeled_custom": np.array(1.0),
            }
        ]
        _write_raw_lmdb(self.lmdb_path, frames, ["H"])
        with self.assertRaisesRegex(LMDBFrameError, "required by LabeledSystem"):
            dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        self.assertNotIn(
            "unlabeled_custom",
            [dt.name for dt in dpdata.System.DTYPES],
        )


class TestLMDBDumpSystems(unittest.TestCase):
    """dump_systems keeps every source system as a distinct frame_system_id."""

    def setUp(self):
        self.lmdb_path = "tmp_dump_systems.lmdb"
        self._remove_temporary_directories()
        # A and C share a formula (water); a plain MultiSystems would merge them.
        self.A = _make_labeled_system([0, 1, 1], ["O", "H"], 2)
        self.B = _make_labeled_system([0, 1, 1, 1, 1], ["C", "H"], 1)
        self.C = _make_labeled_system([0, 1, 1], ["O", "H"], 3)

    def tearDown(self):
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)
        self._remove_temporary_directories()

    def _remove_temporary_directories(self):
        for path in Path(".").glob(f".{self.lmdb_path}.tmp-*"):
            shutil.rmtree(path)
        for path in Path(".").glob(f".{self.lmdb_path}.backup-*"):
            shutil.rmtree(path)

    def test_system_ids_preserved(self):
        from dpdata.formats.lmdb import dump_systems

        dump_systems([self.A, self.B, self.C], self.lmdb_path)
        with lmdb.open(self.lmdb_path, readonly=True, lock=False) as env:
            with env.begin() as txn:
                meta = msgpack.unpackb(txn.get(b"__metadata__"), raw=False)
        self.assertEqual(meta["nframes"], 6)
        # three distinct source systems, in order, none merged
        self.assertEqual(meta["frame_system_ids"], [0, 0, 1, 2, 2, 2])
        self.assertEqual(meta["frame_nlocs"], [3, 3, 5, 3, 3, 3])
        # union type_map in first-appearance order
        self.assertEqual(meta["type_map"], ["O", "H", "C"])

    def test_contrast_with_multisystems_merge(self):
        # the default MultiSystems path merges A and C by formula.
        ms = dpdata.MultiSystems(self.A, self.B, self.C)
        ms.to("lmdb", self.lmdb_path)
        with lmdb.open(self.lmdb_path, readonly=True, lock=False) as env:
            with env.begin() as txn:
                meta = msgpack.unpackb(txn.get(b"__metadata__"), raw=False)
        # only two systems remain after formula merging
        self.assertEqual(sorted(set(meta["frame_system_ids"])), [0, 1])

    def test_roundtrip_total_frames(self):
        from dpdata.formats.lmdb import dump_systems

        dump_systems([self.A, self.B, self.C], self.lmdb_path)
        ms = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")
        self.assertEqual(ms.get_nframes(), 6)

    def test_generator_with_type_map(self):
        from dpdata.formats.lmdb import dump_systems

        def gen():
            yield self.A
            yield self.B
            yield self.C

        dump_systems(gen(), self.lmdb_path, type_map=["H", "C", "N", "O"])
        with lmdb.open(self.lmdb_path, readonly=True, lock=False) as env:
            with env.begin() as txn:
                meta = msgpack.unpackb(txn.get(b"__metadata__"), raw=False)
        self.assertEqual(meta["type_map"], ["H", "C", "N", "O"])
        self.assertEqual(meta["frame_system_ids"], [0, 0, 1, 2, 2, 2])

    def test_empty_input_rejected(self):
        from dpdata.formats.lmdb import dump_systems

        with self.assertRaisesRegex(LMDBError, "empty"):
            dump_systems([], self.lmdb_path)
        self.assertFalse(os.path.exists(self.lmdb_path))

    def test_empty_multisystems_rejected(self):
        with self.assertRaisesRegex(LMDBError, "empty"):
            dpdata.MultiSystems().to("lmdb", self.lmdb_path)
        self.assertFalse(os.path.exists(self.lmdb_path))

    def test_multisystems_input_rejected(self):
        from dpdata.formats.lmdb import dump_systems

        with self.assertRaisesRegex(TypeError, "original ordered systems"):
            dump_systems(dpdata.MultiSystems(self.A, self.B), self.lmdb_path)

    def test_negative_standard_atom_type_rejected(self):
        from dpdata.formats.lmdb import dump_systems

        data = self.A.data.copy()
        data["atom_types"] = np.array([-1, 1, 1])
        with self.assertRaisesRegex(LMDBError, "indices in"):
            dump_systems([data], self.lmdb_path)
        self.assertFalse(os.path.exists(self.lmdb_path))

    def test_floating_atom_type_rejected(self):
        from dpdata.formats.lmdb import dump_systems

        data = self.A.data.copy()
        data["atom_types"] = np.array([0.0, 1.0, 1.0])
        with self.assertRaisesRegex(LMDBError, "integer dtype"):
            dump_systems([data], self.lmdb_path)

    def test_uint64_overflow_not_treated_as_virtual_atom(self):
        from dpdata.formats.lmdb import dump_systems

        data = {
            "atom_numbs": [1],
            "atom_names": ["MIXED_TOKEN"],
            "atom_types": np.array([0]),
            "real_atom_names": ["H"],
            "real_atom_types": np.array([[np.iinfo(np.uint64).max]], dtype=np.uint64),
            "orig": np.zeros(3),
            "cells": np.eye(3)[None],
            "coords": np.zeros((1, 1, 3)),
        }
        with self.assertRaisesRegex(LMDBError, "indices in"):
            dump_systems([data], self.lmdb_path)

    def test_inconsistent_atom_numbs_rejected(self):
        from dpdata.formats.lmdb import dump_systems

        data = self.A.data.copy()
        data["atom_numbs"] = [2, 1]
        with self.assertRaisesRegex(LMDBError, "inconsistent"):
            dump_systems([data], self.lmdb_path)

    def test_nonportable_array_dtypes_rejected(self):
        from dpdata.formats.lmdb import dump_systems

        for dtype in (
            object,
            np.dtype([("x", np.float64)]),
        ):
            with self.subTest(dtype=dtype):
                data = self.A.data.copy()
                data["coords"] = np.zeros(self.A["coords"].shape, dtype=dtype)
                with self.assertRaisesRegex(LMDBError, "not a portable raw-byte dtype"):
                    dump_systems([data], self.lmdb_path)
                self.assertFalse(os.path.exists(self.lmdb_path))

    def test_mixed_raw_virtual_atoms_are_removed(self):
        from dpdata.formats.lmdb import dump_systems

        data = {
            "atom_numbs": [3],
            "atom_names": ["MIXED_TOKEN"],
            "atom_types": np.zeros(3, dtype=int),
            "real_atom_names": ["H", "O"],
            "real_atom_types": np.array([[1, 0, -1]], dtype=int),
            "orig": np.zeros(3),
            "cells": np.eye(3)[None],
            "coords": np.arange(9, dtype=float).reshape(1, 3, 3),
            "energies": np.array([-1.0]),
            "forces": np.arange(9, dtype=float).reshape(1, 3, 3),
        }
        dump_systems([data], self.lmdb_path)
        with lmdb.open(self.lmdb_path, readonly=True, lock=False) as env:
            with env.begin() as txn:
                meta = msgpack.unpackb(txn.get(b"__metadata__"), raw=False)
                frame = msgpack.unpackb(txn.get(b"000000000000"), raw=False)
        self.assertEqual(meta["type_map"], ["H", "O"])
        self.assertEqual(meta["frame_nlocs"], [2])
        self.assertEqual(frame["atom_types"]["shape"], [2])
        self.assertEqual(frame["coords"]["shape"], [2, 3])
        self.assertEqual(frame["forces"]["shape"], [2, 3])

    def test_ambiguous_mixed_custom_field_rejected(self):
        from dpdata.formats.lmdb import dump_systems

        data = {
            "atom_numbs": [3],
            "atom_names": ["MIXED_TOKEN"],
            "atom_types": np.zeros(3, dtype=int),
            "real_atom_names": ["H", "O"],
            "real_atom_types": np.array([[1, 0, -1]], dtype=int),
            "orig": np.zeros(3),
            "cells": np.eye(3)[None],
            "coords": np.zeros((1, 3, 3)),
            "mystery": np.zeros((1, 3, 2)),
        }
        with self.assertRaisesRegex(LMDBError, "shape is ambiguous"):
            dump_systems([data], self.lmdb_path)

    def test_existing_destination_requires_overwrite(self):
        from dpdata.formats.lmdb import dump_systems

        dump_systems([self.A], self.lmdb_path)
        with self.assertRaises(FileExistsError):
            dump_systems([self.B], self.lmdb_path)

    def test_failed_overwrite_preserves_old_database(self):
        from dpdata.formats.lmdb import dump_systems

        old = _make_labeled_system([0], ["H"], 1)
        old.data["energies"][:] = 10.0
        dump_systems([old], self.lmdb_path)

        replacement = _make_labeled_system([0], ["H"], 1)
        replacement.data["energies"][:] = 20.0
        invalid = _make_labeled_system([0], ["H"], 1).data.copy()
        invalid["atom_types"] = np.array([-1])
        with self.assertRaises(LMDBError):
            dump_systems(
                [replacement, invalid],
                self.lmdb_path,
                type_map=["H"],
                overwrite=True,
                write_batch_size=1,
            )
        loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        np.testing.assert_array_equal(loaded["energies"], [10.0])
        temporary = list(
            Path(self.lmdb_path).parent.glob(f".{Path(self.lmdb_path).name}.tmp-*")
        )
        self.assertEqual(temporary, [])

    def test_successful_overwrite_replaces_database(self):
        from dpdata.formats.lmdb import dump_systems

        dump_systems([self.A], self.lmdb_path)
        dump_systems(
            [self.B],
            self.lmdb_path,
            overwrite=True,
            write_batch_size=1,
        )
        loaded = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")
        self.assertEqual(loaded.get_nframes(), 1)
        self.assertEqual(loaded[0].get_natoms(), 5)
        with lmdb.open(self.lmdb_path, readonly=True, lock=False) as env:
            with env.begin() as txn:
                self.assertEqual(txn.stat()["entries"], 2)
                self.assertIsNone(txn.get(b"000000000001"))

    def test_duplicate_type_map_rejected(self):
        from dpdata.formats.lmdb import dump_systems

        with self.assertRaisesRegex(LMDBError, "duplicate"):
            dump_systems(
                [self.A],
                self.lmdb_path,
                type_map=["H", "H", "O"],
            )

    def test_invalid_write_batch_size_rejected_without_artifact(self):
        from dpdata.formats.lmdb import dump_systems

        with self.assertRaisesRegex(ValueError, "write_batch_size"):
            dump_systems(
                [self.A],
                self.lmdb_path,
                write_batch_size=0,
            )
        self.assertFalse(os.path.exists(self.lmdb_path))

    def test_map_full_grows_and_retries_batch(self):
        from dpdata.formats.lmdb import dump_systems

        large = _make_labeled_system([0] * 128, ["H"], 8)
        dump_systems(
            [large],
            self.lmdb_path,
            map_size=4096,
            write_batch_size=3,
        )
        loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        self.assertEqual(loaded.get_nframes(), 8)
        with lmdb.open(self.lmdb_path, readonly=True, lock=False) as env:
            self.assertGreater(env.info()["map_size"], 4096)

    def test_active_dpdata_reader_blocks_overwrite_for_path_alias(self):
        from dpdata.formats.lmdb import dump_systems

        dump_systems([self.A], self.lmdb_path)
        resolved, _ = _open_read_env(str(Path(self.lmdb_path).resolve()))
        try:
            with self.assertRaisesRegex(LMDBError, "reader is active"):
                dump_systems(
                    [self.B],
                    self.lmdb_path,
                    overwrite=True,
                )
        finally:
            _close_read_env(resolved)
        loaded = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")
        self.assertEqual(loaded.get_nframes(), self.A.get_nframes())

    def test_external_lmdb_reader_blocks_overwrite(self):
        from dpdata.formats.lmdb import dump_systems

        dump_systems([self.A], self.lmdb_path)
        external_env = lmdb.open(
            self.lmdb_path,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
        )
        try:
            with self.assertRaisesRegex(
                LMDBError, "Close DeePMD-kit and other readers"
            ):
                dump_systems(
                    [self.B],
                    self.lmdb_path,
                    overwrite=True,
                )
        finally:
            external_env.close()
        loaded = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")
        self.assertEqual(loaded.get_nframes(), self.A.get_nframes())

    def test_publish_guard_closes_external_reader_race(self):
        from dpdata.formats.lmdb import dump_systems

        dump_systems([self.A], self.lmdb_path)
        replace_entered = threading.Event()
        reader_finished = threading.Event()
        reader_was_blocked: list[bool] = []
        real_replace = os.replace

        def concurrent_reader():
            replace_entered.wait(timeout=5)
            try:
                env = lmdb.open(
                    self.lmdb_path,
                    readonly=True,
                    lock=False,
                    readahead=False,
                    meminit=False,
                )
            except lmdb.Error:
                reader_was_blocked.append(True)
            else:
                reader_was_blocked.append(False)
                env.close()
            finally:
                reader_finished.set()

        def synchronized_replace(source, destination):
            replace_entered.set()
            if not reader_finished.wait(timeout=5):
                raise TimeoutError("concurrent reader did not finish")
            return real_replace(source, destination)

        thread = threading.Thread(target=concurrent_reader)
        thread.start()
        try:
            with mock.patch(
                "dpdata.formats.lmdb.format.os.replace",
                side_effect=synchronized_replace,
            ):
                dump_systems(
                    [self.B],
                    self.lmdb_path,
                    overwrite=True,
                )
        finally:
            thread.join(timeout=5)
        self.assertEqual(reader_was_blocked, [True])
        loaded = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")
        self.assertEqual(loaded[0].get_natoms(), self.B.get_natoms())

    def test_staged_validation_failure_preserves_old_database(self):
        from dpdata.formats.lmdb import dump_systems
        from dpdata.formats.lmdb.format import _LMDBWriter

        dump_systems([self.A], self.lmdb_path)
        with (
            mock.patch.object(
                _LMDBWriter,
                "_validate_staged_database",
                side_effect=LMDBError("injected validation failure"),
            ),
            self.assertRaisesRegex(LMDBError, "injected"),
        ):
            dump_systems(
                [self.B],
                self.lmdb_path,
                overwrite=True,
            )
        loaded = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")
        self.assertEqual(loaded.get_nframes(), self.A.get_nframes())

    @unittest.skipUnless(hasattr(os, "fork"), "requires os.fork")
    def test_read_cache_resets_after_fork(self):
        from dpdata.formats.lmdb import dump_systems

        dump_systems([self.A], self.lmdb_path)
        resolved, _ = _open_read_env(self.lmdb_path)
        pid = os.fork()
        if pid == 0:
            try:
                child_resolved, _ = _open_read_env(self.lmdb_path)
                _close_read_env(child_resolved)
            except BaseException:
                os._exit(1)
            os._exit(0)
        _, status = os.waitpid(pid, 0)
        _close_read_env(resolved)
        self.assertEqual(os.waitstatus_to_exitcode(status), 0)


class TestLMDBErrorHandling(unittest.TestCase):
    def setUp(self):
        self.lmdb_missing_meta = "tmp_missing_meta.lmdb"
        self.lmdb_missing_frame = "tmp_missing_frame.lmdb"
        for p in (self.lmdb_missing_meta, self.lmdb_missing_frame):
            if os.path.exists(p):
                shutil.rmtree(p)

        env = lmdb.open(self.lmdb_missing_frame, map_size=1 << 30)
        with env.begin(write=True) as txn:
            meta = {
                "nframes": 1,
                "frame_idx_fmt": "012d",
                "type_map": ["H", "O"],
                "frame_nlocs": [3],
                "frame_system_ids": [0],
            }
            txn.put(b"__metadata__", _reference_packb(meta))
        env.close()

    def tearDown(self):
        for p in (self.lmdb_missing_meta, self.lmdb_missing_frame):
            if os.path.exists(p):
                shutil.rmtree(p)

    def test_missing_metadata(self):
        lmdb.open(self.lmdb_missing_meta, map_size=1 << 30).close()
        with self.assertRaises(LMDBMetadataError):
            dpdata.MultiSystems.from_file(self.lmdb_missing_meta, fmt="lmdb")

    def test_missing_frame(self):
        with self.assertRaises(LMDBFrameError):
            dpdata.MultiSystems.from_file(self.lmdb_missing_frame, fmt="lmdb")

    def test_non_mapping_metadata_rejected(self):
        path = "tmp_bad_metadata.lmdb"
        try:
            env = lmdb.open(path, map_size=1 << 30)
            with env.begin(write=True) as txn:
                txn.put(b"__metadata__", _reference_packb(["not", "a", "mapping"]))
            env.close()
            with self.assertRaisesRegex(LMDBMetadataError, "must contain a mapping"):
                dpdata.MultiSystems.from_file(path, fmt="lmdb")
        finally:
            if os.path.exists(path):
                shutil.rmtree(path)

    def test_invalid_msgpack_metadata_wrapped(self):
        path = "tmp_invalid_msgpack_metadata.lmdb"
        try:
            env = lmdb.open(path, map_size=1 << 30)
            with env.begin(write=True) as txn:
                txn.put(b"__metadata__", b"\xc1")
            env.close()
            with self.assertRaisesRegex(
                LMDBMetadataError, "Cannot decode __metadata__"
            ):
                dpdata.MultiSystems.from_file(path, fmt="lmdb")
        finally:
            if os.path.exists(path):
                shutil.rmtree(path)

    def test_fractional_nframes_rejected(self):
        path = "tmp_fractional_nframes.lmdb"
        try:
            env = lmdb.open(path, map_size=1 << 30)
            with env.begin(write=True) as txn:
                txn.put(
                    b"__metadata__",
                    _reference_packb(
                        {
                            "nframes": 1.9,
                            "frame_idx_fmt": "012d",
                            "type_map": ["H"],
                        }
                    ),
                )
            env.close()
            with self.assertRaisesRegex(
                LMDBMetadataError, "nframes must be an integer"
            ):
                dpdata.MultiSystems.from_file(path, fmt="lmdb")
        finally:
            if os.path.exists(path):
                shutil.rmtree(path)

    def test_non_integer_metadata_vectors_rejected(self):
        for key, value in (
            ("frame_nlocs", [1.5]),
            ("frame_system_ids", [False]),
        ):
            with self.subTest(key=key):
                path = f"tmp_invalid_{key}.lmdb"
                try:
                    frames = [
                        {
                            "atom_types": np.array([0], dtype=np.int32),
                            "coords": np.zeros((1, 3)),
                        }
                    ]
                    _write_raw_lmdb(
                        path,
                        frames,
                        ["H"],
                        metadata_extra={key: value},
                    )
                    with self.assertRaisesRegex(
                        LMDBMetadataError, "must be an integer"
                    ):
                        dpdata.MultiSystems.from_file(path, fmt="lmdb")
                finally:
                    if os.path.exists(path):
                        shutil.rmtree(path)

    def test_malformed_array_payload_rejected(self):
        path = "tmp_bad_array.lmdb"
        try:
            env = lmdb.open(path, map_size=1 << 30)
            with env.begin(write=True) as txn:
                frame = {
                    "atom_types": {
                        "type": "int32",
                        "shape": [2],
                        "data": b"\x00",
                    },
                    "atom_numbs": [2],
                }
                txn.put(b"000000000000", _reference_packb(frame))
                txn.put(
                    b"__metadata__",
                    _reference_packb(
                        {
                            "nframes": 1,
                            "frame_idx_fmt": "012d",
                            "type_map": ["H"],
                            "frame_nlocs": [2],
                        }
                    ),
                )
            env.close()
            with self.assertRaisesRegex(LMDBFrameError, "Cannot decode array"):
                dpdata.MultiSystems.from_file(path, fmt="lmdb")
        finally:
            if os.path.exists(path):
                shutil.rmtree(path)

    def test_invalid_frame_system_ids_length(self):
        path = "tmp_invalid_system_ids.lmdb"
        try:
            system = _make_labeled_system([0], ["H"], 1)
            system.to("lmdb", path)
            env = lmdb.open(path, map_size=1 << 30)
            with env.begin(write=True) as txn:
                metadata = msgpack.unpackb(txn.get(b"__metadata__"), raw=False)
                metadata["frame_system_ids"] = []
                txn.put(
                    b"__metadata__",
                    _reference_packb(metadata),
                )
            env.close()
            with self.assertRaisesRegex(LMDBMetadataError, "frame_system_ids length"):
                dpdata.MultiSystems.from_file(path, fmt="lmdb")
        finally:
            if os.path.exists(path):
                shutil.rmtree(path)

    def test_existing_external_reader_has_actionable_error(self):
        path = "tmp_external_reader.lmdb"
        try:
            _make_labeled_system([0], ["H"], 1).to("lmdb", path)
            env = lmdb.open(path, readonly=True, lock=False)
            try:
                with self.assertRaisesRegex(
                    LMDBError, "another library has already opened"
                ):
                    dpdata.MultiSystems.from_file(path, fmt="lmdb")
            finally:
                env.close()
        finally:
            if os.path.exists(path):
                shutil.rmtree(path)


class TestLMDBConfig(unittest.TestCase):
    def setUp(self):
        self.lmdb_path = "tmp_config.lmdb"
        self.system = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")

    def tearDown(self):
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)

    def test_custom_frame_idx_fmt(self):
        ms = dpdata.MultiSystems(self.system)
        ms.to("lmdb", self.lmdb_path, frame_idx_fmt="06d")
        with lmdb.open(self.lmdb_path, readonly=True, lock=False) as env:
            with env.begin() as txn:
                self.assertIsNotNone(txn.get(b"000000"))
                self.assertIsNone(txn.get(b"000000000000"))
        loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        self.assertEqual(len(loaded), len(self.system))
        np.testing.assert_allclose(loaded.data["coords"], self.system.data["coords"])

    def test_format_annotations_resolve_at_runtime(self):
        from dpdata.format import Format

        to_hints = typing.get_type_hints(Format.to_multi_systems)
        from_hints = typing.get_type_hints(Format.from_multi_systems)
        self.assertIn("return", to_hints)
        self.assertIn("return", from_hints)


if __name__ == "__main__":
    unittest.main()
