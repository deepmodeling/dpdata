from __future__ import annotations

import datetime
import importlib
import os
import shutil
import sys
import tempfile
import unittest
import uuid
from enum import Enum
from pathlib import Path
from unittest import mock

import numpy as np
from context import dpdata  # noqa: F401

from dpdata.serialization import (
    _detect_format,
    _open_binary,
    _open_text,
    dumpfn,
    loadfn,
    process_decoded,
    to_serializable,
)


class Color(Enum):
    RED = "red"
    BLUE = "blue"


class TestDetectFormat(unittest.TestCase):
    def test_explicit_format_wins(self):
        self.assertEqual(_detect_format("data.json", fmt="yaml"), "yaml")

    def test_suffix_dispatch(self):
        for name, expected in (
            ("data.json", "json"),
            ("data.mpk", "mpk"),
            ("data.yaml", "yaml"),
            ("data.yml", "yaml"),
            ("data", "json"),
            ("data.unknown", "json"),
        ):
            with self.subTest(name=name):
                self.assertEqual(_detect_format(name), expected)

    def test_compression_suffix_is_stripped_first(self):
        for name, expected in (
            ("data.yaml.gz", "yaml"),
            ("data.mpk.bz2", "mpk"),
            ("data.json.z", "json"),
            ("data.gz", "json"),
        ):
            with self.subTest(name=name):
                self.assertEqual(_detect_format(name), expected)


class TestCompressedStreams(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def _path(self, name):
        return os.path.join(self.tmp_dir, name)

    def test_text_roundtrip_through_each_codec(self):
        for name in ("plain.txt", "gzipped.gz", "packed.z", "squeezed.bz2"):
            with self.subTest(name=name):
                path = self._path(name)
                with _open_text(path, "wt") as fp:
                    fp.write("hello")
                with _open_text(path, "rt") as fp:
                    self.assertEqual(fp.read(), "hello")

    def test_binary_roundtrip_through_each_codec(self):
        for name in ("plain.bin", "gzipped.gz", "packed.z", "squeezed.bz2"):
            with self.subTest(name=name):
                path = self._path(name)
                with _open_binary(path, "wb") as fp:
                    fp.write(b"\x00\x01")
                with _open_binary(path, "rb") as fp:
                    self.assertEqual(fp.read(), b"\x00\x01")

    def test_compressed_files_are_not_plain_text(self):
        path = self._path("compressed.json.gz")
        dumpfn({"a": 1}, path)
        with open(path, "rb") as fp:
            self.assertEqual(fp.read(2), b"\x1f\x8b")
        self.assertEqual(loadfn(path), {"a": 1})


class TestEncoding(unittest.TestCase):
    def test_complex_array_roundtrip(self):
        arr = np.array([1 + 2j, 3 - 4j])
        encoded = to_serializable(arr)
        self.assertEqual(encoded["dtype"], str(arr.dtype))
        self.assertEqual(len(encoded["data"]), 2)  # real and imaginary parts
        np.testing.assert_allclose(process_decoded(encoded), arr)

    def test_real_array_roundtrip(self):
        arr = np.arange(6, dtype=np.float32).reshape(2, 3)
        np.testing.assert_allclose(process_decoded(to_serializable(arr)), arr)

    def test_numpy_scalar_becomes_python_scalar(self):
        value = to_serializable(np.float64(1.5))
        self.assertIsInstance(value, float)
        self.assertNotIsInstance(value, np.generic)

    def test_uuid_roundtrip(self):
        value = uuid.uuid4()
        self.assertEqual(process_decoded(to_serializable(value)), value)

    def test_path_roundtrip(self):
        value = Path("a") / "b.txt"
        self.assertEqual(process_decoded(to_serializable(value)), value)

    def test_enum_roundtrip(self):
        encoded = to_serializable(Color.BLUE)
        self.assertEqual(encoded["@class"], "Color")
        self.assertIs(process_decoded(encoded), Color.BLUE)

    def test_unserializable_object_passes_through(self):
        marker = object()
        self.assertIs(to_serializable(marker), marker)

    def test_dict_keys_are_converted_too(self):
        encoded = to_serializable({np.int64(3): np.int64(4)})
        self.assertEqual(encoded, {3: 4})


class TestDatetimeDecoding(unittest.TestCase):
    def test_aware_datetime_roundtrip(self):
        value = datetime.datetime(2026, 7, 27, 12, 30, tzinfo=datetime.timezone.utc)
        self.assertEqual(process_decoded(to_serializable(value)), value)

    def test_negative_offset_roundtrip(self):
        value = datetime.datetime(
            2026, 7, 27, 12, 30, tzinfo=datetime.timezone(-datetime.timedelta(hours=5))
        )
        self.assertEqual(process_decoded(to_serializable(value)), value)

    def test_non_iso_payload_with_microseconds(self):
        # `fromisoformat` demands zero-padded fields; the strptime fallbacks
        # keep hand-edited or legacy payloads readable.
        decoded = process_decoded(
            {
                "@module": "datetime",
                "@class": "datetime",
                "string": "2020-1-2 3:4:5.678900",
            }
        )
        self.assertEqual(decoded, datetime.datetime(2020, 1, 2, 3, 4, 5, 678900))

    def test_non_iso_payload_without_microseconds(self):
        decoded = process_decoded(
            {
                "@module": "datetime",
                "@class": "datetime",
                "string": "2020-1-2 3:4:5",
            }
        )
        self.assertEqual(decoded, datetime.datetime(2020, 1, 2, 3, 4, 5))

    def test_non_iso_payload_with_offset_drops_the_offset(self):
        # The fallback splits on "+", so an offset it cannot parse is lost
        # rather than raising.
        decoded = process_decoded(
            {
                "@module": "datetime",
                "@class": "datetime",
                "string": "2020-1-2 3:4:5+00:00",
            }
        )
        self.assertEqual(decoded, datetime.datetime(2020, 1, 2, 3, 4, 5))

    def test_iso_payload_with_space_separator(self):
        decoded = process_decoded(
            {
                "@module": "datetime",
                "@class": "datetime",
                "string": "2020-01-02 03:04:05.678900",
            }
        )
        self.assertEqual(decoded, datetime.datetime(2020, 1, 2, 3, 4, 5, 678900))


class TestUnknownMarkers(unittest.TestCase):
    def test_unimportable_module_is_left_as_a_dict(self):
        payload = {"@module": "no.such.module", "@class": "Thing", "value": 1}
        self.assertEqual(process_decoded(payload), payload)

    def test_missing_class_is_left_as_a_dict(self):
        payload = {"@module": "datetime", "@class": "NotAClass", "value": 1}
        self.assertEqual(process_decoded(payload), payload)

    def test_class_without_from_dict_is_left_as_a_dict(self):
        payload = {"@module": "builtins", "@class": "object", "value": 1}
        self.assertEqual(process_decoded(payload), payload)

    def test_nested_markers_inside_lists_are_decoded(self):
        arr = np.arange(3)
        decoded = process_decoded([to_serializable(arr), {"k": to_serializable(arr)}])
        np.testing.assert_allclose(decoded[0], arr)
        np.testing.assert_allclose(decoded[1]["k"], arr)


class TestFileFormats(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.payload = {
            "array": np.arange(4, dtype=np.float64),
            "when": datetime.datetime(2026, 7, 27, 8, 0),
            "nested": {"list": [1, 2, 3]},
        }

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def _path(self, name):
        return os.path.join(self.tmp_dir, name)

    def _assert_roundtrip(self, path):
        dumpfn(self.payload, path)
        loaded = loadfn(path)
        np.testing.assert_allclose(loaded["array"], self.payload["array"])
        self.assertEqual(loaded["when"], self.payload["when"])
        self.assertEqual(loaded["nested"], self.payload["nested"])

    def test_json_roundtrip(self):
        self._assert_roundtrip(self._path("data.json"))

    def test_yaml_roundtrip(self):
        self._assert_roundtrip(self._path("data.yaml"))

    def test_msgpack_roundtrip(self):
        self._assert_roundtrip(self._path("data.mpk"))

    def test_compressed_msgpack_roundtrip(self):
        self._assert_roundtrip(self._path("data.mpk.gz"))

    def test_invalid_format_on_dump(self):
        with self.assertRaisesRegex(TypeError, "Invalid format: toml"):
            dumpfn({"a": 1}, self._path("data.toml"), fmt="toml")

    def test_invalid_format_on_load(self):
        path = self._path("data.json")
        dumpfn({"a": 1}, path)
        with self.assertRaisesRegex(TypeError, "Invalid format: toml"):
            loadfn(path, fmt="toml")


def _import_without(*blocked):
    """Return an ``import_module`` that pretends ``blocked`` is not installed."""
    real = importlib.import_module

    def fake(name, *args, **kwargs):
        if name in blocked:
            raise ModuleNotFoundError(f"No module named {name!r}")
        return real(name, *args, **kwargs)

    return fake


class TestOptionalDependencies(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def _path(self, name):
        return os.path.join(self.tmp_dir, name)

    def test_yaml_falls_back_to_ruamel(self):
        path = self._path("data.yaml")
        with mock.patch("importlib.import_module", _import_without("yaml")):
            dumpfn({"a": [1, 2]}, path, indent=2)
            with mock.patch("importlib.import_module", _import_without("yaml")):
                loaded = loadfn(path)
        self.assertEqual(loaded, {"a": [1, 2]})

    def test_yaml_without_any_backend_explains_itself(self):
        blocked = _import_without("yaml", "ruamel.yaml")
        with mock.patch("importlib.import_module", blocked):
            with self.assertRaisesRegex(RuntimeError, "requires PyYAML or ruamel"):
                dumpfn({"a": 1}, self._path("data.yaml"))

        path = self._path("readable.yaml")
        dumpfn({"a": 1}, path)
        with mock.patch("importlib.import_module", blocked):
            with self.assertRaisesRegex(RuntimeError, "requires PyYAML or ruamel"):
                loadfn(path)

    def test_msgpack_absence_explains_itself(self):
        path = self._path("data.mpk")
        with mock.patch.dict(sys.modules, {"msgpack": None}):
            with self.assertRaisesRegex(RuntimeError, "requires msgpack"):
                dumpfn({"a": 1}, path)

        dumpfn({"a": 1}, path)
        with mock.patch.dict(sys.modules, {"msgpack": None}):
            with self.assertRaisesRegex(RuntimeError, "requires msgpack"):
                loadfn(path)


if __name__ == "__main__":
    unittest.main()
