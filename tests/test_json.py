from __future__ import annotations

import datetime
import os
import tempfile
import unittest

import numpy as np
from comp_sys import CompLabeledSys, IsPBC
from context import dpdata

from dpdata.serialization import (
    _detect_format,
    dumpfn,
    loadfn,
    process_decoded,
    to_serializable,
)


class NestedSerializable:
    """Test helper that records values passed through ``from_dict``."""

    def __init__(self, value):
        self.value = value

    def as_dict(self):
        """Expose nested values through the same serialization hook as dpdata objects."""
        return {"value": self.value}

    @classmethod
    def from_dict(cls, data):
        """Construct the helper from decoded serialized data."""
        return cls(data["value"])


class TestJsonLoad(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")
        self.system_2 = dpdata.LabeledSystem.load("poscars/h2o.md.json")
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


class TestAsDict(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")
        self.system_2 = dpdata.LabeledSystem.from_dict(self.system_1.as_dict())
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


class TestJsonDumpLoad(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")
        self.tmpdir = tempfile.TemporaryDirectory()
        self.filename = os.path.join(self.tmpdir.name, "h2o.md.json")
        self.system_1.dump(self.filename)
        self.system_2 = dpdata.LabeledSystem.load(self.filename)
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

    def tearDown(self):
        self.tmpdir.cleanup()


class TestSerialization(unittest.TestCase):
    def test_detect_format_uses_final_meaningful_suffix(self):
        self.assertEqual(_detect_format("data.mpk.gz"), "mpk")
        self.assertEqual(_detect_format("data.yaml.bz2"), "yaml")
        self.assertEqual(_detect_format("data.yaml.json"), "json")
        self.assertEqual(_detect_format("data.mpk.backup"), "json")

    def test_from_dict_receives_decoded_nested_values(self):
        value = NestedSerializable(
            {
                "array": np.array([1, 2, 3]),
                "datetime": datetime.datetime(2026, 7, 18, 12, 0),
            }
        )

        decoded = process_decoded(to_serializable(value))

        np.testing.assert_array_equal(decoded.value["array"], value.value["array"])
        self.assertEqual(decoded.value["datetime"], value.value["datetime"])

    def test_timezone_aware_datetime_round_trip(self):
        for value in (
            datetime.datetime(
                2026,
                6,
                19,
                12,
                0,
                0,
                123456,
                tzinfo=datetime.timezone(datetime.timedelta(hours=8)),
            ),
            datetime.datetime(
                2026,
                6,
                19,
                12,
                0,
                0,
                tzinfo=datetime.timezone(datetime.timedelta(hours=-5)),
            ),
        ):
            self.assertEqual(process_decoded(to_serializable(value)), value)

    def test_yaml_dump_load(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(tmpdir, "data.yaml")
            value = {"numbers": [1, 2, 3]}
            dumpfn(value, filename)
            self.assertEqual(loadfn(filename), value)


if __name__ == "__main__":
    unittest.main()
