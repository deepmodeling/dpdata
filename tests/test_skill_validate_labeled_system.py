from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np

_SCRIPT = (
    Path(__file__).parents[1]
    / "skills"
    / "dpdata-inspect"
    / "scripts"
    / "validate_labeled_system.py"
)
_SPEC = importlib.util.spec_from_file_location(
    "dpdata_validate_labeled_system", _SCRIPT
)
assert _SPEC is not None and _SPEC.loader is not None
_VALIDATOR = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _VALIDATOR
_SPEC.loader.exec_module(_VALIDATOR)


class TestValidateLabeledSystem(unittest.TestCase):
    def _data(self):
        return {
            "atom_numbs": [1, 1],
            "atom_names": ["H", "O"],
            "atom_types": np.array([0, 1], dtype=int),
            "orig": np.zeros(3),
            "cells": np.repeat(np.eye(3)[None, ...], 2, axis=0),
            "coords": np.zeros((2, 2, 3)),
            "energies": np.array([-1.0, -0.9]),
            "forces": np.zeros((2, 2, 3)),
            "virials": np.zeros((2, 3, 3)),
        }

    def test_valid_labeled_data(self):
        report = _VALIDATOR.validate_data(
            self._data(),
            expected_type_map=["H", "O"],
            energy_unit="eV",
            force_unit="eV/A",
            virial_unit="eV",
            stress_sign="virial",
            provenance=__file__,
            producer="test",
            producer_version="1",
            parser_version="1",
            strict_metadata=True,
        )
        self.assertTrue(report.ok, report.errors)
        self.assertEqual(report.summary["frames"], 2)

    def test_missing_energy_label_is_an_error(self):
        data = self._data()
        del data["energies"]
        report = _VALIDATOR.validate_data(data)
        self.assertFalse(report.ok)
        self.assertIn("energies: missing required energy labels", report.errors)

    def test_alignment_and_atom_type_errors_are_reported(self):
        data = self._data()
        data["forces"] = np.zeros((1, 2, 3))
        data["atom_types"] = np.array([0, 2], dtype=int)
        report = _VALIDATOR.validate_data(data)
        self.assertFalse(report.ok)
        self.assertTrue(any("forces: shape" in error for error in report.errors))
        self.assertIn("atom_types: index is outside atom_names/type_map", report.errors)

    def test_provenance_check_is_independent_of_parser_version(self):
        report = _VALIDATOR.validate_data(
            self._data(), provenance=__file__, parser_version=None
        )
        self.assertTrue(report.ok, report.errors)

        report = _VALIDATOR.validate_data(
            self._data(), provenance=None, parser_version="1"
        )
        self.assertTrue(report.ok, report.errors)

    def test_empty_cells_are_rejected(self):
        data = self._data()
        data["cells"] = np.empty((0,))
        report = _VALIDATOR.validate_data(data, require_pbc=True)
        self.assertFalse(report.ok)
        self.assertTrue(any("cells: shape" in error for error in report.errors))

    def test_multi_system_type_map_must_match(self):
        first = type("System", (), {"data": self._data()})()
        second_data = self._data()
        second_data["atom_names"] = ["O", "H"]
        second = type("System", (), {"data": second_data})()
        reports = _VALIDATOR.validate_systems([first, second])
        self.assertFalse(reports[1].ok)
        self.assertTrue(
            any("multi-system dataset" in error for error in reports[1].errors)
        )


if __name__ == "__main__":
    unittest.main()
