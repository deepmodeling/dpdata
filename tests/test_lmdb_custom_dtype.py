from __future__ import annotations

import os
import shutil
import unittest

import numpy as np

import dpdata
from dpdata.data_type import Axis, DataType


class TestLMDBCustomDType(unittest.TestCase):
    def setUp(self):
        self.original_system_dtypes = dpdata.System.DTYPES
        self.original_labeled_system_dtypes = dpdata.LabeledSystem.DTYPES

        # Register custom data types as optional
        self.dt_frame = DataType(
            "frame_data", np.ndarray, shape=(Axis.NFRAMES, 2), required=False
        )
        self.dt_static = DataType("static_data", np.ndarray, shape=(2,), required=False)

        dpdata.System.register_data_type(self.dt_frame, self.dt_static)
        dpdata.LabeledSystem.register_data_type(self.dt_frame, self.dt_static)

        self.lmdb_path = "tmp_custom_dtype.lmdb"

        # Create a system with custom data
        # Assuming running from tests/ directory
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
        self.system.data["static_data"] = np.array([1.0, 2.0])
        self.system.check_data()

    def tearDown(self):
        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)

    def test_custom_dtype_preservation(self):
        self.system.to("lmdb", self.lmdb_path)
        system_loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")

        np.testing.assert_allclose(
            system_loaded.data["frame_data"], self.system.data["frame_data"]
        )
        np.testing.assert_allclose(
            system_loaded.data["static_data"], self.system.data["static_data"]
        )

    def test_multi_systems_custom_dtype(self):
        ms = dpdata.MultiSystems(self.system)
        ms.to("lmdb", self.lmdb_path)
        ms_loaded = dpdata.MultiSystems.from_file(self.lmdb_path, fmt="lmdb")

        system_loaded = list(ms_loaded.systems.values())[0]
        np.testing.assert_allclose(
            system_loaded.data["frame_data"], self.system.data["frame_data"]
        )
        np.testing.assert_allclose(
            system_loaded.data["static_data"], self.system.data["static_data"]
        )

    def test_custom_dtype_auto_registration(self):
        # Save with custom data types registered
        self.system.to("lmdb", self.lmdb_path)

        # Simulate a clean session by unregistering the custom types
        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes

        # Verify they are currently missing
        self.assertNotIn("frame_data", [dt.name for dt in dpdata.LabeledSystem.DTYPES])

        # Load from LMDB - should trigger auto-registration
        system_loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")

        # Verify data is loaded and types are registered
        self.assertIn("frame_data", [dt.name for dt in dpdata.LabeledSystem.DTYPES])
        self.assertIn("static_data", [dt.name for dt in dpdata.LabeledSystem.DTYPES])

        np.testing.assert_allclose(
            system_loaded.data["frame_data"], self.system.data["frame_data"]
        )
        np.testing.assert_allclose(
            system_loaded.data["static_data"], self.system.data["static_data"]
        )


class TestLMDBFparamAparam(unittest.TestCase):
    def setUp(self):
        self.original_system_dtypes = dpdata.System.DTYPES
        self.original_labeled_system_dtypes = dpdata.LabeledSystem.DTYPES

        new_datatypes = [
            DataType(
                "fparam",
                np.ndarray,
                shape=(Axis.NFRAMES, 2),
                required=False,
            ),
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
        system_loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")

        np.testing.assert_allclose(
            system_loaded.data["fparam"], self.system.data["fparam"]
        )
        np.testing.assert_allclose(
            system_loaded.data["aparam"], self.system.data["aparam"]
        )

    def test_fparam_aparam_auto_registration(self):
        # Save with fparam/aparam registered
        self.system.to("lmdb", self.lmdb_path)

        # Simulate a clean session by restoring original DTYPES
        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes

        # Load from LMDB
        system_loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")

        # Verify auto-registration and data correctness
        self.assertIn("fparam", [dt.name for dt in dpdata.LabeledSystem.DTYPES])
        self.assertIn("aparam", [dt.name for dt in dpdata.LabeledSystem.DTYPES])

        np.testing.assert_allclose(
            system_loaded.data["fparam"], self.system.data["fparam"]
        )
        np.testing.assert_allclose(
            system_loaded.data["aparam"], self.system.data["aparam"]
        )

    def test_symbolic_axis_natoms_preservation(self):
        # 1. Save system with aparam (which uses Axis.NATOMS)
        self.system.to("lmdb", self.lmdb_path)

        # 2. Simulate new session
        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes

        # 3. Load triggers auto-registration
        dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")

        # 4. Find the newly registered DataType for 'aparam'
        aparam_dt = next(
            dt for dt in dpdata.LabeledSystem.DTYPES if dt.name == "aparam"
        )

        # 5. Assert that it contains Axis.NATOMS, not a fixed integer
        self.assertIn(Axis.NATOMS, aparam_dt.shape)

        # 6. Functional verification
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


if __name__ == "__main__":
    unittest.main()
