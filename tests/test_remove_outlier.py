import os
import unittest

import numpy as np
from comp_sys import CompLabeledSys
from context import dpdata


class TestRemoveOutlier(unittest.TestCase, CompLabeledSys):
    @classmethod
    def setUpClass(cls):
        system = dpdata.LabeledSystem(
            data={
                "atom_names": ["H"],
                "atom_numbs": [1],
                "atom_types": np.zeros((1,), dtype=int),
                "coords": np.zeros((100, 1, 3), dtype=np.float32),
                "cells": np.zeros((100, 3, 3), dtype=np.float32),
                "orig": np.zeros(3, dtype=np.float32),
                "nopbc": True,
                "energies": np.zeros((100,), dtype=np.float32),
                "forces": np.zeros((100, 1, 3), dtype=np.float32),
            }
        )
        system.data["energies"][0] = 100.0
        system.data["energies"][1] = -100.0
        cls.system_1 = system.remove_outlier()
        cls.system_2 = system[2:]
        cls.places = 6
        cls.e_places = 6
        cls.f_places = 6
        cls.v_places = 6
