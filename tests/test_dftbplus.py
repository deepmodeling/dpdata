from __future__ import annotations

import unittest

import numpy as np
from comp_sys import CompLabeledSys, IsNoPBC
from context import dpdata


class TestDeepmdLoadAmmonia(unittest.TestCase, CompLabeledSys, IsNoPBC):
    def setUp(self):
        energy_convert = dpdata.unit.EnergyConversion("hartree", "eV").value()
        force_convert = dpdata.unit.ForceConversion(
            "hartree/bohr", "eV/angstrom"
        ).value()

        self.system_1 = dpdata.LabeledSystem(
            ("dftbplus/dftb_pin.hsd", "dftbplus/detailed.out"), fmt="dftbplus"
        )

        self.system_2 = dpdata.LabeledSystem(
            data={
                "atom_types": np.array([0, 1, 1, 1]),
                "atom_names": ["N", "H"],
                "atom_numbs": [1, 3],
                "coords": np.array(
                    [
                        [
                            [1.014150, 0.112320, 0.047370],
                            [3.909390, 0.037985, -0.101159],
                            [0.702550, -0.851820, -0.060860],
                            [0.702550, 0.603740, -0.789160],
                        ]
                    ]
                ),
                "energies": np.array([-3.2963983884]) * energy_convert,
                "forces": np.array(
                    [
                        [
                            [0.016567056203, 0.002817951422, 0.005634574270],
                            [-0.018803818530, -0.000002880649, -0.000006015442],
                            [0.001118562874, -0.005291070259, -0.000870711110],
                            [0.001118199454, 0.002475999486, -0.004757847718],
                        ]
                    ]
                )
                * force_convert,
                "cells": np.zeros((1, 3, 3)),
                "orig": np.zeros(3),
                "nopbc": True,
            }
        )
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


if __name__ == "__main__":
    unittest.main()
