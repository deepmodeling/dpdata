import unittest

import numpy as np
from comp_sys import CompLabeledSys, IsNoPBC
from context import dpdata


class TestDeepmdLoadDumpHDF5(unittest.TestCase, CompLabeledSys, IsNoPBC):
    def setUp(self):
        length_convert = dpdata.unit.LengthConversion("bohr", "angstrom").value()
        energy_convert = dpdata.unit.EnergyConversion("hartree", "eV").value()
        force_convert = dpdata.unit.ForceConversion(
            "hartree/bohr", "eV/angstrom"
        ).value()

        self.system_1 = dpdata.LabeledSystem("psi4/psi4.out", fmt="psi4/out")

        self.system_2 = dpdata.LabeledSystem(
            data={
                "atom_types": np.array([0, 0, 1, 1, 1, 1, 1, 1]),
                "atom_names": ["C", "H"],
                "atom_numbs": [2, 6],
                "coords": np.array(
                    [
                        [
                            [1.309059187335, -0.530960676560, 0.283395850372],
                            [-1.305263812665, 0.530120323440, -0.297504149628],
                            [2.346254187335, -1.337363676560, -1.334794149628],
                            [0.931624187335, -2.053154676560, 1.656377850372],
                            [2.513533187335, 0.904343323440, 1.196493850372],
                            [-2.757838812665, -0.942676676560, -0.553430149628],
                            [-1.178330812665, 1.650444323440, -2.050622149628],
                            [-1.900432812665, 1.788413323440, 1.253959850372],
                        ]
                    ]
                )
                * length_convert,
                "energies": np.array([-79.8685493165660603]) * energy_convert,
                "forces": -np.array(
                    [
                        [
                            [0.00189577378438, 0.01178186689846, -0.01517052269765],
                            [-0.00054675643432, -0.01239517767892, 0.01285520444405],
                            [0.00862382497882, -0.00438603405641, -0.00576291289370],
                            [-0.01373063001962, -0.00368703316336, 0.00313307980576],
                            [0.00439957658795, 0.00725213801722, 0.00516826201141],
                            [-0.00831692511314, -0.00614283210761, -0.00048696830158],
                            [0.00755493258543, 0.00167237971637, -0.00777559049988],
                            [0.00011879620295, 0.00590450771644, 0.00804206420271],
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
