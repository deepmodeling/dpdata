from __future__ import annotations

import unittest

import numpy as np
from context import dpdata

class TestABACUSDeltaSpin(unittest.TestCase):
    """Make a test for the DeltaSpin module in ABACUS"""

    def setUp(self):
        
        self.system_1 = dpdata.LabeledSystem(
            ("./abacus.scf.deltaspin"), fmt="abacus/scf"
        )

        self.system_2 = dpdata.LabeledSystem(
            data={
                "atom_types": np.array([0, 1, 1, 1]),
                "atom_names": ["Fe"],
                "atom_numbs": [2],
                "coords": np.array(
                    [
                        [
                            [0.       , 0.       , 0.       ],
                            [1.4349999, 1.4349999, 1.4349999],
                        ]
                    ]
                ),
                "energies": np.array([-6437.2128405]),
                "forces": np.array(
                    [
                        [
                            [0., 0., 0.,],
                            [0., 0., 0.,],
                        ]
                    ]
                ),
                "spin": np.array
                (
                    [
                        [
                            [0.,         0.,         2.19999997],
                            [0.,         0.,         2.20000001],
                        ]
                    ]
                ),
                "cells": np.array
                (
                    [
                        [
                            [2.86999979, 0.        ,0.        ],
                            [0.,         2.86999979, 0.        ],
                            [0.,         0.        , 2.86999979],
                        ]
                    ]
                ),
                "virial": np.array(
                    [
                        [
                            [-3.95996034e-01, -3.25934940e-09, -4.71565445e-09],
                            [-3.25934940e-09, -3.95996551e-01,  1.92993618e-09],
                            [-4.71565445e-09,  1.92993618e-09, -3.95998805e-01],
                        ]
                    ]
                ),
                "mag_forces": np.array(
                    [
                        [
                            [ 0.00000000e+00, -1.12095051e-09,  1.80736102e-09],      
                            [ 1.00000000e+00,  3.54604018e-09, -2.33608701e-10],
                        ]
                    ]
                ),
                "coors_deltaspin": np.array(
                    [
                        [
                            [ 0.       ,   0.       ,   0.       ,   0.       ,   0.       , -0.44295301 ],
                            [ 1.4349999,   1.4349999,   1.4349999,   1.4349999,   1.4349999, 0.99204688  ],
                        ]
                    ]
                )
            }
        )
        

if __name__ == "__main__":
    unittest.main()