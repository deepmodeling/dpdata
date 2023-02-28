# The index should map to that in the dump file

import os
import unittest

import numpy as np
from context import dpdata


class TestLmpDumpIdx(unittest.TestCase):
    def setUp(self):
        self.system = dpdata.System(os.path.join("poscars", "conf2.dump"))

    def test_coords(self):
        np.testing.assert_allclose(
            self.system["coords"],
            np.array([[[0.0, 0.0, 0.0], [1.2621856, 0.7018028, 0.5513885]]]),
        )

    def test_type(self):
        np.testing.assert_allclose(
            self.system.get_atom_types(),
            np.array(
                [1, 0],
                dtype=int,
            ),
        )
