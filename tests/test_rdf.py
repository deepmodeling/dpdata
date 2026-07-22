from __future__ import annotations

import unittest

import numpy as np

from dpdata.md.rdf import compute_rdf

try:
    import ase  # noqa: F401
except ModuleNotFoundError:
    skip_ase = True
else:
    skip_ase = False


@unittest.skipIf(skip_ase, "RDF calculation requires ASE")
class TestRDFSelectionOwnership(unittest.TestCase):
    def test_default_selection_is_recomputed_per_system(self):
        box = np.eye(3).reshape(1, 3, 3) * 10
        positions = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]])
        rdf = compute_rdf

        rdf(box, positions, np.array([0, 1]))
        self.assertEqual(rdf.__defaults__[0], None)
        _, values, _ = rdf(box, positions, np.array([2, 2]))
        self.assertTrue(np.isfinite(values).all())

    def test_caller_selection_is_not_mutated(self):
        box = np.eye(3).reshape(1, 3, 3) * 10
        positions = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]])
        selection = [None, None]
        compute_rdf(box, positions, np.array([0, 1]), selection)
        self.assertEqual(selection, [None, None])
