from dpdata.format import Format
import numpy as np

@Format.register("random")
class RandomFormat(Format):
    def from_system(self, N, **kwargs):
        return {
            "atom_numbs": [20],
            "atom_names": ['X'],
            "atom_types": [0] * 20,
            "cells": np.repeat(np.diag(np.diag(np.ones((3, 3))))[np.newaxis,...], N, axis=0) * 100.,
            "coords": np.random.rand(N, 20, 3) * 100.,
        }

    def from_labeled_system(self, N, **kwargs):
        return {
            "atom_numbs": [20],
            "atom_names": ['X'],
            "atom_types": [0] * 20,
            "cells": np.repeat(np.diag(np.diag(np.ones((3, 3))))[np.newaxis,...], N, axis=0) * 100.,
            "coords": np.random.rand(N, 20, 3) * 100.,
            "energies": np.random.rand(N) * 100.,
            "forces": np.random.rand(N, 20, 3) * 100.,
        }