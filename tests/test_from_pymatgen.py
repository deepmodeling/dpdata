import os
import unittest

from comp_sys import CompSys, IsPBC
from context import dpdata

try:
    from pymatgen.core import Structure  # noqa: F401

    exist_module = True
except Exception:
    exist_module = False


@unittest.skipIf(not exist_module, "skip pymatgen")
class TestFormPytmatgen(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        structure = Structure.from_file(os.path.join("poscars", "POSCAR.P42nmc"))
        self.system_1 = dpdata.System(structure, fmt="pymatgen/structure")
        self.system_2 = dpdata.System(
            os.path.join("poscars", "POSCAR.P42nmc"), fmt="poscar"
        )
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


if __name__ == "__main__":
    unittest.main()
