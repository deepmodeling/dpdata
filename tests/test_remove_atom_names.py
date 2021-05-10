import unittest
from context import dpdata
from comp_sys import CompLabeledSys
from comp_sys import IsNoPBC

class TestRemove(unittest.TestCase, CompLabeledSys, IsNoPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem("amber/corr/dp_amber_mask/C6EP0H11HW192O6OW96P1", fmt="deepmd/npy").remove_atom_names('EP')
        self.system_2 = dpdata.LabeledSystem("amber/corr/dataset/C6H11HW192O6OW96P1", fmt="deepmd/npy")
        self.places = 5
        self.e_places = 4
        self.f_places = 6
        self.v_places = 6

if __name__ == '__main__':
    unittest.main()