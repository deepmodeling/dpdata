import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys, IsNoPBC
try:
   import parmed
   exist_module=True
except:
   exist_module=False

class TestPickAtomIdx(unittest.TestCase, CompSys, IsNoPBC):
    
    def setUp(self): 
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6
        self.system_1 = dpdata.LabeledSystem('gaussian/methane_reordered.gaussianlog', fmt='gaussian/log').pick_atom_idx(slice(4))
        self.system_2 = dpdata.LabeledSystem('gaussian/methane_sub.gaussianlog', fmt='gaussian/log')

@unittest.skipIf(not exist_module,"skip")
class TestPickByAmberMask(unittest.TestCase, CompSys, IsNoPBC):
    
    def setUp(self): 
        parmfile="amber/corr/qmmm.parm7"
        ep = r'@%EP'
        target = ":1"
        cutoff = 6.
        interactwith = "(%s)<:%f&!%s" % (target, cutoff, ep)
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6
        self.system_1 = dpdata.LabeledSystem("amber/corr/dp_corr", fmt="deepmd/npy").pick_by_amber_mask(
                            parmfile, interactwith, pass_coords=True, nopbc=True)['C6EP0H11HW192O6OW96P1']
        self.system_2 = dpdata.LabeledSystem("amber/corr/dp_amber_mask/C6EP0H11HW192O6OW96P1", fmt="deepmd/npy")
    

if __name__ == '__main__':
    unittest.main()
    
