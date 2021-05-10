import unittest
from context import dpdata
from comp_sys import CompLabeledSys
from comp_sys import IsPBC

class TestCorr(unittest.TestCase, CompLabeledSys, IsPBC):
    """Make a test to get a correction of two systems.

    Reference
    ---------
    https://doi.org/10.26434/chemrxiv.14120447
    """
    def setUp(self):
        ll="amber/corr/low_level"
        hl="amber/corr/high_level"
        ncfile="amber/corr/rc.nc"
        parmfile="amber/corr/qmmm.parm7"
        ep = r'@%EP'
        target = ":1"
        cutoff = 6.
        interactwith = "(%s)<:%f&!%s" % (target, cutoff, ep)
        s_ll = dpdata.LabeledSystem("amber/corr/dp_ll", fmt="deepmd/npy")
        s_hl = dpdata.LabeledSystem("amber/corr/dp_hl", fmt="deepmd/npy")
        self.system_1 = s_ll.correction(s_hl)
        self.system_2 = dpdata.LabeledSystem("amber/corr/dp_corr" ,fmt="deepmd/npy")
        self.places = 5
        self.e_places = 4
        self.f_places = 6
        self.v_places = 6

if __name__ == '__main__':
    unittest.main()