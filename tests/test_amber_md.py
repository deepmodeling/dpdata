import os
import unittest
import shutil
from context import dpdata
from comp_sys import CompLabeledSys, IsPBC
from comp_sys import MultiSystems
from comp_sys import IsNoPBC


class TestAmberMD(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp (self) :
        self.system_1 = dpdata.LabeledSystem('amber/02_Heat', fmt = 'amber/md')
        self.system_1.to('deepmd/npy','tmp.deepmd.npy')
        self.system_2 = dpdata.LabeledSystem('tmp.deepmd.npy', fmt = 'deepmd/npy')
        self.places = 5
        self.e_places = 4
        self.f_places = 6
        self.v_places = 6

    def tearDown(self) :
        if os.path.exists('tmp.deepmd.npy'):
            shutil.rmtree('tmp.deepmd.npy')



class TestAmberCorr(unittest.TestCase, MultiSystems, CompLabeledSys, IsNoPBC):
    def setUp(self):
        self.systems = self.get_correction("amber/corr/low_level",
                                           "amber/corr/high_level",
                                           "amber/corr/rc.nc",
                                           "amber/corr/qmmm.parm7",
                                           6)
        self.system_1 = self.systems['C6H11HW192O6OW96P1']
        self.system_2 = dpdata.LabeledSystem('amber/corr/dataset/C6H11HW192O6OW96P1', fmt = 'deepmd/npy')
        self.atom_names = ['C', 'H', 'HW', 'O', 'OW', 'P']
        self.system_names = ['C6H11HW192O6OW96P1']
        self.system_sizes = {'C6H11HW192O6OW96P1': 1}
        self.places = 5
        self.e_places = 4
        self.f_places = 6
        self.v_places = 6
    
    def get_correction(self, ll, hl, ncfile, parmfile, cutoff, target=":1"):
        """This function is proposed to get a correction energy and 
        forces from two Amber QM/MM systems, which will test the
        following functions:
        * load amber traj with use_element_symbols
        * load amber traj from mdout
        * LabeledSystem.correction
        * pick_by_amber_mask
        * remove_atom_names
        * set nopbc
        Reference: https://doi.org/10.26434/chemrxiv.14120447
        """
        ms = dpdata.MultiSystems()
        ep = r'@%EP'
        interactwith = "(%s)<:%f&!%s" % (target, cutoff, ep)

        s_ll = dpdata.LabeledSystem(
            ll, nc_file=ncfile, parm7_file=parmfile, fmt='amber/md', use_element_symbols=target)
        s_hl = dpdata.LabeledSystem(
            hl, nc_file=ncfile, parm7_file=parmfile, fmt='amber/md', use_element_symbols=target)
        s_corr = s_ll.correction(s_hl)
        s_corr = s_corr.pick_by_amber_mask(
            parmfile, interactwith, pass_coords=True, nopbc=True)
        for ss in s_corr:
            ss = ss.remove_atom_names('EP')
            ms.append(ss)
        return ms


if __name__ == '__main__':
    unittest.main()