# %%
import os
import unittest

import numpy as np
from comp_sys import CompLabeledSys
from context import dpdata


# %%
class TestCp2kAimdOutput(unittest.TestCase, CompLabeledSys):
    def setUp(self):
<<<<<<< HEAD
        self.system_1 = dpdata.LabeledSystem('cp2k/aimd',fmt='cp2k/aimd_output')
        self.system_2 = dpdata.LabeledSystem('cp2k/aimd/deepmd', fmt='deepmd/npy')
=======
        self.system_1 = dpdata.LabeledSystem("cp2k/aimd", fmt="cp2k/aimd_output")
        self.system_2 = dpdata.LabeledSystem("cp2k/aimd/deepmd", fmt="deepmd/npy")
>>>>>>> upstream/devel
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4

<<<<<<< HEAD
#class TestCp2kAimdRestartOutput(unittest.TestCase, CompLabeledSys):
=======

class TestCp2kAimdStressOutput(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem("cp2k/aimd_stress", fmt="cp2k/aimd_output")
        self.system_2 = dpdata.LabeledSystem(
            "cp2k/aimd_stress/deepmd", fmt="deepmd/raw"
        )
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


# class TestCp2kAimdRestartOutput(unittest.TestCase, CompLabeledSys):
>>>>>>> upstream/devel
#    def setUp(self):
#        self.system_1 = dpdata.LabeledSystem('cp2k/restart_aimd',fmt='cp2k/aimd_output', restart=True)
#        self.system_2 = dpdata.LabeledSystem('cp2k/restart_aimd/deepmd', fmt='deepmd/raw')
#        self.places = 6
#        self.e_places = 6
#        self.f_places = 6
#        self.v_places = 4
#
<<<<<<< HEAD
#class TestCp2kAimdOutputError(unittest.TestCase):
=======
# class TestCp2kAimdOutputError(unittest.TestCase):
>>>>>>> upstream/devel
#    def setUp(self):
#        pass
#
#    def restart_error(self):
#        with self.assertRaises(AssertionError):
#            dpdata.LabeledSystem('cp2k/restart_aimd', fmt='cp2k/aimd_output', restart=False)

if __name__ == "__main__":
    unittest.main()


# %%
# print(1)
# system_1 = dpda.La
# system_1 = dpdata.LabeledSystem('cp2k/restart_aimd',fmt='cp2k/aimd_output', restart=True)

# %%
