from __future__ import annotations

import unittest

from comp_sys import CompLabeledSys, CompSys, IsPBC
from context import dpdata

try:
    import ase  # noqa: F401
except ModuleNotFoundError:
    skip_ase = True
else:
    skip_ase = False


@unittest.skipIf(skip_ase, "skip ase related test. install ase to fix")
@unittest.skip("This test is broken")
class TestASEtraj1(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.multi_systems = dpdata.MultiSystems.from_file(
            "ase_traj/HeAlO.traj", fmt="ase/structure"
        )
        self.system_1 = self.multi_systems.systems["Al0He4O0"]
        self.system_2 = dpdata.LabeledSystem("ase_traj/Al0He4O0", fmt="deepmd")
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


@unittest.skipIf(skip_ase, "skip ase related test. install ase to fix")
class TestASEtraj2(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_temp0 = dpdata.MultiSystems.from_file(
            file_name="ase_traj/HeAlO.traj", fmt="ase/structure"
        )
        self.system_1 = self.system_temp0.systems["Al2He1O3"]  # .sort_atom_types()
        self.system_temp1 = dpdata.LabeledSystem("ase_traj/Al2He1O3", fmt="deepmd")
        self.system_temp2 = dpdata.LabeledSystem("ase_traj/Al4He4O6", fmt="deepmd")
        self.system_temp3 = dpdata.MultiSystems(self.system_temp2, self.system_temp1)
        self.system_2 = self.system_temp3.systems["Al2He1O3"]
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


@unittest.skipIf(skip_ase, "skip ase related test. install ase to fix")
class TestASEtraj3(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.System("ase_traj/MoS2", fmt="deepmd")
        self.system_2 = dpdata.System("ase_traj/MoS2.traj", fmt="ase/traj")
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


@unittest.skipIf(skip_ase, "skip ase related test. install ase to fix")
class TestASEtraj3Labeled(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem("ase_traj/MoS2", fmt="deepmd")
        self.system_2 = dpdata.LabeledSystem("ase_traj/MoS2.traj", fmt="ase/traj")
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


@unittest.skipIf(skip_ase, "skip ase related test. install ase to fix")
class TestASEtraj4(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.System("ase_traj/MoS2", fmt="deepmd")
        self.system_1.to_system(file_name="tmp.traj")
        self.system_2 = dpdata.System("tmp.traj", fmt="ase/traj")
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


@unittest.skipIf(skip_ase, "skip ase related test. install ase to fix")
class TestASEtraj4Labeled(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem("ase_traj/MoS2", fmt="deepmd")
        self.system_1.to_labeled_system(file_name="tmp.traj")
        self.system_2 = dpdata.LabeledSystem("tmp.traj", fmt="ase/traj")
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


if __name__ == "__main__":
    unittest.main()
