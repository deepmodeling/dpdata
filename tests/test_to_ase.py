import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys, IsPBC
try:
    from ase import Atoms
    from ase.io import write
except ModuleNotFoundError:
    exist_module=False
else:
    exist_module=True


@unittest.skipIf(not exist_module,"skip test_ase")
class TestASE(unittest.TestCase, CompSys, IsPBC):
    
    def setUp(self): 
        system_1 = dpdata.System()
        system_1.from_lammps_lmp(os.path.join('poscars', 'conf.lmp'), type_map = ['O', 'H'])
        write('tmp.POSCAR',system_1.to_ase_structure()[0],vasp5=True)
        self.system_1=system_1
        self.system_2=dpdata.System('tmp.POSCAR')
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


@unittest.skipIf(not exist_module, "skip test_ase")
class TestFromASE(unittest.TestCase, CompSys, IsPBC):
    """Test ASEStructureFormat.from_system"""
    def setUp(self): 
        system_1 = dpdata.System()
        system_1.from_lammps_lmp(os.path.join('poscars', 'conf.lmp'), type_map = ['O', 'H'])
        atoms = system_1.to_ase_structure()[0]
        self.system_1 = system_1
        self.system_2 = dpdata.System(atoms, fmt="ase/structure")
        # assign the same type_map
        self.system_2.sort_atom_names(type_map=self.system_1.get_atom_names())
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


if __name__ == '__main__':
    unittest.main()
    
