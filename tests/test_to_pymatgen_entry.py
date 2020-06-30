import os
import numpy as np
import unittest
from context import dpdata
from comp_sys import CompSys, IsPBC
from monty.serialization import  loadfn
try:
   from pymatgen.entries.computed_entries import ComputedStructureEntry
   exist_module=True
except:
   exist_module=False

@unittest.skipIf(not exist_module,"skip pymatgen")
class TestPymatgen(unittest.TestCase):
    
    def test(self): 
         ls1= dpdata.LabeledSystem(os.path.join('poscars', 'OUTCAR.ch4.1step'),fmt='OUTCAR')
         entry1=ls1.to_pymatgen_ComputedStructureEntry()
         self.assertEqual(entry1,[])
         ls2= dpdata.LabeledSystem(os.path.join('poscars', 'OUTCAR.h2o.md.10'),fmt='OUTCAR')
         entry2=ls2.to_pymatgen_ComputedStructureEntry()
         self.assertEqual(len(entry2),10)
         last_entry=loadfn("computed_structure_entry.json")
         self.assertEqual(last_entry.as_dict(),entry2[-1].as_dict())
         

if __name__ == '__main__':
    unittest.main()
    
