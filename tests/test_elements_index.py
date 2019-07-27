import os
import numpy as np
import unittest
from dpdata.system import elements_index_map

class ElementIndexMap(unittest.TestCase):
  def test_func1(self):
      element=["C","N","H"] 
      ref={'C': 0, 'N': 1, 'H': 2}
      self.assertEqual(ref,elements_index_map(element))

  def test_func2(self):
      element=["C","N","H"] 
      ref={'H': 0, 'C': 1, 'N': 2}
      self.assertEqual(ref,elements_index_map(element,standard=True))

  def test_func3(self):
      element=["C","N","H"] 
      ref={0: 'H', 1: 'C', 2: 'N'}
      self.assertEqual(ref,elements_index_map(element,standard=True,inverse=True))

  def test_func4(self):
      element=["C","N","H"] 
      ref={0: 'C', 1: 'N', 2: 'H'}
      self.assertEqual(ref,elements_index_map(element,inverse=True))
if __name__ == '__main__':
    unittest.main()
