import unittest
from context import dpdata

data={"name": "Hydrogen",
   "atomic_no": 1,
   "X": 2.2,
   "atomic_mass": 1.00794,
   "radius": 0.25,
   "calculated_radius": 0.53
    }

class TestPeriodicTable(unittest.TestCase):
    def setUp (self) :
        self.H = dpdata.periodic_table.Element("H")

    def test_H(self):
        H=self.H
        self.assertEqual(H.name,data['name'])
        self.assertEqual(H.Z,data['atomic_no'])
        self.assertEqual(H.X,data['X'])
        self.assertEqual(H.mass,data['atomic_mass'])
        self.assertEqual(H.radius,data['radius'])
        self.assertEqual(H.calculated_radius,data['calculated_radius'])
        self.assertEqual(H.X,dpdata.periodic_table.Element.from_Z(1).X)

if __name__ == '__main__':
    unittest.main()
