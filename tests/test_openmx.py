import unittest

import numpy as np
from context import dpdata

bohr2ang = dpdata.unit.LengthConversion("bohr", "angstrom").value()

### testing result ###
# ['C', 'H']
# [1, 4]
# [0 1 1 1 1]

class TestOPENMXProps:
    def test_atom_names(self):
        self.assertEqual(self.system.data["atom_names"], ["C", "H"])

    def test_atom_numbs(self):
        self.assertEqual(self.system.data["atom_numbs"], [1, 4])

    def test_atom_types(self):
        for ii in range(0, 1):
            self.assertEqual(self.system.data["atom_types"][ii], 0)
        for ii in range(1, 5):
            self.assertEqual(self.system.data["atom_types"][ii], 1)

    def test_cell(self):
        ref = 10.0 * np.eye(3)
        self.assertEqual(self.system.get_nframes(), 200)
        for ff in range(self.system.get_nframes()):
            for ii in range(3):
                for jj in range(3):
                    self.assertEqual(self.system["cells"][ff][ii][jj], ref[ii][jj])

    def test_coord(self):
        with open("openmx/Methane.md") as md_file:
            lines = md_file.readlines()
        lines = lines[-5:]
        atom_names=["C","H"]
        coords=[]
        for index, line in enumerate(lines):
            for atom_name in atom_names:
                atom_name += " "
                if atom_name in line:
                    parts = line.split()
                    for_line = [float(parts[1]), float(parts[2]), float(parts[3])]
                    coords.append(for_line)
        coords=np.array(coords)
        celll=10.0
        ## Applying PBC ##
        for ii in range(5):
            for jj in range(3):
                if coords[ii][jj] < 0:
                    coords[ii][jj] += celll
                elif coords[ii][jj] >= celll:
                    coords[ii][jj] -= celll                
                self.assertAlmostEqual(
                    self.system["coords"][-1][ii][jj], coords[ii][jj]
                )

class TestOPENMXTraj(unittest.TestCase, TestOPENMXProps):
    def setUp(self):
        self.system = dpdata.System("openmx/Methane", fmt="openmx")

class TestOPENMXLabeledTraj(unittest.TestCase, TestOPENMXProps):
    def setUp(self):
        self.system = dpdata.LabeledSystem("openmx/Methane", fmt="openmx")

if __name__ == "__main__":
    # pass
    unittest.main()
