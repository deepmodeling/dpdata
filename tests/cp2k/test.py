import dpdata
import numpy as np


# simple test for cp2k
# first case: with force and virial information
print("first case")
print("--------------------------------------------------------------")
cp2k_output = dpdata.LabeledSystem('cp2k_output', fmt = 'cp2k/output')
print("atom name")
print(cp2k_output['atom_names'])
print("atom number")
print(cp2k_output['atom_numbs'])
print("atom type")
print(cp2k_output['atom_types'])
print("atom cell")
print(cp2k_output['cells'])
print("atom coord")
print(cp2k_output['coords'])
print("energyies")
print(cp2k_output['energies'])
print("forces")
print(cp2k_output['forces'])
print("virials")
print(cp2k_output['virials'])
# no virial

cp2k_output.to_deepmd_raw('dpmd_raw')
cp2k_output.to_deepmd_npy('dpmd_npy')

# second case: with force and no! virial information
# double header information is contained, to test robustness of this parser
print("\n")
print("\n")
print("second case")
print("---------------------------------------------------------------")

cp2k_output = dpdata.LabeledSystem('cp2k_output_2', fmt = 'cp2k/output')
print("atom name")
print(cp2k_output['atom_names'])
print("atom number")
print(cp2k_output['atom_numbs'])
print("atom type")
print(cp2k_output['atom_types'])
#np.savetxt("ref_type", cp2k_output['atom_types'])
print("atom cell")
print(cp2k_output['cells'])
#np.savetxt("ref_cell", cp2k_output['cells'][0])
print("atom coord")
print(cp2k_output['coords'])
#np.savetxt("ref_coord", cp2k_output['coords'][0])
print("energyies")
print(cp2k_output['energies'])
print("forces")
print(cp2k_output['forces'])

