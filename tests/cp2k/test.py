
import dpdata
import numpy as np

cp2k_output = dpdata.LabeledSystem('output_40847', fmt = 'cp2k/output')
print("atom name")
print(cp2k_output['atom_names'])
print("atom number")
print(cp2k_output['atom_numbs'])
print("atom type")
print(cp2k_output['atom_types'])
np.savetxt("ref_type", cp2k_output['atom_types'])
print("atom cell")
print(cp2k_output['cells'])
np.savetxt("ref_cell", cp2k_output['cells'][0])
print("atom coord")
print(cp2k_output['coords'])
np.savetxt("ref_coord", cp2k_output['coords'][0])
print("energyies")
print(cp2k_output['energies'])
print("forces")
print(cp2k_output['forces'])
np.savetxt("ref_force", cp2k_output['forces'][0])
print("virials")
print(cp2k_output['virials'])
np.savetxt("ref_virial", cp2k_output['virials'][0])
# no virial
#cp2k_output.to_deepmd_raw('dpmd_raw')
#cp2k_output.to_deepmd_npy('dpmd_npy')
