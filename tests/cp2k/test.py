
import dpdata

cp2k_output = dpdata.LabeledSystem('cp2k_output', fmt = 'cp2k/output')
print(cp2k_output['atom_names'])
print(cp2k_output['atom_numbs'])
print(cp2k_output['atom_types'])
print(cp2k_output['cells'])
print(cp2k_output['coords'])
print(cp2k_output['energies'])
print(cp2k_output['forces'])
# no virial
cp2k_output.to_deepmd_raw('dpmd_raw')
cp2k_output.to_deepmd_npy('dpmd_npy')
