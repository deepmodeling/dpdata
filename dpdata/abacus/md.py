import os,sys
import numpy as np
from .scf import ry2ev, kbar2evperang3, get_block, get_geometry_in, get_cell, get_coords

# Read in geometries from an ABACUS MD trajectory.
# The atomic coordinates are read in from generated files in OUT.XXXX.
# Energies, forces
# IMPORTANT: the program defaultly takes STRU input file as standard cell information,
# therefore the direct and cartesan coordinates read could be different from the ones in 
# the output cif files!!!
# It is highly recommanded to use ORTHOGANAL coordinates in STRU file if you wish to get
# same coordinates in both dpdata and output cif files. 

def get_path_out(fname, inlines):
    # This function is different from the same-name function in scf.py.
    # This function returns OUT.XXXX's base directory.
    path_out = os.path.join(fname, "OUT.ABACUS/")
    for line in inlines:
        if  len(line)>0 and "suffix" in line and "suffix"==line.split()[0]:
           suffix = line.split()[1]
           path_out = os.path.join(fname, "OUT.%s/" % suffix)
           break
    return path_out

def get_coord_dump_freq(inlines):
    for line in inlines:
        if  len(line)>0 and "md_dumpmdfred" in line and "md_dumpmdfred" == line.split()[0]:
            return int(line.split()[1])
    return 1

# set up a cell according to cell info in cif file.
# maybe useful later
'''
def setup_cell(a, b, c, alpha, beta, gamma):
    cell = np.zeros(3, 3)
    cell[0, 0] = a
    cell[1, 0] = b*np.cos(gamma/180*np.pi)
    cell[1, 1] = b*np.sin(gamma/180*np.pi)
    cell[2, 0] = c*np.cos(beta/180*np.pi)
    cell[2, 1] = c*(b*np.cos(alpha/180*np.pi) - cell[1, 0]*np.cos(beta/180*np.pi))/cell[1, 1]
    cell[2, 2] = np.sqrt(c**2 - cell[2, 0]**2 - cell[2, 1]**2)
    return cell
'''

def get_single_coord_from_cif(pos_file, atom_names, natoms, cell):
    assert(len(atom_names) == len(natoms))
    nele = len(atom_names)
    total_natoms = sum(natoms)
    coord = np.zeros([total_natoms, 3])
    a = 0
    b = 0
    c = 0
    alpha = 0
    beta = 0
    gamma = 0
    with open(pos_file, "r") as fp:
        lines = fp.read().split("\n")
    for line in lines:
        if "_cell_length_a" in line:
            a = float(line.split()[1])
        if "_cell_length_b" in line:
            b = float(line.split()[1])
        if "_cell_length_c" in line:
            c = float(line.split()[1])  
        if "_cell_angle_alpha" in line:
            alpha = float(line.split()[1])
        if "_cell_angle_beta" in line:
            beta = float(line.split()[1])
        if "_cell_angle_gamma" in line:
            gamma = float(line.split()[1])
    assert(a > 0 and b > 0 and c > 0 and alpha > 0 and beta > 0 and gamma > 0)
    #cell = setup_cell(a, b, c, alpha, beta, gamma)
    coord_lines = get_block(lines=lines, keyword="_atom_site_fract_z", skip=0, nlines = total_natoms)
    
    ia_idx = 0
    for it in range(nele):
        for ia in range(natoms[it]):
            coord_line = coord_lines[ia_idx].split()
            assert(coord_line[0] == atom_names[it])
            coord[ia_idx, 0] = float(coord_line[1])
            coord[ia_idx, 1] = float(coord_line[2])
            coord[ia_idx, 2] = float(coord_line[3])
            ia_idx+=1
    coord = np.matmul(coord, cell)
    # important! Coordinates are converted to Cartesian coordinate.
    return coord
    
    
def get_coords_from_cif(ndump, dump_freq, atom_names, natoms, types, path_out, cell):
    total_natoms = sum(natoms)
    #cell = np.zeros(ndump, 3, 3)
    coords = np.zeros([ndump, total_natoms, 3])
    pos_file = os.path.join(path_out, "STRU_READIN_ADJUST.cif")
    # frame 0 file is different from any other frames
    coords[0] = get_single_coord_from_cif(pos_file, atom_names, natoms, cell)
    for dump_idx in range(1, ndump):
        pos_file = os.path.join(path_out, "md_pos_%d.cif" %(dump_idx*dump_freq))
        #print("dump_idx = %s" %dump_idx)
        coords[dump_idx] = get_single_coord_from_cif(pos_file, atom_names, natoms, cell)
    return coords

def get_energy_force_stress(outlines, inlines, dump_freq, ndump, natoms, atom_names):
    stress = None
    total_natoms = sum(natoms)
    for line in inlines:
        if len(line)>0 and "stress" in line and "stress" == line.split()[0] and "1" == line.split()[1]:
            stress = np.zeros([ndump, 3, 3])
            break
    if type(stress) != np.ndarray:
        print("The ABACUS program has no stress output. Stress will not be read.")
    nenergy = 0
    nforce = 0
    nstress = 0
    energy = np.zeros(ndump)
    force = np.zeros([ndump, total_natoms, 3])

    for line_idx, line in enumerate(outlines):
        if "final etot is" in line:
            if nenergy%dump_freq == 0:
                energy[int(nenergy/dump_freq)] = float(line.split()[-2])
            nenergy+=1
        if "TOTAL-FORCE (eV/Angstrom)" in line:
            for iatom in range(0, total_natoms):
                force_line = outlines[line_idx+5+iatom]
                atom_force = [float(i) for i in force_line.split()[1:]]
                assert(len(atom_force) == 3)
                atom_force = np.array(atom_force)
                if nforce%dump_freq == 0:
                    force[int(nforce/dump_freq), iatom] = atom_force
            nforce+=1
            assert(nforce==nenergy)
        if "TOTAL-STRESS (KBAR)" in line:
            for idx in range(0, 3):
                stress_line = outlines[line_idx+4+idx]
                single_stress = [float(i) for i in stress_line.split()]
                if len(single_stress) != 3:
                    print(single_stress)
                assert(len(single_stress) == 3)
                single_stress = np.array(single_stress)
                if nstress%dump_freq == 0:
                    stress[int(nstress/dump_freq), idx] = single_stress
            nstress+=1
            assert(nstress==nforce)
    if type(stress) == np.ndarray:
        stress *= kbar2evperang3
    return energy, force, stress


def get_frame (fname):
    if type(fname) == str:
        # if the input parameter is only one string, it is assumed that it is the 
        # base directory containing INPUT file;
        path_in = os.path.join(fname, "INPUT")
    else:
        raise RuntimeError('invalid input')    
    with open(path_in, 'r') as fp:
        inlines = fp.read().split('\n')
    geometry_path_in = get_geometry_in(fname, inlines) # base dir of STRU
    path_out = get_path_out(fname, inlines) 

    with open(geometry_path_in, 'r') as fp:
        geometry_inlines = fp.read().split('\n')
    celldm, cell = get_cell(geometry_inlines) 
    atom_names, natoms, types, coords = get_coords(celldm, cell, geometry_inlines, inlines) 
    # This coords is not to be used.
    dump_freq = get_coord_dump_freq(inlines = inlines)
    ndump = int(os.popen("ls --l %s | grep 'md_pos_' | wc -l" %path_out).readlines()[0])
    # number of dumped geometry files
    coords = get_coords_from_cif(ndump, dump_freq, atom_names, natoms, types, path_out, cell)
    
    # TODO: Read in energies, forces and pressures.
    with open(os.path.join(path_out, "running_md.log"), 'r') as fp:
        outlines = fp.read().split('\n')
    energy, force, stress = get_energy_force_stress(outlines, inlines, dump_freq, ndump, natoms, atom_names)
    if type(stress) == np.ndarray:
        stress *= np.linalg.det(cell)
    data = {}
    data['atom_names'] = atom_names
    data['atom_numbs'] = natoms
    data['atom_types'] = types
    data['cells'] = np.zeros([ndump, 3, 3])
    for idx in range(ndump):
        data['cells'][:, :, :] = cell
    data['coords'] = coords
    data['energies'] = energy
    data['forces'] = force
    data['virials'] = stress
    data['orig'] = np.zeros(3)

    return data
