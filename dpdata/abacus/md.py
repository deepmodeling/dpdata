from ast import dump
import os,sys
import numpy as np
from .scf import ry2ev, bohr2ang, kbar2evperang3, get_block, get_geometry_in, get_cell, get_coords
import re

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
        if  len(line)>0 and "md_dumpfreq" in line and "md_dumpfreq" == line.split()[0]:
            return int(line.split()[1])
    return 1

def get_coords_from_dump(dumplines, natoms):
    nlines = len(dumplines)
    total_natoms = sum(natoms)
    calc_stress = False
    if "VIRIAL" in dumplines[6]:
        calc_stress = True
    else:
        assert("POSITIONS" in dumplines[6] and "FORCE" in dumplines[6]), "keywords 'POSITIONS' and 'FORCE' cannot be found in the 6th line. Please check."
    nframes_dump = -1
    if calc_stress:
        nframes_dump = int(nlines/(total_natoms + 13))
    else:
        nframes_dump = int(nlines/(total_natoms + 9))
    assert(nframes_dump > 0), "Number of lines in MD_dump file = %d. Number of atoms = %d. The MD_dump file is incomplete."%(nlines, total_natoms)
    cells = np.zeros([nframes_dump, 3, 3])
    stresses = np.zeros([nframes_dump, 3, 3])
    forces = np.zeros([nframes_dump, total_natoms, 3])
    coords = np.zeros([nframes_dump, total_natoms, 3])
    iframe = 0
    for iline in range(nlines):
        if "MDSTEP" in dumplines[iline]:
            # read in LATTICE_CONSTANT
            celldm = float(dumplines[iline+1].split(" ")[-1])
            # read in LATTICE_VECTORS
            for ix in range(3):
                cells[iframe, ix] = np.array([float(i) for i in re.split('\s+', dumplines[iline+3+ix])[-3:]]) * celldm
                if calc_stress:
                    stresses[iframe, ix] = np.array([float(i) for i in re.split('\s+', dumplines[iline+7+ix])[-3:]])
            for iat in range(total_natoms):
                if calc_stress:
                    coords[iframe, iat] = np.array([float(i) for i in re.split('\s+', dumplines[iline+11+iat])[-6:-3]])*celldm
                    forces[iframe, iat] = np.array([float(i) for i in re.split('\s+', dumplines[iline+11+iat])[-3:]])
                else:
                    coords[iframe, iat] = np.array([float(i) for i in re.split('\s+', dumplines[iline+7+iat])[-6:-3]])*celldm
                    forces[iframe, iat] = np.array([float(i) for i in re.split('\s+', dumplines[iline+7+iat])[-3:]])
            iframe += 1
    assert(iframe == nframes_dump), "iframe=%d, nframe_dump=%d. Number of frames does not match number of lines in MD_dump."%(iframe, nframes_dump)
    cells *= bohr2ang
    coords *= bohr2ang
    stresses *= kbar2evperang3
    return coords, cells, forces, stresses

def get_energy(outlines, ndump, dump_freq):
    energy = []
    nenergy = 0
    for line_idx, line in enumerate(outlines):
        if "final etot is" in line:
            if nenergy%dump_freq == 0:
                energy.append(float(line.split()[-2]))
            nenergy+=1
    assert(ndump == len(energy)), "Number of total energies in running_md.log = %d. Number of frames in MD_dump = %d. Please check."%(len(energy), ndump)
    energy = np.array(energy)
    return energy


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
    #ndump = int(os.popen("ls -l %s | grep 'md_pos_' | wc -l" %path_out).readlines()[0])
    # number of dumped geometry files
    #coords = get_coords_from_cif(ndump, dump_freq, atom_names, natoms, types, path_out, cell)
    with open(os.path.join(path_out, "MD_dump"), 'r') as fp:
        dumplines = fp.read().split('\n')
    coords, cells, force, stress = get_coords_from_dump(dumplines, natoms)
    ndump = np.shape(coords)[0]
    with open(os.path.join(path_out, "running_md.log"), 'r') as fp:
        outlines = fp.read().split('\n')
    energy = get_energy(outlines, ndump, dump_freq)
    for iframe in range(ndump):
        stress[iframe] *= np.linalg.det(cells[iframe, :, :].reshape([3, 3]))
    if np.sum(np.abs(stress[0])) < 1e-10:
        stress = None
    data = {}
    data['atom_names'] = atom_names
    data['atom_numbs'] = natoms
    data['atom_types'] = types
    data['cells'] = cells
    #for idx in range(ndump):
    #    data['cells'][:, :, :] = cell
    data['coords'] = coords
    data['energies'] = energy
    data['forces'] = force
    data['virials'] = stress
    if type(data['virials']) != np.ndarray:
        del data['virials']
    data['orig'] = np.zeros(3)

    return data
