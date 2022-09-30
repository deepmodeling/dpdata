import os,sys
import numpy as np
from ..unit import EnergyConversion, PressureConversion, LengthConversion
import re
bohr2ang = LengthConversion("bohr", "angstrom").value()
ry2ev = EnergyConversion("rydberg", "eV").value()
kbar2evperang3 = PressureConversion("kbar", "eV/angstrom^3").value()

def get_block (lines, keyword, skip = 0, nlines = None):
    ret = []
    found = False
    if not nlines:
        nlines = 1e6
    for idx,ii in enumerate(lines) :
        if keyword in ii :
            found = True
            blk_idx = idx + 1 + skip
            line_idx = 0
            while len(re.split("\s+", lines[blk_idx])) == 0:
                blk_idx += 1
            while line_idx < nlines and blk_idx != len(lines):
                if len(re.split("\s+", lines[blk_idx])) == 0 or lines[blk_idx] == "":
                    blk_idx+=1
                    continue
                ret.append(lines[blk_idx])
                blk_idx += 1
                line_idx += 1
            break
    if not found:
        return None
    return ret

def get_geometry_in(fname, inlines):
    geometry_path_in = os.path.join(fname, "STRU")
    for line in inlines:
        if "stru_file" in line and "stru_file"==line.split()[0]:
           atom_file = line.split()[1]
           geometry_path_in = os.path.join(fname, atom_file)
           break
    return geometry_path_in

def get_path_out(fname, inlines):
    path_out = os.path.join(fname, "OUT.ABACUS/running_scf.log")
    for line in inlines:
        if "suffix" in line and "suffix"==line.split()[0]:
           suffix = line.split()[1]
           path_out = os.path.join(fname, "OUT.%s/running_scf.log" % suffix)
           break
    return path_out

def get_cell(geometry_inlines):
    cell_lines = get_block(geometry_inlines, "LATTICE_VECTORS", skip = 0, nlines = 3)
    celldm_lines = get_block(geometry_inlines, "LATTICE_CONSTANT", skip=0, nlines=1)

    celldm = float(celldm_lines[0].split()[0]) * bohr2ang # lattice const is in Bohr
    cell = []
    for ii in range(3):
        cell.append([float(jj) for jj in cell_lines[ii].split()[0:3]])
    cell = celldm*np.array(cell)
    return celldm, cell

def get_coords(celldm, cell, geometry_inlines, inlines):
    coords_lines = get_block(geometry_inlines, "ATOMIC_POSITIONS", skip=0)
    # assuming that ATOMIC_POSITIONS is at the bottom of the STRU file
    coord_type = coords_lines[0].split()[0].lower() # cartisan or direct
    atom_names = [] # element abbr in periodic table 
    atom_types = [] # index of atom_names of each atom in the geometry
    atom_numbs = [] # of atoms for each element
    coords = [] # coordinations of atoms
    ntype = 0
    for line in inlines:
        if "ntype" in line and "ntype"==line.split()[0]:
            ntype = int(line.split()[1])
            break
    if ntype <= 0:
        raise RuntimeError('ntype cannot be found in INPUT file.')
    line_idx = 1 # starting line of first element
    for it in range(ntype):
        atom_names.append(coords_lines[line_idx].split()[0])
        line_idx+=2
        atom_numbs.append(int(coords_lines[line_idx].split()[0]))
        line_idx+=1
        for iline in range(atom_numbs[it]):
            xyz = np.array([float(xx) for xx in coords_lines[line_idx].split()[0:3]])
            if coord_type == "cartesian":
                xyz = xyz*celldm
            elif coord_type == "direct":
                tmp = np.matmul(xyz, cell)
                xyz = tmp
            else:
                print("coord_type = %s" % coord_type)
                raise RuntimeError("Input coordination type is invalid.\n Only direct and cartesian are accepted.")
            coords.append(xyz)
            atom_types.append(it)
            line_idx += 1
    coords = np.array(coords) # need transformation!!!
    atom_types = np.array(atom_types)
    return atom_names, atom_numbs, atom_types, coords

def get_energy(outlines):
    Etot = None
    for line in outlines:
        if "!FINAL_ETOT_IS" in line:
            Etot = float(line.split()[1]) # in eV
            break
    if not Etot:
       raise RuntimeError("Final total energy cannot be found in output. Unknown problem.")
    for line in outlines:
        if "convergence has NOT been achieved!" in line:
            return Etot,False
    return Etot,True

def get_force (outlines, natoms):
    force = []
    force_inlines = get_block (outlines, "TOTAL-FORCE (eV/Angstrom)", skip = 4, nlines=np.sum(natoms))
    if force_inlines is None:
        raise RuntimeError("TOTAL-FORCE (eV/Angstrom) is not found in running_scf.log. Please check.")
    for line in force_inlines:
        force.append([float(f) for f in line.split()[1:4]])
    force = np.array(force)
    return force

def get_stress(outlines):
    stress = []
    stress_inlines = get_block(outlines, "TOTAL-STRESS (KBAR)", skip = 3, nlines=3)
    if stress_inlines is None:
        return None
    for line in stress_inlines:
        stress.append([float(f) for f in line.split()])
    stress = np.array(stress) * kbar2evperang3
    return stress



def get_frame (fname):
    if type(fname) == str:
        # if the input parameter is only one string, it is assumed that it is the 
        # base directory containing INPUT file;
        path_in = os.path.join(fname, "INPUT")
    else:
        raise RuntimeError('invalid input')    
    with open(path_in, 'r') as fp:
        inlines = fp.read().split('\n')
    
    geometry_path_in = get_geometry_in(fname, inlines) 
    path_out = get_path_out(fname, inlines) 
    with open(geometry_path_in, 'r') as fp:
        geometry_inlines = fp.read().split('\n')
    with open(path_out, 'r') as fp:
        outlines = fp.read().split('\n')

    celldm, cell = get_cell(geometry_inlines) 
    atom_names, natoms, types, coords = get_coords(celldm, cell, geometry_inlines, inlines) 
    
    energy,converge = get_energy(outlines) 
    if not converge:
        return {'atom_names':atom_names,\
                'atom_numbs':natoms,\
                'atom_types':types,\
                'cells':[],\
                'coords':[],\
                'energies':[],\
                'forces':[]}
    force = get_force (outlines, natoms) 
    stress = get_stress(outlines) 
    if stress is not None:
        stress *= np.abs(np.linalg.det(cell)) 
    
    data = {}
    data['atom_names'] = atom_names
    data['atom_numbs'] = natoms
    data['atom_types'] = types
    data['cells'] = cell[np.newaxis, :, :]
    data['coords'] = coords[np.newaxis, :, :]
    data['energies'] = np.array(energy)[np.newaxis]
    data['forces'] = force[np.newaxis, :, :]
    if stress is not None:
        data['virials'] = stress[np.newaxis, :, :]
    data['orig'] = np.zeros(3)
    # print("atom_names = ", data['atom_names'])
    # print("natoms = ", data['atom_numbs'])
    # print("types = ", data['atom_types'])
    # print("cells = ", data['cells'])
    # print("coords = ", data['coords'])
    # print("energy = ", data['energies'])
    # print("force = ", data['forces'])
    # print("virial = ", data['virials'])
    return data

def get_nele_from_stru(geometry_inlines):
    key_words_list = ["ATOMIC_SPECIES", "NUMERICAL_ORBITAL", "LATTICE_CONSTANT", "LATTICE_VECTORS", "ATOMIC_POSITIONS", "NUMERICAL_DESCRIPTOR"]
    keyword_sequence = []
    keyword_line_index = []
    atom_names = []
    atom_numbs = []
    for iline, line in enumerate(geometry_inlines):
        if line.split() == []:
            continue
        have_key_word = False
        for keyword in key_words_list:
            if keyword in line and keyword == line.split()[0]:
                keyword_sequence.append(keyword)
                keyword_line_index.append(iline)
    assert(len(keyword_line_index) == len(keyword_sequence))
    assert(len(keyword_sequence) > 0)
    keyword_line_index.append(len(geometry_inlines))

    nele = 0
    for idx, keyword in enumerate(keyword_sequence):
        if keyword == "ATOMIC_SPECIES":
            for iline in range(keyword_line_index[idx]+1, keyword_line_index[idx+1]):
                if len(re.split("\s+", geometry_inlines[iline])) >= 3:
                    nele += 1
    return nele

def get_frame_from_stru(fname):
    assert(type(fname) == str)
    with open(fname, 'r') as fp:
        geometry_inlines = fp.read().split('\n')
    nele = get_nele_from_stru(geometry_inlines)
    inlines = ["ntype %d" %nele]
    celldm, cell = get_cell(geometry_inlines)
    atom_names, natoms, types, coords = get_coords(celldm, cell, geometry_inlines, inlines) 
    data = {}
    data['atom_names'] = atom_names
    data['atom_numbs'] = natoms
    data['atom_types'] = types
    data['cells'] = cell[np.newaxis, :, :]
    data['coords'] = coords[np.newaxis, :, :]
    data['orig'] = np.zeros(3)

    return data

def make_unlabeled_stru(data, frame_idx, pp_file=None, numerical_orbital=None, numerical_descriptor=None, mass=None):
    out = "ATOMIC_SPECIES\n"
    for iele in range(len(data['atom_names'])):
        out += data['atom_names'][iele] + " "
        if mass is not None:
            out += "%.3f "%mass[iele]
        else:
            out += "1 "
        if pp_file is not None:
            out += "%s\n"%pp_file[iele]
        else:
            out += "\n"
    out += "\n"

    if numerical_orbital is not None:
        assert(len(numerical_orbital) == len(data['atom_names']))
        out += "NUMERICAL_ORBITAL\n"
        for iele in range(len(numerical_orbital)):
            out += "%s\n"%numerical_orbital[iele]
        out += "\n"

    if numerical_descriptor is not None:
        assert(type(numerical_descriptor) == str)
        out += "NUMERICAL_DESCRIPTOR\n%s\n"%numerical_descriptor
        out += "\n"
    
    out += "LATTICE_CONSTANT\n"
    out += str(1/bohr2ang) + "\n\n"

    out += "LATTICE_VECTORS\n"
    for ix in range(3):
        for iy in range(3):
            out += str(data['cells'][frame_idx][ix][iy]) + " "
        out += "\n"
    out += "\n"

    out += "ATOMIC_POSITIONS\n"
    out += "Cartesian    # Cartesian(Unit is LATTICE_CONSTANT)\n"
    #ret += "\n"
    natom_tot = 0
    for iele in range(len(data['atom_names'])):
        out += data['atom_names'][iele] + "\n"
        out += "0.0\n"
        out += str(data['atom_numbs'][iele]) + "\n"
        for iatom in range(data['atom_numbs'][iele]):
            out += "%.12f %.12f %.12f %d %d %d\n" % (data['coords'][frame_idx][natom_tot, 0], data['coords'][frame_idx][natom_tot, 1], data['coords'][frame_idx][natom_tot, 2], 1, 1, 1)
            natom_tot += 1
    assert(natom_tot == sum(data['atom_numbs']))
    return out

#if __name__ == "__main__":
#    path = "/home/lrx/work/12_ABACUS_dpgen_interface/dpdata/dpdata/tests/abacus.scf"
#    data = get_frame(path)
