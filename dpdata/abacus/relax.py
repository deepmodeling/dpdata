import os,sys
import numpy as np
from .scf import bohr2ang, kbar2evperang3, get_geometry_in, get_cell, get_coords

# Read in geometries from an ABACUS RELAX(CELL-RELAX) trajectory in OUT.XXXX/runnning_relax/cell-relax.log. 

def get_log_file(fname, inlines):
    suffix = "ABACUS"
    calculation = "scf"
    for line in inlines:
        if "suffix" in line and "suffix"==line.split()[0]:
           suffix = line.split()[1]
        elif "calculation" in line and "calculation" == line.split()[0]:
            calculation = line.split()[1]
    logf = os.path.join(fname, "OUT.%s/running_%s.log"%(suffix,calculation))
    return logf

def get_coords_from_log(loglines,natoms):
    '''
    NOTICE: unit of coords and cells is Angstrom
    '''
    natoms_log = 0
    for line in loglines:
        if line[13:41] == "number of atom for this type":
            natoms_log += int(line.split()[-1])

    assert(natoms_log>0 and natoms_log == natoms),"ERROR: detected atom number in log file is %d" % natoms

    energy = []
    cells = []
    coords = []
    force = []
    stress = []

    for i in range(len(loglines)):
        line = loglines[i]
        if line[18:41] == "lattice constant (Bohr)": 
            a0 = float(line.split()[-1])
        elif len(loglines[i].split()) >=2 and loglines[i].split()[1] == 'COORDINATES':
            coords.append([])
            direct_coord = False
            if loglines[i].split()[0] == 'DIRECT':
                direct_coord = True
                for k in range(2,2+natoms):
                    coords[-1].append(list(map(lambda x: float(x),loglines[i+k].split()[1:4])))
            elif loglines[i].split()[0] == 'CARTESIAN':
                for k in range(2,2+natoms):
                    coords[-1].append(list(map(lambda x: float(x)*a0,loglines[i+k].split()[1:4])))  
            else:
                assert(False),"Unrecongnized coordinate type, %s, line:%d" % (loglines[i].split()[0],i)

            converg = True
            for j in range(i):
                if loglines[i-j-1][1:36] == 'Ion relaxation is not converged yet':
                    converg = False
                    break
                elif loglines[i-j-1][1:29] == 'Ion relaxation is converged!':
                    converg = True
                    break

            if converg:
                for j in range(i+1,len(loglines)):
                    if loglines[j][1:56] == "Lattice vectors: (Cartesian coordinate: in unit of a_0)":
                        cells.append([])
                        for k in range(1,4):
                            cells[-1].append(list(map(lambda x:float(x)*a0,loglines[j+k].split()[0:3]))) 
                        break
            else:
                cells.append(cells[-1])

            if direct_coord:
                coords[-1] = coords[-1].dot(cells[-1])

        elif line[4:15] == "TOTAL-FORCE":
            force.append([])
            for j in range(5,5+natoms):
                force[-1].append(list(map(lambda x:float(x),loglines[i+j].split()[1:4])))
        elif line[1:13] == "TOTAL-STRESS":
            stress.append([])
            for j in range(4,7):
                stress[-1].append(list(map(lambda x:float(x),loglines[i+j].split()[0:3])))
        elif line[1:14] == "final etot is":
            energy.append(float(line.split()[-2]))

    assert(len(cells) == len(coords) or len(cells)+1 == len(coords)),"ERROR: detected %d coordinates and %d cells" % (len(coords),len(cells))
    if len(cells)+1 == len(coords): del(coords[-1])

    energy = np.array(energy)
    cells = np.array(cells)
    coords = np.array(coords)
    stress = np.array(stress)
    force = np.array(force)

    cells *= bohr2ang
    coords *= bohr2ang

    virial = np.zeros([len(cells), 3, 3])
    for i in range(len(cells)):
        volume = np.linalg.det(cells[i, :, :].reshape([3, 3]))
        virial[i] = stress[i] * kbar2evperang3 * volume

    return energy,cells,coords,force,stress,virial

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
    with open(geometry_path_in, 'r') as fp:
        geometry_inlines = fp.read().split('\n')
    celldm, cell = get_cell(geometry_inlines) 
    atom_names, natoms, types, coord_tmp = get_coords(celldm, cell, geometry_inlines, inlines) 
    
    logf = get_log_file(fname, inlines) 
    assert(os.path.isfile(logf)),"Error: can not find %s" % logf
    with open(logf) as f1: lines = f1.readlines()

    atomnumber = 0
    for i in natoms: atomnumber += i
    energy,cells,coords,force,stress,virial = get_coords_from_log(lines,atomnumber)

    data = {}
    data['atom_names'] = atom_names
    data['atom_numbs'] = natoms
    data['atom_types'] = types
    data['cells'] = cells
    data['coords'] = coords
    data['energies'] = energy
    data['forces'] = force
    data['virials'] = virial
    data['stress'] = stress
    data['orig'] = np.zeros(3)

    return data
