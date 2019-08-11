import numpy as np

def get_frames (fname) :
    coord_flag = False
    force_flag = False
    eV = 2.72113838565563E+01
    angstrom = 5.29177208590000E-01
    fp = open(fname)
    atom_symbol_list = []
    cell = []
    coord = []
    force = []

    for idx, ii in enumerate(fp) :
        if 'CELL| Vector' in ii :
            cell.append(ii.split()[4:7])
        if 'Atom  Kind  Element' in ii :
            coord_flag = True
            coord_idx = idx
        # get the coord block info
        if coord_flag :
            if (idx > coord_idx + 1) :
                if (ii == '\n') :
                    coord_flag = False
                else :
                    coord.append(ii.split()[4:7])
                    atom_symbol_list.append(ii.split()[2])
        if 'ENERGY|' in ii :
            energy = (ii.split()[8])
        if ' Atom   Kind ' in ii :
            force_flag = True
            force_idx = idx
        if force_flag :
            if (idx > force_idx) :
                if 'SUM OF ATOMIC FORCES' in ii :
                    force_flag = False
                else :
                    force.append(ii.split()[3:6])
    fp.close()
    assert(coord), "cannot find coords"
    assert(energy), "cannot find energies"
    assert(force), "cannot find forces"

    #conver to float array and add extra dimension for nframes
    cell = np.array(cell)
    cell = cell.astype(np.float)
    cell = cell[np.newaxis, :, :]
    coord = np.array(coord)
    coord = coord.astype(np.float)
    coord = coord[np.newaxis, :, :]
    atom_symbol_list = np.array(atom_symbol_list)
    force = np.array(force)
    force = force.astype(np.float)
    force = force[np.newaxis, :, :]
    force = force * eV / angstrom
    energy = float(energy) * eV
    energy = np.array(energy)
    energy = energy[np.newaxis]
    atom_names, atom_types, atom_numbs = np.unique(atom_symbol_list, return_inverse=True, return_counts=True)


    return list(atom_names), list(atom_numbs), atom_types, cell, coord, energy, force


