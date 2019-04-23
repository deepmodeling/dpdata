import numpy as np

def system_info (lines, type_idx_zero = False) :
    atom_names = []
    atom_numbs = None
    for ii in lines: 
        if 'TITEL  =' in ii : 
            # get atom names from POTCAR info, tested only for PAW_PBE ...
            atom_names.append(ii.split()[3])
        elif 'ions per type' in ii :
            atom_numbs_ = [int(s) for s in ii.split()[4:]]
            if atom_numbs is None :                
                atom_numbs = atom_numbs_
            else :
                assert (atom_numbs == atom_numbs_), "in consistent numb atoms in OUTCAR"
            break
    assert(atom_numbs is not None), "cannot find ion type info in OUTCAR"
    atom_names = atom_names[:len(atom_numbs)]
    atom_types = []
    for idx,ii in enumerate(atom_numbs) :
        for jj in range(ii) :
            if type_idx_zero :
                atom_types.append(idx)
            else :
                atom_types.append(idx+1)
    return atom_names, atom_numbs, np.array(atom_types, dtype = int)


def get_outcar_block(fp) :
    blk = []
    for ii in fp :
        if not ii :
            return blk
        blk.append(ii.rstrip('\n'))
        if 'free  energy   TOTEN' in ii:
            return blk
    return blk

# we assume that the force is printed ...
def get_frames (fname, begin = 0, step = 1) :
    fp = open(fname)
    blk = get_outcar_block(fp)

    atom_names, atom_numbs, atom_types = system_info(blk, type_idx_zero = True)
    ntot = sum(atom_numbs)

    all_coords = []
    all_cells = []
    all_energies = []
    all_forces = []
    all_virials = []    

    cc = 0
    while len(blk) > 0 :
        if cc >= begin and (cc - begin) % step == 0 :
            coord, cell, energy, force, virial = analyze_block(blk, ntot)
            if len(coord) == 0:
                break
            all_coords.append(coord)
            all_cells.append(cell)
            all_energies.append(energy)
            all_forces.append(force)
            if virial is not None :
                all_virials.append(virial)
        blk = get_outcar_block(fp)
        cc += 1

    fp.close()
    return atom_names, atom_numbs, atom_types, np.array(all_cells), np.array(all_coords), np.array(all_energies), np.array(all_forces), np.array(all_virials)


def analyze_block(lines, ntot) :
    coord = []
    cell = []
    energy = None
    force = []
    virial = None
    for idx,ii in enumerate(lines) :
        if 'Iteration' in ii:
            pass
        elif 'free  energy   TOTEN' in ii:
            energy = float(ii.split()[4])
            assert((force is not None) and len(coord) > 0 and len(cell) > 0)
            # all_coords.append(coord)
            # all_cells.append(cell)
            # all_energies.append(energy)
            # all_forces.append(force)
            # if virial is not None :
            #     all_virials.append(virial)
            return coord, cell, energy, force, virial
        elif 'VOLUME and BASIS' in ii:
            for dd in range(3) :
                tmp_l = lines[idx+5+dd]
                cell.append([float(ss) 
                             for ss in tmp_l.replace('-',' -').split()[0:3]])
        elif 'in kB' in ii:
            tmp_v = [float(ss) for ss in ii.split()[2:8]]
            virial = np.zeros([3,3])
            virial[0][0] = tmp_v[0]
            virial[1][1] = tmp_v[1]
            virial[2][2] = tmp_v[2]
            virial[0][1] = tmp_v[3]
            virial[1][0] = tmp_v[3]
            virial[1][2] = tmp_v[4]
            virial[2][1] = tmp_v[4]
            virial[0][2] = tmp_v[5]
            virial[2][0] = tmp_v[5]
        elif 'TOTAL-FORCE' in ii:
            for jj in range(idx+2, idx+2+ntot) :
                tmp_l = lines[jj]
                info = [float(ss) for ss in tmp_l.split()]
                coord.append(info[:3])
                force.append(info[3:])
    return coord, cell, energy, force, virial
