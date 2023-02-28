import warnings

import numpy as np

from ..periodic_table import ELEMENTS


def system_info(lines, type_idx_zero=False):
    atom_names = []
    atom_numbs = []
    nelm = 0
    natoms = int(lines[0].split()[0])
    iteration = float(lines[0].split("Etot")[0].split("=")[1].split(",")[0])
    #    print(iteration)
    if iteration > 0:
        nelm = 40
    else:
        nelm = 100
    atomic_number = []
    for idx, ii in enumerate(lines):
        if "Position" in ii:
            for kk in range(idx + 1, idx + 1 + natoms):
                min = kk
                for jj in range(kk + 1, idx + 1 + natoms):
                    if int(lines[jj].split()[0]) < int(lines[min].split()[0]):
                        min = jj
                        lines[min], lines[kk] = lines[kk], lines[min]
            for gg in range(idx + 1, idx + 1 + natoms):
                tmpn = int(lines[gg].split()[0])
                atomic_number.append(tmpn)
    for ii in np.unique(sorted(atomic_number)):
        atom_numbs.append(atomic_number.count(ii))
    atom_types = []
    for idx, ii in enumerate(atom_numbs):
        for jj in range(ii):
            if type_idx_zero:
                atom_types.append(idx)
            else:
                atom_types.append(idx + 1)
    for ii in np.unique(sorted(atomic_number)):
        atom_names.append(ELEMENTS[ii - 1])
    return atom_names, atom_numbs, np.array(atom_types, dtype=int), nelm


def get_movement_block(fp):
    blk = []
    for ii in fp:
        if not ii:
            return blk
        blk.append(ii.rstrip("\n"))
        if "------------" in ii:
            return blk
    return blk


# we assume that the force is printed ...
def get_frames(fname, begin=0, step=1, convergence_check=True):
    fp = open(fname)
    blk = get_movement_block(fp)

    atom_names, atom_numbs, atom_types, nelm = system_info(blk, type_idx_zero=True)
    ntot = sum(atom_numbs)

    all_coords = []
    all_cells = []
    all_energies = []
    all_atomic_energy = []
    all_forces = []
    all_virials = []

    cc = 0
    rec_failed = []
    while len(blk) > 0:
        if cc >= begin and (cc - begin) % step == 0:
            coord, cell, energy, force, virial, is_converge = analyze_block(
                blk, ntot, nelm
            )
            if len(coord) == 0:
                break
            if is_converge or not convergence_check:
                all_coords.append(coord)
                all_cells.append(cell)
                all_energies.append(energy)
                all_forces.append(force)
                if virial is not None:
                    all_virials.append(virial)
            if not is_converge:
                rec_failed.append(cc + 1)

        blk = get_movement_block(fp)
        cc += 1

    if len(rec_failed) > 0:
        prt = (
            "so they are not collected."
            if convergence_check
            else "but they are still collected due to the requirement for ignoring convergence checks."
        )
        warnings.warn(
            f"The following structures were unconverged: {rec_failed}; " + prt
        )

    if len(all_virials) == 0:
        all_virials = None
    else:
        all_virials = np.array(all_virials)
    fp.close()
    return (
        atom_names,
        atom_numbs,
        atom_types,
        np.array(all_cells),
        np.array(all_coords),
        np.array(all_energies),
        np.array(all_forces),
        all_virials,
    )


def analyze_block(lines, ntot, nelm):
    coord = []
    cell = []
    energy = None
    #    atomic_energy = []
    force = []
    virial = None
    is_converge = True
    sc_index = 0
    for idx, ii in enumerate(lines):
        if "Iteration" in ii:
            sc_index = int(ii.split("SCF =")[1])
            if sc_index >= nelm:
                is_converge = False
            energy = float(
                ii.split("Etot,Ep,Ek (eV)")[1].split()[2]
            )  # use Ep, not Etot=Ep+Ek
        elif "----------" in ii:
            assert (force is not None) and len(coord) > 0 and len(cell) > 0
            # all_coords.append(coord)
            # all_cells.append(cell)
            # all_energies.append(energy)
            # all_forces.append(force)
            # if virial is not None :
            #     all_virials.append(virial)
            return coord, cell, energy, force, virial, is_converge
        #        elif 'NPT' in ii:
        #            tmp_v = []
        elif "Lattice vector" in ii:
            if "stress" in lines[idx + 1]:
                tmp_v = []
                for dd in range(3):
                    tmp_l = lines[idx + 1 + dd]
                    cell.append([float(ss) for ss in tmp_l.split()[0:3]])
                    tmp_v.append([float(stress) for stress in tmp_l.split()[5:8]])
                virial = np.zeros([3, 3])
                virial[0][0] = tmp_v[0][0]
                virial[0][1] = tmp_v[0][1]
                virial[0][2] = tmp_v[0][2]
                virial[1][0] = tmp_v[1][0]
                virial[1][1] = tmp_v[1][1]
                virial[1][2] = tmp_v[1][2]
                virial[2][0] = tmp_v[2][0]
                virial[2][1] = tmp_v[2][1]
                virial[2][2] = tmp_v[2][2]
                volume = np.linalg.det(np.array(cell))
                virial = virial * 160.2 * 10.0 / volume
            else:
                for dd in range(3):
                    tmp_l = lines[idx + 1 + dd]
                    cell.append([float(ss) for ss in tmp_l.split()[0:3]])

        #            else :
        #                for dd in range(3) :
        #                    tmp_l = lines[idx+1+dd]
        #                    cell.append([float(ss)
        #                                 for ss in tmp_l.split()[0:3]])
        #                virial = np.zeros([3,3])
        elif "Position" in ii:
            for kk in range(idx + 1, idx + 1 + ntot):
                min = kk
                for jj in range(kk + 1, idx + 1 + ntot):
                    if int(lines[jj].split()[0]) < int(lines[min].split()[0]):
                        min = jj
                        lines[min], lines[kk] = lines[kk], lines[min]
            for gg in range(idx + 1, idx + 1 + ntot):
                info = [float(jj) for jj in lines[gg].split()[1:4]]
                info = np.matmul(np.array(info), np.array(cell))
                coord.append(info)
        elif "Force" in ii:
            for kk in range(idx + 1, idx + 1 + ntot):
                min = kk
                for jj in range(kk + 1, idx + 1 + ntot):
                    if int(lines[jj].split()[0]) < int(lines[min].split()[0]):
                        min = jj
                        lines[min], lines[kk] = lines[kk], lines[min]
            for gg in range(idx + 1, idx + 1 + ntot):
                info = [
                    -float(ss) for ss in lines[gg].split()
                ]  # forces in MOVEMENT file are dE/dR, lacking a minus sign
                force.append(info[1:4])
    #        elif 'Atomic-Energy' in ii:
    #            for jj in range(idx+1, idx+1+ntot) :
    #                tmp_l = lines[jj]
    #                info = [float(ss) for ss in tmp_l.split()]
    #                atomic_energy.append(info[1])
    return coord, cell, energy, force, virial, is_converge
