import re
import warnings

import numpy as np


def system_info(lines, type_idx_zero=False):
    atom_names = []
    atom_numbs = None
    nelm = None
    for ii in lines:
        ii_word_list = ii.split()
        if "TITEL" in ii:
            # get atom names from POTCAR info, tested only for PAW_PBE ...
            _ii = ii.split()[3]
            if "_" in _ii:
                # for case like : TITEL  = PAW_PBE Sn_d 06Sep2000
                atom_names.append(_ii.split("_")[0])
            else:
                atom_names.append(_ii)
        # a stricker check for "NELM"; compatible with distingct formats in different versions(6 and older, newers_expect-to-work) of vasp
        elif nelm is None:
            m = re.search(r"NELM\s*=\s*(\d+)", ii)
            if m:
                nelm = int(m.group(1))
        if "ions per type" in ii:
            atom_numbs_ = [int(s) for s in ii.split()[4:]]
            if atom_numbs is None:
                atom_numbs = atom_numbs_
            else:
                assert atom_numbs == atom_numbs_, "in consistent numb atoms in OUTCAR"
    assert nelm is not None, "cannot find maximum steps for each SC iteration"
    assert atom_numbs is not None, "cannot find ion type info in OUTCAR"
    atom_names = atom_names[: len(atom_numbs)]
    atom_types = []
    for idx, ii in enumerate(atom_numbs):
        for jj in range(ii):
            if type_idx_zero:
                atom_types.append(idx)
            else:
                atom_types.append(idx + 1)
    return atom_names, atom_numbs, np.array(atom_types, dtype=int), nelm


def get_outcar_block(fp, ml=False):
    blk = []
    energy_token = ["free  energy   TOTEN", "free  energy ML TOTEN"]
    ml_index = int(ml)
    for ii in fp:
        if not ii:
            return blk
        blk.append(ii.rstrip("\n"))
        if energy_token[ml_index] in ii:
            return blk
    return blk


# we assume that the force is printed ...
def get_frames(fname, begin=0, step=1, ml=False, convergence_check=True):
    fp = open(fname)
    blk = get_outcar_block(fp)

    atom_names, atom_numbs, atom_types, nelm = system_info(blk, type_idx_zero=True)
    ntot = sum(atom_numbs)

    all_coords = []
    all_cells = []
    all_energies = []
    all_forces = []
    all_virials = []

    cc = 0
    rec_failed = []
    while len(blk) > 0:
        if cc >= begin and (cc - begin) % step == 0:
            coord, cell, energy, force, virial, is_converge = analyze_block(
                blk, ntot, nelm, ml
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

        blk = get_outcar_block(fp, ml)
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


def analyze_block(lines, ntot, nelm, ml=False):
    coord = []
    cell = []
    energy = None
    force = []
    virial = None
    is_converge = True
    sc_index = 0
    # select different searching tokens based on the ml label
    energy_token = ["free  energy   TOTEN", "free  energy ML TOTEN"]
    energy_index = [4, 5]
    virial_token = ["FORCE on cell =-STRESS in cart. coord.  units", "ML FORCE"]
    virial_index = [14, 4]
    cell_token = ["VOLUME and BASIS", "ML FORCE"]
    cell_index = [5, 12]
    ml_index = int(ml)
    for idx, ii in enumerate(lines):
        # if set ml == True, is_converged will always be True
        if ("Iteration" in ii) and (not ml):
            sc_index = int(ii.split()[3][:-1])
            if sc_index >= nelm:
                is_converge = False
        elif energy_token[ml_index] in ii:
            energy = float(ii.split()[energy_index[ml_index]])
            assert (force is not None) and len(coord) > 0 and len(cell) > 0
            return coord, cell, energy, force, virial, is_converge
        elif cell_token[ml_index] in ii:
            for dd in range(3):
                tmp_l = lines[idx + cell_index[ml_index] + dd]
                cell.append([float(ss) for ss in tmp_l.replace("-", " -").split()[0:3]])
        elif virial_token[ml_index] in ii:
            in_kB_index = virial_index[ml_index]
            while idx + in_kB_index < len(lines) and (
                not lines[idx + in_kB_index].split()[0:2] == ["in", "kB"]
            ):
                in_kB_index += 1
            assert idx + in_kB_index < len(
                lines
            ), 'ERROR: "in kB" is not found in OUTCAR. Unable to extract virial.'
            tmp_v = [float(ss) for ss in lines[idx + in_kB_index].split()[2:8]]
            virial = np.zeros([3, 3])
            virial[0][0] = tmp_v[0]
            virial[1][1] = tmp_v[1]
            virial[2][2] = tmp_v[2]
            virial[0][1] = tmp_v[3]
            virial[1][0] = tmp_v[3]
            virial[1][2] = tmp_v[4]
            virial[2][1] = tmp_v[4]
            virial[0][2] = tmp_v[5]
            virial[2][0] = tmp_v[5]
        elif "TOTAL-FORCE" in ii and (("ML" in ii) == ml):
            for jj in range(idx + 2, idx + 2 + ntot):
                tmp_l = lines[jj]
                info = [float(ss) for ss in tmp_l.split()]
                coord.append(info[:3])
                force.append(info[3:6])
    return coord, cell, energy, force, virial, is_converge
