#!/usr/bin/env python3

import dpdata.vasp.outcar as outcar

# Test the ML OUTCAR parsing
fname = "tests/poscars/OUTCAR.ch4.ml"

print("=== Testing ML mode ===")
result_ml = outcar.get_frames(fname, ml=True)
print(f"ML mode frames: {len(result_ml[4])}")  # coords

print("=== Testing non-ML mode ===")
result_nonml = outcar.get_frames(fname, ml=False)
print(f"Non-ML mode frames: {len(result_nonml[4])}")  # coords

# Let's debug the analyze_block function by patching it temporarily
original_analyze_block = outcar.analyze_block

def debug_analyze_block(lines, ntot, nelm, ml=False):
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
    
    print(f"\n--- Debug analyze_block: ml={ml}, ml_index={ml_index} ---")
    print(f"Looking for energy_token: '{energy_token[ml_index]}'")
    print(f"Looking for cell_token: '{cell_token[ml_index]}'")
    print(f"Looking for virial_token: '{virial_token[ml_index]}'")
    
    found_energy = False
    found_cell = False
    found_virial = False
    found_force = False
    
    for idx, ii in enumerate(lines):
        # if set ml == True, is_converged will always be True
        if ("Iteration" in ii) and (not ml):
            sc_index = int(ii.split()[3][:-1])
            if sc_index >= nelm:
                is_converge = False
        elif energy_token[ml_index] in ii:
            energy = float(ii.split()[energy_index[ml_index]])
            found_energy = True
            print(f"Found energy: {energy}")
            return coord, cell, energy, force, virial, is_converge
        elif cell_token[ml_index] in ii:
            found_cell = True
            print(f"Found cell_token at line {idx}: {ii.strip()}")
            for dd in range(3):
                if idx + cell_index[ml_index] + dd < len(lines):
                    tmp_l = lines[idx + cell_index[ml_index] + dd]
                    print(f"  Cell line {dd}: {tmp_l.strip()}")
                    cell.append([float(ss) for ss in tmp_l.replace("-", " -").split()[0:3]])
        elif virial_token[ml_index] in ii:
            found_virial = True
            print(f"Found virial_token at line {idx}: {ii.strip()}")
            in_kB_index = virial_index[ml_index]
            while idx + in_kB_index < len(lines) and (
                not lines[idx + in_kB_index].split()[0:2] == ["in", "kB"]
            ):
                in_kB_index += 1
            if idx + in_kB_index < len(lines):
                tmp_v = [float(ss) for ss in lines[idx + in_kB_index].split()[2:8]]
                virial = [[tmp_v[0], tmp_v[3], tmp_v[5]], 
                         [tmp_v[3], tmp_v[1], tmp_v[4]], 
                         [tmp_v[5], tmp_v[4], tmp_v[2]]]
        elif "TOTAL-FORCE" in ii and (("ML" in ii) == ml):
            found_force = True
            print(f"Found TOTAL-FORCE at line {idx}: {ii.strip()}")
            for jj in range(idx + 2, min(idx + 2 + ntot, len(lines))):
                tmp_l = lines[jj]
                info = [float(ss) for ss in tmp_l.split()]
                coord.append(info[:3])
                force.append(info[3:6])
    
    print(f"Summary: energy={found_energy}, cell={found_cell}, virial={found_virial}, force={found_force}")
    print(f"Final: coord={len(coord)}, cell={len(cell)}, energy={energy}")
    return coord, cell, energy, force, virial, is_converge

# Temporarily replace the function
outcar.analyze_block = debug_analyze_block

print("\n=== Debug ML mode (first block) ===")
with open(fname) as fp:
    blk = outcar.get_outcar_block(fp, ml=True)
    atom_names, atom_numbs, atom_types, nelm, nwrite = outcar.system_info(blk, type_idx_zero=True)
    ntot = sum(atom_numbs)
    print(f"ntot={ntot}, nelm={nelm}, nwrite={nwrite}")
    coord, cell, energy, force, virial, is_converge = debug_analyze_block(blk, ntot, nelm, ml=True)

print("\n=== Debug non-ML mode (first block) ===")
with open(fname) as fp:
    blk = outcar.get_outcar_block(fp, ml=False)
    coord, cell, energy, force, virial, is_converge = debug_analyze_block(blk, ntot, nelm, ml=False)

# Restore original
outcar.analyze_block = original_analyze_block