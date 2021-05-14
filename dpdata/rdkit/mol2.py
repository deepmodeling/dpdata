import os
import shutil
from .utils import system_data_to_mol
import re
import numpy as np

def sep_mol2_file(fname):
    with open(fname) as f:
        content = f.read()
    title = "@<TRIPOS>MOLECULE"
    mol2_blocks = [title + block for block in content.split(title)[1:]]
    return mol2_blocks

def _format_atom_name(atom_name):
    patt = re.compile("[a-zA-Z]*")
    match = re.search(patt, atom_name)
    fmt_name = match.group().capitalize()
    return fmt_name

def _format_atom_type(atom_type):
    if "." in atom_type:
        fmt_name = atom_type.split('.')[0]
    else:
        fmt_name = atom_type
    if fmt_name in ['Du', 'Any', 'Hal', 'Het', 'Hev']:
        raise ValueError(f"Invalid element: {fmt_name}")
    else:
        return fmt_name

def _read_atom_line(line):
    line = line.split()
    atom_id = int(line[0])
    atom_name = line[1]
    coord = np.array([float(x) for x in line[2: 5]])
    atom_type = line[5]
    try:
        atom_symbol = _format_atom_type(atom_type)
    except ValueError:
        atom_symbol = _format_atom_name(atom_name)
    atom = {
        "atom_id": atom_id,
        "atom_name": atom_name,
        "coord": coord,
        "atom_type": atom_type,
        "atom_symbol": atom_symbol
    }
    return atom


def _read_bond_line(line):
    line = line.split()
    begin_atom_idx = int(line[1]) - 1
    end_atom_idx = int(line[2]) - 1
    if line[3] not in ['du', 'un', 'nc']:
        if line[3] == "ar":
            bond_type = 1.5
        elif line[3] == "am":
            bond_type = 1
        else:
            bond_type = float(line[3])
        return np.array([begin_atom_idx, end_atom_idx, bond_type])
    else:
        print(f"Invalid bond type between atom {begin_atom_idx} and atom {end_atom_idx}: {line[3]}")
        return None

def mol2_block_to_data(block):
    lines = block.split("\n")
    data = {}

    READ_ATOM = 1
    READ_BOND = 2

    flag = 0
    ii = 0
    atoms = []
    bonds = []
    while ii < len(lines):
        if lines[ii] == "@<TRIPOS>MOLECULE":
            mol_name = lines[ii + 1]
            n_atoms, n_bonds, n_substructures = map(int, lines[ii + 2].split()[:3])
            if n_substructures > 1:
                print("Multiple substructures are not supported yet.")
        elif lines[ii] == "@<TRIPOS>ATOM":
            flag = READ_ATOM
        elif lines[ii] == "@<TRIPOS>BOND":
            flag = READ_BOND
        elif lines[ii] == "@<TRIPOS>SUBSTRUCTURE":
            flag = 0
        elif flag == READ_ATOM:
            atoms.append(_read_atom_line(lines[ii]))
        elif flag == READ_BOND:
            bond = _read_bond_line(lines[ii])
            if bond is not None:
                bonds.append(bond)
        ii += 1
    
    atom_symbols = [atom['atom_symbol'] for atom in atoms]
    # print(atoms)
    atom_names, atom_types, atom_numbs = np.unique(atom_symbols, return_inverse=True, return_counts=True)
    data['_name'] = mol_name
    data['orig'] = np.array([0., 0., 0.])
    data['atom_names'] = list(atom_names)
    data['atom_types'] = list(atom_types)
    data['atom_numbs'] = list(atom_numbs)
    data['coords'] = np.array([[atom['coord'] for atom in atoms]])
    data['bonds'] = np.array(bonds)
    data['cells'] = np.array([[[100., 0., 0.],
                               [0., 100., 0.],
                               [0., 0., 100.]]])
    return data
        
        
