from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def mol_to_system_data(mol):
    if not isinstance(mol, Chem.rdchem.Mol):
        raise TypeError(f"rdkit.Chem.Mol required, not {type(mol)}")

    num_confs = mol.GetNumConformers()
    if num_confs:
        atom_symbols = [at.GetSymbol() for at in mol.GetAtoms()]
        atom_names, atom_types, atom_numbs = np.unique(atom_symbols, return_inverse=True, return_counts=True)
        coords = np.array([conf.GetPositions() for conf in mol.GetConformers()])
        bonds = np.array([[bond.GetBeginAtomIdx(),
                           bond.GetEndAtomIdx(),
                           bond.GetBondTypeAsDouble()] for bond in mol.GetBonds()])
        formal_charges = np.array([at.GetFormalCharge() for at in mol.GetAtoms()], dtype=np.int32)
        data = {}
        data['atom_numbs'] = list(atom_numbs)
        data['atom_names'] = list(atom_names)
        data['atom_types'] = atom_types
        data['cells'] = np.array([[[100., 0., 0.],
                                   [0., 100., 0.],
                                   [0., 0., 100.]] for _ in range(num_confs)])
        data['coords'] = coords
        data['bonds'] = bonds
        data['formal_charges'] = formal_charges
        data['orig'] = np.array([0., 0., 0.])
        # other properties
        if mol.HasProp("_Name"):
            data['_name'] = mol.GetProp('_Name')
        return data
    else:
        raise ValueError("The moleclue does not contain 3-D conformers")

def system_data_to_mol(data):
    mol_ed = Chem.RWMol()
    atom_symbols = [data['atom_names'][i] for i in data['atom_types']]
    # add atoms
    for atom_type in data['atom_types']:
        symbol = data['atom_names'][atom_type]
        atom = Chem.Atom(symbol)
        mol_ed.AddAtom(atom)
    # add bonds
    for bond_info in data['bonds']:
        if bond_info[2] == 1:
            mol_ed.AddBond(int(bond_info[0]), int(bond_info[1]), Chem.BondType.SINGLE)
        elif bond_info[2] == 2:
            mol_ed.AddBond(int(bond_info[0]), int(bond_info[1]), Chem.BondType.DOUBLE)
        elif bond_info[2] == 3:
            mol_ed.AddBond(int(bond_info[0]), int(bond_info[1]), Chem.BondType.TRIPLE)
        elif bond_info[2] == 1.5:
            mol_ed.AddBond(int(bond_info[0]), int(bond_info[1]), Chem.BondType.AROMATIC)
    # set conformers
    for frame_idx in range(data['coords'].shape[0]):
        conf = Chem.rdchem.Conformer(len(data['atom_types']))
        for atom_idx in range(len(data['atom_types'])):
            conf.SetAtomPosition(atom_idx, data['coords'][frame_idx][atom_idx])
        mol_ed.AddConformer(conf, assignId=True)
    mol = mol_ed.GetMol()
    # set formal charges
    for idx, atom in enumerate(mol.GetAtoms()):
        atom.SetFormalCharge(data['formal_charges'][idx])
    # set mol name
    if '_name' in list(data.keys()):
        mol.SetProp("_Name", data['_name'])
    # sanitize
    Chem.SanitizeMol(mol_ed)
    return mol        


def check_same_atom(atom_1, atom_2):
    if atom_1.GetIdx() != atom_2.GetIdx():
        return False
    elif atom_1.GetSymbol() != atom_2.GetSymbol():
        return False
    else:
        return True

def check_same_molecule(mol_1, mol_2):
    flag = True
    for bond_1, bond_2 in zip(mol_1.GetBonds(), mol_2.GetBonds()):
        begin_atom_1, end_atom_1 = bond_1.GetBeginAtom(), bond_1.GetEndAtom()
        begin_atom_2, end_atom_2 = bond_2.GetBeginAtom(), bond_2.GetEndAtom()
        if not check_same_atom(begin_atom_1, begin_atom_2):
            flag = False
            break
        elif not check_same_atom(end_atom_1, end_atom_2):
            flag = False
            break
    return flag

def check_molecule_list(mols):
    flag = True
    for mol in mols[1:]:
        if not check_same_molecule(mol, mols[0]):
            flag = False
            break
    return flag

def combine_molecules(mols):
    if check_molecule_list(mols):
        for mol in mols[1:]:
            for conf in mol.GetConformers():
                mols[0].AddConformer(conf, assignId=True)
        return mols[0]
    else:
        raise ValueError("molecules are not of the same topology.")