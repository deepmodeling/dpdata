from rdkit import Chem
import numpy as np

def mol_to_system_data(mol):
    mol = regularize_formal_charges(mol)
    if not mol:
        raise RuntimeError("Sanitize Failed, please check your input file")

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
        return data
    else:
        raise ValueError("The moleclue does not contain 3-D conformers")


def regularize_formal_charges(mol):
    """
    Regularize formal charges of atoms
    """
    for atom in mol.GetAtoms():
        assign_formal_charge_for_atom(atom)
    try:
        Chem.SanitizeMol(mol)
        return mol
    except:
        return None


def assign_formal_charge_for_atom(atom):
    """
        assigen formal charge according to 8-electron rule for element B,C,N,O,S,P,As
    """
    if atom.GetSymbol() == "B":
        atom.SetFormalCharge(3 - atom.GetExplicitValence())
    elif atom.GetSymbol() == "C":
        atom.SetFormalCharge(atom.GetExplicitValence() - 4)
        if atom.GetExplicitValence() == 3:
            print(f"Detect a valence of 3 on carbon #{atom.GetIdx()}, the formal charge of this atom will be assigned to -1")
    elif atom.GetSymbol() == "N":
        atom.SetFormalCharge(atom.GetExplicitValence() - 3)
    elif atom.GetSymbol() == "O":
        atom.SetFormalCharge(atom.GetExplicitValence() - 2)
    elif atom.GetSymbol() == "S":
        if atom.GetExplicitValence() == 1:
            atom.SetFormalCharge(-1)
        elif atom.GetExplicitValence() == 3:
            atom.SetFormalCharge(1)
        else:
            atom.SetFormalCharge(0)
    elif atom.GetSymbol() == "P" or atom.GetSymbol() == "As":
        if atom.GetExplicitValence() == 5:
            atom.SetFormalCharge(0)
        else:
            atom.SetFormalCharge(atom.GetExplicitValence() - 3)


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