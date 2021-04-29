from rdkit import Chem
import numpy as np

def mol_to_system_data(mol):
    num_confs = mol.GetNumConformers()
    if num_confs:
        atom_symbols = [at.GetSymbol() for at in mol.GetAtoms()]
        atom_names, atom_types, atom_numbs = np.unique(atom_symbols, return_inverse=True, return_counts=True)
        coords = np.array([conf.GetPositions() for conf in mol.GetConformers()])
        bonds = [[bond.GetBeginAtomIdx(),
                  bond.GetEndAtomIdx(),
                  bond.GetBondTypeAsDouble()] for bond in mol.GetBonds()]
        formal_charges = [at.GetFormalCharge() for at in mol.GetAtoms()]
        data = {}
        data['atom_numbs'] = atom_numbs
        data['atom_names'] = atom_names
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
    Regularize formal charges of N, O
    """
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "N" and atom.GetExplicitValence() == 4:
            atom.SetFormalCharge(1)
        elif atom.GetSymbol() == "N" and atom.GetExplicitValence() == 2:
            atom.SetFormalCharge(-1)
        elif atom.GetSymbol() == "O" and atom.GetExplicitValence() == 3:
            atom.SetFormalCharge(1)
        elif atom.GetSymbol() == "O" and atom.GetExplicitValence() == 1:
            atom.SetFormalCharge(-1)
        elif atom.GetSymbol() == "S" and atom.GetExplicitValence() == 3:
            atom.SetFormalCharge(1)
        elif atom.GetSymbol() == "S" and atom.GetExplicitValence() == 1:
            atom.SetFormalCharge(-1)
    Chem.SanitizeMol(mol)
    return mol

