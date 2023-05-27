import os
import time
from copy import deepcopy

from rdkit import Chem
from rdkit.Chem.rdchem import BondType

# openbabel
try:
    from openbabel import openbabel

    USE_OBABEL = True
except ModuleNotFoundError as e:
    USE_OBABEL = False


def get_explicit_valence(atom, verbose=False):
    exp_val_calculated_from_bonds = int(
        sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
    )
    try:
        exp_val = atom.GetExplicitValence()
        if exp_val != exp_val_calculated_from_bonds:
            if verbose:
                print(
                    f"Explicit valence given by GetExplicitValence() and sum of bond order are inconsistent on {atom.GetSymbol()}{atom.GetIdx() + 1}, using sum of bond order."
                )
        return exp_val_calculated_from_bonds
    except Exception:
        return exp_val_calculated_from_bonds


def regularize_formal_charges(mol, sanitize=True, verbose=False):
    """Regularize formal charges of atoms."""
    assert isinstance(mol, Chem.rdchem.Mol)
    for atom in mol.GetAtoms():
        assign_formal_charge_for_atom(atom, verbose)
    if sanitize:
        try:
            Chem.SanitizeMol(mol)
            return mol
        except Exception:
            return None
    else:
        return mol


def assign_formal_charge_for_atom(atom, verbose=False):
    """Assigen formal charge according to 8-electron rule for element B,C,N,O,S,P,As."""
    assert isinstance(atom, Chem.rdchem.Atom)
    valence = get_explicit_valence(atom, verbose)
    if atom.GetSymbol() == "B":
        atom.SetFormalCharge(3 - valence)
    elif atom.GetSymbol() == "C":
        atom.SetFormalCharge(valence - 4)
        if valence == 3:
            print(
                f"Detect a valence of 3 on #C{atom.GetIdx() + 1}, the formal charge of this atom will be assigned to -1"
            )
        elif valence > 4:
            raise ValueError(f"#C{atom.GetIdx() + 1} has a valence larger than 4")
    elif atom.GetSymbol() == "N":
        if valence > 4:
            raise ValueError(f"#N{atom.GetIdx() + 1} has a valence larger than 4")
        else:
            atom.SetFormalCharge(valence - 3)
    elif atom.GetSymbol() == "O":
        atom.SetFormalCharge(valence - 2)
    elif atom.GetSymbol() == "S":
        if valence == 1:
            atom.SetFormalCharge(-1)
        elif valence == 3:
            atom.SetFormalCharge(1)
        elif valence > 6:
            raise ValueError(f"#S{atom.GetIdx() + 1} has a valence larger than 6")
        else:
            atom.SetFormalCharge(0)
    elif atom.GetSymbol() == "P" or atom.GetSymbol() == "As":
        if valence == 5:
            atom.SetFormalCharge(0)
        elif valence > 5:
            raise ValueError(
                f"#{atom.GetSymbol()}{atom.GetIdx() + 1} has a valence larger than 5"
            )
        else:
            atom.SetFormalCharge(valence - 3)


# print bond and atom information (for debugger)
def print_bonds(mol):
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        print(
            f"{begin_atom.GetSymbol()}{begin_atom.GetIdx() + 1} {end_atom.GetSymbol()}{end_atom.GetIdx() + 1} {bond.GetBondType()}"
        )


def print_atoms(mol):
    for atom in mol.GetAtoms():
        print(
            f"{atom.GetSymbol()}{atom.GetIdx() + 1} {atom.GetFormalCharge()} {get_explicit_valence(atom)}"
        )


def is_terminal_oxygen(O_atom):
    return len(O_atom.GetNeighbors()) == 1


def get_terminal_oxygens(atom):
    terminal_oxygens = []
    for nei in atom.GetNeighbors():
        if nei.GetSymbol() == "O" or nei.GetSymbol() == "S":
            if is_terminal_oxygen(nei):
                terminal_oxygens.append(nei)
    return terminal_oxygens


def is_terminal_NR2(N_atom):
    return len(N_atom.GetNeighbors()) == 3


def get_terminal_NR2s(atom):
    terminal_NR2s = []
    for nei in atom.GetNeighbors():
        if nei.GetSymbol() == "N":
            if is_terminal_NR2(nei):
                terminal_NR2s.append(nei)
    terminal_NR2s.sort(
        key=lambda N_atom: len(
            [atom for atom in N_atom.GetNeighbors() if atom.GetSymbol() == "H"]
        )
    )
    return terminal_NR2s


def sanitize_phosphate_Patom(P_atom, verbose=True):
    if P_atom.GetSymbol() == "P":
        terminal_oxygens = get_terminal_oxygens(P_atom)
        mol = P_atom.GetOwningMol()
        if len(terminal_oxygens) > 1:
            if verbose:
                print("Phospate group detected, sanitizing it...")
            # set one P=O and two P-O
            bond1 = mol.GetBondBetweenAtoms(
                P_atom.GetIdx(), terminal_oxygens[0].GetIdx()
            )
            bond1.SetBondType(Chem.rdchem.BondType.DOUBLE)
            for ii in range(1, len(terminal_oxygens)):
                bond = mol.GetBondBetweenAtoms(
                    P_atom.GetIdx(), terminal_oxygens[ii].GetIdx()
                )
                bond.SetBondType(Chem.rdchem.BondType.SINGLE)
                terminal_oxygens[ii].SetFormalCharge(-1)


def sanitize_phosphate(mol):
    for atom in mol.GetAtoms():
        sanitize_phosphate_Patom(atom)
    return mol


def sanitize_sulfate_Satom(S_atom, verbose=True):
    if S_atom.GetSymbol() == "S":
        terminal_oxygens = get_terminal_oxygens(S_atom)
        mol = S_atom.GetOwningMol()
        if len(terminal_oxygens) == 3:
            if verbose:
                print("Sulfate group detected, sanitizing it...")
            # set one S-O and two S=O
            bond1 = mol.GetBondBetweenAtoms(
                S_atom.GetIdx(), terminal_oxygens[0].GetIdx()
            )
            bond1.SetBondType(Chem.rdchem.BondType.SINGLE)
            terminal_oxygens[0].SetFormalCharge(-1)
            for ii in range(1, len(terminal_oxygens)):
                bond = mol.GetBondBetweenAtoms(
                    S_atom.GetIdx(), terminal_oxygens[ii].GetIdx()
                )
                bond.SetBondType(Chem.rdchem.BondType.DOUBLE)


def sanitize_sulfate(mol):
    for atom in mol.GetAtoms():
        sanitize_sulfate_Satom(atom)
    return mol


def sanitize_carboxyl_Catom(C_atom, verbose=True):
    if C_atom.GetSymbol() == "C":
        terminal_oxygens = get_terminal_oxygens(C_atom)
        mol = C_atom.GetOwningMol()
        if len(terminal_oxygens) == 2:
            if verbose:
                print("Carbonxyl group detected, sanitizing it...")
            # set one C-O and one C=O
            bond1 = mol.GetBondBetweenAtoms(
                C_atom.GetIdx(), terminal_oxygens[0].GetIdx()
            )
            bond1.SetBondType(Chem.rdchem.BondType.SINGLE)
            terminal_oxygens[0].SetFormalCharge(-1)

            bond2 = mol.GetBondBetweenAtoms(
                C_atom.GetIdx(), terminal_oxygens[1].GetIdx()
            )
            bond2.SetBondType(Chem.rdchem.BondType.DOUBLE)
            terminal_oxygens[1].SetFormalCharge(0)


def sanitize_carboxyl(mol):
    for atom in mol.GetAtoms():
        sanitize_carboxyl_Catom(atom)
    return mol


def sanitize_guanidine_Catom(C_atom, verbose=True):
    if C_atom.GetSymbol() == "C":
        terminal_NR2s = get_terminal_NR2s(C_atom)
        mol = C_atom.GetOwningMol()
        if len(terminal_NR2s) == 3:
            if verbose:
                print("Guanidyl group detected, sanitizing it...")
            # set two C-N and one C=N+
            bond1 = mol.GetBondBetweenAtoms(C_atom.GetIdx(), terminal_NR2s[0].GetIdx())
            bond1.SetBondType(Chem.rdchem.BondType.SINGLE)
            terminal_NR2s[0].SetFormalCharge(-1)

            bond2 = mol.GetBondBetweenAtoms(C_atom.GetIdx(), terminal_NR2s[1].GetIdx())
            bond2.SetBondType(Chem.rdchem.BondType.SINGLE)
            terminal_NR2s[1].SetFormalCharge(0)

            bond3 = mol.GetBondBetweenAtoms(C_atom.GetIdx(), terminal_NR2s[2].GetIdx())
            bond3.SetBondType(Chem.rdchem.BondType.DOUBLE)
            terminal_NR2s[2].SetFormalCharge(1)


def sanitize_guanidine(mol):
    for atom in mol.GetAtoms():
        sanitize_guanidine_Catom(atom)
    return mol


def sanitize_nitro_Natom(N_atom, verbose=True):
    if N_atom.GetSymbol() == "N":
        terminal_oxygens = get_terminal_oxygens(N_atom)
        mol = N_atom.GetOwningMol()
        if len(terminal_oxygens) == 2:
            if verbose:
                print("Nitro group detected, sanitizing it...")
            # set one N-O and one N=O
            bond1 = mol.GetBondBetweenAtoms(
                N_atom.GetIdx(), terminal_oxygens[0].GetIdx()
            )
            bond1.SetBondType(Chem.rdchem.BondType.SINGLE)
            terminal_oxygens[0].SetFormalCharge(-1)

            bond2 = mol.GetBondBetweenAtoms(
                N_atom.GetIdx(), terminal_oxygens[1].GetIdx()
            )
            bond2.SetBondType(Chem.rdchem.BondType.DOUBLE)
            terminal_oxygens[1].SetFormalCharge(0)


def sanitize_nitro(mol):
    for atom in mol.GetAtoms():
        sanitize_nitro_Natom(atom)
    return mol


def is_terminal_nitrogen(N_atom):
    if N_atom.GetSymbol() == "N" and len(N_atom.GetNeighbors()) == 1:
        return True
    else:
        return False


def sanitize_nitrine_Natom(atom, verbose=True):
    if atom.GetSymbol() == "N" and len(atom.GetNeighbors()) == 2:
        mol = atom.GetOwningMol()
        nei1, nei2 = atom.GetNeighbors()[0], atom.GetNeighbors()[1]
        if nei1.GetSymbol() == "N" and nei2.GetSymbol() == "N":
            if is_terminal_nitrogen(nei1):
                N_terminal = nei1
                N_non_terminal = nei2
            elif is_terminal_nitrogen(nei2):
                N_terminal = nei2
                N_non_terminal = nei1
            else:
                N_terminal = None
                N_non_terminal = None
            if (N_terminal is not None) and (N_non_terminal is not None):
                # set X-N=[N+]=[N-]
                if verbose:
                    print("Detecting nitrine group, fixing it...")
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), N_terminal.GetIdx())
                bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
                N_terminal.SetFormalCharge(-1)

                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), N_non_terminal.GetIdx())
                bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
                atom.SetFormalCharge(1)


def contain_hetero_aromatic(mol):
    flag = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C" and atom.GetIsAromatic():
            flag = True
            break
    return flag


# for carbon with explicit valence > 4
def regularize_carbon_bond_order(atom, verbose=True):
    if atom.GetSymbol() == "C" and get_explicit_valence(atom) > 4:
        if verbose:
            print("Detecting carbon with explicit valence > 4, fixing it...")
        mol = atom.GetOwningMol()
        double_bond_idx = -1
        for nei in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nei.GetIdx())
            if bond.GetBondTypeAsDouble() == 2:
                double_bond_idx = bond.GetIdx()
                break
        if double_bond_idx != -1:
            for bond in atom.GetBonds():
                if bond.GetIdx() != double_bond_idx:
                    bond.SetBondType(Chem.rdchem.BondType.SINGLE)


# for nitrogen with explicit valence > 4
def regularize_nitrogen_bond_order(atom, verbose=True):
    mol = atom.GetOwningMol()
    if atom.GetSymbol() == "N" and get_explicit_valence(atom) > 4:
        O_atoms = get_terminal_oxygens(atom)
        for O_atom in O_atoms:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), O_atom.GetIdx())
            if bond.GetBondTypeAsDouble() == 2:
                bond.SetBondType(Chem.rdchem.BondType.SINGLE)
                O_atom.SetFormalCharge(-1)


def sanitize_mol(mol, verbose=False):
    for atom in mol.GetAtoms():
        sanitize_carboxyl_Catom(atom, verbose)
        sanitize_guanidine_Catom(atom, verbose)
        sanitize_phosphate_Patom(atom, verbose)
        sanitize_sulfate_Satom(atom, verbose)
        sanitize_nitro_Natom(atom, verbose)
        sanitize_nitrine_Natom(atom, verbose)
        regularize_carbon_bond_order(atom, verbose)
        regularize_nitrogen_bond_order(atom, verbose)
    return mol


# copy from FEprep
def mol_edit_log(mol, i, j):
    if not mol.HasProp("edit"):
        mol.SetProp("edit", "%d_%d" % (i, j))
    else:
        edited = mol.GetProp("edit")
        mol.SetProp("edit", edited + ",%d_%d" % (i, j))


def kekulize_aromatic_heterocycles(mol_in, assign_formal_charge=True, sanitize=True):
    mol = Chem.RWMol(mol_in)
    rings = Chem.rdmolops.GetSymmSSSR(mol)
    rings = [list(i) for i in list(rings)]
    rings.sort(key=lambda r: len(r))

    def search_and_assign_ring(
        mol, ring, hetero, start, forward=True, start_switch=True
    ):
        j = start
        switch = start_switch
        lring = len(ring)
        delta = 1 if forward else -1
        n_edit = 0
        n_double = 0
        while not ((j in hetero) & (not switch)):
            btype = BondType.SINGLE if switch else BondType.DOUBLE
            bond = mol.GetBondBetweenAtoms(ring[j], ring[(j + delta) % lring])
            if bond.GetBondType() == BondType.AROMATIC:
                bond.SetBondType(btype)
                mol_edit_log(mol, ring[j], ring[(j + delta) % lring])
                # print(ring[j], ring[(j + delta) % lring], bond.GetBondType())
                if btype == BondType.DOUBLE:
                    n_double += 1
                n_edit += 1
            else:
                break
            j = (j + delta) % lring
            switch = not switch
        return n_edit, n_double

    def print_bondtypes(mol, rings):
        for ring in rings:
            lring = len(ring)
            btype = []
            for i in range(lring):
                btype.append(
                    mol.GetBondBetweenAtoms(
                        ring[i], ring[(i + 1) % lring]
                    ).GetBondType()
                )
            atoms = [mol.GetAtomWithIdx(i).GetSymbol() for i in ring]
            print(ring)
            print(atoms)
            print(btype)

    def hetero_priority(idx, mol):
        atom = mol.GetAtomWithIdx(idx)
        sym = atom.GetSymbol()
        valence = len(atom.GetBonds())

        if (sym in ["O", "S"]) & (valence == 2):
            return 0
        elif sym in ["N", "P", "As", "B"]:
            if valence == 3:
                return 1
            elif valence == 2:
                return 2

    # save carbon/hetero aromatic rings
    CAr = []
    HAr = []
    for ring in rings:
        lring = len(ring)
        bAllAr = True
        bAllC = True
        for i in range(lring):
            atom = mol.GetAtomWithIdx(ring[i])
            if atom.GetSymbol() != "C":
                bAllC = False

            bond = mol.GetBondBetweenAtoms(ring[i], ring[(i + 1) % lring])
            if bond.GetBondType() != BondType.AROMATIC:
                bAllAr = False
        if bAllAr and bAllC:
            CAr.append(ring)
        elif bAllAr and not bAllC:
            HAr.append(ring)

    if len(HAr) == 0:
        # no hetrerocycles
        return mol_in
    else:
        # edit heterocycles
        for ring in HAr:
            lring = len(ring)
            cring = len(CAr)
            hetero = []
            hasDouble = []
            fuseCAr = []
            fuseDouble = []
            for i in range(lring):
                fuseCAr.append(-1)
                for j in range(cring):
                    if ring[i] in CAr[j]:
                        fuseCAr[i] = j
                        break
                if i > 1:
                    if (fuseCAr[i] == fuseCAr[i - 1]) & (fuseCAr[i] >= 0):
                        fuseDouble.append(i)
                atom = mol.GetAtomWithIdx(ring[i])
                if atom.GetSymbol() != "C":
                    hetero.append(i)
                atom_bonds = atom.GetBonds()
                btype = [bond.GetBondType() for bond in atom_bonds]
                # print(btype)
                if BondType.DOUBLE in btype:
                    hasDouble.append(i)
                bond = mol.GetBondBetweenAtoms(ring[i], ring[(i + 1) % lring])

            if (fuseCAr[0] == fuseCAr[lring - 1]) & (fuseCAr[0] >= 0):
                fuseDouble.append(0)

            if (len(hetero) > 0) | (len(hasDouble) > 0):
                n_targetDouble = lring // 2
                n_targetEdit = lring
                hetero_prior = {i: hetero_priority(ring[i], mol) for i in hetero}
                hetero.sort(key=lambda i: hetero_prior[i])
                for i in hasDouble:
                    d1, e1 = search_and_assign_ring(mol, ring, hetero, i, forward=True)
                    d2, e2 = search_and_assign_ring(mol, ring, hetero, i, forward=False)
                    n_targetDouble -= d1 + d2 + 1
                    n_targetEdit -= e1 + e2
                for i in fuseDouble:
                    bond = mol.GetBondBetweenAtoms(ring[i], ring[(i - 1) % lring])
                    if bond.GetBondType() == BondType.AROMATIC:
                        bond.SetBondType(BondType.DOUBLE)
                        mol_edit_log(mol, ring[i], ring[(i - 1) % lring])
                    d1, e1 = search_and_assign_ring(mol, ring, hetero, i, forward=True)
                    d2, e2 = search_and_assign_ring(
                        mol, ring, hetero, (i - 1) % lring, forward=False
                    )
                    n_targetDouble -= d1 + d2 + 1
                    n_targetEdit -= e1 + e2 + 1
                for i in hetero:
                    atom = mol.GetAtomWithIdx(ring[i])
                    if (hetero_prior[i] == 2) | (n_targetDouble * 2 >= n_targetEdit):
                        forward_btype = mol.GetBondBetweenAtoms(
                            ring[i], ring[(i + 1) % lring]
                        ).GetBondType()
                        backward_btype = mol.GetBondBetweenAtoms(
                            ring[i], ring[(i - 1) % lring]
                        ).GetBondType()
                        if forward_btype != BondType.AROMATIC:
                            switch = forward_btype == BondType.DOUBLE
                            d1, e1 = search_and_assign_ring(
                                mol, ring, hetero, i, forward=False, start_switch=switch
                            )
                            d2 = e2 = 0
                        elif backward_btype != BondType.AROMATIC:
                            switch = backward_btype == BondType.DOUBLE
                            d1, e1 = search_and_assign_ring(
                                mol, ring, hetero, i, forward=True, start_switch=switch
                            )
                            d2 = e2 = 0
                        else:
                            d1, e1 = search_and_assign_ring(
                                mol, ring, hetero, i, forward=True, start_switch=True
                            )
                            d2, e2 = search_and_assign_ring(
                                mol, ring, hetero, i, forward=False, start_switch=False
                            )
                        n_targetDouble -= d1 + d2
                        n_targetEdit -= e1 + e2
                    else:
                        d1, e1 = search_and_assign_ring(
                            mol, ring, hetero, i, forward=True, start_switch=True
                        )
                        d2, e2 = search_and_assign_ring(
                            mol, ring, hetero, i, forward=False, start_switch=True
                        )
                        n_targetDouble -= d1 + d2
                        n_targetEdit -= e1 + e2

        for ring in CAr:
            lring = len(ring)
            for i in range(lring):
                bond = mol.GetBondBetweenAtoms(ring[i], ring[(i + 1) % lring])
                bond.SetBondType(BondType.AROMATIC)
        print("Manual kekulization for aromatic heterocycles:")
        print_bondtypes(mol, rings)

        atoms = mol.GetAtoms()
        for i in range(len(atoms)):
            mol.ReplaceAtom(i, Chem.Atom(atoms[i].GetSymbol()))
        mol_edited = mol.GetMol()
        # charge assignment
        if assign_formal_charge:
            mol_edited = regularize_formal_charges(mol_edited, sanitize=False)
        if not sanitize:
            return mol_edited
        else:
            try:
                Chem.SanitizeMol(mol_edited)
                return mol_edited
            except Exception as e:
                raise RuntimeError(
                    f"Manual kekulization for aromatic heterocycles failed, below are errors:\n\t {e}"
                )


def convert_by_obabel(
    mol, cache_dir=os.path.join(os.getcwd(), ".cache"), obabel_path="obabel"
):
    if not os.path.exists(cache_dir):
        os.mkdir(cache_dir)
    if mol.HasProp("_Name"):
        name = mol.GetProp("_Name")
    else:
        name = f"mol{int(time.time())}"
    mol_file_in = os.path.join(cache_dir, f"{name}.mol")
    mol_file_out = os.path.join(cache_dir, f"{name}_obabel.mol")
    Chem.MolToMolFile(mol, mol_file_in, kekulize=False)
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "mol")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, mol_file_in)
    obConversion.WriteFile(mol, mol_file_out)
    mol_obabel = Chem.MolFromMolFile(mol_file_out, removeHs=False, sanitize=False)
    return mol_obabel


def super_sanitize_mol(mol, name=None, verbose=True):
    if name is None:
        if mol.HasProp("_Name"):
            name = mol.GetProp("_Name")
        else:
            name = "mol"
    try:
        if verbose:
            print("=====Stage 1: use Hermite procedure=====")
        # use our procedure
        mol = sanitize_mol(mol, verbose)
        mol = regularize_formal_charges(mol, sanitize=False)
        mol_copy = deepcopy(mol)
        Chem.SanitizeMol(mol_copy)
        if verbose:
            print(name, "Success.")
        return mol_copy
    except Exception as e:
        try:
            if verbose:
                print(
                    "Hermite procedure failed, maybe due to unsupported representation of hetero aromatic rings, re-try with obabel"
                )
                print("=====Stage 2: re-try with obabel=====")
            mol = convert_by_obabel(mol)
            mol = sanitize_mol(mol, verbose)
            mol = kekulize_aromatic_heterocycles(
                mol, assign_formal_charge=False, sanitize=False
            )  # aromatic heterocycles
            mol = regularize_formal_charges(mol, sanitize=False)
            mol_copy = deepcopy(mol)
            Chem.SanitizeMol(mol_copy)
            if verbose:
                print(name, "Success.")
            return mol_copy
        except Exception as e:
            if verbose:
                print(e)
                print(name, "Failed!")
            return None


class Sanitizer:
    def __init__(self, level="medium", raise_errors=True, verbose=False):
        """Set up sanitizer.
        --------.

        Parameters
        ----------
        level : 'low', 'medium' or 'high'.
            `low`    - use rdkit.Chem.SanitizeMol() to sanitize
            `medium` - before using rdkit, assign formal charges of each atom first, which requires
                        the rightness of bond order information
            `high`   - try to regularize bond order of nitro, phosphate, sulfate, nitrine, guanidine,
                        pyridine-oxide function groups and aromatic heterocycles. If failed, the program
                        will call obabel to pre-process the mol object and re-try the procedure.
        raise_errors : bool, default=True
            If True, raise SanitizeError when failed.
        verbose : bool, default=False
            If True, print error information when failed.
        """
        self._check_level(level)
        self.level = level
        self.raise_errors = raise_errors
        self.verbose = verbose

    def _check_level(self, level):
        if level not in ["low", "medium", "high"]:
            raise ValueError(
                f"Invalid level '{level}', please set to 'low', 'medium' or 'high'"
            )
        else:
            if level == "high" and not USE_OBABEL:
                raise ModuleNotFoundError(
                    "obabel not installed, high level sanitizer cannot work"
                )

    def _handle_exception(self, error_info):
        if self.raise_errors:
            raise SanitizeError(error_info)
        elif self.verbose:
            print(error_info)

    def sanitize(self, mol):
        """Sanitize mol according to `self.level`. If failed, return None."""
        if self.level == "low":
            try:
                Chem.SanitizeMol(mol)
                return mol
            except Exception as e:
                error_info = f"Sanitization Failed, please use more strict sanitizer by setting 'level' to 'medium' or 'high'. The error occurs:\n\t{e}"
                self._handle_exception(error_info)
                return None
        elif self.level == "medium":
            try:
                mol = regularize_formal_charges(mol, sanitize=False)
                Chem.SanitizeMol(mol)
                return mol
            except Exception as e:
                error_info = f"Sanitization Failed, please use more strict sanitizer by setting 'level' to 'high'. The error occurs:\n\t{e}"
                self._handle_exception(error_info)
                return None
        elif self.level == "high":
            mol = super_sanitize_mol(mol, verbose=self.verbose)
            error_info = "Sanitization Failed. Please check your molecule file."
            if mol is None:
                self._handle_exception(error_info)
            return mol


class SanitizeError(Exception):
    def __init__(self, content="Sanitization Failed."):
        self.content = content

    def __str__(self):
        return self.content

    def __repr__(self):
        return self.__str__()
