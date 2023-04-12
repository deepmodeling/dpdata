# %%
# Bond Order System
from copy import deepcopy

import numpy as np
from rdkit.Chem import Conformer

import dpdata.rdkit.utils
from dpdata.rdkit.sanitize import Sanitizer
from dpdata.system import Axis, DataType, System

# import dpdata.rdkit.mol2


class BondOrderSystem(System):
    """The system with chemical bond and formal charges information.

    For example, a labeled methane system named `d_example` has one molecule (5 atoms, 4 bonds) and `n_frames` frames. The bond order and formal charge information can be accessed by
        - `d_example['bonds']` : a numpy array of size 4 x 3, and
                                    the first column represents the index of begin atom,
                                    the second column represents the index of end atom,
                                    the third columen represents the bond order:
                                        1 - single bond, 2 - double bond, 3 - triple bond, 1.5 - aromatic bond
        - `d_example['formal_charges']` : a numpy array of size 5 x 1
    """

    DTYPES = System.DTYPES + (
        DataType("bonds", np.ndarray, (Axis.NBONDS, 3)),
        DataType("formal_charges", np.ndarray, (Axis.NATOMS,)),
    )

    def __init__(
        self,
        file_name=None,
        fmt="auto",
        type_map=None,
        begin=0,
        step=1,
        data=None,
        rdkit_mol=None,
        sanitize_level="medium",
        raise_errors=True,
        verbose=False,
        **kwargs,
    ):
        """Constructor.

        Parameters
        ----------
        file_name : str
            The file to load the system
        fmt : str
            Format of the file, supported formats are
                - ``auto`` : inferred from `file_name`'s extention
                - ``mol`` : .mol file
                - ``sdf`` : .sdf file
        type_map : list of str
            Needed by formats deepmd/raw and deepmd/npy. Maps atom type to name. The atom with type `ii` is mapped to `type_map[ii]`.
            If not provided the atom names are assigned to `'Type_1'`, `'Type_2'`, `'Type_3'`...
        begin : int
            The beginning frame when loading MD trajectory.
        step : int
            The number of skipped frames when loading MD trajectory.
        data : dict
            System data dict.
        rdkit_mol : rdkit.Chem.rdchem.Mol
            If `file_name` is None, you must init with a rdkit Mol type.
        sanitize_level : str
            The level of sanitizer, 'low', 'medium' or 'high'.
        raise_errors : bool
            whether to raise an Exception if sanitization procedure fails.
        verbose : bool
            whether to print information in the sanitization procedure.
        **kwargs : dict
            Additional arguments for the format.
        """
        System.__init__(self)
        self.sanitizer = Sanitizer(sanitize_level, raise_errors, verbose)

        if data:
            mol = dpdata.rdkit.utils.system_data_to_mol(data)
            self.from_rdkit_mol(mol)
        if file_name:
            self.from_fmt(
                file_name, fmt, type_map=type_map, begin=begin, step=step, **kwargs
            )
        elif rdkit_mol:
            self.from_rdkit_mol(rdkit_mol)
        else:
            raise ValueError("Please specify a mol/sdf file or a rdkit Mol object")

        if type_map:
            self.apply_type_map(type_map)
        self.check_data()

    def from_fmt_obj(self, fmtobj, file_name, **kwargs):
        mol = fmtobj.from_bond_order_system(file_name, **kwargs)
        self.from_rdkit_mol(mol)
        if hasattr(fmtobj.from_bond_order_system, "post_func"):
            for post_f in fmtobj.from_bond_order_system.post_func:
                self.post_funcs.get_plugin(post_f)(self)
        return self

    def to_fmt_obj(self, fmtobj, *args, **kwargs):
        self.rdkit_mol.RemoveAllConformers()
        for ii in range(self.get_nframes()):
            conf = Conformer()
            for idx in range(self.get_natoms()):
                conf.SetAtomPosition(idx, self.data["coords"][ii][idx])
            self.rdkit_mol.AddConformer(conf, assignId=True)
        return fmtobj.to_bond_order_system(self.data, self.rdkit_mol, *args, **kwargs)

    def __str__(self):
        """A brief summary of the system."""
        ret = "Data Summary"
        ret += "\nBondOrder System"
        ret += "\n-------------------"
        ret += f"\nFrame Numbers      : {self.get_nframes()}"
        ret += f"\nAtom Numbers       : {self.get_natoms()}"
        ret += f"\nBond Numbers       : {self.get_nbonds()}"
        ret += "\nElement List       :"
        ret += "\n-------------------"
        ret += "\n" + "  ".join(map(str, self.get_atom_names()))
        ret += "\n" + "  ".join(map(str, self.get_atom_numbs()))
        return ret

    def get_nbonds(self):
        """Return the number of bonds."""
        return len(self.data["bonds"])

    def get_charge(self):
        """Return the total formal charge of the moleclue."""
        return sum(self.data["formal_charges"])

    def get_mol(self):
        """Return the rdkit.Mol object."""
        return self.rdkit_mol

    def get_bond_order(self, begin_atom_idx, end_atom_idx):
        """Return the bond order between given atoms."""
        return self.data["bond_dict"][f"{int(begin_atom_idx)}-{int(end_atom_idx)}"]

    def get_formal_charges(self):
        """Return the formal charges on each atom."""
        return self.data["formal_charges"]

    def copy(self):
        new_mol = deepcopy(self.rdkit_mol)
        self.__class__(data=deepcopy(self.data), rdkit_mol=new_mol)

    def __add__(self, other):
        raise NotImplementedError(
            "magic method '+' has not been implemented on BondOrderSystem"
        )

    #     '''
    #         magic method "+" operation
    #     '''
    #     if isinstance(other, BondOrderSystem):
    #         if dpdata.rdkit.utils.check_same_molecule(self.rdkit_mol, other.rdkit_mol):
    #             self.__class__(self, data=other.data)
    #         else:
    #             raise RuntimeError("The two systems are not of the same topology.")
    #     else:
    #         raise RuntimeError(f"Unsupported data structure: {type(other)}")

    def from_rdkit_mol(self, rdkit_mol):
        """Initialize from a rdkit.Chem.rdchem.Mol object."""
        rdkit_mol = self.sanitizer.sanitize(rdkit_mol)
        self.data = dpdata.rdkit.utils.mol_to_system_data(rdkit_mol)
        self.data["bond_dict"] = dict(
            [(f"{int(bond[0])}-{int(bond[1])}", bond[2]) for bond in self.data["bonds"]]
        )
        self.rdkit_mol = rdkit_mol
