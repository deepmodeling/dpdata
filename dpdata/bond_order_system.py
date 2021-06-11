#%%
# Bond Order System
from dpdata.system import Register, System, LabeledSystem, check_System
import rdkit.Chem
import dpdata.rdkit.utils
from dpdata.rdkit.sanitize import Sanitizer, SanitizeError
from copy import deepcopy
# import dpdata.rdkit.mol2

def check_BondOrderSystem(data):
    check_System(data)
    assert ('bonds' in data.keys())
    
class BondOrderSystem(System):
    '''
    The system with chemical bond and formal charges information

    For example, a labeled methane system named `d_example` has one molecule (5 atoms, 4 bonds) and `n_frames` frames. The bond order and formal charge information can be accessed by
        - `d_example['bonds']` : a numpy array of size 4 x 3, and
                                    the first column represents the index of begin atom,
                                    the second column represents the index of end atom, 
                                    the third columen represents the bond order:
                                        1 - single bond, 2 - double bond, 3 - triple bond, 1.5 - aromatic bond
        - `d_example['formal_charges']` : a numpy array of size 5 x 1
    '''
    def __init__(self,
                 file_name = None,
                 fmt = 'auto',
                 type_map = None,
                 begin = 0,
                 step = 1,
                 data = None,
                 rdkit_mol = None,
                 sanitize_level = "medium",
                 raise_errors = True,
                 verbose = False,
                 **kwargs):
        """
        Constructor

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
        """

        System.__init__(self)
        self.sanitizer = Sanitizer(sanitize_level, raise_errors, verbose)

        if data:
            mol = dpdata.rdkit.utils.system_data_to_mol(data)
            self.from_rdkit_mol(mol)
        if file_name:
            self.from_fmt(file_name, 
                          fmt,
                          type_map=type_map,
                          begin=begin, 
                          step=step,
                          **kwargs)
        elif rdkit_mol:
            self.from_rdkit_mol(rdkit_mol)
        else:
            raise ValueError("Please specify a mol/sdf file or a rdkit Mol object")

        if type_map:
            self.apply_type_map(type_map)

    register_from_funcs = Register()
    register_to_funcs = System.register_to_funcs + Register()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        '''
            A brief summary of the system
        '''
        ret = "Data Summary"
        ret += "\nBondOrder System"
        ret += "\n-------------------"
        ret += f"\nFrame Numbers      : {self.get_nframes()}"
        ret += f"\nAtom Numbers       : {self.get_natoms()}"
        ret += f"\nBond Numbers       : {self.get_nbonds()}"
        ret += "\nElement List       :"
        ret += "\n-------------------"
        ret += "\n"+"  ".join(map(str,self.get_atom_names()))
        ret += "\n"+"  ".join(map(str,self.get_atom_numbs()))
        return ret

    def get_nbonds(self):
        '''
            Return the number of bonds
        '''
        return len(self.data['bonds'])
    
    def get_charge(self):
        '''
            Return the total formal charge of the moleclue
        '''
        return sum(self.data['formal_charges'])
    
    def get_mol(self):
        '''
            Return the rdkit.Mol object
        '''
        return self.rdkit_mol
    
    def get_bond_order(self, begin_atom_idx, end_atom_idx):
        '''
            Return the bond order between given atoms
        '''
        return self.data['bond_dict'][f'{int(begin_atom_idx)}-{int(end_atom_idx)}']
    
    def get_formal_charges(self):
        '''
            Return the formal charges on each atom
        '''
        return self.data['formal_charges']
    
    def copy(self):
        new_mol = deepcopy(self.rdkit_mol)
        self.__class__(data=deepcopy(self.data),
                       rdkit_mol=new_mol)
    
    # def __add__(self, other):
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
        '''
            Initialize from a rdkit.Chem.rdchem.Mol object
        '''
        rdkit_mol = self.sanitizer.sanitize(rdkit_mol)
        self.data = dpdata.rdkit.utils.mol_to_system_data(rdkit_mol)
        self.data['bond_dict'] = dict([(f'{int(bond[0])}-{int(bond[1])}', bond[2]) for bond in self.data['bonds']])
        self.rdkit_mol = rdkit_mol

    @register_from_funcs.register_funcs('mol')
    def from_mol_file(self, file_name):
        mol = rdkit.Chem.MolFromMolFile(file_name, sanitize=False, removeHs=False)
        self.from_rdkit_mol(mol)

    @register_to_funcs.register_funcs("mol")
    def to_mol_file(self, file_name, frame_idx=0):
        assert (frame_idx < self.get_nframes())
        rdkit.Chem.MolToMolFile(self.rdkit_mol, file_name, confId=frame_idx)
    
    @register_from_funcs.register_funcs("sdf")
    def from_sdf_file(self, file_name):
        '''
        Note that it requires all molecules in .sdf file must be of the same topology
        '''
        mols = [m for m in rdkit.Chem.SDMolSupplier(file_name, sanitize=False, removeHs=False)]
        if len(mols) > 1:
            mol = dpdata.rdkit.utils.combine_molecules(mols)
        else:
            mol = mols[0]
        self.from_rdkit_mol(mol)
    
    @register_to_funcs.register_funcs("sdf")
    def to_sdf_file(self, file_name, frame_idx=-1):
        sdf_writer = rdkit.Chem.SDWriter(file_name)
        if frame_idx == -1:
            for ii in self.get_nframes():
                sdf_writer.write(self.rdkit_mol, confId=ii)
        else:
            assert (frame_idx < self.get_nframes())
            sdf_writer.write(self.rdkit_mol, confId=frame_idx)
        sdf_writer.close()
