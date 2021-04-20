#%%
import os
import glob
import inspect
import numpy as np
import dpdata.lammps.lmp
import dpdata.lammps.dump
import dpdata.vasp.poscar
import dpdata.vasp.xml
import dpdata.vasp.outcar
import dpdata.deepmd.raw
import dpdata.deepmd.comp
import dpdata.qe.traj
import dpdata.qe.scf
import dpdata.siesta.output
import dpdata.siesta.aiMD_output
import dpdata.md.pbc
import dpdata.gaussian.log
import dpdata.amber.md
import dpdata.cp2k.output
from dpdata.cp2k.output import Cp2kSystems
import dpdata.pwmat.movement
import dpdata.pwmat.atomconfig
import dpdata.fhi_aims.output
import dpdata.gromacs.gro
from copy import deepcopy
from monty.json import MSONable
from monty.serialization import loadfn,dumpfn
from dpdata.periodic_table import Element
from dpdata.xyz.quip_gap_xyz import QuipGapxyzSystems


class Register:
    def __init__(self):
        self.funcs = {}

    def register_funcs(self, fmt):
        def decorator(func):
            self.funcs[fmt] = func
            return func
        return decorator

    def __add__(self, other):
        self.funcs.update(other.funcs)
        return self


class System (MSONable) :
    '''
    The data System

    A data System (a concept used by `deepmd-kit <https://github.com/deepmodeling/deepmd-kit>`_)
    contains frames (e.g. produced by an MD simulation) that has the same number of atoms of the same type.
    The order of the atoms should be consistent among the frames in one System.

    For example, a water system named `d_example` has two molecules. The properties can be accessed by
        - `d_example['atom_numbs']` : [2, 4]
        - `d_example['atom_names']` : ['O', 'H']
        - `d_example['atom_types']` : [0, 1, 1, 0, 1, 1]
        - `d_example['orig']` : [0, 0, 0]
        - `d_example['cells']` : a numpy array of size nframes x 3 x 3
        - `d_example['coords']` : a numpy array of size nframes x natoms x 3

    It is noted that
        - The order of frames stored in `'atom_types'`, `'cells'` and `'coords'` should be consistent.
        - The order of atoms in **all** frames of `'atom_types'` and  `'coords'` should be consistent.

    Restrictions:
        - `d_example['orig']` is always [0, 0, 0]
        - `d_example['cells'][ii]` is always lower triangular (lammps cell tensor convention)
    '''

    def __init__ (self,
                  file_name = None,
                  fmt = 'auto',
                  type_map = None,
                  begin = 0,
                  step = 1,
                  data = None,
                  **kwargs) :
        """
        Constructor

        Parameters
        ----------
        file_name : str
            The file to load the system
        fmt : str
            Format of the file, supported formats are
                - ``auto``: infered from `file_name`'s extension
                - ``lammps/lmp``: Lammps data
                - ``lammps/dump``: Lammps dump
                - ``deepmd/raw``: deepmd-kit raw
                - ``deepmd/npy``: deepmd-kit compressed format (numpy binary)
                - ``vasp/poscar``: vasp POSCAR
                - ``qe/cp/traj``: Quantum Espresso CP trajectory files. should have: file_name+'.in' and file_name+'.pos'
                - ``siesta/output``: siesta SCF output file
                - ``siesta/aimd_output``: siesta aimd output file
                - ``pwmat/atom.config``: pwmat atom.config
        type_map : list of str
            Needed by formats lammps/lmp and lammps/dump. Maps atom type to name. The atom with type `ii` is mapped to `type_map[ii]`.
            If not provided the atom names are assigned to `'Type_1'`, `'Type_2'`, `'Type_3'`...
        begin : int
            The beginning frame when loading MD trajectory.
        step : int
            The number of skipped frames when loading MD trajectory.
        data : dict
             The raw data of System class.
        """
        self.data = {}
        self.data['atom_numbs'] = []
        self.data['atom_names'] = []
        self.data['atom_types'] = []
        self.data['orig'] = np.array([0, 0, 0])
        self.data['cells'] = []
        self.data['coords'] = []

        if data:
            check_System(data)
            self.data=data
            return
        if file_name is None :
            return
        self.from_fmt(file_name, fmt, type_map=type_map, begin= begin, step=step, **kwargs)

        if type_map is not None:
            self.apply_type_map(type_map)

    register_from_funcs = Register()
    register_to_funcs = Register()

    def from_fmt(self, file_name, fmt='auto', **kwargs):
        fmt = fmt.lower()
        if fmt == 'auto':
            fmt = os.path.basename(file_name).split('.')[-1].lower()
        from_funcs = self.register_from_funcs.funcs
        if fmt in from_funcs:
            func = from_funcs[fmt]
            args = inspect.getfullargspec(func).args
            kwargs = {kk: kwargs[kk] for kk in kwargs if kk in args}
            func(self, file_name, **kwargs)
        else :
            raise RuntimeError('unknow data format ' + fmt)

    def to(self, fmt, *args, **kwargs):
        fmt = fmt.lower()
        to_funcs = self.register_to_funcs.funcs
        if fmt in to_funcs:
            func = to_funcs[fmt]
            func_args = inspect.getfullargspec(func).args
            kwargs = {kk: kwargs[kk] for kk in kwargs if kk in func_args}
            func(self, *args, **kwargs)
        else :
            raise RuntimeError('unknow data format %s. Accepted format: %s' % (fmt, " ".join(to_funcs)))

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        ret="Data Summary"
        ret+="\nUnlabeled System"
        ret+="\n-------------------"
        ret+="\nFrame Numbers     : %d"%self.get_nframes()
        ret+="\nAtom Numbers      : %d"%self.get_natoms()
        ret+="\nElement List      :"
        ret+="\n-------------------"
        ret+="\n"+"  ".join(map(str,self.get_atom_names()))
        ret+="\n"+"  ".join(map(str,self.get_atom_numbs()))
        return ret

    def __getitem__(self, key):
        """Returns proerty stored in System by key or by idx"""
        if isinstance(key, (int, slice)):
            return self.sub_system(key)
        return self.data[key]

    def __len__(self) :
        """Returns number of frames in the system"""
        return self.get_nframes()


    def __add__(self,others) :
       """magic method "+" operation """
       self_copy=self.copy()
       if isinstance(others,System):
          other_copy=others.copy()
          self_copy.append(other_copy)
       elif isinstance(others, list):
          for ii in others:
              assert(isinstance(ii,System))
              ii_copy=ii.copy()
              self_copy.append(ii_copy)
       else:
          raise RuntimeError("Unspported data structure")
       return self.__class__.from_dict({'data':self_copy.data})


    def dump(self,filename,indent=4):
        """dump .json or .yaml file """
        dumpfn(self.as_dict(),filename,indent=indent)


    def map_atom_types(self,type_map=None):
        """
        Map the atom types of the system
        Parameters
        ----------
        type_map :
            dict :  {"H":0,"O":1}
            or list  ["H","C","O","N"]
            The map between elements and index
            if no map_dict is given, index will
            be set according to atomic number

        Returns
        -------
        new_atom_types : list
            The mapped atom types
        """
        if isinstance(type_map,dict) or type_map is None:
           pass
        elif isinstance(type_map,list):
           type_map=dict(zip(type_map,range(len(type_map))))
        else:
           raise RuntimeError("Unknown format")

        if type_map is None:
           type_map=elements_index_map(self.get_atom_names().copy(),standard=True)

        _set1=set(self.get_atom_names())
        _set2=set(list(type_map.keys()))
        assert _set1.issubset(_set2)

        atom_types_list=[]
        for name, numb  in  zip(self.get_atom_names(), self.get_atom_numbs()):
            atom_types_list.extend([name]*numb)
        new_atom_types=np.array([type_map[ii] for ii in atom_types_list],dtype=np.int)

        return new_atom_types

    @register_to_funcs.register_funcs("list")
    def to_list(self):
        """
        convert system to list, usefull for data collection
        """
        if len(self)==0:
           return []
        if len(self)==1:
           return [self]
        else:
           systems=[]
           for ii in range(len(self)):
               systems.append(self.sub_system([ii]))
           return systems

    @staticmethod
    def load(filename):
        """rebuild System obj. from .json or .yaml file """
        return loadfn(filename)

    def as_dict(self):
        """Returns data dict of System instance"""
        d={"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "data": self.data
          }
        return d


    def get_atom_names(self):
        """Returns name of atoms """
        return  self.data['atom_names']


    def get_atom_types(self):
        """Returns type of atoms """
        return self.data['atom_types']


    def get_atom_numbs(self):
        """Returns number of atoms """
        return self.data['atom_numbs']


    def get_nframes(self) :
        """Returns number of frames in the system"""
        return len(self.data['cells'])


    def get_natoms(self) :
        """Returns total number of atoms in the system"""
        return len(self.data['atom_types'])


    def copy(self):
        """Returns a copy of the system.  """
        return self.__class__.from_dict({'data':deepcopy(self.data)})


    def sub_system(self, f_idx) :
        """
        Construct a subsystem from the system

        Parameters
        ----------
        f_idx : int or index
            Which frame to use in the subsystem

        Returns
        -------
        sub_system : System
            The subsystem
        """
        tmp = System()
        for ii in ['atom_numbs', 'atom_names', 'atom_types', 'orig'] :
            tmp.data[ii] = self.data[ii]
        
        tmp.data['cells'] = self.data['cells'][f_idx].reshape(-1, 3, 3)
        tmp.data['coords'] = self.data['coords'][f_idx].reshape(-1, self.data['coords'].shape[1], 3)
        tmp.data['nopbc'] = self.nopbc
        
        return tmp


    def append(self, system) :
        """
        Append a system to this system

        Parameters
        ----------
        system : System
            The system to append
        """
        if not len(system.data['atom_numbs']):
            # skip if the system to append is non-converged
            return False
        elif not len(self.data['atom_numbs']):
            # this system is non-converged but the system to append is converged
            self.data = system.data
            return False
        if system.uniq_formula != self.uniq_formula:
            raise RuntimeError('systems with inconsistent formula could not be append: %s v.s. %s' % (self.uniq_formula, system.uniq_formula))
        if system.data['atom_names'] != self.data['atom_names']:
            # allow to append a system with different atom_names order
            system.sort_atom_names()
            self.sort_atom_names()
        if (system.data['atom_types'] != self.data['atom_types']).any():
            # allow to append a system with different atom_types order
            system.sort_atom_types()
            self.sort_atom_types()
        for ii in ['atom_numbs', 'atom_names'] :
            assert(system.data[ii] == self.data[ii])
        for ii in ['atom_types','orig'] :
            eq = [v1==v2 for v1,v2 in zip(system.data[ii], self.data[ii])]
            assert(all(eq))
        for ii in ['coords', 'cells'] :
            self.data[ii] = np.concatenate((self.data[ii], system[ii]), axis = 0)
        if self.nopbc and not system.nopbc:
            # appended system uses PBC, cancel nopbc
            self.data['nopbc'] = False
        return True

    def sort_atom_names(self, type_map=None):
        """
        Sort atom_names of the system and reorder atom_numbs and atom_types accoarding
        to atom_names. If type_map is not given, atom_names will be sorted by
        alphabetical order. If type_map is given, atom_names will be type_map.

        Parameters
        ----------
        type_map : list
            type_map
        """
        if type_map is not None:
            # assign atom_names index to the specify order
            # atom_names must be a subset of type_map
            assert (set(self.data['atom_names']).issubset(set(type_map)))
            # for the condition that type_map is a proper superset of atom_names
            # new_atoms = set(type_map) - set(self.data["atom_names"])
            new_atoms = [e for e in type_map if e not in self.data["atom_names"]]
            if new_atoms:
                self.add_atom_names(new_atoms)
            # index that will sort an array by type_map
            # a[as[a]] == b[as[b]]  as == argsort
            # as[as[b]] == as^{-1}[b]
            # a[as[a][as[as[b]]]] = b[as[b][as^{-1}[b]]] = b[id]
            idx = np.argsort(self.data['atom_names'])[np.argsort(np.argsort(type_map))]
        else:
            # index that will sort an array by alphabetical order
            idx = np.argsort(self.data['atom_names'])
        # sort atom_names, atom_numbs, atom_types by idx
        self.data['atom_names'] = list(np.array(self.data['atom_names'])[idx])
        self.data['atom_numbs'] = list(np.array(self.data['atom_numbs'])[idx])
        self.data['atom_types'] = np.argsort(idx)[self.data['atom_types']]

    def check_type_map(self, type_map):
        """
        Assign atom_names to type_map if type_map is given and different from
        atom_names.

        Parameters
        ----------
        type_map : list
            type_map
        """
        if type_map is not None and type_map != self.data['atom_names']:
            self.sort_atom_names(type_map=type_map)

    def apply_type_map(self, type_map) :
        if type_map is not None and type(type_map) is list:
            self.check_type_map(type_map)
        else:
            raise RuntimeError('invalid type map, cannot be applied')

    def sort_atom_types(self):
        idx = np.argsort(self.data['atom_types'])
        self.data['atom_types'] = self.data['atom_types'][idx]
        self.data['coords'] = self.data['coords'][:, idx]
        return idx

    @property
    def formula(self):
        """
        Return the formula of this system, like C3H5O2
        """
        return ''.join(["{}{}".format(symbol,numb) for symbol,numb in
            zip(self.data['atom_names'], self.data['atom_numbs'])])

    @property
    def uniq_formula(self):
        """
        Return the uniq_formula of this system.
        The uniq_formula sort the elements in formula by names.
        Systems with the same uniq_formula can be append together.
        """
        return ''.join(["{}{}".format(symbol,numb) for symbol,numb in sorted(
            zip(self.data['atom_names'], self.data['atom_numbs']))])


    def extend(self, systems):
        """
        Extend a system list to this system

        Parameters
        ----------
        systems : [System1, System2, System3 ]
            The list to extend
        """

        for system in systems:
            self.append(system.copy())


    def apply_pbc(self) :
        """
        Append periodic boundary condition
        """
        ncoord = dpdata.md.pbc.dir_coord(self.data['coords'], self.data['cells'])
        ncoord = ncoord % 1
        self.data['coords'] = np.matmul(ncoord, self.data['cells'])


    def remove_pbc(self, protect_layer = 9):
        """
        This method does NOT delete the definition of the cells, it
        (1) revises the cell to a cubic cell and ensures that the cell
        boundary to any atom in the system is no less than `protect_layer`
        (2) translates the system such that the center-of-geometry of the system
        locates at the center of the cell.

        Parameters
        ----------
        protect_layer : the protect layer between the atoms and the cell
                        boundary
        """
        nframes = self.get_nframes()
        natoms = self.get_natoms()
        assert(protect_layer >= 0), "the protect_layer should be no less than 0"
        for ff in range(nframes):
            tmpcoord = self.data['coords'][ff]
            cog = np.average(tmpcoord, axis = 0)
            dist = tmpcoord - np.tile(cog, [natoms, 1])
            max_dist = np.max(np.linalg.norm(dist, axis = 1))
            h_cell_size = max_dist + protect_layer
            cell_size = h_cell_size * 2
            shift = np.array([1,1,1]) * h_cell_size - cog
            self.data['coords'][ff] = self.data['coords'][ff] + np.tile(shift, [natoms, 1])
            self.data['cells'][ff] = cell_size * np.eye(3)


    @register_from_funcs.register_funcs("lmp")
    @register_from_funcs.register_funcs("lammps/lmp")
    def from_lammps_lmp (self, file_name, type_map = None) :
        with open(file_name) as fp:
            lines = [line.rstrip('\n') for line in fp]
            self.data = dpdata.lammps.lmp.to_system_data(lines, type_map)
        self._shift_orig_zero()

    @register_to_funcs.register_funcs("pymatgen/structure")
    def to_pymatgen_structure(self):
        '''
        convert System to Pymatgen Structure obj

        '''
        structures=[]
        try:
           from pymatgen import Structure
        except:
           raise ImportError('No module pymatgen.Structure')

        for system in self.to_list():
            species=[]
            for name,numb in zip(system.data['atom_names'],system.data['atom_numbs']):
                species.extend([name]*numb)
            structure=Structure(system.data['cells'][0],species,system.data['coords'][0],coords_are_cartesian=True)
            structures.append(structure)
        return structures


    @register_to_funcs.register_funcs("ase/structure")
    def to_ase_structure(self):
        '''
        convert System to ASE Atom obj

        '''
        from ase import Atoms
        
        structures=[]

        for system in self.to_list():
            species=[system.data['atom_names'][tt] for tt in system.data['atom_types']]
            structure=Atoms(symbols=species,positions=system.data['coords'][0],pbc=True,cell=system.data['cells'][0])
            structures.append(structure)

        return structures

    @register_to_funcs.register_funcs("lammps/lmp")
    def to_lammps_lmp(self, file_name, frame_idx = 0) :
        """
        Dump the system in lammps data format

        Parameters
        ----------
        file_name : str
            The output file name
        frame_idx : int
            The index of the frame to dump
        """
        assert(frame_idx < len(self.data['coords']))
        w_str = dpdata.lammps.lmp.from_system_data(self.data, frame_idx)
        with open(file_name, 'w') as fp:
            fp.write(w_str)

    @register_from_funcs.register_funcs('dump')
    @register_from_funcs.register_funcs('lammps/dump')
    def from_lammps_dump (self,
                          file_name,
                          type_map = None,
                          begin = 0,
                          step = 1) :
        lines = dpdata.lammps.dump.load_file(file_name, begin = begin, step = step)
        self.data = dpdata.lammps.dump.system_data(lines, type_map)
        self._shift_orig_zero()

    @register_from_funcs.register_funcs('poscar')
    @register_from_funcs.register_funcs('contcar')
    @register_from_funcs.register_funcs('vasp/poscar')
    @register_from_funcs.register_funcs('vasp/contcar')
    def from_vasp_poscar(self, file_name) :
        with open(file_name) as fp:
            lines = [line.rstrip('\n') for line in fp]
            self.data = dpdata.vasp.poscar.to_system_data(lines)
        self.rot_lower_triangular()

    @register_to_funcs.register_funcs("vasp/string")
    def to_vasp_string(self, frame_idx=0):
        """
        Dump the system in vasp POSCAR format string

        Parameters
        ----------
        frame_idx : int
            The index of the frame to dump
        """
        assert(frame_idx < len(self.data['coords']))
        w_str = dpdata.vasp.poscar.from_system_data(self.data, frame_idx)
        return w_str

    @register_to_funcs.register_funcs("vasp/poscar")
    def to_vasp_poscar(self, file_name, frame_idx = 0) :
        """
        Dump the system in vasp POSCAR format

        Parameters
        ----------
        file_name : str
            The output file name
        frame_idx : int
            The index of the frame to dump
        """
        w_str=self.to_vasp_string( frame_idx= frame_idx )
        with open(file_name, 'w') as fp:
            fp.write(w_str)

    @register_from_funcs.register_funcs('qe/cp/traj')
    def from_qe_cp_traj(self,
                        prefix,
                        begin = 0,
                        step = 1) :
        self.data, _ = dpdata.qe.traj.to_system_data(prefix + '.in', prefix, begin = begin, step = step)
        self.data['coords'] \
            = dpdata.md.pbc.apply_pbc(self.data['coords'],
                                      self.data['cells'],
            )
        self.rot_lower_triangular()

    @register_from_funcs.register_funcs('deepmd/npy')
    def from_deepmd_comp(self, folder, type_map = None) :
        self.data = dpdata.deepmd.comp.to_system_data(folder, type_map = type_map, labels = False)

    @register_from_funcs.register_funcs('deepmd')
    @register_from_funcs.register_funcs('deepmd/raw')
    def from_deepmd_raw(self, folder, type_map = None) :
        tmp_data = dpdata.deepmd.raw.to_system_data(folder, type_map = type_map, labels = False)
        if tmp_data is not None :
            self.data = tmp_data

    @register_from_funcs.register_funcs("gro")
    @register_from_funcs.register_funcs("gromacs/gro")
    def from_gromacs_gro(self, file_name) :
        """
        Load gromacs .gro file

        Parameters
        ----------
        file_name : str
            The input file name
        """
        self.data = dpdata.gromacs.gro.file_to_system_data(file_name)

    @register_to_funcs.register_funcs("deepmd/npy")
    def to_deepmd_npy(self, folder, set_size = 5000, prec=np.float32) :
        """
        Dump the system in deepmd compressed format (numpy binary) to `folder`.

        The frames are firstly split to sets, then dumped to seperated subfolders named as `folder/set.000`, `folder/set.001`, ....

        Each set contains `set_size` frames.
        The last set may have less frames than `set_size`.

        Parameters
        ----------
        folder : str
            The output folder
        set_size : int
            The size of each set.
        prec: {numpy.float32, numpy.float64}
            The floating point precision of the compressed data
        """
        dpdata.deepmd.comp.dump(folder, self.data,
                                set_size = set_size,
                                comp_prec = prec)

    @register_to_funcs.register_funcs("deepmd/raw")
    def to_deepmd_raw(self, folder) :
        """
        Dump the system in deepmd raw format to `folder`
        """
        dpdata.deepmd.raw.dump(folder, self.data)

    @register_from_funcs.register_funcs('siesta/output')
    def from_siesta_output(self, fname):
        self.data['atom_names'], \
        self.data['atom_numbs'], \
        self.data['atom_types'], \
        self.data['cells'], \
        self.data['coords'], \
        _e, _f, _v \
            = dpdata.siesta.output.obtain_frame(fname)
        # self.rot_lower_triangular()

    @register_from_funcs.register_funcs('siesta/aimd_output')
    def from_siesta_aiMD_output(self, fname):
        self.data['atom_names'], \
        self.data['atom_numbs'], \
        self.data['atom_types'], \
        self.data['cells'], \
        self.data['coords'], \
        _e, _f, _v \
            = dpdata.siesta.aiMD_output.get_aiMD_frame(fname)
    @register_from_funcs.register_funcs('atom.config')
    @register_from_funcs.register_funcs('final.config')
    @register_from_funcs.register_funcs('pwmat/atom.config')
    @register_from_funcs.register_funcs('pwmat/final.config')
    def from_pwmat_atomconfig(self, file_name) :
        with open(file_name) as fp:
            lines = [line.rstrip('\n') for line in fp]
            self.data = dpdata.pwmat.atomconfig.to_system_data(lines)
        self.rot_lower_triangular()

    @register_to_funcs.register_funcs("pwmat/atom.config")
    def to_pwmat_atomconfig(self, file_name, frame_idx = 0) :
        """
        Dump the system in pwmat atom.config format

        Parameters
        ----------
        file_name : str
            The output file name
        frame_idx : int
            The index of the frame to dump
        """
        assert(frame_idx < len(self.data['coords']))
        w_str = dpdata.pwmat.atomconfig.from_system_data(self.data, frame_idx)
        with open(file_name, 'w') as fp:
            fp.write(w_str)


    def affine_map(self, trans, f_idx = 0) :
        assert(np.linalg.det(trans) != 0)
        self.data['cells'][f_idx] = np.matmul(self.data['cells'][f_idx], trans)
        self.data['coords'][f_idx] = np.matmul(self.data['coords'][f_idx], trans)


    def _shift_orig_zero(self) :
        for ff in self.data['coords'] :
            for ii in ff :
                ii = ii - self.data['orig']
        self.data['orig'] = self.data['orig'] - self.data['orig']
        assert((np.zeros([3]) == self.data['orig']).all())


    def rot_lower_triangular(self) :
        for ii in range(self.get_nframes()) :
            self.rot_frame_lower_triangular(ii)


    def rot_frame_lower_triangular(self, f_idx = 0) :
        qq, rr = np.linalg.qr(self.data['cells'][f_idx].T)
        if np.linalg.det(qq) < 0 :
            qq = -qq
            rr = -rr
        self.affine_map(qq, f_idx = f_idx)
        rot = np.eye(3)
        if self.data['cells'][f_idx][0][0] < 0 :
            rot[0][0] = -1
        if self.data['cells'][f_idx][1][1] < 0 :
            rot[1][1] = -1
        if self.data['cells'][f_idx][2][2] < 0 :
            rot[2][2] = -1
        assert(np.linalg.det(rot) == 1)
        self.affine_map(rot, f_idx = f_idx)
        return np.matmul(qq, rot)


    def add_atom_names(self, atom_names):
        """
        Add atom_names that do not exist.
        """
        self.data['atom_names'].extend(atom_names)
        self.data['atom_numbs'].extend([0 for _ in atom_names])

    def replicate(self, ncopy):
        """
        Replicate the each frame  in the system in 3 dimensions.
        Each frame in the system will become a supercell.

        Parameters
        ----------
        ncopy :
            list: [4,2,3]
            or tuple: (4,2,3,)
            make `ncopy[0]` copys in x dimensions,
            make `ncopy[1]` copys in y dimensions,
            make `ncopy[2]` copys in z dimensions.

        Returns
        -------
        tmp : System
            The system after replication.
        """
        if len(ncopy) !=3:
            raise RuntimeError('ncopy must be a list or tuple with 3 int')
        for ii in ncopy:
            if type(ii) is not int:
                raise RuntimeError('ncopy must be a list or tuple must with 3 int')

        tmp = System()
        nframes = self.get_nframes()
        data = self.data
        tmp.data['atom_names'] = list(np.copy(data['atom_names']))
        tmp.data['atom_numbs'] = list(np.array(np.copy(data['atom_numbs'])) * np.prod(ncopy))
        tmp.data['atom_types'] = np.sort(np.tile(np.copy(data['atom_types']),np.prod(ncopy)))
        tmp.data['cells'] = np.copy(data['cells'])
        for ii in range(3):
            tmp.data['cells'][:,ii,:] *= ncopy[ii]
        tmp.data['coords'] = np.tile(np.copy(data['coords']),tuple(ncopy)+(1,1,1))

        for xx in range(ncopy[0]):
            for yy in range(ncopy[1]):
                for zz in range(ncopy[2]):
                    tmp.data['coords'][xx,yy,zz,:,:,:] += xx * np.reshape(data['cells'][:,0,:], [-1,1,3])\
                                                + yy * np.reshape(data['cells'][:,1,:], [-1,1,3])\
                                                + zz * np.reshape(data['cells'][:,2,:], [-1,1,3])
        tmp.data['coords'] = np.reshape(np.transpose(tmp.data['coords'], [3,4,0,1,2,5]), (nframes, -1 , 3))
        return tmp

    def replace(self, initial_atom_type, end_atom_type, replace_num):
        if type(self) is not dpdata.System:
            raise RuntimeError('Must use method replace() of the instance of class dpdata.System')
        if type(replace_num) is not int:
            raise ValueError("replace_num must be a integer. Now is {replace_num}".format(replace_num=replace_num))
        if replace_num <= 0:
            raise ValueError("replace_num must be larger than 0.Now is {replace_num}".format(replace_num=replace_num))

        try:
            initial_atom_index = self.data['atom_names'].index(initial_atom_type)
        except ValueError as e:
            raise ValueError("atom_type  {initial_atom_type}   not in {atom_names}"
                    .format(initial_atom_type=initial_atom_type, atom_names=self.data['atom_names']))
        max_replace_num = self.data['atom_numbs'][initial_atom_index]

        if replace_num > max_replace_num:
            raise RuntimeError("not enough {initial_atom_type} atom, only {max_replace_num} available, less than {replace_num}.Please check."
                    .format(initial_atom_type=initial_atom_type,max_replace_num=max_replace_num, replace_num=replace_num))

        may_replace_indices = [i for i, x in enumerate(self.data['atom_types']) if x == initial_atom_index]
        to_replace_indices = np.random.choice(may_replace_indices, size=replace_num, replace=False)

        if  end_atom_type not in self.data['atom_names']: 
            self.data['atom_names'].append(end_atom_type)
            self.data['atom_numbs'].append(0)

        end_atom_index = self.data['atom_names'].index(end_atom_type)
        for ii in to_replace_indices:
            self.data['atom_types'][ii] = end_atom_index
        self.data['atom_numbs'][initial_atom_index] -= replace_num
        self.data['atom_numbs'][end_atom_index] += replace_num 
        self.sort_atom_types()
        

    def perturb(self,
        pert_num,
        cell_pert_fraction,
        atom_pert_distance,
        atom_pert_style='normal'):
        """
        Perturb each frame in the system randomly.
        The cell will be deformed randomly, and atoms will be displaced by a random distance in random direction.

        Parameters
        ----------
        pert_num : int
            Each frame in the system will make `pert_num` copies,
            and all the copies will be perturbed.
            That means the system to be returned will contain `pert_num` * frame_num of the input system.
        cell_pert_fraction : float
            A fraction determines how much (relatively) will cell deform.
            The cell of each frame is deformed by a symmetric matrix perturbed from identity.
            The perturbation to the diagonal part is subject to a uniform distribution in [-cell_pert_fraction, cell_pert_fraction),
            and the perturbation to the off-diagonal part is subject to a uniform distribution in [-0.5*cell_pert_fraction, 0.5*cell_pert_fraction).
        atom_pert_distance: float
            unit: Angstrom. A distance determines how far atoms will move.
            Atoms will move about `atom_pert_distance` in random direction.
            The distribution of the distance atoms move is determined by atom_pert_style
        atom_pert_style : str
            Determines the distribution of the distance atoms move is subject to.
            Avaliable options are
                - `'normal'`: the `distance` will be object to `chi-square distribution with 3 degrees of freedom` after normalization.
                    The mean value of the distance is `atom_pert_fraction*side_length`
                - `'uniform'`: will generate uniformly random points in a 3D-balls with radius as `atom_pert_distance`.
                    These points are treated as vector used by atoms to move.
                    Obviously, the max length of the distance atoms move is `atom_pert_distance`.
                - `'const'`: The distance atoms move will be a constant `atom_pert_distance`.

        Returns
        -------
        perturbed_system : System
            The perturbed_system. It contains `pert_num` * frame_num of the input system frames.
        """
        perturbed_system = System()
        nframes = self.get_nframes()
        for ii in range(nframes):
            for jj in range(pert_num):
                tmp_system = self[ii].copy()
                cell_perturb_matrix = get_cell_perturb_matrix(cell_pert_fraction)
                tmp_system.data['cells'][0] = np.matmul(tmp_system.data['cells'][0],cell_perturb_matrix)
                tmp_system.data['coords'][0] = np.matmul(tmp_system.data['coords'][0],cell_perturb_matrix)
                for kk in range(len(tmp_system.data['coords'][0])):
                    atom_perturb_vector = get_atom_perturb_vector(atom_pert_distance, atom_pert_style)
                    tmp_system.data['coords'][0][kk] += atom_perturb_vector
                tmp_system.rot_lower_triangular()
                perturbed_system.append(tmp_system)
        return perturbed_system

    @property
    def nopbc(self):
        if self.data.get("nopbc", False):
            return True
        return False

    def shuffle(self):
        """Shuffle frames randomly."""
        idx = np.random.permutation(self.get_nframes())
        for ii in ['cells', 'coords']:
            self.data[ii] = self.data[ii][idx]
        return idx

    def predict(self, dp):
        """
        Predict energies and forces by deepmd-kit.

        Parameters
        ----------
        dp : deepmd.DeepPot or str
            The deepmd-kit potential class or the filename of the model.

        Returns
        -------
        labeled_sys LabeledSystem
            The labeled system.
        """
        import deepmd.DeepPot as DeepPot
        if not isinstance(dp, DeepPot):
            dp = DeepPot(dp)
        type_map = dp.get_type_map()
        ori_sys = self.copy()
        ori_sys.sort_atom_names(type_map=type_map)
        atype = ori_sys['atom_types']

        labeled_sys = LabeledSystem()

        for ss in self:
            coord = ss['coords'].reshape((-1,1))
            if not ss.nopbc:
                cell = ss['cells'].reshape((-1,1))
            else:
                cell = None
            e, f, v = dp.eval(coord, cell, atype)
            data = ss.data
            data['energies'] = e.reshape((1, 1))
            data['forces'] = f.reshape((1, -1, 3))
            data['virials'] = v.reshape((1, 3, 3))
            this_sys = LabeledSystem.from_dict({'data': data})
            labeled_sys.append(this_sys)
        return labeled_sys

def get_cell_perturb_matrix(cell_pert_fraction):
    if cell_pert_fraction<0:
        raise RuntimeError('cell_pert_fraction can not be negative')
    e0 = np.random.rand(6)
    e = e0 * 2 *cell_pert_fraction - cell_pert_fraction
    cell_pert_matrix = np.array(
        [[1+e[0], 0.5 * e[5], 0.5 * e[4]],
         [0.5 * e[5], 1+e[1], 0.5 * e[3]],
         [0.5 * e[4], 0.5 * e[3], 1+e[2]]]
    )
    return cell_pert_matrix

def get_atom_perturb_vector(atom_pert_distance, atom_pert_style='normal'):
    random_vector = None
    if atom_pert_distance < 0:
        raise RuntimeError('atom_pert_distance can not be negative')

    if atom_pert_style == 'normal':
        e = np.random.randn(3)
        random_vector=(atom_pert_distance/np.sqrt(3))*e
    elif atom_pert_style == 'uniform':
        e = np.random.randn(3)
        while np.linalg.norm(e) < 0.1:
            e = np.random.randn(3)
        random_unit_vector = e/np.linalg.norm(e)
        v0 = np.random.rand(1)
        v = np.power(v0,1/3)
        random_vector = atom_pert_distance*v*random_unit_vector
    elif atom_pert_style == 'const' :
        e = np.random.randn(3)
        while np.linalg.norm(e) < 0.1:
            e = np.random.randn(3)
        random_unit_vector = e/np.linalg.norm(e)
        random_vector = atom_pert_distance*random_unit_vector
    else:
        raise RuntimeError('unsupported options atom_pert_style={}'.format(atom_pert_style))
    return random_vector

class LabeledSystem (System):
    '''
    The labeled data System

    For example, a labeled water system named `d_example` has two molecules (6 atoms) and `nframes` frames. The labels can be accessed by
        - `d_example['energies']` : a numpy array of size nframes
        - `d_example['forces']` : a numpy array of size nframes x 6 x 3
        - `d_example['virials']` : optional, a numpy array of size nframes x 3 x 3

    It is noted that
        - The order of frames stored in `'energies'`, `'forces'` and `'virials'` should be consistent with `'atom_types'`, `'cells'` and `'coords'`.
        - The order of atoms in **every** frame of `'forces'` should be consistent with `'coords'` and `'atom_types'`.
    '''

    def __init__ (self,
                  file_name = None,
                  fmt = 'auto',
                  type_map = None,
                  begin = 0,
                  step = 1,
                  data=None,
                  **kwargs) :
        """
        Constructor

        Parameters
        ----------
        file_name : str
            The file to load the system
        fmt : str
            Format of the file, supported formats are
                - ``auto``: infered from `file_name`'s extension
                - ``vasp/xml``: vasp xml
                - ``vasp/outcar``: vasp OUTCAR
                - ``deepmd/raw``: deepmd-kit raw
                - ``deepmd/npy``: deepmd-kit compressed format (numpy binary)
                - ``qe/cp/traj``: Quantum Espresso CP trajectory files. should have: file_name+'.in', file_name+'.pos', file_name+'.evp' and file_name+'.for'
                - ``qe/pw/scf``: Quantum Espresso PW single point calculations. Both input and output files are required. If file_name is a string, it denotes the output file name. Input file name is obtained by replacing 'out' by 'in' from file_name. Or file_name is a list, with the first element being the input file name and the second element being the output filename.
                - ``siesta/output``: siesta SCF output file
                - ``siesta/aimd_output``: siesta aimd output file
                - ``gaussian/log``: gaussian logs
                - ``gaussian/md``: gaussian ab initio molecular dynamics
                - ``cp2k/output``: cp2k output file
                - ``cp2k/aimd_output``: cp2k aimd output  dir(contains *pos*.xyz and *.log file)
                - ``cp2k/restart_aimd_output``: cp2k restart aimd output  dir(contains *pos*.xyz and *.log file)
                - ``pwmat/movement``: pwmat md output file
                - ``pwmat/out.mlmd``: pwmat scf output file

        type_map : list of str
            Needed by formats deepmd/raw and deepmd/npy. Maps atom type to name. The atom with type `ii` is mapped to `type_map[ii]`.
            If not provided the atom names are assigned to `'Type_1'`, `'Type_2'`, `'Type_3'`...
        begin : int
            The beginning frame when loading MD trajectory.
        step : int
            The number of skipped frames when loading MD trajectory.
        """

        System.__init__(self)

        if data:
           check_LabeledSystem(data)
           self.data=data
           return
        if file_name is None :
            return
        self.from_fmt(file_name, fmt, type_map=type_map, begin= begin, step=step, **kwargs)
        if type_map is not None:
            self.apply_type_map(type_map)

    register_from_funcs = Register()
    register_to_funcs = System.register_to_funcs + Register()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        ret="Data Summary"
        ret+="\nLabeled System"
        ret+="\n-------------------"
        ret+="\nFrame Numbers      : %d"%self.get_nframes()
        ret+="\nAtom Numbers       : %d"%self.get_natoms()
        status= "Yes" if self.has_virial() else "No"
        ret+="\nIncluding Virials  : %s"% status
        ret+="\nElement List       :"
        ret+="\n-------------------"
        ret+="\n"+"  ".join(map(str,self.get_atom_names()))
        ret+="\n"+"  ".join(map(str,self.get_atom_numbs()))
        return ret

    def __add__(self,others) :
       """magic method "+" operation """
       self_copy=self.copy()
       if isinstance(others,LabeledSystem):
          other_copy=others.copy()
          self_copy.append(other_copy)
       elif isinstance(others, list):
          for ii in others:
              assert(isinstance(ii,LabeledSystem))
              ii_copy=ii.copy()
              self_copy.append(ii_copy)
       else:
          raise RuntimeError("Unspported data structure")
       return self.__class__.from_dict({'data':self_copy.data})

    def has_virial(self) :
        # return ('virials' in self.data) and (len(self.data['virials']) > 0)
        return ('virials' in self.data)

    @register_from_funcs.register_funcs('cp2k/aimd_output')
    def from_cp2k_aimd_output(self, file_dir, restart=False):
        xyz_file=sorted(glob.glob("{}/*pos*.xyz".format(file_dir)))[0]
        log_file=sorted(glob.glob("{}/*.log".format(file_dir)))[0]
        for info_dict in Cp2kSystems(log_file, xyz_file, restart):
            l = LabeledSystem(data=info_dict)
            self.append(l)

    @register_from_funcs.register_funcs('cp2k/restart_aimd_output')
    def from_cp2k_aimd_output(self, file_dir, restart=True):
        xyz_file = sorted(glob.glob("{}/*pos*.xyz".format(file_dir)))[0]
        log_file = sorted(glob.glob("{}/*.log".format(file_dir)))[0]
        for info_dict in Cp2kSystems(log_file, xyz_file, restart):
            l = LabeledSystem(data=info_dict)
            self.append(l)

    @register_from_funcs.register_funcs('fhi_aims/md')
    def from_fhi_aims_output(self, file_name, md=True, begin=0, step =1):
        self.data['atom_names'], \
            self.data['atom_numbs'], \
            self.data['atom_types'], \
            self.data['cells'], \
            self.data['coords'], \
            self.data['energies'], \
            self.data['forces'], \
            tmp_virial, \
            = dpdata.fhi_aims.output.get_frames(file_name, md = md, begin = begin, step = step)
        if tmp_virial is not None :
            self.data['virials'] = tmp_virial

    @register_from_funcs.register_funcs('fhi_aims/scf')
    def from_fhi_aims_output(self, file_name ):
        self.data['atom_names'], \
            self.data['atom_numbs'], \
            self.data['atom_types'], \
            self.data['cells'], \
            self.data['coords'], \
            self.data['energies'], \
            self.data['forces'], \
            tmp_virial, \
            = dpdata.fhi_aims.output.get_frames(file_name, md = False, begin = 0, step = 1)
        if tmp_virial is not None :
            self.data['virials'] = tmp_virial

    @register_from_funcs.register_funcs('xml')
    @register_from_funcs.register_funcs('vasp/xml')
    def from_vasp_xml(self, file_name, begin = 0, step = 1) :
        self.data['atom_names'], \
            self.data['atom_types'], \
            self.data['cells'], \
            self.data['coords'], \
            self.data['energies'], \
            self.data['forces'], \
            self.data['virials'], \
            = dpdata.vasp.xml.analyze(file_name, type_idx_zero = True, begin = begin, step = step)
        self.data['atom_numbs'] = []
        for ii in range(len(self.data['atom_names'])) :
            self.data['atom_numbs'].append(sum(self.data['atom_types'] == ii))
        # the vasp xml assumes the direct coordinates
        # apply the transform to the cartesan coordinates
        for ii in range(self.get_nframes()) :
            self.data['coords'][ii] = np.matmul(self.data['coords'][ii], self.data['cells'][ii])
        # scale virial to the unit of eV
        v_pref = 1 * 1e3 / 1.602176621e6
        for ii in range (self.get_nframes()) :
            vol = np.linalg.det(np.reshape(self.data['cells'][ii], [3,3]))
            self.data['virials'][ii] *= v_pref * vol
        # rotate the system to lammps convention
        self.rot_lower_triangular()

    @register_from_funcs.register_funcs('outcar')
    @register_from_funcs.register_funcs('vasp/outcar')
    def from_vasp_outcar(self, file_name, begin = 0, step = 1) :
        # with open(file_name) as fp:
        #     lines = [line.rstrip('\n') for line in fp]
        self.data['atom_names'], \
            self.data['atom_numbs'], \
            self.data['atom_types'], \
            self.data['cells'], \
            self.data['coords'], \
            self.data['energies'], \
            self.data['forces'], \
            tmp_virial, \
            = dpdata.vasp.outcar.get_frames(file_name, begin = begin, step = step)
        if tmp_virial is not None :
            self.data['virials'] = tmp_virial
        # scale virial to the unit of eV
        if 'virials' in self.data :
            v_pref = 1 * 1e3 / 1.602176621e6
            for ii in range (self.get_nframes()) :
                vol = np.linalg.det(np.reshape(self.data['cells'][ii], [3,3]))
                self.data['virials'][ii] *= v_pref * vol
        # rotate the system to lammps convention
        self.rot_lower_triangular()


    def affine_map_fv(self, trans, f_idx) :
        assert(np.linalg.det(trans) != 0)
        self.data['forces'][f_idx] = np.matmul(self.data['forces'][f_idx], trans)
        if self.has_virial():
            self.data['virials'][f_idx] = np.matmul(trans.T, np.matmul(self.data['virials'][f_idx], trans))


    def rot_lower_triangular(self) :
        for ii in range(self.get_nframes()) :
            self.rot_frame_lower_triangular(ii)


    def rot_frame_lower_triangular(self, f_idx = 0) :
        trans = System.rot_frame_lower_triangular(self, f_idx = f_idx)
        self.affine_map_fv(trans, f_idx = f_idx)
        return trans

    @register_from_funcs.register_funcs('deepmd/npy')
    def from_deepmd_comp(self, folder, type_map = None) :
        self.data = dpdata.deepmd.comp.to_system_data(folder, type_map = type_map, labels = True)

    @register_from_funcs.register_funcs('deepmd')
    @register_from_funcs.register_funcs('deepmd/raw')
    def from_deepmd_raw(self, folder, type_map = None) :
        tmp_data = dpdata.deepmd.raw.to_system_data(folder, type_map = type_map, labels = True)
        if tmp_data is not None :
            self.data = tmp_data

    @register_from_funcs.register_funcs('qe/cp/traj')
    def from_qe_cp_traj(self, prefix, begin = 0, step = 1) :
        self.data, cs = dpdata.qe.traj.to_system_data(prefix + '.in', prefix, begin = begin, step = step)
        self.data['coords'] \
            = dpdata.md.pbc.apply_pbc(self.data['coords'],
                                      self.data['cells'],
            )
        self.data['energies'], self.data['forces'], es \
            = dpdata.qe.traj.to_system_label(prefix + '.in', prefix, begin = begin, step = step)
        assert(cs == es), "the step key between files are not consistent"
        self.rot_lower_triangular()

    @register_from_funcs.register_funcs('qe/pw/scf')
    def from_qe_pw_scf(self, file_name) :
        self.data['atom_names'], \
            self.data['atom_numbs'], \
            self.data['atom_types'], \
            self.data['cells'], \
            self.data['coords'], \
            self.data['energies'], \
            self.data['forces'], \
            self.data['virials'], \
            = dpdata.qe.scf.get_frame(file_name)
        self.rot_lower_triangular()

    @register_from_funcs.register_funcs('siesta/output')
    def from_siesta_output(self, file_name) :
        self.data['atom_names'], \
        self.data['atom_numbs'], \
        self.data['atom_types'], \
        self.data['cells'], \
        self.data['coords'], \
        self.data['energies'], \
        self.data['forces'], \
        self.data['virials'] \
            = dpdata.siesta.output.obtain_frame(file_name)
        # self.rot_lower_triangular()

    @register_from_funcs.register_funcs('siesta/aimd_output')
    def from_siesta_aiMD_output(self, file_name):
        self.data['atom_names'], \
        self.data['atom_numbs'], \
        self.data['atom_types'], \
        self.data['cells'], \
        self.data['coords'], \
        self.data['energies'], \
        self.data['forces'], \
        self.data['virials'] \
            = dpdata.siesta.aiMD_output.get_aiMD_frame(file_name)

    @register_from_funcs.register_funcs('gaussian/log')
    def from_gaussian_log(self, file_name, md=False):
        try:
            self.data = dpdata.gaussian.log.to_system_data(file_name, md=md)
        except AssertionError:
            self.data['energies'], self.data['forces']= [], []
            self.data['nopbc'] = True

    @register_from_funcs.register_funcs('gaussian/md')
    def from_gaussian_md(self, file_name):
        self.from_gaussian_log(file_name, md=True)

    @register_from_funcs.register_funcs('amber/md')
    def from_amber_md(self, file_name=None, parm7_file=None, nc_file=None, mdfrc_file=None, mden_file=None):
        # assume the prefix is the same if the spefic name is not given
        if parm7_file is None:
            parm7_file = file_name + ".parm7"
        if nc_file is None:
            nc_file = file_name + ".nc"
        if mdfrc_file is None:
            mdfrc_file = file_name + ".mdfrc"
        if mden_file is None:
            mden_file = file_name + ".mden"
        self.data = dpdata.amber.md.read_amber_traj(parm7_file, nc_file, mdfrc_file, mden_file)

    @register_from_funcs.register_funcs('cp2k/output')
    def from_cp2k_output(self, file_name) :
        self.data['atom_names'], \
            self.data['atom_numbs'], \
            self.data['atom_types'], \
            self.data['cells'], \
            self.data['coords'], \
            self.data['energies'], \
            self.data['forces'], \
            tmp_virial \
            = dpdata.cp2k.output.get_frames(file_name)
        if tmp_virial is not None:
            self.data['virials'] = tmp_virial
    @register_from_funcs.register_funcs('movement')
    @register_from_funcs.register_funcs('MOVEMENT')
    @register_from_funcs.register_funcs('mlmd')
    @register_from_funcs.register_funcs('MLMD')
    @register_from_funcs.register_funcs('pwmat/movement')
    @register_from_funcs.register_funcs('pwmat/MOVEMENT')
    @register_from_funcs.register_funcs('pwmat/mlmd')
    @register_from_funcs.register_funcs('pwmat/MLMD')
    def from_pwmat_output(self, file_name, begin = 0, step = 1) :
        self.data['atom_names'], \
            self.data['atom_numbs'], \
            self.data['atom_types'], \
            self.data['cells'], \
            self.data['coords'], \
            self.data['energies'], \
            self.data['forces'], \
            tmp_virial, \
            = dpdata.pwmat.movement.get_frames(file_name, begin = begin, step = step)
        if tmp_virial is not None :
            self.data['virials'] = tmp_virial
        # scale virial to the unit of eV
        if 'virials' in self.data :
            v_pref = 1 * 1e3 / 1.602176621e6
            for ii in range (self.get_nframes()) :
                vol = np.linalg.det(np.reshape(self.data['cells'][ii], [3,3]))
                self.data['virials'][ii] *= v_pref * vol
        # rotate the system to lammps convention
        self.rot_lower_triangular()


    def sub_system(self, f_idx) :
        """
        Construct a subsystem from the system

        Parameters
        ----------
        f_idx : int or index
            Which frame to use in the subsystem

        Returns
        -------
        sub_system : LabeledSystem
            The subsystem
        """
        tmp_sys = LabeledSystem()
        tmp_sys.data = System.sub_system(self, f_idx).data
        tmp_sys.data['energies'] = np.atleast_1d(self.data['energies'][f_idx])
        tmp_sys.data['forces'] = self.data['forces'][f_idx].reshape(-1, self.data['forces'].shape[1], 3)
        if 'virials' in self.data:
            tmp_sys.data['virials'] = self.data['virials'][f_idx].reshape(-1, 3, 3)
        return tmp_sys

    @register_to_funcs.register_funcs("ase/structure")
    def to_ase_structure(self):
        '''Convert System to ASE Atoms object.'''
        from ase import Atoms
        from ase.calculators.singlepoint import SinglePointCalculator
        
        structures = []

        for system in self.to_list():
            species=[system.data['atom_names'][tt] for tt in system.data['atom_types']]
            structure=Atoms(
                symbols=species,
                positions=system.data['coords'][0],
                pbc=True,
                cell=system.data['cells'][0]
            )

            results = {
                'energy': system.data["energies"][0],
                'forces': system.data["forces"][0]
            }
            if "virials" in system.data:
                # convert to GPa as this is ase convention
                v_pref = 1 * 1e4 / 1.602176621e6
                vol = structure.get_volume()
                results['stress'] = system.data["virials"][0] / (v_pref * vol)

            structure.calc = SinglePointCalculator(structure, **results)
            structures.append(structure)

        return structures

    def append(self, system):
        """
        Append a system to this system

        Parameters
        ----------
        system : System
            The system to append
        """
        if not System.append(self, system):
            # skip if this system or the system to append is non-converged
            return
        tgt = ['energies', 'forces']
        for ii in ['atom_pref']:
            if ii in self.data:
                tgt.append(ii)
        if ('virials' in system.data) and ('virials' not in self.data):
            raise RuntimeError('system has virial, but this does not')
        if ('virials' not in system.data) and ('virials' in self.data):
            raise RuntimeError('this has virial, but system does not')
        if 'virials' in system.data :
            tgt.append('virials')
        for ii in tgt:
            self.data[ii] = np.concatenate((self.data[ii], system[ii]), axis = 0)

    def sort_atom_types(self):
        idx = System.sort_atom_types(self)
        self.data['forces'] = self.data['forces'][:, idx]
        for ii in ['atom_pref']:
            if ii in self.data:
                self.data[ii] = self.data[ii][:, idx]

    def shuffle(self):
        """Also shuffle labeled data e.g. energies and forces."""
        idx = System.shuffle(self)
        for ii in ['energies', 'forces', 'virials', 'atom_pref']:
            if ii in self.data:
                self.data[ii] = self.data[ii][idx]
        return idx

    def to_pymatgen_ComputedStructureEntry(self):
        '''
        convert System to Pymagen ComputedStructureEntry obj

        '''
        try:
           from pymatgen.entries.computed_entries import ComputedStructureEntry
        except:
           raise ImportError('No module ComputedStructureEntry in pymatgen.entries.computed_entries')

        entries=[]
        for system in self.to_list():
            structure=system.to_pymatgen_structure()[0]
            energy=system.data['energies'][0]
            data={'forces':system.data['forces'][0],
                  'virials':system.data['virials'][0]}

            entry=ComputedStructureEntry(structure,energy,data=data)
            entries.append(entry)
        return entries

class MultiSystems:
    '''A set containing several systems.'''

    def __init__(self, *systems,type_map=None):
        """
        Parameters
        ----------
        systems : System
            The systems contained
        type_map : list of str
            Maps atom type to name
        """
        self.systems = {}
        if type_map is not None:
            self.atom_names = type_map
        else:
            self.atom_names = []
        self.append(*systems)

    def __getitem__(self, key):
        """Returns proerty stored in System by key or by idx"""
        if isinstance(key, int):
            return list(self.systems.values())[key]
        return self.systems[key]

    def __len__(self):
        return len(self.systems)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return 'MultiSystems ({} systems containing {} frames)'.format(len(self.systems), self.get_nframes())

    def __add__(self, others) :
       """magic method "+" operation """
       self_copy = deepcopy(self)
       if isinstance(others, System) or isinstance(others, MultiSystems):
          return self.__class__(self, others)
       elif isinstance(others, list):
          return self.__class__(self, *others)
       raise RuntimeError("Unspported data structure")

    @classmethod
    def from_file(cls,file_name,fmt):
        multi_systems = cls()
        multi_systems.load_systems_from_file(file_name=file_name,fmt=fmt)
        return multi_systems

    @classmethod
    def from_dir(cls,dir_name, file_name, fmt='auto', type_map=None):
        multi_systems = cls()
        target_file_list = sorted(glob.glob('./{}/**/{}'.format(dir_name, file_name), recursive=True))
        for target_file in target_file_list:
            multi_systems.append(LabeledSystem(file_name=target_file, fmt=fmt, type_map=type_map))
        return multi_systems


    def load_systems_from_file(self, file_name=None, fmt=None):
        if file_name is not None:
            if fmt is None:
                raise RuntimeError("must specify file format for file {}".format(file_name))
            elif fmt == 'quip/gap/xyz' or 'xyz':
                self.from_quip_gap_xyz_file(file_name)
            else:
                raise RuntimeError("unknown file format for file {} format {},now supported 'quip/gap/xyz'".format(file_name, fmt))


    def get_nframes(self) :
        """Returns number of frames in all systems"""
        return sum(len(system) for system in self.systems.values())

    def append(self, *systems) :
        """
        Append systems or MultiSystems to systems

        Parameters
        ----------
        system : System
            The system to append
        """
        for system in systems:
            if isinstance(system, System):
                self.__append(system)
            elif isinstance(system, MultiSystems):
                for sys in system:
                    self.__append(sys)
            else:
                raise RuntimeError("Object must be System or MultiSystems!")

    def __append(self, system):
        if not system.formula:
            return
        self.check_atom_names(system)
        formula = system.formula
        if formula in self.systems:
            self.systems[formula].append(system)
        else:
            self.systems[formula] = system.copy()

    def check_atom_names(self, system):
        """
        Make atom_names in all systems equal, prevent inconsistent atom_types.
        """
        # new_in_system = set(system["atom_names"]) - set(self.atom_names)
        # new_in_self = set(self.atom_names) - set(system["atom_names"])
        new_in_system = [e for e in system["atom_names"] if e not in self.atom_names]
        new_in_self = [e for e in self.atom_names if e not in system["atom_names"]]
        if len(new_in_system):
            # A new atom_name appear, add to self.atom_names
            self.atom_names.extend(new_in_system)
            # Add this atom_name to each system, and change their names
            new_systems = {}
            for each_system in self.systems.values():
                each_system.add_atom_names(new_in_system)
                each_system.sort_atom_names(type_map=self.atom_names)
                new_systems[each_system.formula] = each_system
            self.systems = new_systems
        if len(new_in_self):
            # Previous atom_name not in this system
            system.add_atom_names(new_in_self)
        system.sort_atom_names(type_map=self.atom_names)

    def from_quip_gap_xyz_file(self,file_name):
        # quip_gap_xyz_systems = QuipGapxyzSystems(file_name)
        # print(next(quip_gap_xyz_systems))
        for info_dict in QuipGapxyzSystems(file_name):
            system=LabeledSystem(data=info_dict)
            system.sort_atom_names()
            self.append(system)


    def to_deepmd_raw(self, folder) :
        """
        Dump systems in deepmd raw format to `folder` for each system.
        """
        for system_name, system in self.systems.items():
            system.to_deepmd_raw(os.path.join(folder, system_name))

    def to_deepmd_npy(self, folder, set_size = 5000, prec=np.float32) :
        """
        Dump the system in deepmd compressed format (numpy binary) to `folder` for each system.

        Parameters
        ----------
        folder : str
            The output folder
        set_size : int
            The size of each set.
        prec: {numpy.float32, numpy.float64}
            The floating point precision of the compressed data
        """
        for system_name, system in self.systems.items():
            system.to_deepmd_npy(os.path.join(folder, system_name),
                                 set_size = set_size,
                                 prec = prec)

    def from_deepmd_raw(self, folder):
        for dd in os.listdir(folder):
            self.append(LabeledSystem(os.path.join(folder, dd), fmt='deepmd/raw'))
        return self

    def from_deepmd_npy(self, folder):
        for dd in os.listdir(folder):
            self.append(LabeledSystem(os.path.join(folder, dd), fmt='deepmd/npy'))
        return self

    def predict(self, dp):
        import deepmd.DeepPot as DeepPot
        if not isinstance(dp, DeepPot):
            dp = DeepPot(dp)
        new_multisystems = dpdata.MultiSystems()
        for ss in self:
            new_multisystems.append(ss.predict(dp))
        return new_multisystems


def check_System(data):
    keys={'atom_names','atom_numbs','cells','coords','orig','atom_types'}
    assert( isinstance(data,dict) )
    assert( keys.issubset(set(data.keys())) )
    if len(data['coords']) > 0 :
        assert( len(data['coords'][0])==len(data['atom_types'])==sum(data['atom_numbs']) )
    else :
        assert( len(data['atom_types'])==sum(data['atom_numbs']) )
    assert( len(data['cells']) == len(data['coords']) )
    assert( len(data['atom_names'])==len(data['atom_numbs']) )

def check_LabeledSystem(data):
    keys={'atom_names', 'atom_numbs', 'atom_types', 'cells', 'coords', 'energies',
           'forces', 'orig'}

    assert( keys.issubset(set(data.keys())) )
    assert( isinstance(data,dict) )
    assert( len(data['atom_names'])==len(data['atom_numbs']) )

    if len(data['coords']) > 0 :
        assert( len(data['coords'][0])==len(data['atom_types']) ==sum(data['atom_numbs'])  )
    else:
        assert( len(data['atom_types']) ==sum(data['atom_numbs'])  )
    if 'virials' in data:
        assert( len(data['cells']) == len(data['coords']) == len(data['virials']) == len(data['energies']) )
    else:
        assert( len(data['cells']) == len(data['coords']) == len(data['energies']) )


def elements_index_map(elements,standard=False,inverse=False):
    if standard:
        elements.sort(key=lambda x: Element(x).Z)
    if inverse:
        return dict(zip(range(len(elements)),elements))
    else:
        return dict(zip(elements,range(len(elements))))



# %%
