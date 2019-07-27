import os
import numpy as np
import dpdata.lammps.lmp
import dpdata.lammps.dump
import dpdata.vasp.poscar
import dpdata.vasp.xml
import dpdata.vasp.outcar
import dpdata.deepmd.raw
import dpdata.deepmd.comp
import dpdata.pwscf.traj
import dpdata.md.pbc
import dpdata.gaussian.log
from copy import deepcopy
from monty.json import MSONable
from monty.serialization import loadfn,dumpfn
from mendeleev import element as Element 

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
                  data = None) :
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
                - ``vasp/poscar``: vasp POSCAR
                - ``pwscf/traj``: pwscf trajectory files. should have: file_name+'.in' and file_name+'.pos'

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
        if fmt == 'auto':
            fmt = os.path.basename(file_name).split('.')[-1] 
        if fmt == 'lmp' or fmt == 'lammps/lmp' :
            self.from_lammps_lmp(file_name, type_map = type_map) 
        elif fmt == 'dump' or fmt == 'lammps/dump' :
            self.from_lammps_dump(file_name, type_map = type_map, begin = begin, step = step)
        elif fmt.lower() == 'poscar' or fmt.lower() == 'contcar' or fmt.lower() == 'vasp/poscar' or fmt.lower() == 'vasp/contcar':
            self.from_vasp_poscar(file_name)
        elif fmt == 'pwscf/traj':
            self.from_pwscf_traj(file_name, begin = begin, step = step)
        else :
            raise RuntimeError('unknow data format ' + fmt)

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
        """Returns proerty stored in System by key"""
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

    def set_atom_types(self,type_map=None):
        """
        Reset the type of the system
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
        system : System with specific type order
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
        _system=self.copy()
        _system.data['atom_types']=new_atom_types
        return _system


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
        tmp.data['cells'] = self.data['cells'][f_idx]
        tmp.data['coords'] = self.data['coords'][f_idx]
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
        assert(system.formula == self.formula)
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
            eq = (system.data[ii] == self.data[ii])
            assert(eq.all())
        for ii in ['coords', 'cells'] :
            self.data[ii] = np.concatenate((self.data[ii], system[ii]), axis = 0)
        return True

    def sort_atom_names(self):
        idx = np.argsort(self.data['atom_names'])
        self.data['atom_names'] = list(np.array(self.data['atom_names'])[idx])
        self.data['atom_numbs'] = list(np.array(self.data['atom_numbs'])[idx])
        self.data['atom_types'] = np.argsort(idx)[self.data['atom_types']]
    
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


    def from_lammps_lmp (self, file_name, type_map = None) :
        with open(file_name) as fp:
            lines = [line.rstrip('\n') for line in fp]
            self.data = dpdata.lammps.lmp.to_system_data(lines, type_map)
        self._shift_orig_zero()


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
    

    def from_lammps_dump (self, 
                          file_name, 
                          type_map = None, 
                          begin = 0,
                          step = 1) :
        lines = dpdata.lammps.dump.load_file(file_name, begin = begin, step = step)
        self.data = dpdata.lammps.dump.system_data(lines, type_map)
        self._shift_orig_zero()


    def from_vasp_poscar(self, file_name) :
        with open(file_name) as fp:
            lines = [line.rstrip('\n') for line in fp]
            self.data = dpdata.vasp.poscar.to_system_data(lines)
        self.rot_lower_triangular()

    
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
        assert(frame_idx < len(self.data['coords']))
        w_str = dpdata.vasp.poscar.from_system_data(self.data, frame_idx)
        with open(file_name, 'w') as fp:
            fp.write(w_str)

    
    def from_pwscf_traj(self, 
                        prefix, 
                        begin = 0,
                        step = 1) :
        self.data = dpdata.pwscf.traj.to_system_data(prefix + '.in', prefix, begin = begin, step = step)
        self.data['coords'] \
            = dpdata.md.pbc.apply_pbc(self.data['coords'], 
                                      self.data['cells'], 
            )
        self.rot_lower_triangular()


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
                  data=None) :
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
                - ``pwscf/traj``: pwscf trajectory files. should have: file_name+'.in', file_name+'.pos', file_name+'.evp' and file_name+'.for'
                - ``gaussian/log``: gaussian logs

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
        if fmt == 'auto':
            fmt = os.path.basename(file_name).split('.')[-1] 
        if fmt == 'xml' or fmt == 'XML' or fmt == 'vasp/xml' :
            self.from_vasp_xml(file_name, begin = begin, step = step) 
        elif fmt == 'outcar' or fmt == 'OUTCAR' or fmt == 'vasp/outcar' :
            self.from_vasp_outcar(file_name, begin = begin, step = step)
        elif fmt == 'deepmd' or fmt == 'deepmd/raw':
            self.from_deepmd_raw(file_name, type_map = type_map)
        elif fmt == 'deepmd/npy':
            self.from_deepmd_comp(file_name, type_map = type_map)
        elif fmt == 'pwscf/traj':
            self.from_pwscf_traj(file_name, begin = begin, step = step)
        elif fmt == 'gaussian/log':
            self.from_gaussian_log(file_name)
        else :
            raise RuntimeError('unknow data format ' + fmt)

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
        return ('virials' in self.data) and (len(self.data['virials']) > 0)


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
            self.data['virials'], \
            = dpdata.vasp.outcar.get_frames(file_name, begin = begin, step = step)
        # scale virial to the unit of eV
        if len(self.data['virials']) != 0 :            
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


    def from_deepmd_raw(self, folder, type_map = None) :
        tmp_data = dpdata.deepmd.raw.to_system_data(folder, type_map = type_map)
        if tmp_data is not None :
            self.data = tmp_data
    

    def to_deepmd_raw(self, folder) :
        """
        Dump the system in deepmd raw format to `folder`
        """
        dpdata.deepmd.raw.dump(folder, self.data)


    def from_deepmd_comp(self, folder, type_map = None) :
        self.data = dpdata.deepmd.comp.to_system_data(folder, type_map = type_map)

        
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

    
    def from_pwscf_traj(self, prefix, begin = 0, step = 1) :
        self.data = dpdata.pwscf.traj.to_system_data(prefix + '.in', prefix, begin = begin, step = step)
        self.data['coords'] \
            = dpdata.md.pbc.apply_pbc(self.data['coords'], 
                                      self.data['cells'], 
            )
        self.data['energies'], self.data['forces'] \
            = dpdata.pwscf.traj.to_system_label(prefix + '.in', prefix, begin = begin, step = step)
        self.data['virials'] = []
        self.rot_lower_triangular()


    def from_gaussian_log(self, file_name):
        try:
            self.data = dpdata.gaussian.log.to_system_data(file_name)
            self.data['cells'] = np.array([[[100., 0., 0.], [0., 100., 0.], [0., 0., 100.]]])
        except AssertionError:
            self.data['energies'], self.data['forces']= [], []
        self.data['virials'] = []


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
        tmp_sys.data['energies'] = self.data['energies'][f_idx]
        tmp_sys.data['forces'] = self.data['forces'][f_idx]
        if len(self.data['virials']) != 0 :
            tmp_sys.data['virials'] = self.data['virials'][f_idx]            
        return tmp_sys


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
        if len(system.data['virials']) != 0 and len(self.data['virials']) == 0:
            raise RuntimeError('system has virial, but this does not')
        if len(system.data['virials']) == 0 and len(self.data['virials']) != 0:
            raise RuntimeError('this has virial, but system does not')
        if len(system.data['virials']) != 0 :
            tgt.append('virials')
        for ii in tgt:
            self.data[ii] = np.concatenate((self.data[ii], system[ii]), axis = 0)
    
    def sort_atom_types(self):
        idx = System.sort_atom_types(self)
        self.data['forces'] = self.data['forces'][:, idx]
        for ii in ['atom_pref']:
            if ii in self.data:
                self.data[ii] = self.data[ii][:, idx]


class MultiSystems:
    '''A set containing several systems.'''

    def __init__(self, *systems):
        """
        Parameters
        ----------
        systems : System
            The systems contained
        """
        self.systems = {}
        self.atom_names = []
        self.append(*systems)

    def __getitem__(self, key):
        """Returns proerty stored in System by key"""
        return self.systems[key]

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
                for sys in system.system:
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
            self.systems[formula] = system
    
    def check_atom_names(self, system):
        """
        Make atom_names in all systems equal, prevent inconsistent atom_types.
        """
        new_in_system = set(system["atom_names"]) - set(self.atom_names)
        new_in_self = set(self.atom_names) - set(system["atom_names"])
        if len(new_in_system):
            # A new atom_name appear, add to self.atom_names
            self.atom_names.extend(new_in_system)
            self.atom_names.sort()
            # Add this atom_name to each system, and change their names
            new_systems = {}
            for each_system in self.systems.values():
                each_system.add_atom_names(new_in_system)  
                each_system.sort_atom_names()
                new_systems[each_system.formula] = each_system
            self.systems = new_systems
        if len(new_in_self):
            # Previous atom_name not in this system
            system.add_atom_names(new_in_self)
        system.sort_atom_names()

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

def check_System(data):
    keys={'atom_names','atom_numbs','cells','coords','orig','atom_types'}
    assert( isinstance(data,dict) )
    assert( set(data.keys())==keys )
    assert( len(data['coords'][0])==len(data['atom_types'])==sum(data['atom_numbs']) )
    assert( len(data['cells']) == len(data['coords']) )
    assert( len(data['atom_names'])==len(data['atom_numbs']) )

def check_LabeledSystem(data):
    keys={'atom_names', 'atom_numbs', 'atom_types', 'cells', 'coords', 'energies',
           'forces', 'orig', 'virials'}
    if 'virials' in data.keys():
        pass
    else:
        data['virials']=[]

    assert( set(data.keys())==keys )
    assert( isinstance(data,dict) )
    assert( len(data['atom_names'])==len(data['atom_numbs']) )

    assert( len(data['coords'][0])==len(data['atom_types']) ==sum(data['atom_numbs'])  )
    if len(data['virials'])>0:
       assert( len(data['cells']) == len(data['coords']) == len(data['virials']) == len(data['energies']) )
    else:
       assert( len(data['cells']) == len(data['coords']) == len(data['energies']) )


def elements_index_map(elements,standard=False,inverse=False):
    if standard:
       elements.sort(key=lambda x: Element(x).atomic_number)
    if inverse:
       return dict(zip(range(len(elements)),elements))
    else:
       return dict(zip(elements,range(len(elements))))
                                                         
