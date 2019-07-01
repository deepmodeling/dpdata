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
from monty.json import MSONable

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


    def as_dict(self):
        """Returns data dict of System instance"""
        d={"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "data": self.data
          }
        return d
    

    @classmethod
    def from_dict(cls,d):
        """Contruct System instance from dict object"""
        data=d['data']
        return cls(data=data)
 
 
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
        return self.__class__.from_dict({'data':self.data})


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
        for ii in ['atom_numbs', 'atom_names'] :
            assert(system.data[ii] == self.data[ii])
        for ii in ['atom_types','orig'] :
            eq = (system.data[ii] == self.data[ii])
            assert(eq.all())
        for ii in ['coords', 'cells'] :
            self.data[ii] = np.concatenate((self.data[ii], system[ii]), axis = 0)


    def extend(self, systems):
        """
        Extend a system list to this system

        Parameters
        ----------
        systems : [System1, System2, System3 ]
            The list to extend
        """
        
        for system in systems:
            self.append(system)

            
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
        
def check_System(data):
    keys={'atom_names','atom_numbs','cells','coords','orig','atom_types'}
    assert( isinstance(data,dict) )
    assert( set(data.keys())==keys )
    assert( len(data['coords'][0])==len(data['atom_types'])==sum(data['atom_numbs']) )
    assert( len(data['atom_names'])==len(data['atom_numbs']) )
    assert( len(data['cells']) == len(data['coords']) )


def check_LabeledSystem(data):
    keys={'atom_names', 'atom_numbs', 'atom_types', 'cells', 'coords', 'energies',
           'forces', 'orig', 'virials'}
    assert( isinstance(data,dict) )
    assert( set(data.keys())==keys )
    assert( len(data['atom_names'])==len(data['atom_numbs']) )
    assert( len(data['coords'][0])==len(data['atom_types']) ==sum(data['atom_numbs'])  )
    if len(data['virials'])>0:
       assert( len(data['cells']) == len(data['coords']) == len(data['virials']) == len(data['energies']) )
    else:
       assert( len(data['cells']) == len(data['coords']) == len(data['energies']) )


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
        ret+="\nIncluding Virials  : %s"% "Yes" if self.has_virial() else "No"
        ret+="\nElement List       :"
        ret+="\n-------------------"
        ret+="\n"+"  ".join(map(str,self.get_atom_names()))
        ret+="\n"+"  ".join(map(str,self.get_atom_numbs()))
        return ret
    
    @classmethod
    def from_dict(cls,d):
        """Contruct LabeledSystem instance from dict object"""
        data=d['data']
        return cls(data=data)

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
        self.data = dpdata.gaussian.log.to_system_data(file_name)
        self.data['cells'] = np.array([[[100., 0., 0.], [0., 100., 0.], [0., 0., 100.]]])
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
        System.append(self, system)
        tgt = ['energies', 'forces']
        if len(system.data['virials']) != 0 and len(self.data['virials']) == 0:
            raise RuntimeError('system has virial, but this does not')
        if len(system.data['virials']) == 0 and len(self.data['virials']) != 0:
            raise RuntimeError('this has virial, but system does not')
        if len(system.data['virials']) != 0 :
            tgt.append('virials')
        for ii in tgt:
            self.data[ii] = np.concatenate((self.data[ii], system[ii]), axis = 0)
