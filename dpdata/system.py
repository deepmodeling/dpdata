import numpy as np
import dpdata.lammps.lmp
import dpdata.vasp.poscar

class System :    
    data = {}
    '''
    for example a water system that has two molecules
    data['atom_numbs'] = [2, 4]
    data['atom_names'] = ['O', 'H']
    data['atom_types'] = [0, 1, 1, 0, 1, 1]
    data['orig'] = [0, 0, 0]
    data['cell'] = [ [...], [...] ]    # two frames
    data['frames'] = [ [...], [...] ]  # two frames

    Assumptions : 
    data['orig'] is always [0, 0, 0]
    data['cell'][ii] is always lower triangular (lammps axes convention)
    '''

    def __init__ (self) :
        self.data = {}
        self.data['atom_numbs'] = []
        self.data['atom_names'] = []
        self.data['atom_types'] = []
        self.data['orig'] = [0, 0, 0]
        self.data['cell'] = []
        self.data['frames'] = []


    def get_data(self) :
        return self.data


    def from_lammps_lmp (self, file_name, type_map = None) :
        with open(file_name) as fp:
            lines = [line.rstrip('\n') for line in fp]
            self.data = dpdata.lammps.lmp.to_system_data(lines, type_map)
        self._shift_orig_zero()


    def to_lammps_lmp(self, file_name, frame_idx = 0) :
        assert(frame_idx < len(self.data['frames']))
        w_str = dpdata.lammps.lmp.from_system_data(self.data, frame_idx)
        with open(file_name, 'w') as fp:
            fp.write(w_str)
    

    def from_vasp_poscar(self, file_name) :
        with open(file_name) as fp:
            lines = [line.rstrip('\n') for line in fp]
            self.data = dpdata.vasp.poscar.to_system_data(lines)
        self._rot_lower_triangular()

    
    def to_vasp_poscar(self, file_name, frame_idx = 0) :
        assert(frame_idx < len(self.data['frames']))
        w_str = dpdata.vasp.poscar.from_system_data(self.data, frame_idx)
        with open(file_name, 'w') as fp:
            fp.write(w_str)


    def affine_map(self, trans, f_idx = 0) :
        assert(np.linalg.det(trans) != 0)
        cell = self.data['cell'][f_idx]
        for ii in range(3) :
            cell[ii] = np.matmul(cell[ii], trans)
        self.data['cell'][f_idx] = cell
        for ff in self.data['frames'] :            
            tmp_ff = ff
            for ii in range(tmp_ff.shape[0]): 
                tmp_ff[ii] = np.matmul(tmp_ff[ii], trans)


    def _shift_orig_zero(self) :
        for ff in self.data['frames'] :
            for ii in ff :
                ii = ii - self.data['orig']
        self.data['orig'] = self.data['orig'] - self.data['orig']
        assert((np.zeros([3]) == self.data['orig']).all())


    def _rot_lower_triangular(self, f_idx = 0) :
        qq, rr = np.linalg.qr(self.data['cell'][f_idx].T)
        self.affine_map(qq)
        rot = np.eye(3)
        if self.data['cell'][f_idx][0][0] < 0 :
            rot[0][0] = -1
        if self.data['cell'][f_idx][1][1] < 0 :
            rot[1][1] = -1
        if self.data['cell'][f_idx][2][2] < 0 :
            rot[2][2] = -1
        assert(np.linalg.det(rot) > 0) 
        self.affine_map(rot)
        
