import numpy as np
import dpdata.lammps.lmp
import dpdata.lammps.dump
import dpdata.vasp.poscar
import dpdata.vasp.xml

class System :    
    data = {}
    '''
    for example a water system that has two molecules
    data['atom_numbs'] = [2, 4]
    data['atom_names'] = ['O', 'H']
    data['atom_types'] = [0, 1, 1, 0, 1, 1]
    data['orig'] = [0, 0, 0]
    data['cells'] = [ [...], [...] ]    # two frames
    data['coords'] = [ [...], [...] ]  # two frames

    Assumptions : 
    data['orig'] is always [0, 0, 0]
    data['cells'][ii] is always lower triangular (lammps axes convention)
    '''

    def __init__ (self) :
        self.data = {}
        self.data['atom_numbs'] = []
        self.data['atom_names'] = []
        self.data['atom_types'] = []
        self.data['orig'] = [0, 0, 0]
        self.data['cells'] = []
        self.data['coords'] = []


    def get_data(self) :
        return self.data

    def get_nframes(self) :
        return len(self.data['cells'])

    def sub_system(self, f_idx) :
        tmp = System()
        for ii in ['atom_numbs', 'atom_names', 'atom_types', 'orig'] :
            tmp.data[ii] = self.data[ii]
        tmp.data['cells'] = self.data['cells'][f_idx]
        tmp.data['coords'] = self.data['coords'][f_idx]
        return tmp
        

    def from_lammps_lmp (self, file_name, type_map = None) :
        with open(file_name) as fp:
            lines = [line.rstrip('\n') for line in fp]
            self.data = dpdata.lammps.lmp.to_system_data(lines, type_map)
        self._shift_orig_zero()


    def to_lammps_lmp(self, file_name, frame_idx = 0) :
        assert(frame_idx < len(self.data['coords']))
        w_str = dpdata.lammps.lmp.from_system_data(self.data, frame_idx)
        with open(file_name, 'w') as fp:
            fp.write(w_str)
    

    def from_lammps_dump (self, file_name, type_map = None) :
        with open(file_name) as fp:
            lines = [line.rstrip('\n') for line in fp]
            self.data = dpdata.lammps.dump.system_data(lines, type_map)
        self._shift_orig_zero()


    def from_vasp_poscar(self, file_name) :
        with open(file_name) as fp:
            lines = [line.rstrip('\n') for line in fp]
            self.data = dpdata.vasp.poscar.to_system_data(lines)
        self.rot_lower_triangular()

    
    def to_vasp_poscar(self, file_name, frame_idx = 0) :
        assert(frame_idx < len(self.data['coords']))
        w_str = dpdata.vasp.poscar.from_system_data(self.data, frame_idx)
        with open(file_name, 'w') as fp:
            fp.write(w_str)


    def affine_map(self, trans, f_idx = 0) :
        assert(np.linalg.det(trans) != 0)
        cell = self.data['cells'][f_idx]
        for ii in range(3) :
            cell[ii] = np.matmul(cell[ii], trans)
        self.data['cells'][f_idx] = cell
        for ff in self.data['coords'] :            
            tmp_ff = ff
            for ii in range(tmp_ff.shape[0]): 
                tmp_ff[ii] = np.matmul(tmp_ff[ii], trans)


    def _shift_orig_zero(self) :
        for ff in self.data['coords'] :
            for ii in ff :
                ii = ii - self.data['orig']
        self.data['orig'] = self.data['orig'] - self.data['orig']
        assert((np.zeros([3]) == self.data['orig']).all())

    def rot_lower_triangular(self) :
        for ii in range(self.get_nframes()) :
            self._rot_lower_triangular(ii)

    def _rot_lower_triangular(self, f_idx = 0) :
        qq, rr = np.linalg.qr(self.data['cells'][f_idx].T)
        self.affine_map(qq)
        rot = np.eye(3)
        if self.data['cells'][f_idx][0][0] < 0 :
            rot[0][0] = -1
        if self.data['cells'][f_idx][1][1] < 0 :
            rot[1][1] = -1
        if self.data['cells'][f_idx][2][2] < 0 :
            rot[2][2] = -1
        assert(np.linalg.det(rot) > 0) 
        self.affine_map(rot, f_idx = f_idx)
        

class LabeledSystem (System): 
    def from_vasp_xml(self, file_name) :
        self.data['atom_names'], \
            self.data['atom_types'], \
            self.data['cells'], \
            self.data['coords'], \
            self.data['energies'], \
            self.data['forces'], \
            self.data['virials'], \
            = dpdata.vasp.xml.analyze(file_name, type_idx_zero = True)
        self.data['atom_numbs'] = []
        for ii in range(len(self.data['atom_names'])) :
            self.data['atom_numbs'].append(sum(self.data['atom_types'] == ii))
        # the vasp xml assumes the direct coordinates
        # apply the transform to the cartesan coordinates
        for ii in range(self.get_nframes()) :
            self.data['coords'][ii] = np.matmul(self.data['coords'][ii], self.data['cells'][ii])
        # rotate the system to lammps convention
        self.rot_lower_triangular()
        # scale virial to the unit of eV
        v_pref = 1 * 1e3 / 1.602176621e6
        for ii in range (self.get_nframes()) :
            vol = np.linalg.det(np.reshape(self.data['cells'][ii], [3,3]))
            self.data['virials'][ii] *= v_pref * vol
        # print(self.data)
