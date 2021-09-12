#%%
import os
import glob
import inspect
import random
import numpy as np
import dpdata.md.pbc
from copy import deepcopy
from monty.json import MSONable
from monty.serialization import loadfn,dumpfn
from dpdata.amber.mask import pick_by_amber_mask, load_param_file

from dpdata.system import System, check_System, load_format, get_atom_perturb_vector
from pymatgen.core.structure import Molecule
from pymatgen.core.operations import SymmOp
from collections import Counter
# from dpdata.molecule import Molecule

# ensure all plugins are loaded!
import dpdata.plugins
from dpdata.plugin import Plugin
from dpdata.format import Format


class CustomSystem(System):
    '''
    Custom data System

    For example, a organic-inorganic hybrid peroviskate system named `d_example` has a molecule (8 atoms) and a lattice (4 atoms). 
        For pure inorganic lattice, the labels can be accessed by
            - `d_lattice['atom_numbs']` : [1, 3]
            - `d_lattice['atom_names']` : ['Pb', 'I']
            - `d_lattice['atom_types']` : [0, 1]
            - `d_lattice['orig']` : [0, 0, 0]
            - `d_lattice['cells']` : a numpy array of size nframes x 3 x 3
            - `d_lattice['coords']` : a numpy array of size nframes x natoms x 3
            - `d_lattice['energies']` : a numpy array of size nframes
            - `d_lattice['forces']` : a numpy array of size nframes x 6 x 3
            - `d_lattice['virials']` : optional, a numpy array of size nframes x 3 x 3
        For hybrid system, the labels can be accessed by
            - `d_example['atom_numbs']` : [1, 3, 1, 2, 5]
            - `d_example['atom_names']` : ['Pb', 'I', 'C', 'N', 'H']
            - `d_example['atom_types']` : [0, 1, 2, 3, 4]
            - `d_example['orig']` : [0, 0, 0]
            - `d_example['cells']` : a numpy array of size nframes x 3 x 3
            - `d_example['coords']` : a numpy array of size nframes x natoms x 3
            - `d_example['energies']` : a numpy array of size nframes
            - `d_example['forces']` : a numpy array of size nframes x 6 x 3
            - `d_example['virials']` : optional, a numpy array of size nframes x 3 x 3

    It is noted that
        - The order of frames stored in `'energies'`, `'forces'` and `'virials'` should be consistent with `'atom_types'`, `'cells'` and `'coords'`.
        - The order of atoms in **every** frame of `'forces'` should be consistent with `'coords'` and `'atom_types'`.
    '''

    fmt_BondOrder = ["mol", "sdf"]

    def __init__ (self,
                  file_name = None,
                  fmt = 'auto',
                  type_map = None,
                  begin = 0,
                  step = 1,
                  data=None,
                  mols=None,
                  **kwargs) :

        System.__init__(self)
        self.mols = []

        if data:
           check_System(data)
           self.data=data
           return
        if file_name is None :
            return
        self.from_fmt(file_name, fmt, type_map=type_map, begin= begin, step=step, **kwargs)
        if type_map is not None:
            self.apply_type_map(type_map)
        if mols:
            self.mols = mols

    post_funcs = Plugin() + System.post_funcs

    def from_fmt_obj(self, fmtobj, file_name, **kwargs):
        data = fmtobj.from_system(file_name, **kwargs)
        if data:
            if isinstance(data, (list, tuple)):
                for dd in data:
                    self.append(CustomSystem(data=dd))
            else:
                self.data = {**self.data, **data}
            if hasattr(fmtobj.from_system, 'post_func'):
                for post_f in fmtobj.from_system.post_func:
                    self.post_funcs.get_plugin(post_f)(self)
        return self

    def to(self, fmt, *args, **kwargs):
        return self.to_fmt_obj(load_format(fmt), *args, **kwargs)
    
    def to_fmt_obj(self, fmtobj, *args, **kwargs):
        data_hybrid = deepcopy(self.data)
        for item in self.mols:
            mol = item["mol"]
            at = item["at"]
            at_is_cartesian = item["at_is_cartesian"]
            data_hybrid = self._merge_data(data_hybrid, mol, at, at_is_cartesian)
        return fmtobj.to_system(data_hybrid, *args, **kwargs)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        ret="Data Summary"
        ret+="\nCustom System"
        ret+="\n-------------------"
        ret+="\nFrame Numbers      : %d"%self.get_nframes()
        ret+="\nAtom Numbers       : %d"%self.get_natoms()
        ret+="\nElement List       :"
        ret+="\n-------------------"
        ret+="\n"+"  ".join(map(str,self.get_atom_names()))
        ret+="\n"+"  ".join(map(str,self.get_atom_numbs()))
        if (len(self.mols) > 0):
            ret+="\nMolecule List:\n"
            for item in self.mols:
                mol = item["mol"]
                at = item["at"]
                at_is_cartesian = item["at_is_cartesian"]
                ret+=mol.__str__()
                if at_is_cartesian:
                    ret+="\n"+"@ cartesian coordinate:"
                    ret+="\n"+"  ".join(map(str,at))
                    ret+="\n"
                else:
                    ret+="\n"+"@ direct coordinate:"
                    ret+="\n"+"  ".join(map(str,at))
                    ret+="\n"
        return ret

    def __getitem__(self, key):
        """Returns proerty stored in System by key or by idx"""
        if isinstance(key, (int, slice)):
            return self.sub_system(key)
        data_hybrid = deepcopy(self.data)
        for item in self.mols:
            mol = item["mol"]
            at = item["at"]
            at_is_cartesian = item["at_is_cartesian"]
            data_hybrid = self._merge_data(data_hybrid, mol, at, at_is_cartesian)
        return data_hybrid[key]

    def get_nmols(self) :
        """Returns number of frames in the system"""
        return len(self.mols)


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
            self.data[ii] = np.concatenate((self.data[ii], system.data[ii]), axis = 0)
        if self.nopbc and not system.nopbc:
            # appended system uses PBC, cancel nopbc
            self.data['nopbc'] = False
        for item in system.mols:
            self.mols.append(item)
        return True


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
        tmp_sys = CustomSystem()
        tmp_sys.data = System.sub_system(self, f_idx).data
        tmp_sys.mols = self.mols
        return tmp_sys

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

        tmp = CustomSystem()
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

        for i in range(np.prod(ncopy)):
            for item in self.mols:
                tmp.mols.append(deepcopy(item))
        long_idx = 0
        for xx in range(ncopy[0]):
            for yy in range(ncopy[1]):
                for zz in range(ncopy[2]):
                    for idx in range(len(self.mols)):
                        item = self.mols[idx]
                        at = deepcopy(item["at"])
                        at_is_cartesian = item["at_is_cartesian"]
                        if at_is_cartesian:
                            raise NotImplementedError("at_is_cartesian of mols doesn't support CustomSystem.replicate")
                        else:
                            assert(long_idx < np.prod(ncopy)*len(self.mols))
                            tmp.mols[long_idx]["at"][0] = at[0] / ncopy[0] + xx / ncopy[0]
                            tmp.mols[long_idx]["at"][1] = at[1] / ncopy[1] + yy / ncopy[1]
                            tmp.mols[long_idx]["at"][2] = at[2] / ncopy[2] + zz / ncopy[2]
                            long_idx += 1

        return tmp


    def rotate_mol(self, axis = None, angle = None):
        if axis:
            assert(len(axis) == len(self.mols))
        else:
            axis = [None]*len(self.mols)
        if angle:
            assert(len(angle) == len(self.mols))
        else:
            angle = [None]*len(self.mols)

        for idx in range(len(self.mols)):
            self.rotate_mol_by_idx(idx, axis[idx], angle[idx])

    def rotate_mol_by_idx(self, idx, axis = None, angle = None):
        mol = self.mols[idx]["mol"]
        self.mols[idx]["mol"] = self._rotate_mol(mol, axis, angle)

    def _rotate_mol(self, mol, axis, angle):
        center = mol.center_of_mass
        centered_coords = np.copy(mol.cart_coords) - center

        if not axis:
            axis = np.random.rand(3)
        if not angle:
            angle=random.uniform(-180, 180)
        print("Rotate by angle %f around axis (%f,  %f,  %f)" % (angle, axis[0], axis[1], axis[2]))
        op = SymmOp.from_origin_axis_angle(
            (0, 0, 0),
            axis=np.array(axis),
            angle=angle,
        )

        m = op.rotation_matrix
        new_coords = np.dot(m, centered_coords.T).T
        return Molecule(
            mol.species_and_occu, 
            new_coords,
            charge=mol._charge,
            spin_multiplicity=mol._spin_multiplicity,
            site_properties=mol.site_properties,
            )


    def getLatticeData(self):
        return self.data

    def popmol(self, idx = -1):
        return self.mols.pop(idx)

    def addmol(self, file_name = None, at = [0.0, 0.0, 0.0], at_is_cartesian = False, mol = None):
        if mol:
            if isinstance (mol, Molecule):
                self.mols.append(dict({"mol": mol, "at": at, "at_is_cartesian": at_is_cartesian}))
            else:
                raise TypeError("Unsupported object type %r" , mol)
            return

        if file_name is None:
            return
        mol = Molecule.from_file(file_name)
        self.mols.append(dict({"mol": mol, "at": at, "at_is_cartesian": at_is_cartesian}))

        
    def _merge_data(self, data, mol, at, at_is_cartesian):
        elem_mol = Extract_elem(mol)
        elem_latt = data['atom_names']
        atomNumbs_latt = data['atom_numbs']
        elemIdx_mol, orderedUniqElem, atomNumbs = assignIdxElem(elem_mol, elem_latt, atomNumbs_latt)
        # self.data = {}
        # self.data['atom_numbs'] = []
        # self.data['atom_names'] = []
        # self.data['atom_types'] = []
        # self.data['orig'] = np.array([0, 0, 0])
        # self.data['cells'] = []
        # self.data['coords'] = []
        new_data = {}
        for key in data.keys():
            if key == "atom_names":
                new_data['atom_names'] = orderedUniqElem
                continue
            if key == "atom_types":
                new_data['atom_types'] = np.concatenate((np.copy(data['atom_types']), elemIdx_mol), axis = 0)
                continue
            if key == "atom_numbs":
                new_data['atom_numbs'] =  atomNumbs
                continue
            if key == "coords":
                continue
            new_data[key]= deepcopy(data[key])
        
        
        nframes = self.get_nframes()
        for to_frame in range(nframes):
            if at_is_cartesian:
                cart_at = deepcopy(at)
            else:
                cart_at = np.array([0.0, 0.0, 0.0])
                cart_at[0] = at[0]*data['cells'][to_frame][0][0] + at[1]*data['cells'][to_frame][1][0] + at[2]*data['cells'][to_frame][2][0]
                cart_at[1] = at[0]*data['cells'][to_frame][0][1] + at[1]*data['cells'][to_frame][1][1] + at[2]*data['cells'][to_frame][2][1]
                cart_at[2] = at[0]*data['cells'][to_frame][0][2] + at[1]*data['cells'][to_frame][1][2] + at[2]*data['cells'][to_frame][2][2]
            center = mol.center_of_mass
            newcoords_mol = np.copy(mol.cart_coords) - center + cart_at
            if "coords" in new_data:
                new_data["coords"].append(np.concatenate((np.copy(data["coords"][to_frame]), np.array(newcoords_mol))))
            else:
                new_data["coords"] = np.concatenate((np.copy(data["coords"][to_frame]), np.array(newcoords_mol))).reshape(1, -1, 3)
        return new_data
        
        
        
def Extract_elem(mol):
    return list(str(site.species.elements[0]) for site in mol.sites)

def assignIdxElem(elem_mol, elem_latt, atomNumbs_latt):
    elemIdx_mol = list()
    orderedUniqElem = deepcopy(elem_latt)
    atomNumbs = deepcopy(atomNumbs_latt)
    for item in elem_mol:
        if item in orderedUniqElem:
            idx = orderedUniqElem.index(item)
            atomNumbs[idx] += 1
            elemIdx_mol.append(idx)
        else:
            idx = len(orderedUniqElem)
            elemIdx_mol.append( idx)
            orderedUniqElem.append(item)
            atomNumbs.append(1)
    return list(elemIdx_mol), orderedUniqElem, list(atomNumbs)
    
def removeElem(elem_mol, elem, atomNumbs):
    assert(len(elem) == len(atomNumbs))
    elem_latt = deepcopy(elem)
    atomNumbs_latt = deepcopy(atomNumbs)
    for idx in range(len(elem)):
        item = elem[idx]
        if item in elem_mol:
            atomNumbs_latt[idx] -= 1
    pickidces = np.where(atomNumbs_latt >= 0)
    return elem_latt[pickidces], atomNumbs_latt[pickidces], 


def c2d(cart, cell):
    reccell = np.zeros((3,3))
    Vcell = np.dot(cell[0], np.cross(cell[1], cell[2]))
    reccell[0] = np.cross(cell[1], cell[2])/Vcell
    reccell[1] = np.cross(cell[2], cell[0])/Vcell
    reccell[2] = np.cross(cell[0], cell[1])/Vcell
    direct = np.zeros(3)
    direct[0] = cart[0]*reccell[0][0] + cart[1]*reccell[1][0] + cart[2]*reccell[2][0]
    direct[1] = cart[1]*reccell[0][1] + cart[1]*reccell[1][1] + cart[2]*reccell[2][1]
    direct[2] = cart[2]*reccell[0][2] + cart[1]*reccell[1][2] + cart[2]*reccell[2][2]
    return direct

