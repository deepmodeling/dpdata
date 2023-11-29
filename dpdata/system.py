# %%
import glob
import hashlib
import os
import warnings
from copy import deepcopy
from typing import Any, Dict, Optional, Tuple, Union

import numpy as np
from monty.json import MSONable
from monty.serialization import dumpfn, loadfn

import dpdata
import dpdata.md.pbc

# ensure all plugins are loaded!
import dpdata.plugins
from dpdata.amber.mask import load_param_file, pick_by_amber_mask
from dpdata.data_type import Axis, DataError, DataType, get_data_types
from dpdata.driver import Driver, Minimizer
from dpdata.format import Format
from dpdata.plugin import Plugin
from dpdata.utils import (
    add_atom_names,
    elements_index_map,
    remove_pbc,
    sort_atom_names,
    utf8len,
)


def load_format(fmt):
    fmt = fmt.lower()
    formats = Format.get_formats()
    if fmt in formats:
        return formats[fmt]()
    raise NotImplementedError(
        "Unsupported data format {}. Supported formats: {}".format(
            fmt, " ".join(formats)
        )
    )


class System(MSONable):
    """The data System.

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

    Attributes
    ----------
    DTYPES : tuple[DataType]
        data types of this class
    """

    DTYPES = (
        DataType("atom_numbs", list, (Axis.NTYPES,)),
        DataType("atom_names", list, (Axis.NTYPES,)),
        DataType("atom_types", np.ndarray, (Axis.NATOMS,)),
        DataType("orig", np.ndarray, (3,)),
        DataType("cells", np.ndarray, (Axis.NFRAMES, 3, 3)),
        DataType("coords", np.ndarray, (Axis.NFRAMES, Axis.NATOMS, 3)),
        DataType(
            "real_atom_types", np.ndarray, (Axis.NFRAMES, Axis.NATOMS), required=False
        ),
        DataType("real_atom_names", list, (Axis.NTYPES,), required=False),
        DataType("nopbc", bool, required=False),
    )

    def __init__(
        self,
        file_name=None,
        fmt="auto",
        type_map=None,
        begin=0,
        step=1,
        data=None,
        convergence_check=True,
        **kwargs,
    ):
        """Constructor.

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
                - ``vasp/contcar``: vasp contcar
                - ``vasp/string``: vasp string
                - ``vasp/outcar``: vasp outcar
                - ``vasp/xml``: vasp xml
                - ``qe/cp/traj``: Quantum Espresso CP trajectory files. should have: file_name+'.in' and file_name+'.pos'
                - ``qe/pw/scf``: Quantum Espresso PW single point calculations. Both input and output files are required. If file_name is a string, it denotes the output file name. Input file name is obtained by replacing 'out' by 'in' from file_name. Or file_name is a list, with the first element being the input file name and the second element being the output filename.
                - ``abacus/scf``: ABACUS pw/lcao scf. The directory containing INPUT file is required.
                - ``abacus/md``: ABACUS pw/lcao MD. The directory containing INPUT file is required.
                - ``abacus/relax``: ABACUS pw/lcao relax or cell-relax. The directory containing INPUT file is required.
                - ``abacus/stru``: abacus stru
                - ``abacus/lcao/scf``: abacus lcao scf
                - ``abacus/pw/scf``: abacus pw scf
                - ``abacus/lcao/md``: abacus lcao md
                - ``abacus/pw/md``: abacus pw md
                - ``abacus/lcao/relax``: abacus lcao relax
                - ``abacus/pw/relax``: abacus pw relax
                - ``siesta/output``: siesta SCF output file
                - ``siesta/aimd_output``: siesta aimd output file
                - ``pwmat/atom.config``: pwmat atom.config
                - ``pwmat/movement``: pwmat movement
                - ``pwmat/output``: pwmat output
                - ``pwmat/mlmd``: pwmat mlmd
                - ``pwmat/final.config``: pwmat final.config
                - ``quip/gap/xyz_file``: quip gap xyz_file
                - ``quip/gap/xyz``: quip gap xyz
                - ``fhi_aims/output``: fhi_aims output
                - ``fhi_aims/md``: fhi_aims md
                - ``fhi_aims/scf``: fhi_aims scf
                - ``pymatgen/structure``: pymatgen structure
                - ``pymatgen/molecule``: pymatgen molecule
                - ``pymatgen/computedstructureentry``: pymatgen computedstructureentry
                - ``amber/md``: amber md
                - ``sqm/out``: sqm out
                - ``sqm/in``: sqm in
                - ``ase/structure``: ase structure
                - ``gaussian/log``: gaussian log
                - ``gaussian/md``: gaussian md
                - ``gaussian/gjf``: gaussian gjf
                - ``deepmd/comp``: deepmd comp
                - ``deepmd/hdf5``: deepmd hdf5
                - ``gromacs/gro``: gromacs gro
                - ``cp2k/aimd_output``: cp2k aimd_output
                - ``cp2k/output``: cp2k output
        type_map : list of str
            Needed by formats lammps/lmp and lammps/dump. Maps atom type to name. The atom with type `ii` is mapped to `type_map[ii]`.
            If not provided the atom names are assigned to `'Type_1'`, `'Type_2'`, `'Type_3'`...
        begin : int
            The beginning frame when loading MD trajectory.
        step : int
            The number of skipped frames when loading MD trajectory.
        data : dict
            The raw data of System class.
        convergence_check : boolean
            Whether to request a convergence check.
        **kwargs : dict
            other parameters
        """
        self.data = {}
        self.data["atom_numbs"] = []
        self.data["atom_names"] = []
        self.data["atom_types"] = []
        self.data["orig"] = np.array([0, 0, 0])
        self.data["cells"] = []
        self.data["coords"] = []

        if data:
            self.data = data
            self.check_data()
            return
        if file_name is None:
            return
        self.from_fmt(
            file_name,
            fmt,
            type_map=type_map,
            begin=begin,
            step=step,
            convergence_check=convergence_check,
            **kwargs,
        )

        if type_map is not None:
            self.apply_type_map(type_map)

    def check_data(self):
        """Check if data is correct.

        Raises
        ------
        DataError
            if data is not correct
        """
        if not isinstance(self.data, dict):
            raise DataError("data is not a dict!")
        for dd in self.DTYPES:
            dd.check(self)
        if sum(self.get_atom_numbs()) != self.get_natoms():
            raise DataError(
                "Sum of atom_numbs (%d) is not equal to natoms (%d)."
                % (sum(self.get_atom_numbs()), self.get_natoms())
            )

    post_funcs = Plugin()

    def from_fmt(self, file_name, fmt="auto", **kwargs):
        fmt = fmt.lower()
        if fmt == "auto":
            fmt = os.path.basename(file_name).split(".")[-1].lower()
        return self.from_fmt_obj(load_format(fmt), file_name, **kwargs)

    def from_fmt_obj(self, fmtobj, file_name, **kwargs):
        data = fmtobj.from_system(file_name, **kwargs)
        if data:
            if isinstance(data, (list, tuple)):
                for dd in data:
                    self.append(System(data=dd))
            else:
                self.data = {**self.data, **data}
                self.check_data()
            if hasattr(fmtobj.from_system, "post_func"):
                for post_f in fmtobj.from_system.post_func:
                    self.post_funcs.get_plugin(post_f)(self)
        return self

    def to(self, fmt: str, *args, **kwargs) -> "System":
        """Dump systems to the specific format.

        Parameters
        ----------
        fmt : str
            format
        *args
            arguments
        **kwargs
            keyword arguments

        Returns
        -------
        System
            self
        """
        return self.to_fmt_obj(load_format(fmt), *args, **kwargs)

    def to_fmt_obj(self, fmtobj, *args, **kwargs):
        return fmtobj.to_system(self.data, *args, **kwargs)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        ret = "Data Summary"
        ret += "\nUnlabeled System"
        ret += "\n-------------------"
        ret += "\nFrame Numbers     : %d" % self.get_nframes()
        ret += "\nAtom Numbers      : %d" % self.get_natoms()
        ret += "\nElement List      :"
        ret += "\n-------------------"
        ret += "\n" + "  ".join(map(str, self.get_atom_names()))
        ret += "\n" + "  ".join(map(str, self.get_atom_numbs()))
        return ret

    def __getitem__(self, key):
        """Returns proerty stored in System by key or by idx."""
        if isinstance(key, (int, slice, list, np.ndarray)):
            return self.sub_system(key)
        return self.data[key]

    def __len__(self):
        """Returns number of frames in the system."""
        return self.get_nframes()

    def __add__(self, others):
        """Magic method "+" operation."""
        self_copy = self.copy()
        if isinstance(others, System):
            other_copy = others.copy()
            self_copy.append(other_copy)
        elif isinstance(others, list):
            for ii in others:
                assert isinstance(ii, System)
                ii_copy = ii.copy()
                self_copy.append(ii_copy)
        else:
            raise RuntimeError("Unspported data structure")
        return self.__class__.from_dict({"data": self_copy.data})

    def dump(self, filename, indent=4):
        """Dump .json or .yaml file."""
        dumpfn(self.as_dict(), filename, indent=indent)

    def map_atom_types(self, type_map=None) -> np.ndarray:
        """Map the atom types of the system.

        Parameters
        ----------
        type_map
            dict :  {"H":0,"O":1}
            or list  ["H","C","O","N"]
            The map between elements and index
            if no map_dict is given, index will
            be set according to atomic number

        Returns
        -------
        new_atom_types : np.ndarray
            The mapped atom types
        """
        if isinstance(type_map, dict) or type_map is None:
            pass
        elif isinstance(type_map, list):
            type_map = dict(zip(type_map, range(len(type_map))))
        else:
            raise RuntimeError("Unknown format")

        if type_map is None:
            type_map = elements_index_map(self.get_atom_names().copy(), standard=True)

        _set1 = set(self.get_atom_names())
        _set2 = set(list(type_map.keys()))
        assert _set1.issubset(_set2)

        atom_types_list = []
        for name, numb in zip(self.get_atom_names(), self.get_atom_numbs()):
            atom_types_list.extend([name] * numb)
        new_atom_types = np.array([type_map[ii] for ii in atom_types_list], dtype=int)

        return new_atom_types

    @staticmethod
    def load(filename):
        """Rebuild System obj. from .json or .yaml file."""
        return loadfn(filename)

    def as_dict(self):
        """Returns data dict of System instance."""
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "data": self.data,
        }
        return d

    def get_atom_names(self):
        """Returns name of atoms."""
        return self.data["atom_names"]

    def get_atom_types(self):
        """Returns type of atoms."""
        return self.data["atom_types"]

    def get_atom_numbs(self):
        """Returns number of atoms."""
        return self.data["atom_numbs"]

    def get_nframes(self):
        """Returns number of frames in the system."""
        return len(self.data["cells"])

    def get_natoms(self):
        """Returns total number of atoms in the system."""
        return len(self.data["atom_types"])

    def get_ntypes(self) -> int:
        """Returns total number of atom types in the system."""
        return len(self.data["atom_names"])

    def copy(self):
        """Returns a copy of the system."""
        return self.__class__.from_dict({"data": deepcopy(self.data)})

    def sub_system(self, f_idx):
        """Construct a subsystem from the system.

        Parameters
        ----------
        f_idx : int or index
            Which frame to use in the subsystem

        Returns
        -------
        sub_system : System
            The subsystem
        """
        tmp = self.__class__()
        # convert int to array_like
        if isinstance(f_idx, (int, np.int64)):
            f_idx = np.array([f_idx])
        for tt in self.DTYPES:
            if tt.name not in self.data:
                # skip optional data
                continue
            if tt.shape is not None and Axis.NFRAMES in tt.shape:
                axis_nframes = tt.shape.index(Axis.NFRAMES)
                new_shape = [slice(None) for _ in self.data[tt.name].shape]
                new_shape[axis_nframes] = f_idx
                tmp.data[tt.name] = self.data[tt.name][tuple(new_shape)]
            else:
                # keep the original data
                tmp.data[tt.name] = self.data[tt.name]
        return tmp

    def append(self, system):
        """Append a system to this system.

        Parameters
        ----------
        system : System
            The system to append
        """
        if not len(system.data["atom_numbs"]):
            # skip if the system to append is non-converged
            return False
        elif not len(self.data["atom_numbs"]):
            # this system is non-converged but the system to append is converged
            self.data = system.data.copy()
            return False
        if system.uniq_formula != self.uniq_formula:
            raise RuntimeError(
                f"systems with inconsistent formula could not be append: {self.uniq_formula} v.s. {system.uniq_formula}"
            )
        if system.data["atom_names"] != self.data["atom_names"]:
            # prevent original system to be modified
            system = system.copy()
            # allow to append a system with different atom_names order
            system.sort_atom_names()
            self.sort_atom_names()
        if (system.data["atom_types"] != self.data["atom_types"]).any():
            # prevent original system to be modified
            system = system.copy()
            # allow to append a system with different atom_types order
            system.sort_atom_types()
            self.sort_atom_types()
        for ii in ["atom_numbs", "atom_names"]:
            assert system.data[ii] == self.data[ii]
        for ii in ["atom_types", "orig"]:
            eq = [v1 == v2 for v1, v2 in zip(system.data[ii], self.data[ii])]
            assert all(eq)
        for tt in self.DTYPES:
            # check if the first shape is nframes
            if tt.shape is not None and Axis.NFRAMES in tt.shape:
                if tt.name not in self.data and tt.name in system.data:
                    raise RuntimeError("system has %s, but this does not" % tt.name)
                elif tt.name in self.data and tt.name not in system.data:
                    raise RuntimeError("this has %s, but system does not" % tt.name)
                elif tt.name not in self.data and tt.name not in system.data:
                    # skip if both not exist
                    continue
                # concat any data in nframes axis
                axis_nframes = tt.shape.index(Axis.NFRAMES)
                self.data[tt.name] = np.concatenate(
                    (self.data[tt.name], system[tt.name]), axis=axis_nframes
                )
        if self.nopbc and not system.nopbc:
            # appended system uses PBC, cancel nopbc
            self.data["nopbc"] = False
        return True

    def convert_to_mixed_type(self, type_map=None):
        """Convert the data dict to mixed type format structure, in order to append systems
        with different formula but the same number of atoms. Change the 'atom_names' to
        one placeholder type 'MIXED_TOKEN' and add 'real_atom_types' to store the real type
        vectors according to the given type_map.

        Parameters
        ----------
        type_map : list
            type_map
        """
        if "real_atom_types" in self.data.keys():
            return
        if type_map is None:
            type_map = self.get_atom_names()
        type_index = [type_map.index(i) for i in self.data["atom_names"]]
        frames = self.get_nframes()
        self.data["real_atom_types"] = np.tile(
            np.array([type_index[i] for i in self.data["atom_types"]]), [frames, 1]
        )
        self.data["real_atom_names"] = type_map
        natoms = self.get_natoms()
        self.data["atom_types"] = np.zeros((natoms,), dtype=int)
        self.data["atom_numbs"] = [natoms]
        self.data["atom_names"] = ["MIXED_TOKEN"]

    def sort_atom_names(self, type_map=None):
        """Sort atom_names of the system and reorder atom_numbs and atom_types accoarding
        to atom_names. If type_map is not given, atom_names will be sorted by
        alphabetical order. If type_map is given, atom_names will be type_map.

        Parameters
        ----------
        type_map : list
            type_map
        """
        self.data = sort_atom_names(self.data, type_map=type_map)

    def check_type_map(self, type_map):
        """Assign atom_names to type_map if type_map is given and different from
        atom_names.

        Parameters
        ----------
        type_map : list
            type_map
        """
        if type_map is not None and type_map != self.data["atom_names"]:
            self.sort_atom_names(type_map=type_map)

    def apply_type_map(self, type_map):
        if type_map is not None and isinstance(type_map, list):
            self.check_type_map(type_map)
        else:
            raise RuntimeError("invalid type map, cannot be applied")

    def sort_atom_types(self) -> np.ndarray:
        """Sort atom types.

        Returns
        -------
        idx : np.ndarray
            new atom index in the Axis.NATOMS
        """
        idx = np.argsort(self.data["atom_types"], kind="stable")
        for tt in self.DTYPES:
            if tt.name not in self.data:
                # skip optional data
                continue
            if tt.shape is not None and Axis.NATOMS in tt.shape:
                axis_natoms = tt.shape.index(Axis.NATOMS)
                new_shape = [slice(None) for _ in self.data[tt.name].shape]
                new_shape[axis_natoms] = idx
                self.data[tt.name] = self.data[tt.name][tuple(new_shape)]
        return idx

    @property
    def formula(self):
        """Return the formula of this system, like C3H5O2."""
        return "".join(
            [
                f"{symbol}{numb}"
                for symbol, numb in zip(
                    self.data["atom_names"], self.data["atom_numbs"]
                )
            ]
        )

    @property
    def uniq_formula(self):
        """Return the uniq_formula of this system.
        The uniq_formula sort the elements in formula by names.
        Systems with the same uniq_formula can be append together.
        """
        return "".join(
            [
                f"{symbol}{numb}"
                for symbol, numb in sorted(
                    zip(self.data["atom_names"], self.data["atom_numbs"])
                )
            ]
        )

    @property
    def short_formula(self) -> str:
        """Return the short formula of this system. Elements with zero number
        will be removed.
        """
        return "".join(
            [
                f"{symbol}{numb}"
                for symbol, numb in zip(
                    self.data["atom_names"], self.data["atom_numbs"]
                )
                if numb
            ]
        )

    @property
    def formula_hash(self) -> str:
        """Return the hash of the formula of this system."""
        return hashlib.sha256(self.formula.encode("utf-8")).hexdigest()

    @property
    def short_name(self) -> str:
        """Return the short name of this system (no more than 255 bytes), in
        the following order:
            - formula
            - short_formula
            - formula_hash.
        """
        formula = self.formula
        if utf8len(formula) <= 255:
            return formula
        short_formula = self.short_formula
        if utf8len(short_formula) <= 255:
            return short_formula
        return self.formula_hash

    def extend(self, systems):
        """Extend a system list to this system.

        Parameters
        ----------
        systems : [System1, System2, System3 ]
            The list to extend
        """
        for system in systems:
            self.append(system.copy())

    def apply_pbc(self):
        """Append periodic boundary condition."""
        ncoord = dpdata.md.pbc.dir_coord(self.data["coords"], self.data["cells"])
        ncoord = ncoord % 1
        self.data["coords"] = np.matmul(ncoord, self.data["cells"])

    @post_funcs.register("remove_pbc")
    def remove_pbc(self, protect_layer=9):
        """This method does NOT delete the definition of the cells, it
        (1) revises the cell to a cubic cell and ensures that the cell
        boundary to any atom in the system is no less than `protect_layer`
        (2) translates the system such that the center-of-geometry of the system
        locates at the center of the cell.

        Parameters
        ----------
        protect_layer : the protect layer between the atoms and the cell
            boundary
        """
        assert protect_layer >= 0, "the protect_layer should be no less than 0"
        remove_pbc(self.data, protect_layer)

    def affine_map(self, trans, f_idx=0):
        assert np.linalg.det(trans) != 0
        self.data["cells"][f_idx] = np.matmul(self.data["cells"][f_idx], trans)
        self.data["coords"][f_idx] = np.matmul(self.data["coords"][f_idx], trans)

    @post_funcs.register("shift_orig_zero")
    def _shift_orig_zero(self):
        for ff in self.data["coords"]:
            for ii in ff:
                ii = ii - self.data["orig"]
        self.data["orig"] = self.data["orig"] - self.data["orig"]
        assert (np.zeros([3]) == self.data["orig"]).all()

    @post_funcs.register("rot_lower_triangular")
    def rot_lower_triangular(self):
        for ii in range(self.get_nframes()):
            self.rot_frame_lower_triangular(ii)

    def rot_frame_lower_triangular(self, f_idx=0):
        qq, rr = np.linalg.qr(self.data["cells"][f_idx].T)
        if np.linalg.det(qq) < 0:
            qq = -qq
            rr = -rr
        self.affine_map(qq, f_idx=f_idx)
        rot = np.eye(3)
        if self.data["cells"][f_idx][0][0] < 0:
            rot[0][0] = -1
        if self.data["cells"][f_idx][1][1] < 0:
            rot[1][1] = -1
        if self.data["cells"][f_idx][2][2] < 0:
            rot[2][2] = -1
        assert np.linalg.det(rot) == 1
        self.affine_map(rot, f_idx=f_idx)
        return np.matmul(qq, rot)

    def add_atom_names(self, atom_names):
        """Add atom_names that do not exist."""
        self.data = add_atom_names(self.data, atom_names)

    def replicate(self, ncopy):
        """Replicate the each frame  in the system in 3 dimensions.
        Each frame in the system will become a supercell.

        Parameters
        ----------
        ncopy
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
        if len(ncopy) != 3:
            raise RuntimeError("ncopy must be a list or tuple with 3 int")
        for ii in ncopy:
            if not isinstance(ii, int):
                raise RuntimeError("ncopy must be a list or tuple must with 3 int")

        tmp = System()
        nframes = self.get_nframes()
        data = self.data
        tmp.data["atom_names"] = list(np.copy(data["atom_names"]))
        tmp.data["atom_numbs"] = list(
            np.array(np.copy(data["atom_numbs"])) * np.prod(ncopy)
        )
        tmp.data["atom_types"] = np.sort(
            np.tile(np.copy(data["atom_types"]), np.prod(ncopy)), kind="stable"
        )
        tmp.data["cells"] = np.copy(data["cells"])
        for ii in range(3):
            tmp.data["cells"][:, ii, :] *= ncopy[ii]
        tmp.data["coords"] = np.tile(np.copy(data["coords"]), tuple(ncopy) + (1, 1, 1))

        for xx in range(ncopy[0]):
            for yy in range(ncopy[1]):
                for zz in range(ncopy[2]):
                    tmp.data["coords"][xx, yy, zz, :, :, :] += (
                        xx * np.reshape(data["cells"][:, 0, :], [-1, 1, 3])
                        + yy * np.reshape(data["cells"][:, 1, :], [-1, 1, 3])
                        + zz * np.reshape(data["cells"][:, 2, :], [-1, 1, 3])
                    )
        tmp.data["coords"] = np.reshape(
            np.transpose(tmp.data["coords"], [3, 4, 0, 1, 2, 5]), (nframes, -1, 3)
        )
        return tmp

    def replace(self, initial_atom_type, end_atom_type, replace_num):
        if type(self) is not dpdata.System:
            raise RuntimeError(
                "Must use method replace() of the instance of class dpdata.System"
            )
        if not isinstance(replace_num, int):
            raise ValueError(f"replace_num must be a integer. Now is {replace_num}")
        if replace_num <= 0:
            raise ValueError(f"replace_num must be larger than 0.Now is {replace_num}")

        try:
            initial_atom_index = self.data["atom_names"].index(initial_atom_type)
        except ValueError as e:
            raise ValueError(
                "atom_type  {initial_atom_type}   not in {atom_names}".format(
                    initial_atom_type=initial_atom_type,
                    atom_names=self.data["atom_names"],
                )
            )
        max_replace_num = self.data["atom_numbs"][initial_atom_index]

        if replace_num > max_replace_num:
            raise RuntimeError(
                f"not enough {initial_atom_type} atom, only {max_replace_num} available, less than {replace_num}.Please check."
            )

        may_replace_indices = [
            i for i, x in enumerate(self.data["atom_types"]) if x == initial_atom_index
        ]
        to_replace_indices = np.random.choice(
            may_replace_indices, size=replace_num, replace=False
        )

        if end_atom_type not in self.data["atom_names"]:
            self.data["atom_names"].append(end_atom_type)
            self.data["atom_numbs"].append(0)

        end_atom_index = self.data["atom_names"].index(end_atom_type)
        for ii in to_replace_indices:
            self.data["atom_types"][ii] = end_atom_index
        self.data["atom_numbs"][initial_atom_index] -= replace_num
        self.data["atom_numbs"][end_atom_index] += replace_num
        self.sort_atom_types()

    def perturb(
        self, pert_num, cell_pert_fraction, atom_pert_distance, atom_pert_style="normal"
    ):
        """Perturb each frame in the system randomly.
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
        atom_pert_distance : float
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
        if type(self) is not dpdata.System:
            raise RuntimeError(
                f"Using method perturb() of an instance of {type(self)}. "
                f"Must use method perturb() of the instance of class dpdata.System."
            )
        perturbed_system = System()
        nframes = self.get_nframes()
        for ii in range(nframes):
            for jj in range(pert_num):
                tmp_system = self[ii].copy()
                cell_perturb_matrix = get_cell_perturb_matrix(cell_pert_fraction)
                tmp_system.data["cells"][0] = np.matmul(
                    tmp_system.data["cells"][0], cell_perturb_matrix
                )
                tmp_system.data["coords"][0] = np.matmul(
                    tmp_system.data["coords"][0], cell_perturb_matrix
                )
                for kk in range(len(tmp_system.data["coords"][0])):
                    atom_perturb_vector = get_atom_perturb_vector(
                        atom_pert_distance, atom_pert_style
                    )
                    tmp_system.data["coords"][0][kk] += atom_perturb_vector
                tmp_system.rot_lower_triangular()
                perturbed_system.append(tmp_system)
        return perturbed_system

    @property
    def nopbc(self):
        if self.data.get("nopbc", False):
            return True
        return False

    @nopbc.setter
    def nopbc(self, value):
        self.data["nopbc"] = value

    def shuffle(self):
        """Shuffle frames randomly."""
        idx = np.random.permutation(self.get_nframes())
        self.data = self.sub_system(idx).data
        return idx

    def predict(self, *args: Any, driver: str = "dp", **kwargs: Any) -> "LabeledSystem":
        """Predict energies and forces by a driver.

        Parameters
        ----------
        *args : iterable
            Arguments passing to the driver
        driver : str, default=dp
            The assigned driver. For compatibility, default is dp
        **kwargs : dict
            Other arguments passing to the driver

        Returns
        -------
        labeled_sys : LabeledSystem
            A new labeled system.

        Examples
        --------
        The default driver is DP:

        >>> labeled_sys = ori_sys.predict("frozen_model_compressed.pb")
        """
        if not isinstance(driver, Driver):
            driver = Driver.get_driver(driver)(*args, **kwargs)
        data = driver.label(self.data.copy())
        return LabeledSystem(data=data)

    def minimize(
        self, *args: Any, minimizer: Union[str, Minimizer], **kwargs: Any
    ) -> "LabeledSystem":
        """Minimize the geometry.

        Parameters
        ----------
        *args : iterable
            Arguments passing to the minimizer
        minimizer : str or Minimizer
            The assigned minimizer
        **kwargs : dict
            Other arguments passing to the minimizer

        Returns
        -------
        labeled_sys : LabeledSystem
            A new labeled system.
        """
        if not isinstance(minimizer, Minimizer):
            minimizer = Minimizer.get_minimizer(minimizer)(*args, **kwargs)
        data = minimizer.minimize(self.data.copy())
        return LabeledSystem(data=data)

    def pick_atom_idx(self, idx, nopbc=None):
        """Pick atom index.

        Parameters
        ----------
        idx : int or list or slice
            atom index
        nopbc : Boolen (default: None)
            If nopbc is True or False, set nopbc

        Returns
        -------
        new_sys: System
            new system
        """
        new_sys = self.copy()
        if isinstance(idx, (int, np.int64)):
            idx = np.array([idx])
        for tt in self.DTYPES:
            if tt.name not in self.data:
                # skip optional data
                continue
            if tt.shape is not None and Axis.NATOMS in tt.shape:
                axis_natoms = tt.shape.index(Axis.NATOMS)
                new_shape = [slice(None) for _ in self.data[tt.name].shape]
                new_shape[axis_natoms] = idx
                new_sys.data[tt.name] = self.data[tt.name][tuple(new_shape)]
        # recalculate atom_numbs according to atom_types
        atom_numbs = np.bincount(
            new_sys.data["atom_types"], minlength=len(self.get_atom_names())
        )
        new_sys.data["atom_numbs"] = list(atom_numbs)
        if nopbc is True or nopbc is False:
            new_sys.nopbc = nopbc
        return new_sys

    def remove_atom_names(self, atom_names):
        """Remove atom names and all such atoms.
        For example, you may not remove EP atoms in TIP4P/Ew water, which
        is not a real atom.
        """
        if isinstance(atom_names, str):
            atom_names = [atom_names]
        removed_atom_idx = []
        for an in atom_names:
            # get atom name idx
            idx = self.data["atom_names"].index(an)
            atom_idx = self.data["atom_types"] == idx
            removed_atom_idx.append(atom_idx)
        picked_atom_idx = ~np.any(removed_atom_idx, axis=0)
        new_sys = self.pick_atom_idx(picked_atom_idx)
        # let's remove atom_names
        # firstly, rearrange atom_names and put these atom_names in the end
        new_atom_names = list(
            [xx for xx in new_sys.data["atom_names"] if xx not in atom_names]
        )
        new_sys.sort_atom_names(type_map=new_atom_names + atom_names)
        # remove atom_names and atom_numbs
        new_sys.data["atom_names"] = new_atom_names
        new_sys.data["atom_numbs"] = new_sys.data["atom_numbs"][: len(new_atom_names)]
        return new_sys

    def pick_by_amber_mask(self, param, maskstr, pass_coords=False, nopbc=None):
        """Pick atoms by amber mask.

        Parameters
        ----------
        param : str or parmed.Structure
            filename of Amber param file or parmed.Structure
        maskstr : str
            Amber masks
        pass_coords : Boolen (default: False)
            If pass_coords is true, the function will pass coordinates and
            return a MultiSystem. Otherwise, the result is
            coordinate-independent, and the function will return System or
            LabeledSystem.
        nopbc : Boolen (default: None)
            If nopbc is True or False, set nopbc
        """
        parm = load_param_file(param)
        if pass_coords:
            ms = MultiSystems()
            for sub_s in self:
                # TODO: this can computed in pararrel
                idx = pick_by_amber_mask(parm, maskstr, sub_s["coords"][0])
                ms.append(sub_s.pick_atom_idx(idx, nopbc=nopbc))
            return ms
        else:
            idx = pick_by_amber_mask(parm, maskstr)
            return self.pick_atom_idx(idx, nopbc=nopbc)

    @classmethod
    def register_data_type(cls, *data_type: Tuple[DataType]):
        """Register data type.

        Parameters
        ----------
        *data_type : tuple[DataType]
            data type to be regiestered
        """
        all_dtypes = cls.DTYPES + tuple(data_type)
        dtypes_dict = {}
        for dt in all_dtypes:
            if dt.name in dtypes_dict:
                warnings.warn(
                    f"Data type {dt.name} is registered twice; only the newly registered one will be used.",
                    UserWarning,
                )
            dtypes_dict[dt.name] = dt
        cls.DTYPES = tuple(dtypes_dict.values())


def get_cell_perturb_matrix(cell_pert_fraction):
    if cell_pert_fraction < 0:
        raise RuntimeError("cell_pert_fraction can not be negative")
    e0 = np.random.rand(6)
    e = e0 * 2 * cell_pert_fraction - cell_pert_fraction
    cell_pert_matrix = np.array(
        [
            [1 + e[0], 0.5 * e[5], 0.5 * e[4]],
            [0.5 * e[5], 1 + e[1], 0.5 * e[3]],
            [0.5 * e[4], 0.5 * e[3], 1 + e[2]],
        ]
    )
    return cell_pert_matrix


def get_atom_perturb_vector(atom_pert_distance, atom_pert_style="normal"):
    random_vector = None
    if atom_pert_distance < 0:
        raise RuntimeError("atom_pert_distance can not be negative")

    if atom_pert_style == "normal":
        e = np.random.randn(3)
        random_vector = (atom_pert_distance / np.sqrt(3)) * e
    elif atom_pert_style == "uniform":
        e = np.random.randn(3)
        while np.linalg.norm(e) < 0.1:
            e = np.random.randn(3)
        random_unit_vector = e / np.linalg.norm(e)
        v0 = np.random.rand(1)
        v = np.power(v0, 1 / 3)
        random_vector = atom_pert_distance * v * random_unit_vector
    elif atom_pert_style == "const":
        e = np.random.randn(3)
        while np.linalg.norm(e) < 0.1:
            e = np.random.randn(3)
        random_unit_vector = e / np.linalg.norm(e)
        random_vector = atom_pert_distance * random_unit_vector
    else:
        raise RuntimeError(f"unsupported options atom_pert_style={atom_pert_style}")
    return random_vector


class LabeledSystem(System):
    """The labeled data System.

    For example, a labeled water system named `d_example` has two molecules (6 atoms) and `nframes` frames. The labels can be accessed by
        - `d_example['energies']` : a numpy array of size nframes
        - `d_example['forces']` : a numpy array of size nframes x 6 x 3
        - `d_example['virials']` : optional, a numpy array of size nframes x 3 x 3

    It is noted that
        - The order of frames stored in `'energies'`, `'forces'` and `'virials'` should be consistent with `'atom_types'`, `'cells'` and `'coords'`.
        - The order of atoms in **every** frame of `'forces'` should be consistent with `'coords'` and `'atom_types'`.

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
                - ``cp2k/aimd_output``: cp2k aimd output  dir(contains *pos*.xyz and *.log file); optional `restart=True` if it is a cp2k restarted task.
                - ``pwmat/movement``: pwmat md output file
                - ``pwmat/out.mlmd``: pwmat scf output file

        type_map : list of str
            Maps atom type to name. The atom with type `ii` is mapped to `type_map[ii]`.
            If not provided the atom names are assigned to `'Type_1'`, `'Type_2'`, `'Type_3'`...
        begin : int
            The beginning frame when loading MD trajectory.
        step : int
            The number of skipped frames when loading MD trajectory.
    """

    DTYPES = System.DTYPES + (
        DataType("energies", np.ndarray, (Axis.NFRAMES,)),
        DataType("forces", np.ndarray, (Axis.NFRAMES, Axis.NATOMS, 3)),
        DataType("virials", np.ndarray, (Axis.NFRAMES, 3, 3), required=False),
        DataType("atom_pref", np.ndarray, (Axis.NFRAMES, Axis.NATOMS), required=False),
    )

    post_funcs = Plugin() + System.post_funcs

    def from_fmt_obj(self, fmtobj, file_name, **kwargs):
        data = fmtobj.from_labeled_system(file_name, **kwargs)
        if data:
            if isinstance(data, (list, tuple)):
                for dd in data:
                    self.append(LabeledSystem(data=dd))
            else:
                self.data = {**self.data, **data}
                self.check_data()
            if hasattr(fmtobj.from_labeled_system, "post_func"):
                for post_f in fmtobj.from_labeled_system.post_func:
                    self.post_funcs.get_plugin(post_f)(self)
        return self

    def to_fmt_obj(self, fmtobj, *args, **kwargs):
        return fmtobj.to_labeled_system(self.data, *args, **kwargs)

    def __str__(self):
        ret = "Data Summary"
        ret += "\nLabeled System"
        ret += "\n-------------------"
        ret += "\nFrame Numbers      : %d" % self.get_nframes()
        ret += "\nAtom Numbers       : %d" % self.get_natoms()
        status = "Yes" if self.has_virial() else "No"
        ret += "\nIncluding Virials  : %s" % status
        ret += "\nElement List       :"
        ret += "\n-------------------"
        ret += "\n" + "  ".join(map(str, self.get_atom_names()))
        ret += "\n" + "  ".join(map(str, self.get_atom_numbs()))
        return ret

    def __add__(self, others):
        """Magic method "+" operation."""
        self_copy = self.copy()
        if isinstance(others, LabeledSystem):
            other_copy = others.copy()
            self_copy.append(other_copy)
        elif isinstance(others, list):
            for ii in others:
                assert isinstance(ii, LabeledSystem)
                ii_copy = ii.copy()
                self_copy.append(ii_copy)
        else:
            raise RuntimeError("Unspported data structure")
        return self.__class__.from_dict({"data": self_copy.data})

    def has_virial(self):
        # return ('virials' in self.data) and (len(self.data['virials']) > 0)
        return "virials" in self.data

    def affine_map_fv(self, trans, f_idx):
        assert np.linalg.det(trans) != 0
        self.data["forces"][f_idx] = np.matmul(self.data["forces"][f_idx], trans)
        if self.has_virial():
            self.data["virials"][f_idx] = np.matmul(
                trans.T, np.matmul(self.data["virials"][f_idx], trans)
            )

    def rot_frame_lower_triangular(self, f_idx=0):
        trans = System.rot_frame_lower_triangular(self, f_idx=f_idx)
        self.affine_map_fv(trans, f_idx=f_idx)
        return trans

    def correction(self, hl_sys):
        """Get energy and force correction between self and a high-level LabeledSystem.
        The self's coordinates will be kept, but energy and forces will be replaced by
        the correction between these two systems.

        Note: The function will not check whether coordinates and elements of two systems
        are the same. The user should make sure by itself.

        Parameters
        ----------
        hl_sys : LabeledSystem
            high-level LabeledSystem

        Returns
        -------
        corrected_sys: LabeledSystem
            Corrected LabeledSystem
        """
        if not isinstance(hl_sys, LabeledSystem):
            raise RuntimeError("high_sys should be LabeledSystem")
        corrected_sys = self.copy()
        corrected_sys.data["energies"] = hl_sys.data["energies"] - self.data["energies"]
        corrected_sys.data["forces"] = hl_sys.data["forces"] - self.data["forces"]
        if "virials" in self.data and "virials" in hl_sys.data:
            corrected_sys.data["virials"] = (
                hl_sys.data["virials"] - self.data["virials"]
            )
        return corrected_sys

    def remove_outlier(self, threshold: float = 8.0) -> "LabeledSystem":
        r"""Remove outlier frames from the system.

        Remove the frames whose energies satisfy the condition

        .. math::

            \frac{\left \| E - \bar{E} \right \|}{\sigma(E)} \geq \text{threshold}

        where :math:`\bar{E}` and :math:`\sigma(E)` are the mean and standard deviation
        of the energies in the system.

        Parameters
        ----------
        threshold : float
            The threshold of outlier detection. The default value is 8.0.

        Returns
        -------
        LabeledSystem
            The system without outlier frames.

        References
        ----------
        .. [1] Gao, X.; Ramezanghorbani, F.; Isayev, O.; Smith, J. S.;
           Roitberg, A. E. TorchANI: A Free and Open Source PyTorch-Based
           Deep Learning Implementation of the ANI Neural Network
           Potentials. J. Chem. Inf. Model. 2020, 60, 3408-3415.
        .. [2] Zeng, J.; Tao, Y.; Giese, T. J.; York, D. M.. QDπ: A Quantum
           Deep Potential Interaction Model for Drug Discovery. J. Comput.
           Chem. 2023, 19, 1261-1275.
        """
        energies = self.data["energies"]
        std = np.std(energies)
        if np.isclose(std, 0.0):
            return self.copy()
        idx = np.abs(energies - np.mean(energies)) / std < threshold
        return self.sub_system(idx)


class MultiSystems:
    """A set containing several systems."""

    def __init__(self, *systems, type_map=None):
        """Parameters
        ----------
        *systems : System
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

    def from_fmt_obj(self, fmtobj, directory, labeled=True, **kwargs):
        if not isinstance(fmtobj, dpdata.plugins.deepmd.DeePMDMixedFormat):
            for dd in fmtobj.from_multi_systems(directory, **kwargs):
                if labeled:
                    system = LabeledSystem().from_fmt_obj(fmtobj, dd, **kwargs)
                else:
                    system = System().from_fmt_obj(fmtobj, dd, **kwargs)
                system.sort_atom_names()
                self.append(system)
            return self
        else:
            system_list = []
            for dd in fmtobj.from_multi_systems(directory, **kwargs):
                if labeled:
                    data_list = fmtobj.from_labeled_system_mix(dd, **kwargs)
                    for data_item in data_list:
                        system_list.append(LabeledSystem(data=data_item, **kwargs))
                else:
                    data_list = fmtobj.from_system_mix(dd, **kwargs)
                    for data_item in data_list:
                        system_list.append(System(data=data_item, **kwargs))
            self.append(*system_list)
            return self

    def to_fmt_obj(self, fmtobj, directory, *args, **kwargs):
        if not isinstance(fmtobj, dpdata.plugins.deepmd.DeePMDMixedFormat):
            for fn, ss in zip(
                fmtobj.to_multi_systems(
                    [ss.short_name for ss in self.systems.values()], directory, **kwargs
                ),
                self.systems.values(),
            ):
                ss.to_fmt_obj(fmtobj, fn, *args, **kwargs)
        else:
            mixed_systems = fmtobj.mix_system(
                *list(self.systems.values()), type_map=self.atom_names, **kwargs
            )
            for fn in mixed_systems:
                mixed_systems[fn].to_fmt_obj(
                    fmtobj, os.path.join(directory, fn), *args, **kwargs
                )
        return self

    def to(self, fmt: str, *args, **kwargs) -> "MultiSystems":
        """Dump systems to the specific format.

        Parameters
        ----------
        fmt : str
            format
        *args : list
            arguments
        **kwargs : dict
            keyword arguments

        Returns
        -------
        MultiSystems
            self
        """
        return self.to_fmt_obj(load_format(fmt), *args, **kwargs)

    def __getitem__(self, key):
        """Returns proerty stored in System by key or by idx."""
        if isinstance(key, int):
            return list(self.systems.values())[key]
        return self.systems[key]

    def __len__(self):
        return len(self.systems)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "MultiSystems ({} systems containing {} frames)".format(
            len(self.systems), self.get_nframes()
        )

    def __add__(self, others):
        """Magic method "+" operation."""
        self_copy = deepcopy(self)
        if isinstance(others, System) or isinstance(others, MultiSystems):
            return self.__class__(self, others)
        elif isinstance(others, list):
            return self.__class__(self, *others)
        raise RuntimeError("Unspported data structure")

    @classmethod
    def from_file(cls, file_name, fmt, **kwargs):
        multi_systems = cls()
        multi_systems.load_systems_from_file(file_name=file_name, fmt=fmt, **kwargs)
        return multi_systems

    @classmethod
    def from_dir(cls, dir_name, file_name, fmt="auto", type_map=None):
        multi_systems = cls()
        target_file_list = sorted(
            glob.glob(f"./{dir_name}/**/{file_name}", recursive=True)
        )
        for target_file in target_file_list:
            multi_systems.append(
                LabeledSystem(file_name=target_file, fmt=fmt, type_map=type_map)
            )
        return multi_systems

    def load_systems_from_file(self, file_name=None, fmt=None, **kwargs):
        fmt = fmt.lower()
        return self.from_fmt_obj(load_format(fmt), file_name, **kwargs)

    def get_nframes(self):
        """Returns number of frames in all systems."""
        return sum(len(system) for system in self.systems.values())

    def append(self, *systems):
        """Append systems or MultiSystems to systems.

        Parameters
        ----------
        *systems : System
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
        # prevent changing the original system
        system = system.copy()
        self.check_atom_names(system)
        formula = system.formula
        if formula in self.systems:
            self.systems[formula].append(system)
        else:
            self.systems[formula] = system.copy()

    def check_atom_names(self, system):
        """Make atom_names in all systems equal, prevent inconsistent atom_types."""
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

    def predict(self, *args: Any, driver="dp", **kwargs: Any) -> "MultiSystems":
        """Predict energies and forces by a driver.

        Parameters
        ----------
        *args : iterable
            Arguments passing to the driver
        driver : str, default=dp
            The assigned driver. For compatibility, default is dp
        **kwargs : dict
            Other arguments passing to the driver

        Returns
        -------
        MultiSystems
            A new labeled MultiSystems.
        """
        if not isinstance(driver, Driver):
            driver = Driver.get_driver(driver)(*args, **kwargs)
        new_multisystems = dpdata.MultiSystems(type_map=self.atom_names)
        for ss in self:
            new_multisystems.append(ss.predict(*args, driver=driver, **kwargs))
        return new_multisystems

    def minimize(
        self, *args: Any, minimizer: Union[str, Minimizer], **kwargs: Any
    ) -> "MultiSystems":
        """Minimize geometry by a minimizer.

        Parameters
        ----------
        *args : iterable
            Arguments passing to the minimizer
        minimizer : str or Minimizer
            The assigned minimizer
        **kwargs : dict
            Other arguments passing to the minimizer

        Returns
        -------
        MultiSystems
            A new labeled MultiSystems.

        Examples
        --------
        Minimize a system using ASE BFGS along with a DP driver:

        >>> from dpdata.driver import Driver
        >>> from ase.optimize import BFGS
        >>> driver = Driver.get_driver("dp")("some_model.pb")
        >>> some_system.minimize(minimizer="ase", driver=driver, optimizer=BFGS, fmax=1e-5)
        """
        if not isinstance(minimizer, Minimizer):
            minimizer = Minimizer.get_minimizer(minimizer)(*args, **kwargs)
        new_multisystems = dpdata.MultiSystems(type_map=self.atom_names)
        for ss in self:
            new_multisystems.append(ss.minimize(*args, minimizer=minimizer, **kwargs))
        return new_multisystems

    def pick_atom_idx(self, idx, nopbc=None):
        """Pick atom index.

        Parameters
        ----------
        idx : int or list or slice
            atom index
        nopbc : Boolen (default: None)
            If nopbc is True or False, set nopbc

        Returns
        -------
        new_sys: MultiSystems
            new system
        """
        new_sys = MultiSystems()
        for ss in self:
            new_sys.append(ss.pick_atom_idx(idx, nopbc=nopbc))
        return new_sys

    def correction(self, hl_sys: "MultiSystems"):
        """Get energy and force correction between self (assumed low-level) and a high-level MultiSystems.
        The self's coordinates will be kept, but energy and forces will be replaced by
        the correction between these two systems.

        Notes
        -----
        This method will not check whether coordinates and elements of two systems
        are the same. The user should make sure by itself.

        Parameters
        ----------
        hl_sys : MultiSystems
            high-level MultiSystems

        Returns
        -------
        corrected_sys : MultiSystems
            Corrected MultiSystems

        Examples
        --------
        Get correction between a low-level system and a high-level system:

        >>> low_level = dpdata.MultiSystems().from_deepmd_hdf5("low_level.hdf5")
        >>> high_level = dpdata.MultiSystems().from_deepmd_hdf5("high_level.hdf5")
        >>> corr = low_level.correction(high_lebel)
        >>> corr.to_deepmd_hdf5("corr.hdf5")
        """
        if not isinstance(hl_sys, MultiSystems):
            raise RuntimeError("high_sys should be MultiSystems")
        corrected_sys = MultiSystems(type_map=self.atom_names)
        for nn in self.systems.keys():
            ll_ss = self[nn]
            hl_ss = hl_sys[nn]
            corrected_sys.append(ll_ss.correction(hl_ss))
        return corrected_sys

    def train_test_split(
        self, test_size: Union[float, int], seed: Optional[int] = None
    ) -> Tuple["MultiSystems", "MultiSystems", Dict[str, np.ndarray]]:
        """Split systems into random train and test subsets.

        Parameters
        ----------
        test_size : float or int
            If float, should be between 0.0 and 1.0 and represent the proportion of the dataset to include in the test split.
            If int, represents the absolute number of test samples.
        seed : int, default=None
            Random seed

        Returns
        -------
        MultiSystems
            The training set
        MultiSystems
            The testing set
        Dict[str, np.ndarray]
            The bool array of training and testing sets for each system. False for training set and True for testing set.
        """
        nframes = self.get_nframes()
        if isinstance(test_size, float):
            assert 0 <= test_size <= 1
            test_size = int(np.floor(test_size * nframes))
        elif isinstance(test_size, int):
            assert 0 <= test_size <= nframes
        else:
            raise RuntimeError("test_size should be float or int")
        # get random indices
        rng = np.random.default_rng(seed=seed)
        test_idx = rng.choice(nframes, test_size, replace=False)
        select_test = np.zeros(nframes, dtype=bool)
        select_test[test_idx] = True
        select_train = np.logical_not(select_test)
        # flatten systems dict
        system_names, system_sizes = zip(
            *((kk, len(vv)) for (kk, vv) in self.systems.items())
        )
        system_idx = np.empty(len(system_sizes) + 1, dtype=int)
        system_idx[0] = 0
        np.cumsum(system_sizes, out=system_idx[1:])
        # make new systems
        train_systems = MultiSystems(type_map=self.atom_names)
        test_systems = MultiSystems(type_map=self.atom_names)
        test_system_idx = {}
        for ii, nn in enumerate(system_names):
            sub_train = self[nn][select_train[system_idx[ii] : system_idx[ii + 1]]]
            if len(sub_train):
                train_systems.append(sub_train)
            sub_test = self[nn][select_test[system_idx[ii] : system_idx[ii + 1]]]
            if len(sub_test):
                test_systems.append(sub_test)
            test_system_idx[nn] = select_test[system_idx[ii] : system_idx[ii + 1]]
        return train_systems, test_systems, test_system_idx


def get_cls_name(cls: object) -> str:
    """Returns the fully qualified name of a class, such as `np.ndarray`.

    Parameters
    ----------
    cls : object
        the class

    Returns
    -------
    str
        the fully qualified name of a class
    """
    return ".".join([cls.__module__, cls.__name__])


def add_format_methods():
    """Add format methods to System, LabeledSystem, and MultiSystems; add data types
    to System and LabeledSystem.

    Notes
    -----
    Ensure all plugins have been loaded before execuating this function!
    """
    # automatically register from/to functions for formats
    # for example, deepmd/npy will be registered as from_deepmd_npy and to_deepmd_npy
    for key, formatcls in Format.get_formats().items():
        formattedkey = key.replace("/", "_").replace(".", "")
        from_func_name = "from_" + formattedkey
        to_func_name = "to_" + formattedkey
        Format.register_from(from_func_name)(formatcls)
        Format.register_to(to_func_name)(formatcls)

    for method, formatcls in Format.get_from_methods().items():

        def get_func(ff):
            # ff is not initized when defining from_format so cannot be polluted
            def from_format(self, file_name, **kwargs):
                return self.from_fmt_obj(ff(), file_name, **kwargs)

            from_format.__doc__ = "Read data from :class:`%s` format." % (
                get_cls_name(ff)
            )
            return from_format

        setattr(System, method, get_func(formatcls))
        setattr(LabeledSystem, method, get_func(formatcls))
        setattr(MultiSystems, method, get_func(formatcls))

    for method, formatcls in Format.get_to_methods().items():

        def get_func(ff):
            def to_format(self, *args, **kwargs):
                return self.to_fmt_obj(ff(), *args, **kwargs)

            to_format.__doc__ = "Dump data to :class:`%s` format." % (get_cls_name(ff))
            return to_format

        setattr(System, method, get_func(formatcls))
        setattr(LabeledSystem, method, get_func(formatcls))
        setattr(MultiSystems, method, get_func(formatcls))

    # at this point, System.DTYPES and LabeledSystem.DTYPES has been initialized
    System.register_data_type(*get_data_types(labeled=False))
    LabeledSystem.register_data_type(*get_data_types(labeled=False))
    LabeledSystem.register_data_type(*get_data_types(labeled=True))


add_format_methods()
