import re
import os
from scipy.io import netcdf
import numpy as np
from dpdata.amber.mask import pick_by_amber_mask

kcalmol2eV= 0.04336410390059322
symbols = ['X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

energy_convert = kcalmol2eV
force_convert = energy_convert


def read_amber_traj(parm7_file, nc_file, mdfrc_file=None, mden_file = None, mdout_file = None,
        use_element_symbols=None, labeled=True,
    ):
    """The amber trajectory includes:
    * nc, NetCDF format, stores coordinates
    * mdfrc, NetCDF format, stores forces
    * mden (optional), text format, stores energies
    * mdout (optional), text format, may store energies if there is no mden_file
    * parm7, text format, stores types

    Parameters
    ----------
    parm7_file, nc_file, mdfrc_file, mden_file, mdout_file:
      filenames
    use_element_symbols: None or list or str
      If use_element_symbols is a list of atom indexes, these atoms will use element symbols
      instead of amber types. For example, a ligand will use C, H, O, N, and so on
      instead of h1, hc, o, os, and so on.
      IF use_element_symbols is str, it will be considered as Amber mask.
    """

    flag_atom_type = False
    flag_atom_numb = False
    amber_types = []
    atomic_number = []
    with open(parm7_file) as f:
        for line in f:
            if line.startswith("%FLAG"):
                flag_atom_type = line.startswith("%FLAG AMBER_ATOM_TYPE")
                flag_atom_numb = (use_element_symbols is not None) and line.startswith("%FLAG ATOMIC_NUMBER")
            elif flag_atom_type or flag_atom_numb:
                if line.startswith("%FORMAT"):
                    fmt = re.findall(r'\d+', line)
                    fmt0 = int(fmt[0])
                    fmt1 = int(fmt[1])
                else:
                    for ii in range(fmt0):
                        start_index = ii * fmt1
                        end_index = (ii + 1) * fmt1
                        if end_index >= len(line):
                            continue
                        content = line[start_index:end_index].strip()
                        if flag_atom_type:
                            amber_types.append(content)
                        elif flag_atom_numb:
                            atomic_number.append(int(content))    
    if use_element_symbols is not None:
        if isinstance(use_element_symbols, str):
            use_element_symbols = pick_by_amber_mask(parm7_file, use_element_symbols)
        for ii in use_element_symbols:
            amber_types[ii] = symbols[atomic_number[ii]]

    with netcdf.netcdf_file(nc_file, 'r') as f:
        coords = np.array(f.variables["coordinates"][:])
        cell_lengths = np.array(f.variables["cell_lengths"][:])
        cell_angles = np.array(f.variables["cell_angles"][:])
        if np.all(cell_angles > 89.99 ) and np.all(cell_angles < 90.01):
            # only support 90
            # TODO: support other angles
            shape = cell_lengths.shape
            cells = np.zeros((shape[0], 3, 3))
            for ii in range(3):
                cells[:, ii, ii] = cell_lengths[:, ii]
        else:
            raise RuntimeError("Unsupported cells")

    if labeled:
        with netcdf.netcdf_file(mdfrc_file, 'r') as f:
            forces = np.array(f.variables["forces"][:])

        # load energy from mden_file or mdout_file
        energies = []
        if mden_file is not None and os.path.isfile(mden_file):
            with open(mden_file) as f:
                for line in f:
                    if line.startswith("L6"):
                        s = line.split()
                        if s[2] != "E_pot":
                            energies.append(float(s[2]))
        elif mdout_file is not None and os.path.isfile(mdout_file):
            with open(mdout_file) as f:
                for line in f:
                    if "EPtot" in line:
                        s = line.split()
                        energies.append(float(s[-1]))
        else:
            raise RuntimeError("Please provide one of mden_file and mdout_file")

    atom_names, atom_types, atom_numbs = np.unique(amber_types, return_inverse=True, return_counts=True)

    data = {}
    data['atom_names'] = list(atom_names)
    data['atom_numbs'] = list(atom_numbs)
    data['atom_types'] = atom_types
    if labeled:
        data['forces'] = forces * force_convert
        data['energies'] = np.array(energies) * energy_convert
    data['coords'] = coords
    data['cells'] = cells
    data['orig'] = np.array([0, 0, 0])
    return data

