import os
import warnings

import numpy as np

from .scf import (
    bohr2ang,
    get_cell,
    get_coords,
    get_geometry_in,
    kbar2evperang3,
)

# Read in geometries from an ABACUS MD trajectory.
# The atomic coordinates are read in from generated files in OUT.XXXX.
# Energies, forces
# IMPORTANT: the program defaultly takes STRU input file as standard cell information,
# therefore the direct and cartesan coordinates read could be different from the ones in
# the output cif files!!!
# It is highly recommanded to use ORTHOGANAL coordinates in STRU file if you wish to get
# same coordinates in both dpdata and output cif files.


def get_path_out(fname, inlines):
    # This function is different from the same-name function in scf.py.
    # This function returns OUT.XXXX's base directory.
    path_out = os.path.join(fname, "OUT.ABACUS/")
    for line in inlines:
        if len(line) > 0 and "suffix" in line and "suffix" == line.split()[0]:
            suffix = line.split()[1]
            path_out = os.path.join(fname, "OUT.%s/" % suffix)
            break
    return path_out


def get_coord_dump_freq(inlines):
    for line in inlines:
        if len(line) > 0 and "md_dumpfreq" in line and "md_dumpfreq" == line.split()[0]:
            return int(line.split()[1])
    return 1


def get_coords_from_dump(dumplines, natoms):
    nlines = len(dumplines)
    total_natoms = sum(natoms)
    # The output of VIRIAL, FORCE, and VELOCITY are controlled by INPUT parameters dump_virial, dump_force, and dump_vel, respectively.
    # So the search of keywords can determine whether these datas are printed into MD_dump.
    calc_stress = False
    calc_force = False
    check_line = 6
    if "VIRIAL" in dumplines[6]:
        calc_stress = True
        check_line = 10
    assert (
        "POSITION" in dumplines[check_line]
    ), "keywords 'POSITION' cannot be found in the 6th line. Please check."
    if "FORCE" in dumplines[check_line]:
        calc_force = True

    nframes_dump = -1
    if calc_stress:
        nframes_dump = int(nlines / (total_natoms + 13))
    else:
        nframes_dump = int(nlines / (total_natoms + 9))
    assert nframes_dump > 0, (
        "Number of lines in MD_dump file = %d. Number of atoms = %d. The MD_dump file is incomplete."
        % (nlines, total_natoms)
    )
    cells = np.zeros([nframes_dump, 3, 3])
    stresses = np.zeros([nframes_dump, 3, 3])
    forces = np.zeros([nframes_dump, total_natoms, 3])
    coords = np.zeros([nframes_dump, total_natoms, 3])
    iframe = 0
    for iline in range(nlines):
        if "MDSTEP" in dumplines[iline]:
            # read in LATTICE_CONSTANT
            # for abacus version >= v3.1.4, the unit is angstrom, and "ANGSTROM" is added at the end
            # for abacus version <  v3.1.4, the unit is bohr
            celldm = float(dumplines[iline + 1].split()[1])
            newversion = True
            if "Angstrom" not in dumplines[iline + 1]:
                celldm *= bohr2ang  # transfer unit to ANGSTROM
                newversion = False

            # read in LATTICE_VECTORS
            for ix in range(3):
                cells[iframe, ix] = (
                    np.array([float(i) for i in dumplines[iline + 3 + ix].split()[0:3]])
                    * celldm
                )
                if calc_stress:
                    stresses[iframe, ix] = np.array(
                        [float(i) for i in dumplines[iline + 7 + ix].split()[0:3]]
                    )

            if calc_stress:
                skipline = 11
            else:
                skipline = 7

            for iat in range(total_natoms):
                # INDEX    LABEL    POSITION (Angstrom)    FORCE (eV/Angstrom)    VELOCITY (Angstrom/fs)
                # 0  Sn  0.000000000000  0.000000000000  0.000000000000  -0.000000000000  -0.000000000001  -0.000000000001  0.001244557166  -0.000346684288  0.000768457739
                # 1  Sn  0.000000000000  3.102800034079  3.102800034079  -0.000186795145  -0.000453823768  -0.000453823768  0.000550996187  -0.000886442775  0.001579501983
                # for abacus version >= v3.1.4, the value of POSITION is the real cartessian position, and unit is angstrom, and if cal_force the VELOCITY is added at the end.
                # for abacus version < v3.1.4, the real position = POSITION * celldm
                coords[iframe, iat] = np.array(
                    [float(i) for i in dumplines[iline + skipline + iat].split()[2:5]]
                )

                if not newversion:
                    coords[iframe, iat] *= celldm

                if calc_force:
                    forces[iframe, iat] = np.array(
                        [
                            float(i)
                            for i in dumplines[iline + skipline + iat].split()[5:8]
                        ]
                    )
            iframe += 1
    assert iframe == nframes_dump, (
        "iframe=%d, nframe_dump=%d. Number of frames does not match number of lines in MD_dump."
        % (iframe, nframes_dump)
    )
    stresses *= kbar2evperang3
    return coords, cells, forces, stresses


def get_energy(outlines, ndump, dump_freq):
    energy = []
    nenergy = 0
    for line_idx, line in enumerate(outlines):
        if "final etot is" in line:
            if nenergy % dump_freq == 0:
                energy.append(float(line.split()[-2]))
            nenergy += 1
        elif "!! convergence has not been achieved" in line:
            if nenergy % dump_freq == 0:
                energy.append(np.nan)
            nenergy += 1
    assert ndump == len(energy), (
        "Number of total energies in running_md.log = %d. Number of frames in MD_dump = %d. Please check."
        % (len(energy), ndump)
    )
    energy = np.array(energy)
    return energy


def get_frame(fname):
    if isinstance(fname, str):
        # if the input parameter is only one string, it is assumed that it is the
        # base directory containing INPUT file;
        path_in = os.path.join(fname, "INPUT")
    else:
        raise RuntimeError("invalid input")
    with open(path_in) as fp:
        inlines = fp.read().split("\n")
    geometry_path_in = get_geometry_in(fname, inlines)  # base dir of STRU
    path_out = get_path_out(fname, inlines)

    with open(geometry_path_in) as fp:
        geometry_inlines = fp.read().split("\n")
    celldm, cell = get_cell(geometry_inlines)
    atom_names, natoms, types, coords = get_coords(
        celldm, cell, geometry_inlines, inlines
    )
    # This coords is not to be used.
    dump_freq = get_coord_dump_freq(inlines=inlines)
    # ndump = int(os.popen("ls -l %s | grep 'md_pos_' | wc -l" %path_out).readlines()[0])
    # number of dumped geometry files
    # coords = get_coords_from_cif(ndump, dump_freq, atom_names, natoms, types, path_out, cell)
    with open(os.path.join(path_out, "MD_dump")) as fp:
        dumplines = fp.read().split("\n")
    coords, cells, force, stress = get_coords_from_dump(dumplines, natoms)
    ndump = np.shape(coords)[0]
    with open(os.path.join(path_out, "running_md.log")) as fp:
        outlines = fp.read().split("\n")
    energy = get_energy(outlines, ndump, dump_freq)

    unconv_stru = ""
    for i, iene in enumerate(energy):
        if np.isnan(iene):
            coords = np.delete(coords, i - ndump, axis=0)
            cells = np.delete(cells, i - ndump, axis=0)
            force = np.delete(force, i - ndump, axis=0)
            stress = np.delete(stress, i - ndump, axis=0)
            energy = np.delete(energy, i - ndump, axis=0)
            unconv_stru += "%d " % i
    ndump = len(energy)
    if unconv_stru != "":
        warnings.warn("Structure %s are unconverged and not collected!" % unconv_stru)

    for iframe in range(ndump):
        stress[iframe] *= np.linalg.det(cells[iframe, :, :].reshape([3, 3]))
    if np.sum(np.abs(stress[0])) < 1e-10:
        stress = None
    data = {}
    data["atom_names"] = atom_names
    data["atom_numbs"] = natoms
    data["atom_types"] = types
    data["cells"] = cells
    # for idx in range(ndump):
    #    data['cells'][:, :, :] = cell
    data["coords"] = coords
    data["energies"] = energy
    data["forces"] = force
    data["virials"] = stress
    if not isinstance(data["virials"], np.ndarray):
        del data["virials"]
    data["orig"] = np.zeros(3)

    return data
