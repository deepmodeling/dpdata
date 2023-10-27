import os

import numpy as np

from .scf import (
    bohr2ang,
    collect_force,
    collect_stress,
    get_cell,
    get_coords,
    get_geometry_in,
    kbar2evperang3,
)

# Read in geometries from an ABACUS RELAX(CELL-RELAX) trajectory in OUT.XXXX/runnning_relax/cell-relax.log.


def get_log_file(fname, inlines):
    suffix = "ABACUS"
    calculation = "scf"
    for line in inlines:
        if "suffix" in line and "suffix" == line.split()[0]:
            suffix = line.split()[1]
        elif "calculation" in line and "calculation" == line.split()[0]:
            calculation = line.split()[1]
    logf = os.path.join(fname, f"OUT.{suffix}/running_{calculation}.log")
    return logf


def get_coords_from_log(loglines, natoms):
    """NOTICE: unit of coords and cells is Angstrom
    order:
        coordinate
        cell (no output if cell is not changed)
        energy (no output, if SCF is not converged)
        force (no output, if cal_force is not setted or abnormal ending)
        stress (no output, if set cal_stress is not setted or abnormal ending).
    """
    natoms_log = 0
    for line in loglines:
        if line[13:41] == "number of atom for this type":
            natoms_log += int(line.split()[-1])

    assert natoms_log > 0 and natoms_log == natoms, (
        "ERROR: detected atom number in log file is %d" % natoms
    )

    energy = []
    cells = []
    coords = []
    coord_direct = []  # if the coordinate is direct type or not

    for i in range(len(loglines)):
        line = loglines[i]
        if line[18:41] == "lattice constant (Bohr)":
            a0 = float(line.split()[-1])
        elif len(loglines[i].split()) >= 2 and loglines[i].split()[1] == "COORDINATES":
            # read coordinate information
            coords.append([])
            direct_coord = False
            if loglines[i].split()[0] == "DIRECT":
                coord_direct.append(True)
                for k in range(2, 2 + natoms):
                    coords[-1].append(
                        list(map(lambda x: float(x), loglines[i + k].split()[1:4]))
                    )
            elif loglines[i].split()[0] == "CARTESIAN":
                coord_direct.append(False)
                for k in range(2, 2 + natoms):
                    coords[-1].append(
                        list(map(lambda x: float(x) * a0, loglines[i + k].split()[1:4]))
                    )
            else:
                assert False, "Unrecongnized coordinate type, %s, line:%d" % (
                    loglines[i].split()[0],
                    i,
                )

        elif (
            loglines[i][1:56]
            == "Lattice vectors: (Cartesian coordinate: in unit of a_0)"
        ):
            # add the cell information for previous structures
            while len(cells) < len(coords) - 1:
                cells.append(cells[-1])
            # get current cell information
            cells.append([])
            for k in range(1, 4):
                cells[-1].append(
                    list(map(lambda x: float(x) * a0, loglines[i + k].split()[0:3]))
                )

        elif line[1:14] == "final etot is":
            # add the energy for previous structures whose SCF is not converged
            while len(energy) < len(coords) - 1:
                energy.append(np.nan)
            # get the energy of current structure
            energy.append(float(line.split()[-2]))

    force = collect_force(loglines)
    stress = collect_stress(loglines)

    # delete last structures which has no energy
    while len(energy) < len(coords):
        del coords[-1]
        del coord_direct[-1]

    # add cells for last structures whose cell is not changed
    while len(cells) < len(coords):
        cells.append(cells[-1])

    # only keep structures that have all of coord, force and stress
    if len(stress) == 0 and len(force) == 0:
        minl = len(coords)
    elif len(stress) == 0:
        minl = min(len(coords), len(force))
        force = force[:minl]
    elif len(force) == 0:
        minl = min(len(coords), len(stress))
        stress = stress[:minl]
    else:
        minl = min(len(coords), len(force), len(stress))
        force = force[:minl]
        stress = stress[:minl]

    coords = coords[:minl]
    energy = energy[:minl]
    cells = cells[:minl]

    # delete structures whose energy is np.nan
    for i in range(minl):
        if np.isnan(energy[i - minl]):
            del energy[i - minl]
            del coords[i - minl]
            del cells[i - minl]
            del coord_direct[i - minl]
            if len(force) > 0:
                del force[i - minl]
            if len(stress) > 0:
                del stress[i - minl]

    energy = np.array(energy)
    cells = np.array(cells)
    coords = np.array(coords)
    stress = np.array(stress)
    force = np.array(force)

    # transfer direct coordinate to cartessian type
    for i in range(len(coords)):
        if coord_direct[i]:
            coords[i] = coords[i].dot(cells[i])

    # transfer bohrium to angstrom
    cells *= bohr2ang
    coords *= bohr2ang

    if len(stress) > 0:
        virial = np.zeros([len(cells), 3, 3])
        for i in range(len(cells)):
            volume = np.linalg.det(cells[i, :, :].reshape([3, 3]))
            virial[i] = stress[i] * kbar2evperang3 * volume
    else:
        virial = None

    return energy, cells, coords, force, stress, virial


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
    with open(geometry_path_in) as fp:
        geometry_inlines = fp.read().split("\n")
    celldm, cell = get_cell(geometry_inlines)
    atom_names, natoms, types, coord_tmp = get_coords(
        celldm, cell, geometry_inlines, inlines
    )

    logf = get_log_file(fname, inlines)
    assert os.path.isfile(logf), "Error: can not find %s" % logf
    with open(logf) as f1:
        lines = f1.readlines()

    atomnumber = 0
    for i in natoms:
        atomnumber += i
    energy, cells, coords, force, stress, virial = get_coords_from_log(
        lines, atomnumber
    )

    data = {}
    data["atom_names"] = atom_names
    data["atom_numbs"] = natoms
    data["atom_types"] = types
    data["cells"] = cells
    data["coords"] = coords
    data["energies"] = energy
    data["forces"] = force
    if isinstance(virial, np.ndarray):
        data["virials"] = virial
    data["stress"] = stress
    data["orig"] = np.zeros(3)

    return data
