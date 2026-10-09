#!/usr/bin/env python3
from __future__ import annotations

import os
import re

import numpy as np

from dpdata.utils import open_file

from .traj import (
    convert_celldm,
    kbar2evperang3,
    ry2ev,
)
from .traj import (
    length_convert as bohr2ang,
)

_QE_BLOCK_KEYWORDS = [
    "ATOMIC_SPECIES",
    "ATOMIC_POSITIONS",
    "K_POINTS",
    "ADDITIONAL_K_POINTS",
    "CELL_PARAMETERS",
    "CONSTRAINTS",
    "OCCUPATIONS",
    "ATOMIC_VELOCITIES",
    "ATOMIC_FORCES",
    "SOLVENTS",
    "HUBBARD",
]


def get_block(lines, keyword, skip=0):
    ret = []
    for idx, ii in enumerate(lines):
        if keyword in ii:
            blk_idx = idx + 1 + skip
            while len(lines[blk_idx].split()) == 0:
                blk_idx += 1
            while (
                len(lines[blk_idx].split()) != 0
                and (lines[blk_idx].split()[0] not in _QE_BLOCK_KEYWORDS)
            ) and blk_idx != len(lines):
                ret.append(lines[blk_idx])
                blk_idx += 1
            break
    return ret


def get_cell(lines):
    ret = []
    idx = None
    match = None
    for line_idx, ii in enumerate(lines):
        match = re.search(r"ibrav\s*=\s*([+-]?\d+)", ii, re.IGNORECASE)
        if match is not None:
            idx = line_idx
            break
    if idx is None or match is None:
        raise RuntimeError("ibrav cannot be found in QE input")
    ibrav = int(match.group(1))
    if ibrav == 0:
        for line_idx, line in enumerate(lines):
            if "CELL_PARAMETERS" not in line.upper():
                continue
            block = lines[line_idx + 1 : line_idx + 4]
            if len(block) != 3:
                raise RuntimeError("CELL_PARAMETERS block must contain three rows")
            ret = np.array(
                [
                    [
                        float(value.replace("D", "E").replace("d", "e"))
                        for value in row.split()[:3]
                    ]
                    for row in block
                ]
            )
            unit = _cell_unit(line)
            if unit == "bohr":
                ret *= bohr2ang
            elif unit == "alat":
                ret *= (_card_alat(line) or _input_alat(lines)) * bohr2ang
            elif unit not in (None, "angstrom"):
                raise RuntimeError(
                    "CELL_PARAMETERS must be written in Angstrom, bohr, or alat."
                )
            break
        else:
            raise RuntimeError("CELL_PARAMETERS block cannot be found")
    elif ibrav == 1:
        a = None
        for iline in lines:
            line = iline.replace("=", " ").replace(",", "").split()
            if len(line) >= 2 and "a" == line[0]:
                # print("line = ", line)
                a = float(line[1].replace("D", "E").replace("d", "e"))
            if len(line) >= 2 and "celldm(1)" == line[0]:
                a = float(line[1].replace("D", "E").replace("d", "e")) * bohr2ang
        # print("a = ", a)
        if not a:
            raise RuntimeError("parameter 'a' or 'celldm(1)' cannot be found.")
        ret = np.array([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]])
    elif ibrav in (2, 3, -3):
        ret = convert_celldm(ibrav, [_input_alat(lines)]) * bohr2ang
    else:
        raise RuntimeError("ibrav > 1 not supported yet.")
    return ret


def _cell_unit(line):
    match = re.search(r"[({]\s*([^})]+?)\s*[})]", line, re.IGNORECASE)
    if match is None:
        return None
    unit = match.group(1).strip().lower()
    if unit.startswith("alat"):
        return "alat"
    if unit.startswith("angstrom"):
        return "angstrom"
    if unit.startswith("bohr"):
        return "bohr"
    return unit


def _input_alat(lines):
    for line in lines:
        words = line.replace("=", " ").replace(",", " ").split()
        if len(words) >= 2 and words[0].lower() == "celldm(1)":
            return float(words[1].replace("D", "E").replace("d", "e"))
        if len(words) >= 2 and words[0].lower() == "a":
            return float(words[1].replace("D", "E").replace("d", "e"))
    raise RuntimeError("parameter 'a' or 'celldm(1)' cannot be found")


def get_coords(lines, cell):
    """Read the first ``ATOMIC_POSITIONS`` block from a QE input file.

    PWscf prints an ``ATOMIC_POSITIONS`` block for every MD step.  The old
    implementation walked over every matching line and converted the list to
    an array inside the loop, so a second block attempted ``append`` on an
    ndarray.  Reading one block here also makes this helper useful for the
    static input metadata needed by the AIMD reader.
    """
    for line_idx, line in enumerate(lines):
        if "ATOMIC_POSITIONS" not in line.upper():
            continue
        unit = _position_unit(line)
        block = _read_atom_block(lines, line_idx + 1)
        if not block:
            continue
        atom_symbol_list = np.array([parts[0] for parts in block])
        coord = np.array(
            [
                [
                    float(value.replace("D", "E").replace("d", "e"))
                    for value in parts[1:4]
                ]
                for parts in block
            ]
        )
        coord = _convert_positions(coord, unit, cell, alat=_card_alat(line))
        break
    else:
        raise RuntimeError("ATOMIC_POSITIONS block cannot be found")

    tmp_names, symbol_idx = np.unique(atom_symbol_list, return_index=True)
    atom_types = []
    atom_numbs = []
    # preserve the atom_name order
    atom_names = atom_symbol_list[np.sort(symbol_idx, kind="stable")]
    for jj in atom_symbol_list:
        for idx, ii in enumerate(atom_names):
            if jj == ii:
                atom_types.append(idx)
    for idx in range(len(atom_names)):
        atom_numbs.append(atom_types.count(idx))
    atom_types = np.array(atom_types)

    return list(atom_names), atom_numbs, atom_types, coord


def _position_unit(line):
    """Return the coordinate unit in an ``ATOMIC_POSITIONS`` header."""
    match = re.search(r"[({]\s*([^})]+?)\s*[})]", line, re.IGNORECASE)
    if match is None:
        return None
    unit = match.group(1).strip().lower()
    if unit.startswith("alat"):
        return "alat"
    if unit.startswith("crystal"):
        return "crystal"
    if unit.startswith("angstrom"):
        return "angstrom"
    if unit.startswith("bohr"):
        return "bohr"
    return unit


def _card_alat(line):
    """Return an explicit ``alat`` value (in bohr) from a QE card header."""
    match = re.search(r"alat\s*=\s*([-+0-9.eEdD]+)", line, re.IGNORECASE)
    if match is None:
        return None
    return float(match.group(1).replace("D", "E").replace("d", "e"))


def _read_atom_block(lines, start, natoms=None):
    """Read atom rows following an ``ATOMIC_POSITIONS`` header."""
    block = []
    idx = start
    while idx < len(lines) and (natoms is None or len(block) < natoms):
        words = lines[idx].split()
        if not words:
            if block:
                break
            idx += 1
            continue
        if words[0].upper() in _QE_BLOCK_KEYWORDS:
            break
        if len(words) < 4:
            break
        try:
            [float(value.replace("D", "E").replace("d", "e")) for value in words[1:4]]
        except ValueError:
            break
        block.append(words)
        idx += 1
    return block


def _convert_positions(coord, unit, cell, alat=None):
    """Convert QE positions to Cartesian Angstrom coordinates."""
    if unit in (None, "angstrom"):
        return coord
    if unit == "bohr":
        return coord * bohr2ang
    if unit == "crystal":
        return np.matmul(coord, cell)
    if unit == "alat":
        if alat is None:
            alat = np.linalg.norm(cell[0]) / bohr2ang
        return coord * alat * bohr2ang
    raise RuntimeError(
        "ATOMIC_POSITIONS must be written in Angstrom, bohr, alat, or crystal. "
        f"Unit {unit!r} is not supported."
    )


def get_energy(lines):
    energy = None
    for ii in lines:
        if "!    total energy" in ii:
            energy = ry2ev * float(ii.split("=")[1].split()[0])
    return energy


def get_force(lines, natoms):
    blk = get_block(lines, "Forces acting on atoms", skip=1)
    ret = []
    blk = blk[0 : sum(natoms)]
    for ii in blk:
        ret.append([float(jj) for jj in ii.split("=")[1].split()])
    ret = np.array(ret)
    ret *= ry2ev / bohr2ang
    return ret


def get_stress(lines):
    blk = get_block(lines, "total   stress")
    if len(blk) == 0:
        return None
    ret = []
    for ii in blk:
        ret.append([float(jj) for jj in ii.split()[3:6]])
    ret = np.array(ret)
    ret *= kbar2evperang3
    return ret


def _output_alat(lines, fallback_cell):
    """Get the QE ``alat`` value in bohr from an output header."""
    for line in lines:
        if "lattice parameter (alat)" not in line.lower():
            continue
        match = re.search(r"=\s*([-+0-9.eEdD]+)", line)
        if match:
            return float(match.group(1).replace("D", "E").replace("d", "e"))
    return np.linalg.norm(fallback_cell[0]) / bohr2ang


def _output_position_frames(lines, natoms, cell):
    """Read all Cartesian position frames printed by ``pw.x``."""
    frames = []
    units = []
    symbols = []
    for line_idx, line in enumerate(lines):
        if "ATOMIC_POSITIONS" not in line.upper():
            continue
        block = _read_atom_block(lines, line_idx + 1, natoms)
        if len(block) != natoms:
            continue
        unit = _position_unit(line)
        symbols_in_block = [parts[0] for parts in block]
        coords = np.array(
            [
                [
                    float(value.replace("D", "E").replace("d", "e"))
                    for value in parts[1:4]
                ]
                for parts in block
            ]
        )
        frames.append(coords)
        units.append((unit, _card_alat(line)))
        if not symbols:
            symbols = symbols_in_block
    return np.array(frames), symbols, units


def _output_cell_frames(lines, nframes, fallback_cell):
    """Read all ``CELL_PARAMETERS`` blocks from a ``pw.x`` output."""
    cells = []
    alat = _output_alat(lines, fallback_cell)
    for line_idx, line in enumerate(lines):
        if "CELL_PARAMETERS" not in line.upper():
            continue
        block = lines[line_idx + 1 : line_idx + 4]
        if len(block) != 3:
            continue
        try:
            current = np.array(
                [
                    [
                        float(value.replace("D", "E").replace("d", "e"))
                        for value in row.split()[:3]
                    ]
                    for row in block
                ]
            )
        except (ValueError, IndexError):
            continue
        unit = _cell_unit(line)
        card_alat = _card_alat(line)
        if unit == "bohr":
            current *= bohr2ang
        elif unit == "alat":
            current *= (card_alat or alat) * bohr2ang
        elif unit not in (None, "angstrom"):
            continue
        cells.append(current)
    if not cells:
        return np.tile(fallback_cell, (nframes, 1, 1))
    if len(cells) == 1:
        return np.tile(cells[0], (nframes, 1, 1))
    if len(cells) != nframes:
        raise RuntimeError(
            "the number of CELL_PARAMETERS blocks does not match "
            f"the number of coordinate frames ({len(cells)} != {nframes})"
        )
    return np.array(cells)


def _output_energies(lines):
    energies = []
    for line in lines:
        if re.search(r"!\s*total energy", line, re.IGNORECASE) is None:
            continue
        try:
            value = line.split("=", 1)[1].split()[0]
            energies.append(ry2ev * float(value.replace("D", "E").replace("d", "e")))
        except (IndexError, ValueError):
            continue
    return np.array(energies)


def _output_forces(lines, natoms):
    forces = []
    for line_idx, line in enumerate(lines):
        if "forces acting on atoms" not in line.lower():
            continue
        frame = []
        for candidate in lines[line_idx + 1 :]:
            if "force =" not in candidate.lower():
                if frame and (
                    "total force" in candidate.lower()
                    or "atomic positions" in candidate.lower()
                ):
                    break
                continue
            try:
                values = candidate.lower().split("force =", 1)[1].split()[:3]
                frame.append([float(value.replace("d", "e")) for value in values])
            except (IndexError, ValueError):
                continue
            if len(frame) == natoms:
                break
        if len(frame) == natoms:
            forces.append(frame)
    if not forces:
        return np.empty((0, natoms, 3))
    return np.asarray(forces) * (ry2ev / bohr2ang)


def _output_stresses(lines):
    stresses = []
    for line_idx, line in enumerate(lines):
        if re.search(r"total\s+stress", line, re.IGNORECASE) is None:
            continue
        frame = []
        for candidate in lines[line_idx + 1 :]:
            words = candidate.split()
            if len(words) < 6:
                if frame:
                    break
                continue
            try:
                # QE prints the kbar values in columns 4--6.
                frame.append(
                    [
                        float(value.replace("D", "E").replace("d", "e"))
                        for value in words[3:6]
                    ]
                )
            except ValueError:
                continue
            if len(frame) == 3:
                break
        if len(frame) == 3:
            stresses.append(frame)
    if not stresses:
        return None
    return np.asarray(stresses) * kbar2evperang3


def _frame_count(*arrays):
    counts = [len(array) for array in arrays if len(array)]
    return min(counts) if counts else 0


def get_frames(fname, begin=0, step=1):
    """Read one or more PWscf single-point or AIMD frames.

    ``pw.x`` writes an ``ATOMIC_POSITIONS`` block for each MD step while the
    historical ``qe/pw/scf`` reader assumed a single block.  This reader keeps
    the single-point behavior and additionally aligns repeated positions,
    energies, forces, cells, and stresses into frame-major arrays.
    """
    if begin < 0 or step <= 0:
        raise ValueError("begin must be non-negative and step must be positive")
    if isinstance(fname, str):
        path_out = fname
        outname = os.path.basename(path_out)
        # The input file is conventionally obtained by replacing ``out`` by ``in``.
        inname = outname.replace("out", "in")
        path_in = os.path.join(os.path.dirname(path_out), inname)
    elif isinstance(fname, list) and len(fname) == 2:
        path_in, path_out = fname
    else:
        raise RuntimeError("invalid input")
    with open_file(path_out) as fp:
        outlines = fp.read().splitlines()
    # Restarted pw.x jobs may append a second run to the same output file.
    # Parse only the latest run, matching the trajectory files written by QE.
    program_headers = [
        idx for idx, line in enumerate(outlines) if "program pwscf" in line.lower()
    ]
    if program_headers:
        outlines = outlines[program_headers[-1] :]
    with open_file(path_in) as fp:
        inlines = fp.read().splitlines()

    cell = get_cell(inlines)
    atom_names, natoms, types, input_coords = get_coords(inlines, cell)
    natom = sum(natoms)
    output_coords, output_symbols, output_units = _output_position_frames(
        outlines, natom, cell
    )
    energies = _output_energies(outlines)
    forces = _output_forces(outlines, natom)
    stresses = _output_stresses(outlines)

    # A regular SCF output has no positions block.  Use the input structure in
    # that case; this preserves the existing one-frame ``qe/pw/scf`` behavior.
    if len(output_coords) == 0:
        output_coords = input_coords[np.newaxis, :, :]
    # A static pw.x run prints one ``total energy`` line for every SCF
    # iteration but only one force frame.  Keep the converged (last) energy in
    # that case; AIMD output normally has one final energy per position frame.
    if len(energies) > len(output_coords):
        energies = energies[-len(output_coords) :]
    nframes = _frame_count(output_coords, energies, forces)
    if nframes == 0:
        raise RuntimeError("QE output does not contain a complete labeled frame")
    if output_symbols:
        output_names = list(dict.fromkeys(output_symbols))
        if output_names != atom_names:
            raise RuntimeError("atom order in QE output does not match the input")
    output_coords = output_coords[:nframes]
    energies = energies[:nframes]
    forces = forces[:nframes]
    cell_frame_count = len(output_coords) if output_units else nframes
    cells = _output_cell_frames(outlines, cell_frame_count, cell)[:nframes]
    if output_units:
        alat = _output_alat(outlines, cell)
        output_coords = np.array(
            [
                _convert_positions(coords, unit, frame_cell, alat=card_alat or alat)
                for coords, (unit, card_alat), frame_cell in zip(
                    output_coords, output_units[:nframes], cells
                )
            ]
        )
        # In a pw.x MD output, the energy and forces are printed for the
        # current geometry before the updated positions for the next step.
        # Associate the first labels with the input geometry and discard the
        # final unlabelled update.
        output_coords = np.concatenate(
            (input_coords[np.newaxis, :, :], output_coords[:-1]), axis=0
        )
        cells = np.concatenate((cell[np.newaxis, :, :], cells[:-1]), axis=0)
    virials = None
    if stresses is not None:
        if len(stresses) < nframes:
            raise RuntimeError(
                "the number of stress blocks does not match the number of frames"
            )
        virials = stresses[:nframes] * np.linalg.det(cells)[:, None, None]
    selection = slice(begin, None, step)
    output_coords = output_coords[selection]
    energies = energies[selection]
    forces = forces[selection]
    cells = cells[selection]
    if virials is not None:
        virials = virials[selection]
    return atom_names, natoms, types, cells, output_coords, energies, forces, virials


def get_frame(fname):
    """Read the first frame from a PWscf output for backward compatibility."""
    data = get_frames(fname)
    atom_names, atom_numbs, atom_types, cells, coords, energies, forces, virials = data
    if virials is not None:
        virials = virials[:1]
    return (
        atom_names,
        atom_numbs,
        atom_types,
        cells[:1],
        coords[:1],
        energies[:1],
        forces[:1],
        virials,
    )
