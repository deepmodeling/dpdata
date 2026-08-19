from __future__ import annotations

import re
import warnings

import numpy as np


def atom_name_from_potcar_string(instr: str) -> str:
    """Get atom name from a potcar element name.

    e.g. Sn_d -> Sn

    Parameters
    ----------
    instr : str
        input potcar elemenet name

    Returns
    -------
    name: str
        name of atoms
    """
    if "_" in instr:
        # for case like : TITEL  = PAW_PBE Sn_d 06Sep2000
        return instr.split("_")[0]
    else:
        return instr


# VASP echoes the POSCAR title into a fixed-width field, so a long formula is
# cut off; anything reaching the full width may have lost its last species.
_POSCAR_TITLE_WIDTH = 40
_FORMULA_TOKEN = re.compile(r"^([A-Z][a-z]?)(\d*)$")


def composition_from_poscar_title(title: str) -> dict[str, int] | None:
    """Read a composition from a POSCAR title when every count is explicit.

    The title is free-form text, so this returns ``None`` unless every token
    reads as an element symbol with a count -- ``Li3 F39 K3`` yields
    ``{"Li": 3, "F": 39, "K": 3}``, while both ``H C`` and ``POSCAR file
    written by OVITO`` yield ``None``. Requiring counts avoids interpreting a
    comment that merely lists elements as a declaration of their order. A
    title filling the whole field is treated as truncated and its final token
    is dropped.

    Parameters
    ----------
    title : str
        the text following ``POSCAR =`` in the OUTCAR

    Returns
    -------
    Optional[dict[str, int]]
        element counts, or None if the title is not an explicit composition
    """
    from dpdata.periodic_table import ELEMENTS

    truncated = len(title.rstrip()) >= _POSCAR_TITLE_WIDTH
    tokens = title.split()
    if truncated:
        tokens = tokens[:-1]
    if len(tokens) < 2:
        # A single symbol cannot disagree on ordering, and a one-word title is
        # far more likely to be prose than a formula.
        return None
    composition = {}
    for token in tokens:
        matched = _FORMULA_TOKEN.match(token)
        if matched is None or matched.group(1) not in ELEMENTS or not matched.group(2):
            return None
        name = matched.group(1)
        composition[name] = composition.get(name, 0) + int(matched.group(2))
    return composition


def check_potcar_poscar_order(
    atom_names: list[str], atom_numbs: list[int], poscar_title: str | None
) -> None:
    """Warn when POTCAR-paired counts contradict a POSCAR-title composition.

    VASP pairs the ``ions per type`` counts, which come from the POSCAR, with
    the species order of the POTCAR. When the two files disagree, VASP neither
    reorders nor complains, so counts can silently land on the wrong elements
    and dpdata faithfully reports the mislabeled system. The POSCAR title is
    only a comment, however, so compare compositions rather than token order.
    """
    if poscar_title is None:
        return
    title_composition = composition_from_poscar_title(poscar_title)
    if title_composition is None:
        return

    potcar_composition = {}
    for name, count in zip(atom_names, atom_numbs):
        potcar_composition[name] = potcar_composition.get(name, 0) + count

    truncated = len(poscar_title.rstrip()) >= _POSCAR_TITLE_WIDTH
    if truncated:
        matches = all(
            potcar_composition.get(name) == count
            for name, count in title_composition.items()
        )
    else:
        matches = title_composition == potcar_composition
    if matches:
        return

    def format_composition(composition: dict[str, int]) -> str:
        return " ".join(f"{name}{count}" for name, count in composition.items())

    warnings.warn(
        "the composition produced by pairing the POTCAR species with 'ions per "
        f"type' in this OUTCAR ({format_composition(potcar_composition)}) does "
        "not match the explicit composition in the POSCAR title "
        f"({format_composition(title_composition)}). The atom names reported "
        "here are what VASP actually computed; if the title describes the "
        "intended structure, the POTCAR may have been concatenated in a "
        "different order and the calculation used the wrong potentials."
    )


def system_info(
    lines: list[str],
    type_idx_zero: bool = False,
) -> tuple[list[str], list[int], np.ndarray, int | None, int | None]:
    """Get system information from lines of an OUTCAR file.

    Parameters
    ----------
    lines : list[str]
        the lines of the OUTCAR file
    type_idx_zero : bool
        if true atom types starts from 0 otherwise from 1.

    Returns
    -------
    atom_names: list[str]
        name of atoms
    atom_numbs: list[int]
        number of atoms that have a certain name. same length as atom_names
    atom_types: np.ndarray
        type of each atom, the array has same lenght as number of atoms
    nelm: optional[int]
        the value of NELM parameter
    nwrite: optional[int]
        the value of NWRITE parameter
    """
    atom_names = []
    atom_names_potcar = []
    atom_numbs = None
    nelm = None
    nwrite = None
    poscar_title = None
    for ii in lines:
        if "TITEL" in ii:
            # get atom names from POTCAR info, tested only for PAW_PBE ...
            # for case like : TITEL  = PAW_PBE Sn_d 06Sep2000
            _ii = ii.split()[3]
            atom_names.append(atom_name_from_potcar_string(_ii))
        elif "POTCAR:" in ii:
            # get atom names from POTCAR info, tested only for PAW_PBE ...
            # for case like : POTCAR:  PAW_PBE Ti 08Apr2002
            _ii = ii.split()[2]
            atom_names_potcar.append(atom_name_from_potcar_string(_ii))
        # a stricker check for "NELM"; compatible with distingct formats in different versions(6 and older, newers_expect-to-work) of vasp
        elif nelm is None:
            m = re.search(r"NELM\s*=\s*(\d+)", ii)
            if m:
                nelm = int(m.group(1))
        elif nwrite is None:
            m = re.search(r"NWRITE\s*=\s*(\d+)", ii)
            if m:
                nwrite = int(m.group(1))
        if poscar_title is None:
            # the POSCAR title echoed among the start parameters, e.g.
            # POSCAR =  Li3 F39 K3 Mg3 Ca3 Na3 Al10 O6
            m = re.match(r"\s*POSCAR\s*=\s*(.*)$", ii)
            if m:
                poscar_title = m.group(1)
        if "ions per type" in ii:
            atom_numbs_ = [int(s) for s in ii.split()[4:]]
            if atom_numbs is None:
                atom_numbs = atom_numbs_
            else:
                assert atom_numbs == atom_numbs_, "in consistent numb atoms in OUTCAR"
    if len(atom_names) == 0:
        # try to use atom_names_potcar
        if len(atom_names_potcar) == 0:
            raise ValueError("cannot get atom names from potcar")
        nnames = len(atom_names_potcar)
        # the names are repeated. check if it is the case
        assert atom_names_potcar[: nnames // 2] == atom_names_potcar[nnames // 2 :]
        atom_names = atom_names_potcar[: nnames // 2]
    assert nelm is not None, "cannot find maximum steps for each SC iteration"
    assert atom_numbs is not None, "cannot find ion type info in OUTCAR"
    if len(atom_numbs) != len(atom_names):
        raise RuntimeError(
            f"The number of the atom numbers per each type ({len(atom_numbs)}) "
            f"does not match that of the atom types ({len(atom_names)}) detected "
            f"from the OUTCAR. This issue may be cause by a bug in vasp <= 6.3. "
            f"Please try to convert data from vasprun.xml instead."
        )
    atom_names = atom_names[: len(atom_numbs)]
    check_potcar_poscar_order(atom_names, atom_numbs, poscar_title)
    atom_types = []
    for idx, ii in enumerate(atom_numbs):
        for jj in range(ii):
            if type_idx_zero:
                atom_types.append(idx)
            else:
                atom_types.append(idx + 1)
    return atom_names, atom_numbs, np.array(atom_types, dtype=int), nelm, nwrite


def get_outcar_block(fp, ml=False):
    blk = []
    energy_token = ["free  energy   TOTEN", "free  energy ML TOTEN"]
    ml_index = int(ml)
    for ii in fp:
        if not ii:
            return blk
        blk.append(ii.rstrip("\n"))
        if energy_token[ml_index] in ii:
            return blk
    return blk


class _IncompleteForceTableError(ValueError):
    """Indicate that a ``TOTAL-FORCE`` table ended before all atom rows."""


def check_outputs(coord, cell, force):
    if len(force) == 0:
        raise ValueError("cannot find forces in OUTCAR block")
    if len(coord) == 0:
        raise ValueError("cannot find coordinates in OUTCAR block")
    if len(cell) == 0:
        raise ValueError("cannot find cell in OUTCAR block")
    return True


# we assume that the force is printed ...
def get_frames(fname, begin=0, step=1, ml=False, convergence_check=True):
    with open(fname) as fp:
        return _get_frames_lower(
            fp,
            fname,
            begin=begin,
            step=step,
            ml=ml,
            convergence_check=convergence_check,
        )


def _get_frames_lower(fp, fname, begin=0, step=1, ml=False, convergence_check=True):
    # Split on the same token as every later block. An ML_ISTART=2 run writes
    # only "free  energy ML TOTEN", so asking for the ab initio token here
    # swallowed the whole file into one block and yielded a single frame.
    blk = get_outcar_block(fp, ml)

    atom_names, atom_numbs, atom_types, nelm, nwrite = system_info(
        blk, type_idx_zero=True
    )
    ntot = sum(atom_numbs)

    all_coords = []
    all_cells = []
    all_energies = []
    all_forces = []
    all_virials = []

    cc = 0
    rec_failed = []
    while len(blk) > 0:
        if cc >= begin and (cc - begin) % step == 0:
            try:
                coord, cell, energy, force, virial, is_converge = analyze_block(
                    blk, ntot, nelm, ml
                )
            except _IncompleteForceTableError as exc:
                # An interrupted VASP process can leave one damaged ionic
                # step between otherwise usable frames.  Its atom arrays are
                # incomplete, so skip precisely that step instead of aborting
                # the complete trajectory.
                warnings.warn(
                    f"incomplete labels in frame {cc + 1}; it is ignored: {exc}"
                )
                blk = get_outcar_block(fp, ml)
                cc += 1
                continue
            if energy is None:
                # Reaching here normally just means the trailing timing report
                # after the last ionic step. Only an interrupted run leaves a
                # block that had started emitting step data, and nothing after
                # it can be read either, since the block splitter keys on the
                # energy token.
                if len(coord) > 0 or len(cell) > 0 or len(force) > 0:
                    warnings.warn(
                        f"no energy found in frame {cc + 1}; it and any frame "
                        "after it are not read. This usually means the "
                        "calculation was interrupted before it finished writing."
                    )
                break
            if nwrite == 0:
                has_label = len(force) > 0 and len(coord) > 0 and len(cell) > 0
                if not has_label:
                    warnings.warn("cannot find labels in the frame, ingore")
            else:
                has_label = check_outputs(coord, cell, force)
            if (is_converge or not convergence_check) and has_label:
                all_coords.append(coord)
                all_cells.append(cell)
                all_energies.append(energy)
                all_forces.append(force)
                if virial is not None:
                    all_virials.append(virial)
            if not is_converge:
                rec_failed.append(cc + 1)

        blk = get_outcar_block(fp, ml)
        cc += 1

    if len(rec_failed) > 0:
        prt = (
            "so they are not collected."
            if convergence_check
            else "but they are still collected due to the requirement for ignoring convergence checks."
        )
        warnings.warn(
            f"The following structures were unconverged: {rec_failed}; " + prt
        )

    if len(all_virials) == 0:
        all_virials = None
    else:
        all_virials = np.array(all_virials)
    return (
        atom_names,
        atom_numbs,
        atom_types,
        np.array(all_cells),
        np.array(all_coords),
        np.array(all_energies),
        np.array(all_forces),
        all_virials,
    )


def analyze_block(lines, ntot, nelm, ml=False):
    coord = []
    cell = []
    energy = None
    force = []
    virial = None
    is_converge = True
    sc_index = 0
    # select different searching tokens based on the ml label
    energy_token = ["free  energy   TOTEN", "free  energy ML TOTEN"]
    energy_index = [4, 5]
    virial_token = ["FORCE on cell =-STRESS in cart. coord.  units", "ML FORCE"]
    virial_index = [14, 4]
    cell_token = ["VOLUME and BASIS", "ML FORCE"]
    cell_index = [5, 12]
    ml_index = int(ml)
    for idx, ii in enumerate(lines):
        # if set ml == True, is_converged will always be True
        if ("Iteration" in ii) and (not ml):
            sc_index = int(ii.split()[3][:-1])
            if sc_index >= nelm:
                is_converge = False
        elif energy_token[ml_index] in ii:
            energy = float(ii.split()[energy_index[ml_index]])
            return coord, cell, energy, force, virial, is_converge
        elif cell_token[ml_index] in ii:
            for dd in range(3):
                tmp_l = lines[idx + cell_index[ml_index] + dd]
                cell.append([float(ss) for ss in tmp_l.replace("-", " -").split()[0:3]])
        elif virial_token[ml_index] in ii:
            in_kB_index = virial_index[ml_index]
            while idx + in_kB_index < len(lines) and (
                not lines[idx + in_kB_index].split()[0:2] == ["in", "kB"]
            ):
                in_kB_index += 1
            assert idx + in_kB_index < len(lines), (
                'ERROR: "in kB" is not found in OUTCAR. Unable to extract virial.'
            )
            tmp_v = [float(ss) for ss in lines[idx + in_kB_index].split()[2:8]]
            virial = np.zeros([3, 3])
            virial[0][0] = tmp_v[0]
            virial[1][1] = tmp_v[1]
            virial[2][2] = tmp_v[2]
            virial[0][1] = tmp_v[3]
            virial[1][0] = tmp_v[3]
            virial[1][2] = tmp_v[4]
            virial[2][1] = tmp_v[4]
            virial[0][2] = tmp_v[5]
            virial[2][0] = tmp_v[5]
        elif "TOTAL-FORCE" in ii and (("ML" in ii) == ml):
            block_coord = []
            block_force = []
            for jj in range(idx + 2, idx + 2 + ntot):
                if jj >= len(lines):
                    # VASP output can end while a force table is being
                    # written.  Do not index beyond the block or retain the
                    # preceding subset of atoms as a complete frame.
                    raise _IncompleteForceTableError(
                        f"expected {ntot} atom rows, found {len(block_force)}"
                    )
                fields = lines[jj].split()
                if len(fields) < 6:
                    raise _IncompleteForceTableError(
                        f"expected 6 numeric columns, got {len(fields)} in {lines[jj]!r}"
                    )
                try:
                    info = [float(ss) for ss in fields[:6]]
                except ValueError as exc:
                    # A non-numeric record with at least six fields, such as
                    # the following energy record, means the table ended
                    # before ``ntot`` atom rows.  Parse into temporary lists
                    # so none of the incomplete table leaks into the frame.
                    raise _IncompleteForceTableError(
                        f"non-numeric atom row {lines[jj]!r}"
                    ) from exc
                block_coord.append(info[:3])
                block_force.append(info[3:6])
            coord = block_coord
            force = block_force
    return coord, cell, energy, force, virial, is_converge
