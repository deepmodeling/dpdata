#!/usr/bin/python3
from __future__ import annotations

import numpy as np


def _safe_float(value: str, context: str = "") -> float:
    """Safely convert string to float with informative error message.
    
    Parameters
    ----------
    value : str
        String value to convert to float
    context : str, optional
        Context information for error message
        
    Returns
    -------
    float
        Converted float value
        
    Raises
    ------
    ValueError
        If conversion fails, with informative error message
    """
    try:
        return float(value)
    except ValueError as e:
        if context:
            raise ValueError(f"Failed to parse {context}: {e}") from e
        else:
            raise ValueError(f"Failed to convert to float: {e}") from e


def _to_system_data_lower(lines, cartesian=True, selective_dynamics=False):
    def move_flag_mapper(flag):
        if flag == "T":
            return True
        elif flag == "F":
            return False
        else:
            raise RuntimeError(f"Invalid move flag: {flag}")

    """Treat as cartesian poscar."""
    system = {}
    system["atom_names"] = [str(ii) for ii in lines[5].split()]
    system["atom_numbs"] = [int(ii) for ii in lines[6].split()]
    scale = _safe_float(lines[1].strip(), "scale factor from POSCAR line 2")
    cell = []
    move_flags = []
    for ii in range(2, 5):
        line_content = lines[ii].split()
        try:
            boxv = [_safe_float(jj, f"cell vector component on line {ii+1}") for jj in line_content]
        except ValueError as e:
            raise ValueError(f"Failed to parse cell vectors in POSCAR: {e}") from e
        boxv = np.array(boxv) * scale
        cell.append(boxv)
    system["cells"] = [np.array(cell)]
    natoms = sum(system["atom_numbs"])
    coord = []
    for ii in range(8, 8 + natoms):
        tmp = lines[ii].split()
        try:
            tmpv = [_safe_float(jj, f"coordinate component on line {ii+1}") for jj in tmp[:3]]
        except ValueError as e:
            raise ValueError(f"Failed to parse coordinates in POSCAR: {e}") from e
        if cartesian:
            tmpv = np.array(tmpv) * scale
        else:
            tmpv = np.matmul(np.array(tmpv), system["cells"][0])
        coord.append(tmpv)
        if selective_dynamics:
            if len(tmp) == 6:
                move_flags.append(list(map(move_flag_mapper, tmp[3:])))
            else:
                raise RuntimeError(
                    f"Invalid move flags, should be 6 columns, got {tmp}"
                )

    system["coords"] = [np.array(coord)]
    system["orig"] = np.zeros(3)
    atom_types = []
    for idx, ii in enumerate(system["atom_numbs"]):
        for jj in range(ii):
            atom_types.append(idx)
    system["atom_types"] = np.array(atom_types, dtype=int)
    system["cells"] = np.array(system["cells"])
    system["coords"] = np.array(system["coords"])
    if move_flags:
        move_flags = np.array(move_flags, dtype=bool)
        move_flags = move_flags.reshape((1, natoms, 3))
        system["move"] = np.array(move_flags, dtype=bool)
    return system


def to_system_data(lines):
    # remove the line that has 'selective dynamics'
    selective_dynamics = False
    if lines[7][0] == "S" or lines[7][0] == "s":
        selective_dynamics = True
        lines.pop(7)
    is_cartesian = lines[7][0] in ["C", "c", "K", "k"]
    if not is_cartesian:
        if lines[7][0] not in ["d", "D"]:
            raise RuntimeError(
                "seem not to be a valid POSCAR of vasp 5.x, may be a POSCAR of vasp 4.x?"
            )
    return _to_system_data_lower(lines, is_cartesian, selective_dynamics)


def from_system_data(system, f_idx=0, skip_zeros=True):
    ret = ""
    for ii, name in zip(system["atom_numbs"], system["atom_names"]):
        if ii == 0:
            continue
        ret += "%s%d " % (name, ii)  # noqa: UP031
    ret += "\n"
    ret += "1.0\n"
    for ii in system["cells"][f_idx]:
        for jj in ii:
            ret += f"{jj:.16e} "
        ret += "\n"
    for idx, ii in enumerate(system["atom_names"]):
        if system["atom_numbs"][idx] == 0:
            continue
        ret += f"{ii} "
    ret += "\n"
    for ii in system["atom_numbs"]:
        if ii == 0:
            continue
        ret += "%d " % ii  # noqa: UP031
    ret += "\n"
    move = system.get("move", None)
    if move is not None and len(move) > 0:
        ret += "Selective Dynamics\n"

    # should use Cartesian for VESTA software
    ret += "Cartesian\n"
    atype = system["atom_types"]
    posis = system["coords"][f_idx]
    # atype_idx = [[idx,tt] for idx,tt in enumerate(atype)]
    # sort_idx = np.argsort(atype, kind = 'mergesort')
    sort_idx = np.lexsort((np.arange(len(atype)), atype))
    atype = atype[sort_idx]
    posis = posis[sort_idx]
    if move is not None and len(move) > 0:
        move = move[f_idx][sort_idx]

    if isinstance(move, np.ndarray):
        move = move.tolist()

    posi_list = []
    for idx in range(len(posis)):
        ii_posi = posis[idx]
        line = f"{ii_posi[0]:15.10f} {ii_posi[1]:15.10f} {ii_posi[2]:15.10f}"
        if move is not None and len(move) > 0:
            move_flags = move[idx]
            if not isinstance(move_flags, list) or len(move_flags) != 3:
                raise RuntimeError(
                    f"Invalid move flags: {move_flags}, should be a list of 3 bools"
                )
            line += " " + " ".join("T" if flag else "F" for flag in move_flags)

        posi_list.append(line)

    posi_list.append("")
    ret += "\n".join(posi_list)
    return ret
