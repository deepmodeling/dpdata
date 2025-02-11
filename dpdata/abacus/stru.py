import os, re
import warnings
import numpy as np
from ..unit import LengthConversion

bohr2ang = LengthConversion("bohr", "angstrom").value()


def split_stru_block(lines):
    """Split the ABACUS STRU file into blocks by keyword.
    
    Args:
        lines (list): list of lines in the ABACUS STRU file.
        
    Returns:
        dict: dictionary of blocks.
    """
    def clean_comment(line):
        return re.split("[#]", line)[0]

    ABACUS_STRU_KEYS = [
        "ATOMIC_SPECIES",
        "NUMERICAL_ORBITAL",
        "LATTICE_CONSTANT",
        "LATTICE_VECTORS",
        "ATOMIC_POSITIONS",
        "NUMERICAL_DESCRIPTOR",
        "PAW_FILES",
    ]
    blocks = {i:[] for i in ABACUS_STRU_KEYS}
    i = 0
    while i < len(lines):
        line = clean_comment(lines[i]).strip()
        if line in ABACUS_STRU_KEYS:
            key = line
            for j in range(i + 1, len(lines)):
                if clean_comment(lines[j]).strip() == "":
                    continue
                elif clean_comment(lines[j]).strip() in ABACUS_STRU_KEYS:
                    break
                else:
                    blocks[key].append(clean_comment(lines[j]))
            i = j
        else:
            i += 1
    
    return blocks

def parse_atomic_species_block(lines):
    """Parse the ATOMIC_SPECIES block.
    
    Args:
        lines (list): list of lines in the ATOMIC_SPECIES block.
        
    Returns:
        tuple: tuple of atom_names, masses, and pp_files.
        
    """
    atom_names, masses, pp_files = [], [], []
    for line in lines:
        line = line.split()
        atom_names.append(line[0])
        masses.append(float(line[1]))
        
        # for standard STRU, the pseudo potential file is required, 
        # but it is not required for dpdata.
        if len(line) > 2:
            pp_files.append(line[2])
        else:
            pp_files.append(None) 
    
    return atom_names, masses, pp_files

def parse_numerical_orbital_block(lines):
    """Parse the NUMERICAL_ORBITAL block.
    
    Args:
        lines (list): list of lines in the NUMERICAL_ORBITAL block.
        
    Returns:
        list: list of orbital files.
    """
    return [line.strip() for line in lines] 

def parse_lattice_constant_block(lines):
    """Parse the LATTICE_CONSTANT block.
    
    Args:
        lines (list): list of lines in the LATTICE_CONSTANT block.
        
    Returns:
        float: the lattice constant.
    """
    return float(lines[0])

def parse_lattice_vectors_block(lines):
    """Parse the LATTICE_VECTORS block.
    
    Args:
        lines (list): list of lines in the LATTICE_VECTORS block.
        
    Returns:
        np.ndarray: the cell vectors.
    """
    cell = np.zeros((3, 3))
    for i, line in enumerate(lines):
        cell[i] = [float(x) for x in line.split()]
    return cell

def parse_pos_oneline(pos_line):
    """Parses a line from the atom position block in a structure file.

    The content in atom position block can include:
    - `m` or NO key word: Three numbers (0 or 1) controlling atom movement in geometry relaxation calculations.
    - `v`, `vel`, or `velocity`: Three components of initial velocity of atoms in geometry relaxation calculations.
    - `mag` or `magmom`: Start magnetization for each atom. Can be one number (colinear) or three numbers (non-colinear).
    - `angle1`: In non-colinear case, angle between c-axis and real spin (in degrees).
    - `angle2`: In non-colinear case, angle between a-axis and real spin projection in ab-plane (in degrees).
    - `cs` or `constrain`: Three numbers (0 or 1) controlling the spin constraint of the atom.
    - `lambda`: Three numbers controlling the lambda of the atom.

    Parameters
    ----------
    pos_line : A line from the atom position block.

    Returns
    -------
    tuple: A tuple containing:
          - pos (list of float): The position coordinates.
          - move (list of int or None): Movement control values.
          - velocity (list of float or None): Initial velocity components.
          - magmom (float, list of float, or None): Magnetization values.
          - angle1 (float or None): Angle1 value.
          - angle2 (float or None): Angle2 value.
          - constrain (list of bool or None): Spin constraint values.
          - lambda1 (float, list of float, or None): Lambda values.

        e.g.:
        ```
        Fe
        1.0
        2
        0.0 0.0 0.0 m 0 0 0 mag 1.0 angle1 90 angle2 0 cs 0 0 0
        0.5 0.5 0.5 m 1 1 1 mag 1.0 angle1 90 angle2 180
        ```
    """
    pos_line = pos_line.split("#")[0]  # remove comments
    sline = pos_line.split()
    pos = [float(i) for i in sline[:3]]
    move = None
    velocity = None
    magmom = None
    angle1 = None
    angle2 = None
    constrain = None
    lambda1 = None
    if len(sline) > 3:
        mag_list = None
        velocity_list = None
        move_list = []
        angle1_list = None
        angle2_list = None
        constrain_list = None
        lambda_list = None
        label = "move"
        for i in range(3, len(sline)):
            # firstly read the label
            if sline[i] == "m":
                label = "move"
            elif sline[i] in ["v", "vel", "velocity"]:
                label = "velocity"
                velocity_list = []
            elif sline[i] in ["mag", "magmom"]:
                label = "magmom"
                mag_list = []
            elif sline[i] == "angle1":
                label = "angle1"
                angle1_list = []
            elif sline[i] == "angle2":
                label = "angle2"
                angle2_list = []
            elif sline[i] in ["constrain", "sc"]:
                label = "constrain"
                constrain_list = []
            elif sline[i] in ["lambda"]:
                label = "lambda"
                lambda_list = []

            # the read the value to the list
            elif label == "move":
                move_list.append(int(sline[i]))
            elif label == "velocity":
                velocity_list.append(float(sline[i]))
            elif label == "magmom":
                mag_list.append(float(sline[i]))
            elif label == "angle1":
                angle1_list.append(float(sline[i]))
            elif label == "angle2":
                angle2_list.append(float(sline[i]))
            elif label == "constrain":
                constrain_list.append(bool(int(sline[i])))
            elif label == "lambda":
                lambda_list.append(float(sline[i]))

        if move_list is not None and len(move_list) > 0:
            if len(move_list) == 3:
                move = move_list
            else:
                raise RuntimeError(f"Invalid setting of move: {pos_line}")

        if velocity_list is not None:
            if len(velocity_list) == 3:
                velocity = velocity_list
            else:
                raise RuntimeError(f"Invalid setting of velocity: {pos_line}")

        if mag_list is not None:
            if len(mag_list) == 3:
                magmom = mag_list
            elif len(mag_list) == 1:
                magmom = mag_list[0]
            else:
                raise RuntimeError(f"Invalid magnetic moment {pos_line}")

        if angle1_list is not None:
            if len(angle1_list) == 1:
                angle1 = angle1_list[0]
            else:
                raise RuntimeError(f"Invalid angle1 {pos_line}")

        if angle2_list is not None:
            if len(angle2_list) == 1:
                angle2 = angle2_list[0]
            else:
                raise RuntimeError(f"Invalid angle2 {pos_line}")

        if constrain_list is not None:
            if len(constrain_list) == 3:
                constrain = constrain_list
            elif len(constrain_list) == 1:
                constrain = constrain_list[0]
            else:
                raise RuntimeError(f"Invalid constrain {pos_line}")

        if lambda_list is not None:
            if len(lambda_list) == 3:
                lambda1 = lambda_list
            elif len(lambda_list) == 1:
                lambda1 = lambda_list[0]
            else:
                raise RuntimeError(f"Invalid lambda {pos_line}")

    return pos, move, velocity, magmom, angle1, angle2, constrain, lambda1

def get_atom_mag_cartesian(atommag, angle1, angle2):
    """Transform atommag, angle1, angle2 to magmom in cartesian coordinates.

    Parameters
    ----------
    atommag : float/list of float/None
        Atom magnetic moment.
    angle1 : float/None
        value of angle1.
    angle2 : float/None
        value of angle2.
    ABACUS support defining mag, angle1, angle2 at the same time.
    angle1 is the angle between z-axis and real spin (in degrees).
    angle2 is the angle between x-axis and real spin projection in xy-plane (in degrees).
    If only mag is defined, then transfer it to magmom directly.
    And if mag, angle1, angle2 are defined, then mag is only the norm of magmom, and the direction is defined by angle1 and angle2.
    """
    if atommag is None:
        return None
    if not (isinstance(atommag, list) or isinstance(atommag, float)):
        raise RuntimeError(f"Invalid atommag: {atommag}")

    if angle1 is None and angle2 is None:
        if isinstance(atommag, list):
            return atommag
        else:
            return [0, 0, atommag]
    else:
        a1 = 0
        a2 = 0
        if angle1 is not None:
            a1 = angle1
        if angle2 is not None:
            a2 = angle2
        if isinstance(atommag, list):
            mag_norm = np.linalg.norm(atommag)
        else:
            mag_norm = atommag
        return [
            mag_norm * np.sin(np.radians(a1)) * np.cos(np.radians(a2)),
            mag_norm * np.sin(np.radians(a1)) * np.sin(np.radians(a2)),
            mag_norm * np.cos(np.radians(a1)),
        ]

def get_carteisan_coords(coords, coord_type, celldm, cell):
    """Transform the atomic coordinates to cartesian coordinates.
    
    Args:
        coords (np.ndarray): atomic coordinates read from the STRU file.
        coord_type (str): the coordination type, either "cartesian" or "direct".
        celldm (float): the lattice constant.
        cell (np.ndarray): the cell vectors in angstrom.
        
    Returns:
        np.ndarray: the cartesian coordinates in angstrom.
    """
    if coord_type == "cartesian":
        return coords * celldm * bohr2ang
    elif coord_type == "direct":
        return np.matmul(coords, cell)
    else:
        raise RuntimeError(f"Invalid coordination type: {coord_type}")


def parse_pos(coords_lines, atom_names, celldm, cell):
    """Read the atomic positions block in the ABACUS STRU file.
    
    Args:
        coords_lines (list): list of lines in the atomic positions block.
        atom_names (list): list of atom names.
        celldm (float): the lattice constant.
        cell (np.ndarray): the cell vectors in angstrom, and has multipy celldm.
        
    Returns:
        tuple: tuple of atom_numbs, coords, move, mags, velocity, sc, lambda_
        
    Note: for atomic magnetic moment, we finnaly transform it to non-collinear magnetic moment in cartesian coordinates,
          and do not return the angle1 and angle2, and the magnetic moment of each atom type.
        
    """ 
    coord_type = coords_lines[0].split()[0].lower()  # cartisan or direct
    atom_numbs = []  # the number of each atom type
    coords = []  # coordinations of atoms
    move = []  # move flag of each atom
    velocity = []  # velocity of each atom
    mags = []  # magnetic moment of each atom
    sc = []  # spin constraint flag of each atom
    lambda_ = []  # lambda of each atom

    ntype = len(atom_names)
    line_idx = 1  # starting line of first element
    for it in range(ntype):
        atom_name = coords_lines[line_idx].split()[0]
        if atom_name != atom_names[it]:
            raise RuntimeError(f"Read atom name '{atom_name}' is not equal to the expected atom name '{atom_names[it]}'")
        atom_type_mag = float(coords_lines[line_idx + 1].split()[0])
        line_idx += 2
        atom_numbs.append(int(coords_lines[line_idx].split()[0]))
        line_idx += 1
        for iline in range(atom_numbs[it]):
            pos, imove, ivelocity, imagmom, iangle1, iangle2, iconstrain, ilambda1 = (
                parse_pos_oneline(coords_lines[line_idx])
            )
            
            coords.append(get_carteisan_coords(np.array(pos), coord_type, celldm, cell))

            move.append(imove)
            velocity.append(ivelocity)
            sc.append(iconstrain)
            lambda_.append(ilambda1)

            # calculate the magnetic moment in cartesian coordinates
            mag = get_atom_mag_cartesian(imagmom, iangle1, iangle2)
            if mag is None:
                mag = [0, 0, atom_type_mag]
            mags.append(mag)

            line_idx += 1
    coords = np.array(coords)  # need transformation!!!
    
    if all([i is None for i in move]):
        move = []
    else:
        move = np.array(move, dtype=bool)
    
    if all([i is None for i in velocity]):
        velocity = []
    else:
        velocity = np.array(velocity)
        
    if all([i is None for i in sc]):
        sc = []
    
    if all([i is None for i in lambda_]):
        lambda_ = []
        
    mags = np.array(mags)
    
    return atom_numbs, coords, move, mags, velocity, sc, lambda_

def get_frame_from_stru(stru):
    """Read the ABACUS STRU file and return the dpdata frame.
    
    The description of ABACUS STRU can be found in https://abacus.deepmodeling.com/en/latest/advanced/input_files/stru.html
    
    Args:
        stru (str): path to the ABACUS STRU file.
        
    Returns:
        data: the parsed stru information in dictionary.
        {
            "atom_names": list of atom names,
            "atom_numbs": list of atom numbers,
            "atom_types": list of atom types,
            "masses": list of atomic masses,
            "pp_files", list of pseudo potential files,
            "orb_files", list of orbital files,
            "dpks_descriptor": the deepks descriptor file,
            
            # below are the information in each frame
            
            "cells": list of cell vectors,
            "coords": list of atomic coordinates,
            "spins": list of magnetic moments,
            "moves": list of move flags,  
        } 
        For some keys, if the information is not provided in the STRU file, then it will not be included in the dictionary.
    """
    if not os.path.isfile(stru):
        raise FileNotFoundError(f"ABACUS STRU file {stru} not found!!!")
    
    # 1. read the file and split the lines to blocks
    with open(stru) as f:
        lines = f.readlines()
    blocks = split_stru_block(lines)
        
    # 2. parse the blocks
    atom_names, masses, pp_files = parse_atomic_species_block(blocks["ATOMIC_SPECIES"])
    orb_files = parse_numerical_orbital_block(blocks.get("NUMERICAL_ORBITAL", []))
    dpks_descriptor = blocks.get("NUMERICAL_DESCRIPTOR", [])
    celldm = parse_lattice_constant_block(blocks["LATTICE_CONSTANT"])
    cell = parse_lattice_vectors_block(blocks["LATTICE_VECTORS"])
    cell = np.array(cell) * celldm * bohr2ang
    atom_numbs, coords, move, mags, velocity, sc, lambda_ = parse_pos(
        blocks["ATOMIC_POSITIONS"], atom_names, celldm, cell
    )
    
    data = {
        "atom_names": atom_names,
        "atom_numbs": atom_numbs,
        "atom_types": np.array([i for i in range(len(atom_numbs)) for j in range(atom_numbs[i])]),
        "masses": np.array(masses),
        "pp_files": pp_files,
        "cells": np.array([cell]),
        "coords": np.array([coords]),
        "spins": np.array([mags]),
    }
    if len(orb_files) > 0:
        data["orb_files"] = orb_files
    if len(dpks_descriptor) > 0:
        data["dpks_descriptor"] = dpks_descriptor[0].strip()
    if len(move) > 0:
        data["move"] = np.array([move])
    
    return data

def make_unlabeled_stru(
    data,
    frame_idx,
    pp_file=None,
    numerical_orbital=None,
    numerical_descriptor=None,
    mass=None,
    move=None,
    velocity=None,
    mag=None,
    angle1=None,
    angle2=None,
    sc=None,
    lambda_=None,
    link_file=False,
    dest_dir=None,
    **kwargs,
):
    """Make an unlabeled STRU file from a dictionary.

    Parameters
    ----------
    data : dict
        System data
    frame_idx : int
        The index of the frame to dump
    pp_file : list of string or dict
        List of pseudo potential files, or a dictionary of pseudo potential files for each atomnames
    numerical_orbital : list of string or dict, optional
        List of orbital files, or a dictionary of orbital files for each atomnames
    numerical_descriptor : str, optional
        numerical descriptor file
    mass : list of float, optional
        List of atomic masses
    move : list of (list of list of bool), optional
        List of the move flag of each xyz direction of each atom for each frame
    velocity : list of list of float, optional
        List of the velocity of each xyz direction of each atom
    mag : list of (list of float or float), optional
        List of the magnetic moment of each atom, can be a list of three floats or one float
        For noncollinear, three floats are the xyz component of the magnetic moment.
        For collinear, one float is the norm of the magnetic moment.
    angle1 : list of float, optional
        List of the angle1 of each atom. For noncollinear calculation, it is the angle between the magnetic moment and the z-axis.
    angle2 : list of float, optional
        List of the angle2 of each atom. For noncollinear calculation, it is the angle between the projection of magnetic moment on xy plane and the x-axis.
    sc : list of (bool or list of 3 bool), optional
        List of the spin constraint flag of each atom. Each element can be a bool or a list of three bools or None.
    lambda_ : list of (float or list of 3 float), optional
        List of the lambda of each atom. Each element can be a float or a list of three floats.
    link_file : bool, optional
        Whether to link the pseudo potential files and orbital files in the STRU file.
        If True, then only filename will be written in the STRU file, and make a soft link to the real file.
    dest_dir : str, optional
        The destination directory to make the soft link of the pseudo potential files and orbital files.
    For velocity, mag, angle1, angle2, sc, and lambda_, if the value is None, then the corresponding information will not be written.
    ABACUS support defining "mag" and "angle1"/"angle2" at the same time, and in this case, the "mag" only define the norm of the magnetic moment, and "angle1" and "angle2" define the direction of the magnetic moment.
    If data has spins, then it will be written as mag to STRU file; while if mag is passed at the same time, then mag will be used.
    """

    def _link_file(dest_dir, src_file):
        if not os.path.isfile(src_file):
            print(f"ERROR: link_file: {src_file} is not a file.")
            return False
        src_file = os.path.abspath(src_file)
        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)
        dest_file = os.path.join(dest_dir, os.path.basename(src_file))
        if os.path.isfile(dest_file):
            if os.path.samefile(src_file, dest_file):
                return True
            else:
                os.remove(dest_file)
        os.symlink(src_file, dest_file)
        return True

    def ndarray2list(i):
        if isinstance(i, np.ndarray):
            return i.tolist()
        else:
            return i

    def process_file_input(file_input, atom_names, input_name):
        # For pp_file and numerical_orbital, process the file input, and return a list of file names
        # file_input can be a list of file names, or a dictionary of file names for each atom names
        if isinstance(file_input, (list, tuple)):
            if len(file_input) != len(atom_names):
                raise ValueError(
                    f"{input_name} length is not equal to the number of atom types"
                )
            return file_input
        elif isinstance(file_input, dict):
            for element in atom_names:
                if element not in file_input:
                    raise KeyError(f"{input_name} does not contain {element}")
            return [file_input[element] for element in atom_names]
        else:
            raise ValueError(f"Invalid {input_name}: {file_input}")

    if link_file and dest_dir is None:
        print(
            "WARNING: make_unlabeled_stru: link_file is True, but dest_dir is None. Will write the filename to STRU but not making soft link."
        )
    if dest_dir is not None and dest_dir.strip() == "":
        dest_dir = "."
    
    # check the input data
    if mass is None and data.get("masses") is not None and len(data["masses"]) > 0:
        mass = data["masses"]
        
    if pp_file is None and data.get("pp_files") is not None and len(data["pp_files"]) > 0:
        pp_file = data["pp_files"]
        
    if numerical_orbital is None and data.get("orb_files") is not None and len(data["orb_files"]) > 0:
        numerical_orbital = data["orb_files"]
    
    if numerical_descriptor is None and data.get("dpks_descriptor") is not None:
        numerical_descriptor = data["dpks_descriptor"]

    if mag is None and data.get("spins") is not None and len(data["spins"]) > 0:
        mag = data["spins"][frame_idx]

    if move is None and data.get("move", None) is not None and len(data["move"]) > 0:
        move = data["move"][frame_idx]

    # check the length of the input data
    atom_numbs = sum(data["atom_numbs"])
    for key in [move, velocity, mag, angle1, angle2, sc, lambda_]:
        if key is not None:
            if (
                not isinstance(ndarray2list(key), (list, tuple))
                and len(key) != atom_numbs
            ):
                key_name = [name for name, value in locals().items() if value is key][0]
                print(
                    f"ERROR: make_unlabeled_stru: the length of '{key_name}' ({len(key)}) should be equal to the number of atom number ({atom_numbs})."
                )
                return ""

    # ATOMIC_SPECIES block
    out = "ATOMIC_SPECIES\n"
    if pp_file is not None:
        ppfiles = process_file_input(
            ndarray2list(pp_file), data["atom_names"], "pp_file"
        )
    else:
        warnings.warn(
            "pp_file is not provided, will use empty string for pseudo potential file."
        )
        ppfiles = [""] * len(data["atom_names"])

    for iele in range(len(data["atom_names"])):
        if data["atom_numbs"][iele] == 0:
            continue
        out += data["atom_names"][iele] + " "
        if mass is not None:
            out += f"{mass[iele]:.3f} "
        else:
            out += "1 "

        ipp_file = ppfiles[iele]
        if ipp_file != "":
            if not link_file:
                out += ipp_file
            else:
                out += os.path.basename(ipp_file.rstrip("/"))
                if dest_dir is not None:
                    _link_file(dest_dir, ipp_file)
        out += "\n"
    out += "\n"

    # NUMERICAL_ORBITAL block
    if numerical_orbital is not None:
        numerical_orbital = ndarray2list(numerical_orbital)
        orbfiles = process_file_input(
            numerical_orbital, data["atom_names"], "numerical_orbital"
        )
        orbfiles = [
            orbfiles[i]
            for i in range(len(data["atom_names"]))
            if data["atom_numbs"][i] != 0
        ]
        out += "NUMERICAL_ORBITAL\n"
        for iorb in orbfiles:
            if not link_file:
                out += iorb
            else:
                out += os.path.basename(iorb.rstrip("/"))
                if dest_dir is not None:
                    _link_file(dest_dir, iorb)
            out += "\n"
        out += "\n"

    # deepks block
    if numerical_descriptor is not None:
        assert isinstance(numerical_descriptor, str)
        if not link_file:
            out += f"NUMERICAL_DESCRIPTOR\n{numerical_descriptor}\n"
        else:
            out += f"NUMERICAL_DESCRIPTOR\n{os.path.basename(numerical_descriptor)}\n"
            if dest_dir is not None:
                _link_file(dest_dir, numerical_descriptor)
        out += "\n"

    # LATTICE_CONSTANT and LATTICE_VECTORS block
    out += "LATTICE_CONSTANT\n"
    out += str(1 / bohr2ang) + "\n\n"

    out += "LATTICE_VECTORS\n"
    for ix in range(3):
        for iy in range(3):
            out += str(data["cells"][frame_idx][ix][iy]) + " "
        out += "\n"
    out += "\n"

    # ATOMIC_POSITIONS block
    out += "ATOMIC_POSITIONS\n"
    out += "Cartesian    # Cartesian(Unit is LATTICE_CONSTANT)\n"
    # ret += "\n"
    natom_tot = 0  # in for loop, it is also the atom index
    for iele in range(len(data["atom_names"])):
        if data["atom_numbs"][iele] == 0:
            continue
        out += data["atom_names"][iele] + "\n"
        out += "0.0\n"
        out += str(data["atom_numbs"][iele]) + "\n"
        for iatom in range(data["atom_numbs"][iele]):
            iatomtype = np.nonzero(data["atom_types"] == iele)[0][iatom]
            iout = f"{data['coords'][frame_idx][iatomtype, 0]:.12f} {data['coords'][frame_idx][iatomtype, 1]:.12f} {data['coords'][frame_idx][iatomtype, 2]:.12f}"
            # add flags for move, velocity, mag, angle1, angle2, and sc
            if move is not None:
                if (
                    isinstance(ndarray2list(move[natom_tot]), (list, tuple))
                    and len(move[natom_tot]) == 3
                ):
                    iout += " " + " ".join(
                        ["1" if ii else "0" for ii in move[natom_tot]]
                    )
                elif isinstance(ndarray2list(move[natom_tot]), (int, float, bool)):
                    iout += " 1 1 1" if move[natom_tot] else " 0 0 0"
            else:
                iout += " 1 1 1"

            if (
                velocity is not None
                and isinstance(ndarray2list(velocity[natom_tot]), (list, tuple))
                and len(velocity[natom_tot]) == 3
            ):
                iout += " v " + " ".join([f"{ii:.12f}" for ii in velocity[natom_tot]])

            if mag is not None:
                if isinstance(ndarray2list(mag[natom_tot]), (list, tuple)) and len(
                    mag[natom_tot]
                ) in [1, 3]:
                    iout += " mag " + " ".join([f"{ii:.12f}" for ii in mag[natom_tot]])
                elif isinstance(ndarray2list(mag[natom_tot]), (int, float)):
                    iout += " mag " + f"{mag[natom_tot]:.12f}"

            if angle1 is not None and isinstance(
                ndarray2list(angle1[natom_tot]), (int, float)
            ):
                iout += " angle1 " + f"{angle1[natom_tot]:.12f}"

            if angle2 is not None and isinstance(
                ndarray2list(angle2[natom_tot]), (int, float)
            ):
                iout += " angle2 " + f"{angle2[natom_tot]:.12f}"

            if sc is not None:
                if isinstance(ndarray2list(sc[natom_tot]), (list, tuple)) and len(
                    sc[natom_tot]
                ) in [1, 3]:
                    iout += " sc " + " ".join(
                        ["1" if ii else "0" for ii in sc[natom_tot]]
                    )
                elif isinstance(ndarray2list(sc[natom_tot]), (int, float, bool)):
                    iout += " sc " + "1" if sc[natom_tot] else "0"

            if lambda_ is not None:
                if isinstance(ndarray2list(lambda_[natom_tot]), (list, tuple)) and len(
                    lambda_[natom_tot]
                ) in [1, 3]:
                    iout += " lambda " + " ".join(
                        [f"{ii:.12f}" for ii in lambda_[natom_tot]]
                    )
                elif isinstance(ndarray2list(lambda_[natom_tot]), (int, float)):
                    iout += " lambda " + f"{lambda_[natom_tot]:.12f}"

            out += iout + "\n"
            natom_tot += 1
    assert natom_tot == sum(data["atom_numbs"])
    return out