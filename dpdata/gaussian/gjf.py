# The initial code of this file is based on
# https://github.com/deepmodeling/dpgen/blob/0767dce7cad29367edb2e4a55fd0d8724dbda642/dpgen/generator/lib/gaussian.py#L1-L190
# under LGPL 3.0 license
"""Generate Gaussian input file."""

import itertools
import re
import uuid
import warnings
from typing import List, Optional, Tuple, Union

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

try:
    from openbabel import openbabel
except ImportError:
    try:
        import openbabel
    except ImportError:
        openbabel = None
from dpdata.periodic_table import Element


def _crd2frag(symbols: List[str], crds: np.ndarray) -> Tuple[int, List[int]]:
    """Detect fragments from coordinates.

    Parameters
    ----------
    symbols : list[str]
        element symbols; virtual elements are not supported
    crds : np.ndarray
        atomic coordinates, shape: (N, 3)

    Returns
    -------
    frag_numb : int
        number of fragments
    frag_index : list[int]
        frament index that each atom belongs to

    Notes
    -----
    In this method, Open Babel is used to detect bond connectivity. The threshold
    is the sum of covalent radii with a slight tolerance (0.45 A). Note that
    this threshold has errors.

    PBC support is removed from this method as Gaussian does not support PBC calculation.

    Raises
    ------
    ImportError
        if Open Babel is not installed
    """
    if openbabel is None:
        raise ImportError(
            "Open Babel (Python interface) should be installed to detect fragmentation!"
        )
    atomnumber = len(symbols)
    # Use openbabel to connect atoms
    mol = openbabel.OBMol()
    mol.BeginModify()
    for idx, (symbol, position) in enumerate(zip(symbols, crds.astype(np.float64))):
        num = Element(symbol).Z
        atom = mol.NewAtom(idx)
        atom.SetAtomicNum(int(num))
        atom.SetVector(*position)
    mol.ConnectTheDots()
    mol.PerceiveBondOrders()
    mol.EndModify()
    bonds = []
    for ii in range(mol.NumBonds()):
        bond = mol.GetBond(ii)
        a = bond.GetBeginAtom().GetId()
        b = bond.GetEndAtom().GetId()
        bo = bond.GetBondOrder()
        bonds.extend([[a, b, bo], [b, a, bo]])
    bonds = np.array(bonds, ndmin=2).reshape((-1, 3))
    graph = csr_matrix(
        (bonds[:, 2], (bonds[:, 0], bonds[:, 1])), shape=(atomnumber, atomnumber)
    )
    frag_numb, frag_index = connected_components(graph, 0)
    return frag_numb, frag_index


def detect_multiplicity(symbols: np.ndarray) -> int:
    """Find the minimal multiplicity of the given molecules.

    Parameters
    ----------
    symbols : np.ndarray
        element symbols; virtual elements are not supported

    Returns
    -------
    int
        spin multiplicity
    """
    # currently only support charge=0
    # oxygen -> 3
    if np.count_nonzero(symbols == ["O"]) == 2 and symbols.size == 2:
        return 3
    # calculates the total number of electrons, assumes they are paired as much as possible
    n_total = sum([Element(s).Z for s in symbols])
    return n_total % 2 + 1


def make_gaussian_input(
    sys_data: dict,
    keywords: Union[str, List[str]],
    multiplicity: Union[str, int] = "auto",
    charge: int = 0,
    fragment_guesses: bool = False,
    basis_set: Optional[str] = None,
    keywords_high_multiplicity: Optional[str] = None,
    nproc: int = 1,
) -> str:
    """Make gaussian input file.

    Parameters
    ----------
    sys_data : dict
        system data
    keywords : str or list[str]
        Gaussian keywords, e.g. force b3lyp/6-31g**. If a list,
        run multiple steps
    multiplicity : str or int, default=auto
        spin multiplicity state. It can be a number. If auto,
        multiplicity will be detected automatically, with the
        following rules:
            fragment_guesses=True
                multiplicity will +1 for each radical, and +2
                for each oxygen molecule
            fragment_guesses=False
                multiplicity will be 1 or 2, but +2 for each
                oxygen molecule
    charge : int, default=0
        molecule charge. Only used when charge is not provided
        by the system
    fragment_guesses : bool, default=False
        initial guess generated from fragment guesses. If True,
        multiplicity should be auto
    basis_set : str, default=None
        custom basis set
    keywords_high_multiplicity : str, default=None
        keywords for points with multiple raicals. multiplicity
        should be auto. If not set, fallback to normal keywords
    nproc : int, default=1
        Number of CPUs to use

    Returns
    -------
    str
        gjf output string
    """
    coordinates = sys_data["coords"][0]
    atom_names = sys_data["atom_names"]
    atom_numbs = sys_data["atom_numbs"]
    atom_types = sys_data["atom_types"]
    # get atom symbols list
    symbols = [atom_names[atom_type] for atom_type in atom_types]

    # assume default charge is zero and default spin multiplicity is 1
    if "charge" in sys_data.keys():
        charge = sys_data["charge"]

    use_fragment_guesses = False
    if isinstance(multiplicity, int):
        mult_auto = False
    elif multiplicity == "auto":
        mult_auto = True
    else:
        raise RuntimeError('The keyword "multiplicity" is illegal.')

    if fragment_guesses:
        # Initial guess generated from fragment guesses
        # New feature of Gaussian 16
        use_fragment_guesses = True
        if not mult_auto:
            warnings.warn("Automatically set multiplicity to auto!")
            mult_auto = True

    if mult_auto:
        frag_numb, frag_index = _crd2frag(symbols, coordinates)
        if frag_numb == 1:
            use_fragment_guesses = False
        mult_frags = []
        for i in range(frag_numb):
            idx = frag_index == i
            mult_frags.append(detect_multiplicity(np.array(symbols)[idx]))
        if use_fragment_guesses:
            multiplicity = sum(mult_frags) - frag_numb + 1 - charge % 2
            chargekeywords_frag = "%d %d" % (charge, multiplicity) + "".join(
                [" %d %d" % (charge, mult_frag) for mult_frag in mult_frags]
            )
        else:
            multi_frags = np.array(mult_frags)
            multiplicity = (
                1
                + np.count_nonzero(multi_frags == 2) % 2
                + np.count_nonzero(multi_frags == 3) * 2
                - charge % 2
            )

        if (
            keywords_high_multiplicity is not None
            and np.count_nonzero(multi_frags == 2) >= 2
        ):
            # at least 2 radicals
            keywords = keywords_high_multiplicity

    if isinstance(keywords, str):
        keywords = [keywords]
    else:
        keywords = keywords.copy()

    buff = []
    # keywords, e.g., force b3lyp/6-31g**
    if use_fragment_guesses:
        keywords[0] = f"{keywords[0]} guess=fragment={frag_numb}"

    chkkeywords = []
    if len(keywords) > 1:
        chkkeywords.append(f"%chk={str(uuid.uuid1())}.chk")

    nprockeywords = f"%nproc={nproc:d}"
    # use formula as title
    titlekeywords = "".join(
        [f"{symbol}{numb}" for symbol, numb in zip(atom_names, atom_numbs)]
    )
    chargekeywords = f"{charge} {multiplicity}"

    buff = [
        *chkkeywords,
        nprockeywords,
        f"#{keywords[0]}",
        "",
        titlekeywords,
        "",
        (chargekeywords_frag if use_fragment_guesses else chargekeywords),
    ]

    for ii, (symbol, coordinate) in enumerate(zip(symbols, coordinates)):
        if use_fragment_guesses:
            buff.append(
                "%s(Fragment=%d) %f %f %f" % (symbol, frag_index[ii] + 1, *coordinate)
            )
        else:
            buff.append("{} {:f} {:f} {:f}".format(symbol, *coordinate))
    if not sys_data.get("nopbc", False):
        # PBC condition
        cell = sys_data["cells"][0]
        for ii in range(3):
            # use TV as atomic symbol, see https://gaussian.com/pbc/
            buff.append("TV {:f} {:f} {:f}".format(*cell[ii]))
    if basis_set is not None:
        # custom basis set
        buff.extend(["", basis_set, ""])
    for kw in itertools.islice(keywords, 1, None):
        buff.extend(
            [
                "\n--link1--",
                *chkkeywords,
                nprockeywords,
                f"#{kw}",
                "",
                titlekeywords,
                "",
                chargekeywords,
                "",
            ]
        )
    buff.append("\n")
    return "\n".join(buff)


def read_gaussian_input(inp: str):
    """Read Gaussian input.

    Parameters
    ----------
    inp : str
        Gaussian input str

    Returns
    -------
    dict
        system data
    """
    flag = 0
    coords = []
    elements = []
    cells = []
    for line in inp.split("\n"):
        if not line.strip():
            # empty line
            flag += 1
        elif flag == 0:
            # keywords
            if line.startswith("#"):
                # setting
                keywords = line.split()
            elif line.startswith("%"):
                pass
        elif flag == 1:
            # title
            pass
        elif flag == 2:
            # multi and coords
            s = line.split()
            if len(s) == 2:
                pass
            elif len(s) == 4:
                if s[0] == "TV":
                    cells.append(list(map(float, s[1:4])))
                else:
                    # element
                    elements.append(re.sub("\\(.*?\\)|\\{.*?}|\\[.*?]", "", s[0]))
                    coords.append(list(map(float, s[1:4])))
        elif flag == 3:
            # end
            break
    atom_names, atom_types, atom_numbs = np.unique(
        elements, return_inverse=True, return_counts=True
    )
    if len(cells):
        nopbc = False
    else:
        nopbc = True
        cells = np.array([np.eye(3)]) * 100
    return {
        "atom_names": list(atom_names),
        "atom_numbs": list(atom_numbs),
        "atom_types": atom_types,
        "cells": np.array(cells).reshape(1, 3, 3),
        "nopbc": nopbc,
        "coords": np.array(coords).reshape(1, -1, 3),
        "orig": np.zeros(3),
    }
