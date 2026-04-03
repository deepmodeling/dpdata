from __future__ import annotations

import io
import os
from contextlib import contextmanager
from typing import TYPE_CHECKING, Literal, overload

import numpy as np

from dpdata.periodic_table import Element


@overload
def elements_index_map(
    elements: list[str], standard: bool, inverse: Literal[True]
) -> dict[int, str]: ...
@overload
def elements_index_map(
    elements: list[str], standard: bool, inverse: Literal[False] = ...
) -> dict[str, int]: ...
@overload
def elements_index_map(
    elements: list[str], standard: bool, inverse: bool = False
) -> dict[str, int] | dict[int, str]: ...


def elements_index_map(
    elements: list[str], standard: bool = False, inverse: bool = False
) -> dict:
    if standard:
        elements.sort(key=lambda x: Element(x).Z)
    if inverse:
        return dict(zip(range(len(elements)), elements))
    else:
        return dict(zip(elements, range(len(elements))))


# %%


def remove_pbc(system, protect_layer=9):
    nframes = len(system["coords"])
    natoms = len(system["coords"][0])
    for ff in range(nframes):
        tmpcoord = system["coords"][ff]
        cog = np.average(tmpcoord, axis=0)
        dist = tmpcoord - np.tile(cog, [natoms, 1])
        max_dist = np.max(np.linalg.norm(dist, axis=1))
        h_cell_size = max_dist + protect_layer
        cell_size = h_cell_size * 2
        shift = np.array([1, 1, 1]) * h_cell_size - cog
        system["coords"][ff] = system["coords"][ff] + np.tile(shift, [natoms, 1])
        system["cells"][ff] = cell_size * np.eye(3)
    return system


def add_atom_names(data, atom_names):
    """Add atom_names that do not exist."""
    data["atom_names"].extend(atom_names)
    data["atom_numbs"].extend([0 for _ in atom_names])
    return data


def sort_atom_names(data, type_map=None):
    """Sort atom_names of the system and reorder atom_numbs and atom_types according
    to atom_names. If type_map is not given, atom_names will be sorted by
    alphabetical order. If type_map is given, atom_names will be set to type_map,
    and zero-count elements are kept.

    Parameters
    ----------
    data : dict
        system data
    type_map : list
        type_map
    """
    if type_map is not None:
        # assign atom_names index to the specified order
        # only active (numb > 0) atom names must be in type_map
        orig_names = data["atom_names"]
        orig_numbs = data["atom_numbs"]
        active_names = {name for name, numb in zip(orig_names, orig_numbs) if numb > 0}
        type_map_set = set(type_map)
        if not active_names.issubset(type_map_set):
            missing = active_names - type_map_set
            raise ValueError(f"Active atom types {missing} not in provided type_map.")

        # for the condition that type_map is a proper superset of atom_names,
        # we allow new elements with atom_numb = 0
        new_names = list(type_map)
        new_numbs = []
        name_to_old_idx = {name: i for i, name in enumerate(orig_names)}

        for name in new_names:
            if name in name_to_old_idx:
                new_numbs.append(orig_numbs[name_to_old_idx[name]])
            else:
                new_numbs.append(0)

        # build mapping from old atom type index to new one
        # old_types[i] = j  -->  new_types[i] = type_map.index(atom_names[j])
        old_to_new_index = {}
        for old_idx, name in enumerate(orig_names):
            if name in type_map_set:
                new_idx = type_map.index(name)
                old_to_new_index[old_idx] = new_idx

        # remap atom_types using the index mapping
        old_types = np.array(data["atom_types"])
        new_types = np.empty_like(old_types)
        for old_idx, new_idx in old_to_new_index.items():
            new_types[old_types == old_idx] = new_idx

        # update data in-place
        data["atom_names"] = new_names
        data["atom_numbs"] = new_numbs
        data["atom_types"] = new_types

    else:
        # index that will sort an array by alphabetical order
        # idx = argsort(atom_names)  -->  atom_names[idx] is sorted
        idx = np.argsort(data["atom_names"], kind="stable")
        # sort atom_names and atom_numbs by idx
        data["atom_names"] = list(np.array(data["atom_names"])[idx])
        data["atom_numbs"] = list(np.array(data["atom_numbs"])[idx])
        # to update atom_types: we need the inverse permutation of idx
        # because if old_type = i, and atom_names[i] moves to position j,
        # then the new type should be j.
        # inv_idx = argsort(idx) satisfies: inv_idx[idx[i]] = i
        inv_idx = np.argsort(idx, kind="stable")
        data["atom_types"] = inv_idx[data["atom_types"]]

    return data


def uniq_atom_names(data):
    """Make the atom names uniq. For example
    ['O', 'H', 'O', 'H', 'O'] -> ['O', 'H'].

    Parameters
    ----------
    data : dict
        data dict of `System`, `LabeledSystem`

    """
    unames = []
    uidxmap = []
    for idx, ii in enumerate(data["atom_names"]):
        if ii not in unames:
            unames.append(ii)
        uidxmap.append(unames.index(ii))
    data["atom_names"] = unames
    tmp_type = list(data["atom_types"]).copy()
    data["atom_types"] = np.array([uidxmap[jj] for jj in tmp_type], dtype=int)
    data["atom_numbs"] = [
        sum(ii == data["atom_types"]) for ii in range(len(data["atom_names"]))
    ]
    return data


def utf8len(s: str) -> int:
    """Return the byte length of a string."""
    return len(s.encode("utf-8"))


if TYPE_CHECKING:
    from collections.abc import Generator

    FileType = io.IOBase | str | os.PathLike


@contextmanager
def open_file(file: FileType, *args, **kwargs) -> Generator[io.IOBase, None, None]:
    """A context manager that yields a file object.

    Parameters
    ----------
    file : file object or file path
        A file object or a file path.

    Yields
    ------
    file : io.IOBase
        A file object.
    *args
        parameters to open
    **kwargs
        other parameters

    Raises
    ------
    ValueError
        If file is not a file object or a file

    Examples
    --------
    >>> with open_file("file.txt") as file:
    ...     print(file.read())
    """
    if isinstance(file, io.IOBase):
        yield file
    elif isinstance(file, (str, os.PathLike)):
        with open(file, *args, **kwargs) as f:
            yield f
    else:
        raise ValueError("file must be a file object or a file path.")
