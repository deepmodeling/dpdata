import numpy as np

from dpdata.periodic_table import Element


def elements_index_map(elements, standard=False, inverse=False):
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
    """Sort atom_names of the system and reorder atom_numbs and atom_types accoarding
    to atom_names. If type_map is not given, atom_names will be sorted by
    alphabetical order. If type_map is given, atom_names will be type_map.

    Parameters
    ----------
    data : dict
        system data
    type_map : list
        type_map
    """
    if type_map is not None:
        # assign atom_names index to the specify order
        # atom_names must be a subset of type_map
        assert set(data["atom_names"]).issubset(set(type_map))
        # for the condition that type_map is a proper superset of atom_names
        # new_atoms = set(type_map) - set(data["atom_names"])
        new_atoms = [e for e in type_map if e not in data["atom_names"]]
        if new_atoms:
            data = add_atom_names(data, new_atoms)
        # index that will sort an array by type_map
        # a[as[a]] == b[as[b]]  as == argsort
        # as[as[b]] == as^{-1}[b]
        # a[as[a][as[as[b]]]] = b[as[b][as^{-1}[b]]] = b[id]
        idx = np.argsort(data["atom_names"], kind="stable")[
            np.argsort(np.argsort(type_map, kind="stable"), kind="stable")
        ]
    else:
        # index that will sort an array by alphabetical order
        idx = np.argsort(data["atom_names"], kind="stable")
    # sort atom_names, atom_numbs, atom_types by idx
    data["atom_names"] = list(np.array(data["atom_names"])[idx])
    data["atom_numbs"] = list(np.array(data["atom_numbs"])[idx])
    data["atom_types"] = np.argsort(idx, kind="stable")[data["atom_types"]]
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
