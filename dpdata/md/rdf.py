import numpy as np


def rdf(sys, sel_type=[None, None], max_r=5, nbins=100):
    """Compute the rdf of a system.

    Parameters
    ----------
    sys : System or LabeledSystem
        The dpdata system
    sel_type : list
        List of size 2. The first element specifies the type of the first atom,
        while the second element specifies the type of the second atom.
        Both elements can be ints or list of ints.
        If the element is None, all types are specified.
        Examples are sel_type = [0, 0], sel_type = [0, [0, 1]] or sel_type = [0, None]
    max_r : float
        Maximal range of rdf calculation
    nbins : int
        Number of bins for rdf calculation

    Returns
    -------
    xx: np.array
        The lattice of r
    rdf: np.array
        The value of rdf at r
    coord: np.array
        The coordination number up to r
    """
    return compute_rdf(
        sys["cells"],
        sys["coords"],
        sys["atom_types"],
        sel_type=sel_type,
        max_r=max_r,
        nbins=nbins,
    )


def compute_rdf(box, posis, atype, sel_type=[None, None], max_r=5, nbins=100):
    nframes = box.shape[0]
    xx = None
    all_rdf = []
    all_cod = []
    for ii in range(nframes):
        xx, rdf, cod = _compute_rdf_1frame(
            box[ii], posis[ii], atype, sel_type, max_r, nbins
        )
        all_rdf.append(rdf)
        all_cod.append(cod)
    all_rdf = np.array(all_rdf).reshape([nframes, -1])
    all_cod = np.array(all_cod).reshape([nframes, -1])
    all_rdf = np.average(all_rdf, axis=0)
    all_cod = np.average(all_cod, axis=0)
    return xx, all_rdf, all_cod


def _compute_rdf_1frame(box, posis, atype, sel_type=[None, None], max_r=5, nbins=100):
    all_types = list(set(list(np.sort(atype, kind="stable"))))
    if sel_type[0] is None:
        sel_type[0] = all_types
    if sel_type[1] is None:
        sel_type[1] = all_types
    if not isinstance(sel_type[0], list):
        sel_type[0] = [sel_type[0]]
    if not isinstance(sel_type[1], list):
        sel_type[1] = [sel_type[1]]
    natoms = len(posis)
    import ase.neighborlist
    from ase import Atoms

    atoms = Atoms(positions=posis, cell=box, pbc=[1, 1, 1])
    nlist = ase.neighborlist.NeighborList(
        max_r,
        self_interaction=False,
        bothways=True,
        primitive=ase.neighborlist.NewPrimitiveNeighborList,
    )
    nlist.update(atoms)
    stat = np.zeros(nbins)
    stat_acc = np.zeros(nbins)
    hh = max_r / float(nbins)
    for ii in range(natoms):
        # atom "0"
        if atype[ii] in sel_type[0]:
            indices, offsets = nlist.get_neighbors(ii)
            for jj, os in zip(indices, offsets):
                # atom "1"
                if atype[jj] in sel_type[1]:
                    posi_jj = atoms.positions[jj] + np.dot(os, atoms.get_cell())
                    diff = posi_jj - atoms.positions[ii]
                    dr = np.linalg.norm(diff)
                    # if (np.linalg.norm(diff- diff_1)) > 1e-12 :
                    #     raise RuntimeError
                    si = int(dr / hh)
                    if si < nbins:
                        stat[si] += 1
    # count the number of atom1
    c0 = 0
    for ii in sel_type[0]:
        c0 += np.sum(atype == ii)
    # count the number of atom1
    c1 = 0
    for ii in sel_type[1]:
        c1 += np.sum(atype == ii)
    rho1 = c1 / np.linalg.det(box)
    # compute coordination number
    for ii in range(1, nbins):
        stat_acc[ii] = stat_acc[ii - 1] + stat[ii - 1]
    stat_acc = stat_acc / c0
    # compute rdf
    for ii in range(nbins):
        vol = 4.0 / 3.0 * np.pi * (((ii + 1) * hh) ** 3 - ((ii) * hh) ** 3)
        rho = stat[ii] / vol
        stat[ii] = rho / rho1 / c0
    xx = np.arange(0, max_r - 1e-12, hh)
    return xx, stat, stat_acc


if __name__ == "__main__":
    import dpdata

    sys = dpdata.System("out.lmp")
    xx, stat = rdf(sys, sel_type=[[0], None], max_r=8, nbins=100)
    res = np.concatenate([xx, stat]).reshape([2, -1])
    np.savetxt("rdf.out", res.T)
