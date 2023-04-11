"""Amber mask."""
try:
    import parmed
except ImportError:
    pass


def pick_by_amber_mask(param, maskstr, coords=None):
    """Pick atoms by amber masks.

    Parameters
    ----------
    param : str or parmed.Structure
        filename of Amber param file or parmed.Structure
    maskstr : str
        Amber masks
    coords : np.ndarray (optional)
        frame coordinates, shape: N*3
    """
    parm = load_param_file(param)
    if coords is not None:
        parm.initialize_topology(xyz=coords)
    sele = []
    if len(maskstr) > 0:
        newmaskstr = maskstr.replace("@0", "!@*")
        sele = [
            parm.atoms[i].idx
            for i in parmed.amber.mask.AmberMask(parm, newmaskstr).Selected()
        ]
    return sele


def load_param_file(param_file):
    if isinstance(param_file, str):
        return parmed.load_file(param_file)
    elif isinstance(param_file, parmed.Structure):
        return param_file
    else:
        raise RuntimeError("Unsupported structure")
