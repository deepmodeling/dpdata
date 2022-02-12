import numpy as np

def coord_to_xyz(coord: np.ndarray, types: list)->str:
    """Convert coordinates and types to xyz format.
    
    Parameters
    ----------
    coord: np.ndarray
        coordinates, Nx3 array
    types: list
        list of types
    
    Returns
    -------
    str
        xyz format string
    
    Examples
    --------
    >>> coord_to_xyz(np.ones((1,3)), ["C"])
    1

    C 1.000000 1.000000 1.000000
    """
    buff = [str(len(types)), '']
    for at, cc in zip(types, coord):
        buff.append("{} {:.6f} {:.6f} {:.6f}".format(at, *cc))
    return "\n".join(buff)
