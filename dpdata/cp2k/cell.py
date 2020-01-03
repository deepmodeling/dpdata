
#%%
import numpy as np
from collections import OrderedDict
import re

def cell_to_low_triangle(A,B,C,alpha,beta,gamma):
    """
        Convert cell to low triangle matrix.

        Parameters
        ----------
        A : float
            cell length A
        B : float
            cell length B
        C : float
            cell length C
        alpha : float
            radian. The angle between vector B and  vector C.
        beta : float
            radian. The angle between vector A and  vector C.
        gamma : float
            radian. The angle between vector B and  vector C.
        
        Returns
        -------
        cell : list
            The cell matrix used by dpdata in low triangle form.
    """
    if not np.pi*5/180<alpha< np.pi*175/180:
        raise RuntimeError("alpha=={}: must be a radian, and \
            must be in np.pi*5/180 < alpha < np.pi*175/180".format(alpha))
    if not np.pi*5/180<beta< np.pi*175/180:
        raise RuntimeError("beta=={}: must be a radian, and \
            must be in np.pi*5/180 < beta < np.pi*175/180".format(beta))
    if not np.pi*5/180<gamma< np.pi*175/180:
            raise RuntimeError("gamma=={}: must be a radian, and \
                must be in np.pi*5/180 < gamma < np.pi*175/180".format(gamma))
    if not A > 0.2:
        raise RuntimeError("A=={}, must be greater than 0.2".format(A))
    if not B > 0.2:
        raise RuntimeError("B=={}, must be greater than 0.2".format(B))
    if not C > 0.2:
        raise RuntimeError("C=={}, must be greater than 0.2".format(C))

    lx = A
    xy = B * np.cos(gamma)
    xz = C * np.cos(beta)
    ly = B* np.sin(gamma)
    if not ly > 0.1:
        raise RuntimeError("ly:=B* np.sin(gamma)=={}, must be greater than 0.1",format(ly))
    yz = (B*C*np.cos(alpha)-xy*xz)/ly
    if not C**2-xz**2-yz**2 > 0.01:
        raise RuntimeError("lz^2:=C**2-xz**2-yz**2=={}, must be greater than 0.01",format(C**2-xz**2-yz**2))
    lz = np.sqrt(C**2-xz**2-yz**2)
    cell = np.asarray([[lx, 0 , 0],
            [xy, ly, 0 ],
            [xz, yz, lz]]).astype('float32')
    return cell

