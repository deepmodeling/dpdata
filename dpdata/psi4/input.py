from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

# Angston is used in Psi4 by default
template = """molecule {{
{atoms:s}
{charge:d} {multiplicity:d}
}}
set basis {basis:s}
set gradient_write on
G, wfn = gradient("WB97M-D3BJ", return_wfn=True)
wfn.energy()
wfn.gradient().print_out()
"""


def write_psi4_input(
    types: np.ndarray,
    coords: np.ndarray,
    method: str,
    basis: str,
    charge: int = 0,
    multiplicity: int = 1,
) -> str:
    """Write Psi4 input file.

    Parameters
    ----------
    types : np.ndarray
        atomic symbols
    coords : np.ndarray
        atomic coordinates
    method : str
        computational method
    basis : str
        basis set; see https://psicode.org/psi4manual/master/basissets_tables.html
    charge : int, default=0
        charge of system
    multiplicity : int, default=1
        multiplicity of system

    Returns
    -------
    str
        content of Psi4 input file
    """
    return template.format(
        atoms="\n".join(
            [
                "{:s} {:16.9f} {:16.9f} {:16.9f}".format(*ii)
                for ii in zip(types, *coords.T)
            ]
        ),
        charge=charge,
        multiplicity=multiplicity,
        method=method,
        basis=basis,
    )
