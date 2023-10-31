# monty needs lzma
# See https://github.com/pandas-dev/pandas/pull/27882
try:
    import lzma  # noqa: F401
except ImportError:

    class fakemodule:
        pass

    import sys

    sys.modules["lzma"] = fakemodule

from . import lammps, md, vasp
from .system import LabeledSystem, MultiSystems, System

try:
    from ._version import version as __version__
except ImportError:
    from .__about__ import __version__

# BondOrder System has dependency on rdkit
try:
    # prevent conflict with dpdata.rdkit
    import rdkit as _  # noqa: F401

    USE_RDKIT = True
except ModuleNotFoundError:
    USE_RDKIT = False

if USE_RDKIT:
    from .bond_order_system import BondOrderSystem

__all__ = [
    "__version__",
    "lammps",
    "md",
    "vasp",
    "System",
    "LabeledSystem",
    "MultiSystems",
    "BondOrderSystem",
]
