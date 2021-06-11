from . import vasp
from . import lammps
from . import md
from .system import System
from .system import LabeledSystem
from .system import MultiSystems

try:
    from ._version import version as __version__
except ImportError:
    from .__about__ import __version__

# BondOrder System has dependency on rdkit
try:
    import rdkit
    USE_RDKIT = True
except ModuleNotFoundError:
    USE_RDKIT = False

if USE_RDKIT:
    from .bond_order_system import BondOrderSystem


