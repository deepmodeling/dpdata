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
