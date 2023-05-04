import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import dpdata
import dpdata.gaussian.gjf  # noqa: F401
import dpdata.md.msd  # noqa: F401
import dpdata.md.water  # noqa: F401
import dpdata.stat  # noqa: F401
import dpdata.system  # noqa: F401

__all__ = [
    "dpdata",
]
