"""Quantum ESPRESSO ``pw.x`` molecular-dynamics output reader.

The implementation is shared with the ``qe/pw/scf`` reader because a PWscf
output can contain either one labeled structure or a sequence of MD frames.
"""

from __future__ import annotations

from .scf import get_frames


def to_system_data(fname, begin=0, step=1):
    """Read PWscf AIMD frames from an output file and its matching input."""
    return get_frames(fname, begin=begin, step=step)
