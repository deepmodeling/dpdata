.. dpdata documentation master file, created by
   sphinx-quickstart on Fri Apr 12 14:24:37 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dpdata's documentation!
==================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


.. mdinclude:: ../README.md


API documentation
=================

.. automodule:: dpdata

.. autoclass:: System
    :members: __init__, __getitem__, get_nframes, get_natoms, sub_system, append, apply_pbc, to_lammps_lmp, to_vasp_poscar

.. autoclass:: LabeledSystem
    :members: __init__, sub_system, to_deepmd_raw, to_deepmd_npy


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
