Supported Formats
=================

dpdata uses *format aliases* such as ``vasp/outcar`` or ``deepmd/npy`` to
select a reader or writer. Several aliases may point to the same
implementation; each reference page lists all equivalent names, supported
data models, conversion parameters, and copyable examples.

Choosing a data model
---------------------

Use :class:`dpdata.System` for cells, coordinates, and atom types. Use
:class:`dpdata.LabeledSystem` when the input also contains energies, forces,
virials, or other labels. :class:`dpdata.MultiSystems` groups systems with
different formulas, while :class:`dpdata.BondOrderSystem` preserves molecular
bonds and formal charges.

The supported-conversions column below is authoritative: a format may be
read-only, write-only, or support only one data model.

Common conversion patterns
--------------------------

Load a file explicitly and convert it to another format:

.. code-block:: python

   import dpdata

   labeled = dpdata.LabeledSystem("OUTCAR", fmt="vasp/outcar")
   labeled.to("deepmd/npy", "training_data")

Pass format-specific parameters to the constructor or :meth:`~dpdata.System.to`:

.. code-block:: python

   system = dpdata.System(
       "dump.lammpstrj",
       fmt="lammps/dump",
       type_map=["O", "H"],
       begin=100,
       step=10,
       unwrap=True,
   )
   system.to("vasp/poscar", "POSCAR", frame_idx=-1)

For files or directories containing multiple formulas, use
:meth:`dpdata.MultiSystems.from_file`:

.. code-block:: python

   systems = dpdata.MultiSystems.from_file("dataset.xyz", fmt="extxyz")
   systems.to("deepmd/hdf5", "dataset.hdf5")

The command-line interface supports the same basic conversion flow:

.. code-block:: bash

   dpdata OUTCAR -i vasp/outcar -o deepmd/npy -O training_data
   dpdata POSCAR -i vasp/poscar -n -o lammps/lmp -O data.lmp

Format reference
----------------

.. csv-table:: Supported Formats
   :file: formats.csv
   :header-rows: 1


.. toctree::
   :maxdepth: 1
   :glob:

   formats/*
