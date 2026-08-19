.. image:: _static/logo.svg
   :alt: dpdata logo
   :width: 600px
   :align: center

======
dpdata
======

.. rst-class:: lead

   Turn atomistic simulation outputs into interoperable, machine-learning-ready
   datasets.

.. important::

   **One data model, many atomistic formats.** Load structures, trajectories,
   energies, forces, and virials from simulation codes, manipulate them through
   :class:`dpdata.System`, :class:`dpdata.LabeledSystem`, and
   :class:`dpdata.MultiSystems`, then export the data in the format your next
   tool expects.

dpdata connects electronic-structure codes, molecular-dynamics engines, and
atomistic machine-learning workflows. Use the command line for one-off
conversion, or the Python API to build reproducible data-processing pipelines
that preserve atomistic structures and labels.

Choose your path
================

.. grid:: 1 2 3 3
   :gutter: 3

   .. grid-item-card:: 🔄 Convert a file
      :link: cli
      :link-type: doc
      :shadow: md

      Move structures or labeled trajectories between registered formats from
      the command line.

   .. grid-item-card:: 🏷️ Build labeled datasets
      :link: systems/system
      :link-type: doc
      :shadow: md

      Load coordinates, cells, energies, forces, virials, and other registered
      fields into a common Python data model.

   .. grid-item-card:: 🧩 Combine many systems
      :link: systems/multi
      :link-type: doc
      :shadow: md

      Organize data with different compositions or atom counts using
      ``MultiSystems``.

   .. grid-item-card:: 💾 Scale with LMDB
      :link: systems/lmdb
      :link-type: doc
      :shadow: md

      Store frames from one or more systems in a single DeePMD-compatible LMDB
      database.

   .. grid-item-card:: 🗂️ Browse formats
      :link: formats
      :link-type: doc
      :shadow: md

      Find the registered readers and writers for electronic structure,
      molecular dynamics, atomistic ML, chemistry, and analysis tools.

   .. grid-item-card:: 🚀 Try dpdata online
      :link: try_dpdata
      :link-type: doc
      :shadow: md

      Open an interactive browser session without installing dpdata locally.

Why dpdata
==========

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: One data model, many formats
      :link: formats
      :link-type: doc
      :shadow: sm

      Read and write atomistic formats through the same ``System`` and
      ``LabeledSystem`` interfaces instead of maintaining one converter per
      code pair.

   .. grid-item-card:: ML-ready labels
      :link: systems/system
      :link-type: doc
      :shadow: sm

      Keep structures together with energies, forces, virials, and registered
      extra fields while converting or selecting frames.

   .. grid-item-card:: Heterogeneous datasets
      :link: systems/multi
      :link-type: doc
      :shadow: sm

      Work with multiple compositions and atom counts as a collection rather
      than forcing every frame into one homogeneous system.

   .. grid-item-card:: Dataset-scale storage
      :link: systems/lmdb
      :link-type: doc
      :shadow: sm

      Export DeePMD NumPy data or use LMDB when a single database is a better
      fit for heterogeneous atomistic data.

   .. grid-item-card:: Structure operations
      :link: systems/system
      :link-type: doc
      :shadow: sm

      Select frames, build supercells, perturb structures, replace species,
      and compose transformations in Python.

   .. grid-item-card:: Open plugin ecosystem
      :link: plugin
      :link-type: doc
      :shadow: sm

      Add new formats as installable Python packages through the
      ``dpdata.plugins`` entry point.

From calculation output to training data
========================================

A typical dpdata workflow has four steps:

1. **Load** structures or labeled trajectories from a supported simulation
   format.
2. **Manipulate** frames, structures, species, or dataset partitions in Python.
3. **Organize** homogeneous data in ``System`` / ``LabeledSystem`` or
   heterogeneous data in ``MultiSystems``.
4. **Export** to the format expected by the next simulation, analysis, or
   machine-learning tool.

The same registry powers both the command-line interface and Python API, so a
one-line conversion can grow into a larger data pipeline without changing the
underlying format model.

Start in minutes
================

dpdata requires Python 3.10 or later. The fastest installation path is:

.. code-block:: bash

   python -m pip install dpdata
   dpdata --version

Conda-forge and source installation are covered in the
:doc:`installation guide <installation>`.

Convert from the command line
-----------------------------

Convert a VASP ``OUTCAR`` directly to a DeePMD NumPy dataset:

.. code-block:: bash

   dpdata OUTCAR -i vasp/outcar -o deepmd/npy -O deepmd_data

See the full :doc:`command-line reference <cli>` for input and output options.

Build a labeled dataset in Python
---------------------------------

.. code-block:: python

   import dpdata

   # OUTCAR is recognized as a labeled VASP trajectory.
   data = dpdata.LabeledSystem("OUTCAR")

   # Keep selected frames and write a DeePMD NumPy dataset.
   data.sub_system([0, -1]).to("deepmd/npy", "deepmd_data")

:class:`dpdata.LabeledSystem` keeps structures together with energies, forces,
and virials when they are available. Continue with :doc:`System and
LabeledSystem <systems/system>` for data access, frame selection, replication,
perturbation, and species replacement.

Combine many systems and scale out
----------------------------------

.. code-block:: python

   import dpdata

   systems = dpdata.MultiSystems.from_dir(
       "./calculations",
       file_name="OUTCAR",
       fmt="vasp/outcar",
   )
   systems.to("deepmd/lmdb", "training.lmdb")

``MultiSystems`` groups heterogeneous structures by composition, while
``deepmd/lmdb`` stores frames from one or more systems in a single database.
See :doc:`MultiSystems <systems/multi>` and :doc:`LMDB datasets <systems/lmdb>`.

Core data model
===============

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Object
     - Use it for
   * - :class:`dpdata.System`
     - Structures and trajectories: atom types, coordinates, cells, and other
       non-label fields.
   * - :class:`dpdata.LabeledSystem`
     - Reference data for atomistic ML: a ``System`` plus energies, forces,
       virials, and other registered labels.
   * - :class:`dpdata.MultiSystems`
     - Collections containing multiple systems, compositions, or atom counts.

Specialized representations, including bond-order and mixed-type systems, are
available in the :doc:`Systems guide <systems/index>`.

Scientific ecosystem
====================

dpdata is designed to sit between the tools already used in computational
chemistry and materials science.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Area
     - Representative built-in integrations
   * - Electronic structure and quantum chemistry
     - VASP, ABACUS, Quantum ESPRESSO, Gaussian, CP2K, ORCA, FHI-aims,
       SIESTA, OpenMX, and DFTB+
   * - Molecular dynamics
     - LAMMPS and GROMACS
   * - Atomistic ML and data
     - DeePMD-kit formats, LMDB datasets, ASE, and pymatgen-compatible
       structures
   * - Chemistry and visualization
     - RDKit and 3Dmol.js

The :doc:`supported-formats table <formats>` is generated from dpdata's format
registry and is the source of truth for available readers and writers.

Documentation map
=================

* **Get started:** :doc:`Installation <installation>` · :doc:`Try dpdata online
  <try_dpdata>` · :doc:`Command line <cli>`.
* **Work with data:** :doc:`Systems <systems/index>` · :doc:`Supported formats
  <formats>` · :doc:`Drivers <drivers>` · :doc:`Minimizers <minimizers>`.
* **Extend dpdata:** :doc:`Plugin guide <plugin>` · :doc:`Python API <api/api>`.
* **Project:** :doc:`Credits <credits>` · `GitHub
  <https://github.com/deepmodeling/dpdata>`_.

Citation
========

If dpdata contributes to published work, please cite:

Jinzhe Zeng, Xingliang Peng, Yong-Bin Zhuang, Haidi Wang, Fengbo Yuan, Duo
Zhang, Renxi Liu, Yingze Wang, Ping Tuo, Yuzhi Zhang, Yixiao Chen, Yifan Li,
Cao Thang Nguyen, Jiameng Huang, Anyang Peng, Marián Rynik, Wei-Hong Xu,
Zezhong Zhang, Xu-Yuan Zhou, Tao Chen, Jiahao Fan, Wanrun Jiang, Bowen Li,
Denan Li, Haoxi Li, Wenshuo Liang, Ruihao Liao, Liping Liu, Chenxing Luo,
Logan Ward, Kaiwei Wan, Junjie Wang, Pan Xiang, Chengqian Zhang, Jinchao
Zhang, Rui Zhou, Jia-Xin Zhu, Linfeng Zhang, and Han Wang. “dpdata: A Scalable
Python Toolkit for Atomistic Machine Learning Data Sets.” *Journal of Chemical
Information and Modeling* **65** (21), 11497–11504 (2025). `DOI:
10.1021/acs.jcim.5c01767 <https://doi.org/10.1021/acs.jcim.5c01767>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents
   :hidden:

   installation
   systems/index
   try_dpdata
   cli
   formats
   drivers
   minimizers
   plugin
   api/api
   credits

Indices
=======

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
