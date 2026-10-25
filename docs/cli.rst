Command line interface
======================

.. argparse::
   :module: dpdata.cli
   :func: dpdata_parser
   :prog: dpdata

Offline workflow guides
-----------------------

The installed ``dpdata`` package includes the public workflow guides from
its ``skills/`` directory. Read the instructions matching your installed
version without installing a separate agent integration:

.. code-block:: console

   dpdata skills list
   dpdata skills read dpdata-cli
   dpdata skills read dpdata-cli README.md

The optional file argument must be a relative path inside the selected guide;
without one, ``skills read`` displays ``SKILL.md``. These commands work
offline and only print documentation: they do not execute workflow examples
or convert data. To convert a file literally named ``skills``, pass
``./skills`` or use ``--`` to separate the filename from CLI options.
