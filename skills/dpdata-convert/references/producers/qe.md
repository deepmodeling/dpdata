# Quantum ESPRESSO producer contract

Use this reference when converting Quantum ESPRESSO PWscf or CP trajectory
output (`qe/pw/scf` or `qe/cp/traj`). Record the QE release, input file, output
file, and any trajectory sidecar files used by the parser.

Before accepting the converted data, verify:

- the input atom order is the order used in every force and trajectory frame;
- the cell and periodic boundary conditions are retained;
- energies and forces have the same frame count; and
- stress/virial values, units, and sign conventions are checked rather than
  inferred from a label name.

QE output may omit labels when a run stops early. Preserve the raw files and
fail the conversion validation when a required label is missing.
