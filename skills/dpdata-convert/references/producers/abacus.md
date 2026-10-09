# ABACUS producer contract

Use this reference when converting ABACUS SCF, MD, relax, or `STRU` output.
Record the ABACUS version, calculation mode (PW or LCAO), and the input/output
paths used by the parser.

Before accepting the converted data, verify:

- the `STRU` atom type labels and atom order match `atom_names` and
  `atom_types`;
- SCF, MD, or relax frames have matching coordinates, cells, energies, and
  forces;
- cell-relax changes are preserved as per-frame cells; and
- any stress or force sign and unit conversion is explicit.

ABACUS output formats evolve. Check the producer release and parser contract
for the exact output layout; do not assume a successful parse proves that all
optional labels were present.
