# VASP producer contract

Use this reference when converting POSCAR/CONTCAR, OUTCAR, or vasprun XML.
Choose the corresponding dpdata format and record the VASP release and whether
the source is a structure-only file or a labeled calculation.

Before accepting the converted data, verify:

- atom names, counts, and the coordinate order agree with the selected
  `type_map`;
- selective dynamics, Cartesian/direct coordinates, and cell vectors were
  interpreted as intended;
- energies and forces are frame-aligned, and stress/virial conversion has
  been checked against the downstream convention; and
- an incomplete final ionic step was not retained as a labeled frame.

Use the VASP documentation and the exact output version for units and stress
semantics. Keep the original POSCAR/OUTCAR/vasprun XML with the converted
system.
