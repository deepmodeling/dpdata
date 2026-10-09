# Semantic validation contract

A parser can succeed while producing a dataset that is unusable for its next
scientific step. Validate the following contract after conversion. Checks that
can be made mechanically should be delegated to
`dpdata-inspect/scripts/validate_labeled_system.py`; conventions that are not
encoded in dpdata must be declared by the producer contract.

| Invariant       | Check                                                                                                                                                                                                         |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Labels          | Decide whether energies are required. If forces or virials are present, they must be complete for the same frames. Missing labels are reported, not filled with zeros or silently dropped.                    |
| Frame alignment | `coords` and `cells` have `nframes`; required labels cover every frame; optional labels, when present, have matching frame dimensions. Forces have `(nframes, natoms, 3)` and virials have `(nframes, 3, 3)`. |
| Atoms and types | `atom_types` has `natoms`, indices are in `atom_names`, and `sum(atom_numbs) == natoms`. A multi-system dataset uses one documented type-map order.                                                           |
| Values          | Coordinates, cells, energies, forces, and virials contain finite values.                                                                                                                                      |
| Cell and PBC    | Decide whether the source is periodic. For periodic data, cells are present and nonsingular; for molecules, preserve the explicit non-periodic (`nopbc`) choice.                                              |
| Units           | Record energy and force units. Record the virial/stress unit and the sign convention separately; they are not inferred from array names.                                                                      |
| Provenance      | Preserve the raw input, producer and version, parser version, conversion command, selected format, type map, and any filtering or dropped labels.                                                             |

`dpdata`'s standard labeled fields are `energies`, `forces`, and `virials`.
Some producers call the last field stress or pressure. Confirm whether the
producer reports stress or virial and whether its sign is the negative of the
convention expected by the consumer.

Producer references describe where to look for a contract. Project-specific
requirements such as a required print level, convergence threshold, or
scheduler setup belong in the project record, not in this reusable skill.

## Deterministic validator

For one labeled input, run:

```bash
python skills/dpdata-inspect/scripts/validate_labeled_system.py \
    converted -f deepmd/raw \
    --energy-unit eV --force-unit eV/A \
    --producer VASP --producer-version 6.4 --parser-version 0.1 \
    --provenance raw/OUTCAR
```

Use `--multi` for a directory of systems, `--type-map C H O` to assert atom
ordering, and `--require-pbc` when periodic cells are required. Add
`--strict-metadata` when undeclared units, stress sign, or provenance should
fail the check. `--json` emits a machine-readable result for CI.
