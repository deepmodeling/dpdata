---
name: dpdata-inspect
description: Inspect and validate dpdata systems against labeled-data, shape, atom-order, cell/PBC, finite-value, units, and provenance contracts. Use before training, merging, filtering, or publishing converted data.
compatibility: Requires dpdata for loading inputs; the deterministic validator requires NumPy.
metadata:
  repository: https://github.com/deepmodeling/dpdata
---

# Inspect atomistic data

Use this workflow when the user asks whether a dataset is valid, complete, or
ready for a downstream model. Start by identifying whether the input is one
system or a collection and whether labels are required. Then load the
producer-specific contract only if the source producer is known.

## Procedure

1. Preserve the raw source and record its producer/version, parser version,
   format, and conversion command.
2. Load one `LabeledSystem` when energies are required. Load `System` only for
   an intentionally unlabeled structure. Use `MultiSystems` for independent
   systems and check that their atom/type order is compatible.
3. Run the deterministic validator below. It checks array shapes, frame
   alignment, atom counts and type indices, finite values, labels, and (when
   requested) periodic-cell determinants.
4. Review the remaining contract manually: units, stress versus virial and
   sign, periodicity, atom ordering, missing-label policy, producer/version,
   and provenance. These cannot be inferred safely from NumPy arrays.
5. Stop with an actionable error when a required label or contract declaration
   is absent. Do not repair missing labels by inserting zeros.

## Deterministic validator

For a labeled system:

```bash
python skills/dpdata-inspect/scripts/validate_labeled_system.py \
  converted -f deepmd/raw \
  --energy-unit eV --force-unit eV/A \
  --producer VASP --producer-version 6.4 --parser-version 0.1 \
  --provenance raw/source.out
```

For multiple systems or an explicitly unlabeled structure:

```bash
python skills/dpdata-inspect/scripts/validate_labeled_system.py \
  dataset --format deepmd/comp --multi --type-map C H O
python skills/dpdata-inspect/scripts/validate_labeled_system.py \
  POSCAR --format vasp/poscar --allow-unlabeled --require-pbc
```

Add `--strict-metadata` to make undeclared units, stress sign, producer/version,
or provenance fail. Add `--json` for CI or another machine consumer. The script exits zero
only when every loaded system satisfies the selected contract.

The validator is intentionally small and deterministic. It does not infer
units, transform forces, choose a stress sign, or decide whether a project
should accept partial labels. Those are producer and project contracts and
must be recorded by the caller.
