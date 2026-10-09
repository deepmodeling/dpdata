---
name: dpdata-convert
description: Convert atomistic producer outputs with dpdata while preserving and validating the labeled-data contract. Use for format conversion, dataset preparation, and exports from VASP, CP2K, ABACUS, Quantum ESPRESSO, LAMMPS, or other supported producers.
compatibility: Requires dpdata and the dependencies of the selected format.
metadata:
  repository: https://github.com/deepmodeling/dpdata
---

# Convert atomistic data

Use this workflow when a user starts with a producer output and wants a
`dpdata` system, a training dataset, or another file format. The command is
only one part of the task: the resulting arrays must retain their scientific
meaning.

## Procedure

1. **Identify the producer contract.** Record the producer and version, input
   and output paths, parser/`dpdata` version, and whether the source contains
   energies, forces, and virials/stress. Load the matching reference only when
   the producer is known:
   [`CP2K`](references/producers/cp2k.md),
   [`VASP`](references/producers/vasp.md),
   [`ABACUS`](references/producers/abacus.md), or
   [`Quantum ESPRESSO`](references/producers/qe.md).
2. **Choose the data shape.** Use `LabeledSystem` for data that must carry
   energies, and `System` for structures without labels. Use `MultiSystems`
   or `--multi` only when directories contain independent systems. Do not
   silently turn partial labels into a complete dataset; state which labels
   are required and what happens when one is absent.
3. **Convert without discarding the source.** Keep the raw files and write the
   converted output to a new path. Include an explicit `--from_format` when
   auto-detection could be ambiguous. Supply `--type-map` when the producer's
   atom ordering is not self-describing.
4. **Validate the output.** Check frame alignment, atom ordering, type-map
   consistency, finite values, cell/PBC assumptions, units, and virial/stress
   conventions. The reusable checks and their limits are in
   [`semantic-validation.md`](references/semantic-validation.md); the
   deterministic validator is part of [`dpdata-inspect`](../dpdata-inspect/SKILL.md).
5. **Record provenance.** Keep the raw path, producer/version, parser version,
   command, format, type map, units, and any dropped or missing labels next to
   the converted data. A successful parser run is not evidence that a dataset
   is suitable for training.

## Command reference

The CLI uses the same format names as the Python API:

```bash
uvx dpdata INPUT -i PRODUCER_FORMAT -O OUTPUT -o TARGET_FORMAT
```

Examples:

```bash
# A labeled VASP trajectory to DeePMD raw data
uvx dpdata OUTCAR -i vasp/outcar -O deepmd_data -o deepmd/raw

# A directory containing independent systems
uvx dpdata data_dir -i vasp/outcar -O output_dir -o deepmd/comp --multi

# Make the atom order explicit when reading a LAMMPS dump
uvx dpdata dump.lammps -i lammps/dump -O POSCAR -o vasp/poscar -t C H O
```

For all options and format aliases, load
[`cli-reference.md`](references/cli-reference.md). The reference is kept
separate so a conversion request does not require loading every format at
once.
