---
name: dpdata-cli
description: Convert and manipulate atomic simulation data formats using dpdata CLI. Use when converting between DFT/MD output formats (VASP, LAMMPS, QE, CP2K, Gaussian, ABACUS, etc.), preparing training data for DeePMD-kit, or working with DeePMD formats. Supports 50+ formats including deepmd/raw, deepmd/comp, deepmd/npy, deepmd/hdf5.
compatibility: Requires uvx (uv) for running dpdata
metadata:
  author: njzjz-bot
  version: '1.0'
  repository: https://github.com/deepmodeling/dpdata
---

# dpdata CLI

dpdata is a tool for manipulating multiple atomic simulation data formats. This skill enables format conversion between various DFT/MD software outputs via command line.

## Quick Start

Run dpdata via uvx:

```bash
uvx dpdata <from_file> [options]
```

## Command Line Usage

```text
dpdata: Manipulating multiple atomic simulation data formats
usage: dpdata [-h] [--to_file TO_FILE] [--from_format FROM_FORMAT]
              [--to_format TO_FORMAT] [--no-labeled] [--multi]
              [--type-map TYPE_MAP [TYPE_MAP ...]] [--version]
              from_file
```

### Arguments

| Argument              | Description                                           |
| --------------------- | ----------------------------------------------------- |
| `from_file`           | Read data from a file (positional)                    |
| `--to_file`, `-O`     | Dump data to a file                                   |
| `--from_format`, `-i` | Format of from_file (default: "auto")                 |
| `--to_format`, `-o`   | Format of to_file                                     |
| `--no-labeled`, `-n`  | Labels aren't provided (default: False)               |
| `--multi`, `-m`       | System contains multiple directories (default: False) |
| `--type-map`, `-t`    | Type map for atom types                               |
| `--version`           | Show dpdata version and exit                          |

## Common Examples

### Convert VASP OUTCAR to deepmd format

```bash
uvx dpdata OUTCAR -i vasp/outcar -O deepmd_data -o deepmd/raw
```

### Convert LAMMPS dump to VASP POSCAR

```bash
uvx dpdata dump.lammps -i lammps/dump -O POSCAR -o vasp/poscar
```

### Convert with type map

```bash
uvx dpdata OUTCAR -i vasp/outcar -O deepmd_data -o deepmd/raw -t C H O N
```

### Convert multiple systems

```bash
uvx dpdata data_dir -i vasp/outcar -O output_dir -o deepmd/comp --multi
```

### Convert to deepmd/npy (compressed format)

```bash
uvx dpdata OUTCAR -i vasp/outcar -O deepmd_npy -o deepmd/npy
```

### Convert to deepmd/hdf5

```bash
uvx dpdata OUTCAR -i vasp/outcar -O data.h5 -o deepmd/hdf5
```

## Supported Formats

Formats may be updated. For the complete and latest list, see:

- [Formats Reference (stable)](https://docs.deepmodeling.com/projects/dpdata/en/stable/formats.html)

### DeePMD-kit Formats

| Format Name                  | Description                        |
| ---------------------------- | ---------------------------------- |
| `deepmd/raw`                 | DeePMD-kit raw text format         |
| `deepmd/comp` / `deepmd/npy` | DeePMD-kit compressed numpy format |
| `deepmd/npy/mixed`           | DeePMD-kit mixed type format       |
| `deepmd/hdf5`                | DeePMD-kit HDF5 format             |

### VASP Formats

| Format Name                                           | Description          |
| ----------------------------------------------------- | -------------------- |
| `vasp/poscar` / `vasp/contcar` / `poscar` / `contcar` | VASP structure files |
| `vasp/outcar` / `outcar`                              | VASP OUTCAR output   |
| `vasp/xml` / `xml`                                    | VASP XML output      |
| `vasp/string`                                         | VASP string format   |

### LAMMPS Formats

| Format Name            | Description      |
| ---------------------- | ---------------- |
| `lammps/lmp` / `lmp`   | LAMMPS data file |
| `lammps/dump` / `dump` | LAMMPS dump file |

### ABACUS Formats

| Format Name                                              | Description           |
| -------------------------------------------------------- | --------------------- |
| `stru` / `abacus/stru`                                   | ABACUS structure file |
| `abacus/lcao/scf` / `abacus/pw/scf` / `abacus/scf`       | ABACUS SCF output     |
| `abacus/lcao/md` / `abacus/pw/md` / `abacus/md`          | ABACUS MD output      |
| `abacus/lcao/relax` / `abacus/pw/relax` / `abacus/relax` | ABACUS relax output   |

### Quantum ESPRESSO Formats

| Format Name  | Description      |
| ------------ | ---------------- |
| `qe/cp/traj` | QE CP trajectory |
| `qe/pw/scf`  | QE PWscf output  |

### CP2K Formats

| Format Name        | Description      |
| ------------------ | ---------------- |
| `cp2k/output`      | CP2K output      |
| `cp2k/aimd_output` | CP2K AIMD output |

### Gaussian Formats

| Format Name     | Description                   |
| --------------- | ----------------------------- |
| `gaussian/log`  | Gaussian log file             |
| `gaussian/fchk` | Gaussian formatted checkpoint |
| `gaussian/md`   | Gaussian MD output            |
| `gaussian/gjf`  | Gaussian input file           |

### Other Formats

| Format Name                                                         | Description           |
| ------------------------------------------------------------------- | --------------------- |
| `xyz`                                                               | XYZ format            |
| `mace/xyz` / `nequip/xyz` / `gpumd/xyz` / `extxyz` / `quip/gap/xyz` | Extended XYZ variants |
| `ase/structure`                                                     | ASE structure format  |
| `ase/traj`                                                          | ASE trajectory        |
| `pymatgen/structure`                                                | pymatgen structure    |
| `pymatgen/molecule`                                                 | pymatgen molecule     |
| `gromacs/gro` / `gro`                                               | GROMACS gro file      |
| `siesta/output`                                                     | SIESTA output         |
| `siesta/aimd_output`                                                | SIESTA AIMD output    |
| `pwmat/output` / `pwmat/mlmd` / `pwmat/movement`                    | PWmat output          |
| `pwmat/final.config` / `pwmat/atom.config`                          | PWmat config          |
| `orca/spout`                                                        | ORCA output           |
| `psi4/out`                                                          | PSI4 output           |
| `dftbplus`                                                          | DFTB+ output          |
| `fhi_aims/output` / `fhi_aims/md`                                   | FHI-aims output       |
| `amber/md`                                                          | AMBER MD              |
| `n2p2`                                                              | n2p2 format           |
| `mol_file` / `mol`                                                  | MOL file              |
| `sdf_file` / `sdf`                                                  | SDF file              |
| `openmx/md`                                                         | OpenMX MD             |
| `sqm/out`                                                           | SQM output            |
| `sqm/in`                                                            | SQM input             |
| `list`                                                              | List format           |
| `3dmol`                                                             | 3Dmol visualization   |

## Tips

1. **Auto-detection**: Use `-i auto` (default) to let dpdata detect format automatically
1. **Type mapping**: Use `-t` to specify atom type order for deepmd formats
1. **Multi-system**: Use `--multi` for directories containing multiple systems
1. **Compressed output**: Use `deepmd/npy` or `deepmd/hdf5` for smaller file sizes

## References

- [dpdata Documentation](https://docs.deepmodeling.com/projects/dpdata/)
- [CLI Reference](https://docs.deepmodeling.com/projects/dpdata/en/stable/cli.html)
- [Formats Reference](https://docs.deepmodeling.com/projects/dpdata/en/stable/formats.html)
- [GitHub Repository](https://github.com/deepmodeling/dpdata)
