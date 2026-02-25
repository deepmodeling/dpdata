# dpdata CLI Agent Skill

An Agent Skill for using [dpdata](https://github.com/deepmodeling/dpdata) CLI to convert and manipulate atomic simulation data formats.

## Installation

To install this skill, provide the skill's GitHub URL to your agent (e.g., OpenClaw):

```text
Install the dpdata-cli skill from https://github.com/deepmodeling/dpdata/tree/devel/skills/dpdata-cli
```

Or manually copy to your agent's skills directory:

```bash
cp -r skills/dpdata-cli /path/to/your/agent/skills/
```

## Usage

Once installed, ask your AI agent to work with dpdata:

```text
Convert my VASP OUTCAR to DeePMD-kit format
```

```text
Convert LAMMPS dump file to VASP POSCAR
```

## Features

- **Format Conversion**: Convert between 50+ DFT/MD formats
- **DeePMD-kit Support**: Prepare training data for machine learning potentials
- **Auto-detection**: Automatic format detection for common file types
- **Multi-system Support**: Handle directories with multiple systems

## Supported Formats

Formats may be updated. See [Formats Reference (stable)](https://docs.deepmodeling.com/projects/dpdata/en/stable/formats.html) for the latest list.

### DeePMD-kit

`deepmd/raw`, `deepmd/comp`, `deepmd/npy`, `deepmd/npy/mixed`, `deepmd/hdf5`

### VASP

`vasp/poscar`, `vasp/contcar`, `vasp/outcar`, `vasp/xml`

### LAMMPS

`lammps/lmp`, `lammps/dump`

### ABACUS

`stru`, `abacus/scf`, `abacus/md`, `abacus/relax`

### And many more...

QE, CP2K, Gaussian, ORCA, PSI4, FHI-aims, SIESTA, PWmat, AMBER, GROMACS, ASE, pymatgen, XYZ, etc.

## Requirements

- [uv](https://docs.astral.sh/uv/) for running dpdata via `uvx`

## References

- [dpdata Documentation](https://docs.deepmodeling.com/projects/dpdata/)
- [CLI Reference](https://docs.deepmodeling.com/projects/dpdata/en/stable/cli.html)
- [Formats Reference](https://docs.deepmodeling.com/projects/dpdata/en/stable/formats.html)
- [GitHub Repository](https://github.com/deepmodeling/dpdata)
