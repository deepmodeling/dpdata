# dpdata

[![conda-forge](https://img.shields.io/conda/dn/conda-forge/dpdata?color=red&label=conda-forge&logo=conda-forge)](https://anaconda.org/conda-forge/dpdata)
[![pip install](https://img.shields.io/pypi/dm/dpdata?label=pip%20install&logo=pypi)](https://pypi.org/project/dpdata)
[![Documentation Status](https://readthedocs.org/projects/dpdata/badge/)](https://dpdata.readthedocs.io/)

**dpdata** is a Python package for manipulating atomistic data of software in computational science.

## Installation

dpdata only supports Python 3.8 and above. You can [setup a conda/pip environment](https://docs.deepmodeling.com/faq/conda.html), and then use one of the following methods to install dpdata:

- Install via pip: `pip install dpdata`
- Install via conda: `conda install -c conda-forge dpdata`
- Install from source code: `git clone https://github.com/deepmodeling/dpdata && pip install ./dpdata`

To test if the installation is successful, you may execute

```bash
dpdata --version
```

## Supported packages

`dpdata` is aimmed to support different kinds of atomistic packages:

- Atomistic machine learning packages, such as [DeePMD-kit](https://github.com/deepmodeling/deepmd-kit);
- Molecular dynamics packages, such as [LAMMPS](https://github.com/lammps/lammps) and [GROMACS](https://gitlab.com/gromacs/gromacs);
- Quantum chemistry packages, such as [VASP](https://www.vasp.at/), [Gaussian](https://gaussian.com), and [ABACUS](https://github.com/deepmodeling/abacus-develop);
- Atomistic visualization packages, such as [3Dmol.js](https://3dmol.csb.pitt.edu/).
- Other atomistic tools, such as [ASE](https://gitlab.com/ase/ase).
- Common formats such as `xyz`.

All supported formats are listed [here](https://docs.deepmodeling.com/projects/dpdata/en/master/formats.html).

## Quick start

The quickest way to convert a simple file from one format to another one is to use the [command line](https://docs.deepmodeling.com/projects/dpdata/en/master/cli.html).

```sh
dpdata OUTCAR -i vasp/outcar -o deepmd/npy -O deepmd_data
```

For advanced usage with Python APIs, [read dpdata documentation](https://docs.deepmodeling.com/projects/dpdata/).

## Plugins

- [cp2kdata](https://github.com/robinzyb/cp2kdata) adds the latest CP2K support for dpdata.

For how to create your own plugin packages, [read dpdata documentation](https://docs.deepmodeling.com/projects/dpdata/).
