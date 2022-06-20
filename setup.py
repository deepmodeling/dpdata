# -*- coding: utf-8 -*-

from os import path
import setuptools

readme_file = path.join(path.dirname(path.abspath(__file__)), 'README.md')
try:
    from m2r import parse_from_file
    readme = parse_from_file(readme_file)
except ImportError:
    with open(readme_file) as f:
        readme = f.read()

# install_requires = ['xml']
install_requires=['numpy>=1.14.3', 'monty', 'scipy', 'h5py', 'wcmatch', 'importlib_metadata>=1.4; python_version < "3.8"']

setuptools.setup(
    name="dpdata",
    use_scm_version={'write_to': 'dpdata/_version.py'},
    setup_requires=['setuptools_scm'],
    author="Han Wang",
    author_email="wang_han@iapcm.ac.cn",
    description="Manipulating data formats of DeePMD-kit, VASP, QE, PWmat, and LAMMPS, etc.",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/deepmodeling/dpdata",
    packages=['dpdata',
              'dpdata/vasp', 
              'dpdata/lammps', 
              'dpdata/md', 
              'dpdata/deepmd', 
              'dpdata/qe', 
              'dpdata/siesta', 
              'dpdata/gaussian', 
              'dpdata/cp2k',
              'dpdata/xyz',
              'dpdata/pwmat', 
              'dpdata/amber',
              'dpdata/fhi_aims',
              'dpdata/gromacs',
              'dpdata/abacus',
              'dpdata/rdkit',
              'dpdata/plugins',
              'dpdata/pymatgen',
    ],
    package_data={'dpdata':['*.json']},
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
    ],
    keywords='lammps vasp deepmd-kit',
    install_requires=install_requires,
    extras_require={
        'ase': ['ase'],
        'amber': ['parmed'],
        'pymatgen': ['pymatgen'],
        'docs': ['sphinx', 'recommonmark', 'sphinx_rtd_theme>=1.0.0rc1', 'numpydoc', 'm2r2', 'deepmodeling-sphinx', 'sphinx-argparse'],
    },
    entry_points={"console_scripts": ["dpdata = dpdata.cli:dpdata_cli"]},
)

