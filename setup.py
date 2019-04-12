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
install_requires=[]

setuptools.setup(
    name="dpdata",
    version_format='{tag}.dev{commitcount}+{gitsha}',
    setup_requires=['setuptools-git-version'],
    # version="0.0.2",
    author="Han Wang",
    author_email="wang_han@iapcm.ac.cn",
    description="Manipulating DeePMD-kit, VASP and LAMMPS data formats",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/deepmodeling/dpdata",
    packages=['dpdata', 'dpdata/vasp', 'dpdata/lammps', 'dpdata/md', 'dpdata/deepmd'],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
    ],
    keywords='lammps vasp deepmd-kit',
    install_requires=install_requires,    
)

