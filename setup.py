# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='dpfu',
    version='0.0.0',
    description='Manipulating DeePMD-kit data, VASP, LAMMPS file formats',
    long_description=readme,
    author='Deep modeling developers',
    author_email='wang_han@iapcm.ac.cn',
    url='https://github.com/deepmodeling/dpfu',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)
