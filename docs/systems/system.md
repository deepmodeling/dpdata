# `System` and `LabeledSystem`

This section gives some examples on how dpdata works. Firstly one needs to import the module in a python 3.x compatible code.
```python
import dpdata
```
The typicall workflow of `dpdata` is

1. Load data from vasp or lammps or deepmd-kit data files.
2. Manipulate data
3. Dump data to in a desired format


### Load data
```python
d_poscar = dpdata.System("POSCAR", fmt="vasp/poscar")
```
or let dpdata infer the format (`vasp/poscar`) of the file from the file name extension
```python
d_poscar = dpdata.System("my.POSCAR")
```
The number of atoms, atom types, coordinates are loaded from the `POSCAR` and stored to a data {class}`System` called `d_poscar`.
A data {class}`System` (a concept used by [deepmd-kit](https://github.com/deepmodeling/deepmd-kit)) contains frames that has the same number of atoms of the same type. The order of the atoms should be consistent among the frames in one {class}`System`.
It is noted that `POSCAR` only contains one frame.
If the multiple frames stored in, for example, a `OUTCAR` is wanted,
```python
d_outcar = dpdata.LabeledSystem("OUTCAR")
```
The labels provided in the `OUTCAR`, i.e. energies, forces and virials (if any), are loaded by {class}`LabeledSystem`. It is noted that the forces of atoms are always assumed to exist. {class}`LabeledSystem` is a derived class of {class}`System`.

The {class}`System` or {class}`LabeledSystem` can be constructed from the following file formats with the `format key` in the table passed to argument `fmt`:



### Access data
These properties stored in {class}`System` and {class}`LabeledSystem` can be accessed by operator `[]` with the key of the property supplied, for example
```python
coords = d_outcar["coords"]
```
Available properties are (nframe: number of frames in the system, natoms: total number of atoms in the system)

| key		|  type		| dimension		| are labels	| description
| ---		| ---		| ---			| ---		| ---
| 'atom_names'	| list of str	| ntypes		| False		| The name of each atom type
| 'atom_numbs'	| list of int	| ntypes		| False		| The number of atoms of each atom type
| 'atom_types'	| np.ndarray	| natoms		| False		| Array assigning type to each atom
| 'cells'	| np.ndarray	| nframes x 3 x 3	| False		| The cell tensor of each frame
| 'coords'	| np.ndarray	| nframes x natoms x 3	| False		| The atom coordinates
| 'energies'	| np.ndarray	| nframes		| True		| The frame energies
| 'forces'	| np.ndarray	| nframes x natoms x 3	| True		| The atom forces
| 'virials'	| np.ndarray	| nframes x 3 x 3	| True		| The virial tensor of each frame


### Dump data
The data stored in {class}`System` or {class}`LabeledSystem` can be dumped in 'lammps/lmp' or 'vasp/poscar' format, for example:
```python
d_outcar.to("lammps/lmp", "conf.lmp", frame_idx=0)
```
The first frames of `d_outcar` will be dumped to 'conf.lmp'
```python
d_outcar.to("vasp/poscar", "POSCAR", frame_idx=-1)
```
The last frames of `d_outcar` will be dumped to 'POSCAR'.

The data stored in `LabeledSystem` can be dumped to deepmd-kit raw format, for example
```python
d_outcar.to("deepmd/raw", "dpmd_raw")
```
Or a simpler command:
```python
dpdata.LabeledSystem("OUTCAR").to("deepmd/raw", "dpmd_raw")
```
Frame selection can be implemented by
```python
dpdata.LabeledSystem("OUTCAR").sub_system([0, -1]).to("deepmd/raw", "dpmd_raw")
```
by which only the first and last frames are dumped to `dpmd_raw`.


### replicate
dpdata will create a super cell of the current atom configuration.
```python
dpdata.System("./POSCAR").replicate(
    (
        1,
        2,
        3,
    )
)
```
tuple(1,2,3) means don't copy atom configuration in x direction, make 2 copys in y direction, make 3 copys in z direction.


### perturb
By the following example, each frame of the original system (`dpdata.System('./POSCAR')`) is perturbed to generate three new frames. For each frame, the cell is perturbed by 5% and the atom positions are perturbed by 0.6 Angstrom. `atom_pert_style` indicates that the perturbation to the atom positions is subject to normal distribution. Other available options to `atom_pert_style` are`uniform` (uniform in a ball), and `const` (uniform on a sphere).
```python
perturbed_system = dpdata.System("./POSCAR").perturb(
    pert_num=3,
    cell_pert_fraction=0.05,
    atom_pert_distance=0.6,
    atom_pert_style="normal",
)
print(perturbed_system.data)
```

### replace
By the following example, Random 8 Hf atoms in the system will be replaced by Zr atoms with the atom postion unchanged.
```python
s = dpdata.System("tests/poscars/POSCAR.P42nmc", fmt="vasp/poscar")
s.replace("Hf", "Zr", 8)
s.to_vasp_poscar("POSCAR.P42nmc.replace")
```
