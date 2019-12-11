**dpdata** is a python package for manipulating DeePMD-kit, VASP, LAMMPS data formats.
dpdata only works with python 3.x.


# Installation
One can download the source code of dpdata by 
```bash
git clone https://github.com/deepmodeling/dpdata.git dpdata
```
then use `setup.py` to install the module
```bash
cd dpdata
python setup.py install
```

`dpdata` can also by install via pip
```bash
pip3 install dpdata
```


# Quick start

This section gives some examples on how dpdata works. Firstly one needs to import the module in a python 3.x compatible code.
```python
import dpdata
```
The typicall workflow of `dpdata` is 

1. Load data from vasp or lammps or deepmd-kit data files.
2. Manipulate data 
3. Dump data to in a desired format


## Load data
```python
d_poscar = dpdata.System('POSCAR', fmt = 'vasp/poscar')
```
or let dpdata infer the format (`vasp/poscar`) of the file from the file name extension
```python
d_poscar = dpdata.System('my.POSCAR')
```
The number of atoms, atom types, coordinates are loaded from the `POSCAR` and stored to a data `System` called `d_poscar`.
A data `System` (a concept used by [deepmd-kit](https://github.com/deepmodeling/deepmd-kit)) contains frames that has the same number of atoms of the same type. The order of the atoms should be consistent among the frames in one `System`. 
It is noted that `POSCAR` only contains one frame.
If the multiple frames stored in, for example, a `OUTCAR` is wanted, 
```python
d_outcar = dpdata.LabeledSystem('OUTCAR')
```
The labels provided in the `OUTCAR`, i.e. energies, forces and virials (if any), are loaded by `LabeledSystem`. It is noted that the forces of atoms are always assumed to exist. `LabeledSystem` is a derived class of `System`.

The `System` or `LabeledSystem` can be constructed from the following file formats with the `format key` in the table passed to argument `fmt`:

The Class `dpdata.MultiSystems`  can read data  from a dir which may contains many files of different systems, or from single xyz file which contains different systems.

Use `dpdata.MultiSystems.from_dir` to read from a  directory, `dpdata.MultiSystems` will walk in the directory 
Recursively  and  find all file with specific file_name. Supports all the file formats that `dpdata.LabeledSystem` supports.

Use  `dpdata.MultiSystems.from_file` to read from single file. Now only support quip/gap/xyz  format file.

For example, for `quip/gap xyz` files, single .xyz file may contain many different configurations with different atom numbers and atom type.

The following commands relating to `Class dpdata.MultiSystems` may be useful.
```python
# load data

xyz_multi_systems = dpdata.MultiSystems.from_file(file_name='tests/xyz/xyz_unittest.xyz',fmt='quip/gap/xyz')
vasp_multi_systems = dpdata.MultiSystems.from_dir(dir_name='./mgal_outcar', file_name='OUTCAR', fmt='vasp/outcar')

# use wildcard
vasp_multi_systems = dpdata.MultiSystems.from_dir(dir_name='./mgal_outcar', file_name='*OUTCAR', fmt='vasp/outcar')

# print the multi_system infomation
print(xyz_multi_systems)
print(xyz_multi_systems.systems) # return a dictionaries

# print the system infomation
print(xyz_multi_systems.systems['B1C9'].data)

# dump a system's data to ./my_work_dir/B1C9_raw folder
xyz_multi_systems.systems['B1C9'].to_deepmd_raw('./my_work_dir/B1C9_raw')

# dump all systems
xyz_multi_systems.to_deepmd_raw('./my_deepmd_data/')


```

| Software| format | multi frames | labeled | class	    | format key    |
| ------- | :---   | :---:        | :---:   | :---          | :---          |
| vasp	  | poscar | False        | False   | System	    | 'vasp/poscar' | 
| vasp    | outcar | True         | True    | LabeledSystem | 'vasp/outcar' |	
| vasp    | xml    | True         | True    | LabeledSystem | 'vasp/xml'    |	
| lammps  | lmp    | False        | False   | System        | 'lammps/lmp'  |
| lammps  | dump   | True         | False   | System        | 'lammps/dump' |
| deepmd  | raw    | True         | False   | System	    | 'deepmd/raw'  |
| deepmd  | npy    | True         | False   | System        | 'deepmd/npy'  |
| deepmd  | raw    | True         | True    | LabeledSystem | 'deepmd/raw'  |
| deepmd  | npy    | True         | True    | LabeledSystem | 'deepmd/npy'  |
| gaussian| log    | False        | True    | LabeledSystem | 'gaussian/log'|
| gaussian| log    | True         | True    | LabeledSystem | 'gaussian/md' |
| siesta| output    | False        | True    | LabeledSystem | 'siesta/output'|
| siesta| aimd_output    | True         | True    | LabeledSystem | 'siesta/aimd_output' |
| cp2k    | output | False        | True    | LabeledSystem | 'cp2k/output' |
| QE      | log    | False        | True    | LabeledSystem | 'qe/pw/scf'   |
| QE      | log    | True         | False   | System        | 'qe/cp/traj'  |
| QE      | log    | True         | True    | LabeledSystem | 'qe/cp/traj'  |
|quip/gap|xyz|True|True|MultiSystems|'quip/gap/xyz'|

## Access data
These properties stored in `System` and `LabeledSystem` can be accessed by operator `[]` with the key of the property supplied, for example
```python
coords = d_outcar['coords']
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



## Dump data
The data stored in `System` or `LabeledSystem` can be dumped in 'lammps/lmp' or 'vasp/poscar' format, for example:
```python
d_outcar.to_lammps_lmp('conf.lmp', frame_idx=0)
```
The first frames of `d_outcar` will be dumped to 'conf.lmp'
```python
d_outcar.to_vasp_poscar('POSCAR', frame_idx=-1)
```
The last frames of `d_outcar` will be dumped to 'POSCAR'.


The data stored in `LabeledSystem` can be dumped to deepmd-kit raw format, for example
```python
d_outcar.to_deepmd_raw('dpmd_raw')
```
Or a simpler command:
```python
dpdata.LabeledSystem('OUTCAR').to_deepmd_raw('dpmd_raw')
```
Frame selection can be implemented by
```python
dpdata.LabeledSystem('OUTCAR').sub_system([0,-1]).to_deepmd_raw('dpmd_raw')
```
by which only the first and last frames are dumped to `dpmd_raw`.

# replicate 
dpdata will create a super cell of the current atom configuration.

```python
dpdata.System('./POSCAR').replicate((1,2,3,) )
```
tuple(1,2,3) means don't copy atom configuration in x direction, make 2 copys in y direction, make 3 copys in z direction.

# perturb
dpdata will disturb the box size and shape and change atom coordinates  randomly.
```python3
dpdata.System('./POSCAR').perturb(pert_num=3, 
    box_pert_fraction=0.02, 
    atom_pert_fraction=0.01, 
    atom_pert_style='normal')
```

### `pert_num`

Each frame in the input system will generate `pert_num` frames.
 That means the command will return a system containing `frames of input system * pert_num` frames.
 
### `box_pert_fraction`
A relative length that determines the length the box size change in each frame. It is just a fraction and doesn't have unit. Typlicaly, for cubic box with side length `a` ,  `a` will increase or decrease a random value from the intervel `[-a*box_pert_fraction, a*box_pert_fraction]`.

It will also change the shape of the box. That means an orthogonal box will become a non-orthogonal box after perturbing. The angle of inclination of the box is a random variable.
 `box_pert_fraction` is also relating to the probability distribution function of the angle.
 
See more details about how it will change the box below.

### `atom_pert_fraction`
A relative length that determines the length atom moves in each frame. It is just a fraction and doesn't have unit. Typlicaly, for a cubic box with side length `a` , the mean value of the distance that atom  moves is approximate `a*atom_pert_fraction`.

### `atom_pert_style`
The probability distribution function used to change atom coordinates.
available options:`'uniform' 'normal' 'const'`. 

`uniform` means that if the box is a cube with side length `a`, how far atoms in the cube move is a random vector with max length `a*atom_pert_fraction `

`normal` means the squares of the distance atoms move are subject to  a chi-square distribution (chi-square distribution can be seen as the sum of squares of normal distributed random variable. This is why we  name this option as 'normal'.).
If the box is a cube with side length `a`, the mean value of the distance atom moves is `a*atom_pert_fraction `

`const` means that if the box is a cube with side length `a`, the distance atom moves is always `a*atom_pert_fraction `(For triclinic box, the distances are not equal.)

The direction atoms move and box deformation is random.

See more details about how atoms will move below.

## The perturb details
---
For each frame in the input system,dpdata will repeat the following steps `pert_num` times. That means the command will return a system containing `frames of input system * pert_num` frames.

### first step: box deform, and atom moves correspondingly

#### 1. generate a perturb matrix for box and atoms
---
dpdata will generate a perturb matrix

$L_1=\begin{pmatrix}
n_1+1 & n_2/2 & n_4/2 \\
n_2/2 & n_3+1 & n_5/2 \\
n_4/2 & n_5/2 & n_6+1
\end{pmatrix}$

$n_i,i\in\{1,2,3,4,5,6\}$ are independent variables which are subject to uniform
distrubution  on interval $[-box\_pert\_fraction，box\_pert\_fraction]$. 
That is $n_i \sim U(-box\_pert\_fraction，box\_pert\_fraction)$

#### 2.box deforms
---
The origin box matrix of this frame is defined as $B_o$, usually with lower lower triangular matrix form. That is:

$origin\_box\_matrix \equiv B_o=\begin{pmatrix}
xx &  0 & 0 \\
xy & yy &  0\\
xz & yz & zz
\end{pmatrix}$

The box matrix after deformation will be
$deformed\_box\_matrix \equiv B_d = B_o L_1$

#### 3.atoms relocate in the new box.
---
The atom coordinates in this frame will change correspondingly. 
That is,
for atom with index $j, j \in\{0,1,2,3,...,atom\_num-1\}$ in the frame with coordinate $\vec r_j=(x_j, y_j, z_j)$,
the coordinate after deformation $\vec r_{jd}$ will become the matrix multiplication of $r_j$ and $L_1$.
That is:
$\vec r_{jd} = \vec r_j L_1$

### second step: perturb atoms randomly.

#### 1. generate a random vector $\vec l_{j}$ for each atom 
---

For each atom with index $j, j \in\{0,1,2,3,...,atom\_num-1\}$ in the frame with coordinate $\vec r_{jd}=(x_{jd}, y_{jd}, z_{jd})$, dpdata will generate a random vector $\vec l_{j}=(l_{jx}, l_{jy}, l_{jz})$ .

The method used to generate $\vec l_j$ is described below.

:tada: if pert_style is 'normal':


$l_{jx},l_{jy},l_{jz}$ are independent random variables which are subject to normal distribution with mean $\mu=0$ and variance $\sigma^2 = atom\_pert\_fraction^2/3$ . 
That is:
$l_{jx},l_{jy},l_{jz} \sim N(0,atom\_pert\_fraction^2/3)$

$\Vert \vec l_{j} \Vert_2$ is the distance the atom moves, and it satisfies the following equation   $\Vert \vec l_{j}  \Vert_2^2=l_{jx}^2 + l_{jy}^2+l_{jz}^2$. 
$3\Vert \vec l_{j}  \Vert_2^2/atom\_pert\_fraction^2$ will be object to the chi-square distribution with 3 degrees of freedom. That is:
$\frac{3\Vert \vec l_{j}  \Vert_2^2}{atom\_pert\_fraction^2} \sim \chi(3)$
The expectation value of $\Vert \vec l_{j}  \Vert_2$ will be `atom_pert_fraction` .That is:
$E(\Vert \vec l_{j}  \Vert_2)=atom\_pert\_fraction$

:tada: if pert_style is 'uniform':

$\vec l_{j}$ will be a vector point to a random point inside a 3D unit sphere and its internal space.

That is:
$\{ \vec l_{j} \in \mathbb{R^3} |  \Vert \vec l_{j} \Vert_2 \leq atom\_pert\_fraction \}$.

The point is chosen with equal probability. That means for arbitrary two subsets of the 3D sphere and its internal space $\{(x,y,z)|x,y,z\in\mathbb{R},x^2+y^2+z^2\leq atom\_pert\_fraction^2\}$, if the 'volume' of the subsets are equal, the probability that the random point is in them are equal.

The following method is used to generate such $\vec l_{j}$.

random direction  of equal probability
>let $\vec x$ become a 3 dimension standard normal distribution. That is
$\vec x=(x_1,x_2,x_3), x_i(i\in\{1,2,3\})\ are\ independent.$
$x_i \sim N(0,1)$
and 
$\vec l_{j,unit\_surface}=\vec x /\Vert \vec x \Vert_2$
Now $\vec l_{j,unit\_surface}$ point to a point located at the surface of the 3D unit sphere. 

random point of equal probability
> let $u$ become a random variable which is object to the uniform on interval [0,1]. That is:
$u \sim U(0,1)$
and define $v$ as the 3th root of $u$ (because it is 3 dimension), That is:
$v\equiv u^{1/3}$
Then the target vector $\vec l_{j}$ is
$\vec l_{j}=atom\_pert\_fraction \cdot v \cdot \vec l_{j,unit\_surface}$

:tada: if pert_style is 'const':

$\vec l_j \equiv atom\_pert\_fraction \cdot \vec l_{j,unit\_surface}$

The definition of $\vec l_{j,unit\_surface}$  is described in detail above.
and it is obvious that
$\Vert \vec l_{j}  \Vert_2=atom\_pert\_fraction$
####  2. atom moves according to $\vec l_{j}$
---

After the $\vec l_j$ generated, the atom coordinates disturbed $\vec {r_{jd,disturbed}}$ will become
$\vec {r_{jd,disturbed}}=\vec r_{jd}+\vec l_{j} B_d$


$B_d$ is  the box matrix after deformatiion, described in first step
$\vec r_{jd}$ is the coordinate of atom j after deformatiion, described in first step.

### the results the method System.perturb() return
---

At the beginning, dpdata will create an empty system.
```python
perturbed_system = System()
```
After every single perturbed frame is generated, it will be appended to the `perturbed_system`.

perturbed_system will contain `pert_num * frames of the input system` frames.

$B_d$ will be used as the final box matrix in the perturbed frame(that is `System.data['cells'][0]`)
$\vec {r_{jd,disturbed}}$ will be used as the final coordinates in the perturbed frame (that is `System.data['coords'][0][j]`)
> note: 0 means the index of the frame, j means atom index, 

 Finally dpdata will return `perturbed_system` as results.
 
