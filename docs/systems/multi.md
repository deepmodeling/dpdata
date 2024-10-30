# `MultiSystems`

The Class {class}`dpdata.MultiSystems`  can read data  from a dir which may contains many files of different systems, or from single xyz file which contains different systems.

Use {meth}`dpdata.MultiSystems.from_dir` to read from a  directory, {class}`dpdata.MultiSystems` will walk in the directory
Recursively  and  find all file with specific file_name. Supports all the file formats that {class}`dpdata.LabeledSystem` supports.

Use {meth}`dpdata.MultiSystems.from_file` to read from single file. Single-file support is available for the `quip/gap/xyz` and `ase/structure` formats.

For example, for `quip/gap xyz` files, single .xyz file may contain many different configurations with different atom numbers and atom type.

The following commands relating to {class}`dpdata.MultiSystems` may be useful.
```python
# load data

xyz_multi_systems = dpdata.MultiSystems.from_file(
    file_name="tests/xyz/xyz_unittest.xyz", fmt="quip/gap/xyz"
)
vasp_multi_systems = dpdata.MultiSystems.from_dir(
    dir_name="./mgal_outcar", file_name="OUTCAR", fmt="vasp/outcar"
)

# use wildcard
vasp_multi_systems = dpdata.MultiSystems.from_dir(
    dir_name="./mgal_outcar", file_name="*OUTCAR", fmt="vasp/outcar"
)

# print the multi_system infomation
print(xyz_multi_systems)
print(xyz_multi_systems.systems)  # return a dictionaries

# print the system infomation
print(xyz_multi_systems.systems["B1C9"].data)

# dump a system's data to ./my_work_dir/B1C9_raw folder
xyz_multi_systems.systems["B1C9"].to_deepmd_raw("./my_work_dir/B1C9_raw")

# dump all systems
xyz_multi_systems.to_deepmd_raw("./my_deepmd_data/")
```

You may also use the following code to parse muti-system:
```python
from dpdata import LabeledSystem, MultiSystems
from glob import glob

"""
process multi systems
"""
fs = glob("./*/OUTCAR")  # remeber to change here !!!
ms = MultiSystems()
for f in fs:
    try:
        ls = LabeledSystem(f)
    except:
        print(f)
    if len(ls) > 0:
        ms.append(ls)

ms.to_deepmd_raw("deepmd")
ms.to_deepmd_npy("deepmd")
```
