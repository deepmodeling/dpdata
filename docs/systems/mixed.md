# Mixed Type Format

The format `deepmd/npy/mixed` is the mixed type numpy format for DeePMD-kit, and can be loaded or dumped through class {class}`dpdata.MultiSystems`.

Under this format, systems with the same number of atoms but different formula can be put together
for a larger system, especially when the frame numbers in systems are sparse.

This also helps to mixture the type information together for model training with type embedding network.

Here are examples using `deepmd/npy/mixed` format:

- Dump a MultiSystems into a mixed type numpy directory:
```python
import dpdata

dpdata.MultiSystems(*systems).to_deepmd_npy_mixed("mixed_dir")
```

- Load a mixed type data into a MultiSystems:
```python
import dpdata

dpdata.MultiSystems().load_systems_from_file("mixed_dir", fmt="deepmd/npy/mixed")
```
