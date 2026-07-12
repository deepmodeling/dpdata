# LMDB Format

The format `lmdb` stores the frames of one or more systems in a single [LMDB](http://www.lmdb.tech/doc/) database, and can be loaded or dumped through {class}`dpdata.System`, {class}`dpdata.LabeledSystem`, and {class}`dpdata.MultiSystems`. The on-disk layout coincides with the LMDB datasets read by the DeePMD-kit data loader. Core fields and registered additional fields are therefore available to DeePMD-kit under their `deepmd_name`.

In contrast to the directory-based `deepmd/npy` format, every frame is stored as an independent record indexed by a global frame number. Frames within one database may therefore differ in the number of atoms and in chemical composition, which is suited to data sets in which the number of frames per system is small.

## Layout

A database contains one metadata record together with one record per frame.

The metadata record is stored under the key `__metadata__` and holds the following entries.

| key                | description                                                        |
| ------------------ | ------------------------------------------------------------------ |
| `nframes`          | total number of frames                                             |
| `frame_idx_fmt`    | format string of the integer frame key (default `012d`)            |
| `type_map`         | the global element table                                           |
| `frame_nlocs`      | the number of atoms of each frame                                  |
| `frame_system_ids` | the index of the source system of each frame                       |
| `dp_data_shapes`   | optional symbolic shapes of additional dpdata fields               |
| `dp_data_names`    | optional mapping from DeePMD-kit field names to dpdata field names |

Each frame is stored under its zero-padded global index (`000000000000`, `000000000001`, ...). The per-frame arrays `coords`, `cells`, `energies`, `forces`, `virials`, and `atom_types` are stored together with their data type and shape; `atom_types` records global indices into `type_map`, and `atom_numbs` is a list of per-element counts over `type_map`. The numerical dtype of each source array is retained, except that `atom_types` is stored as `int32`.

## Dump data

A single system or a collection of systems is written to one database.

```python
import dpdata

dpdata.LabeledSystem("OUTCAR", fmt="vasp/outcar").to("lmdb", "data.lmdb")

dpdata.MultiSystems(*systems).to("lmdb", "data.lmdb")
```

The element table recorded in `type_map` defaults to the union of the elements present in the data. An explicit table may be supplied through the `type_map` argument, for example the full periodic table.

```python
from dpdata.periodic_table import ELEMENTS

dpdata.MultiSystems(*systems).to("lmdb", "data.lmdb", type_map=list(ELEMENTS))
```

Frames are committed in batches of 1,000 by default. The `write_batch_size` argument changes the transaction size. If a transaction exceeds `map_size`, the map is enlarged and the same encoded batch is retried before frame counters are advanced. The destination must not exist unless `overwrite=True` is supplied.

Data are first written to a temporary sibling database. After closing it, dpdata reopens the staged database and validates the metadata, entry count, and continuous frame-key sequence. A new destination is published by one directory rename. On POSIX systems, `overwrite=True` replaces the existing `data.mdb` with the staged `data.mdb` in one file rename, while the LMDB runtime lock file remains in place. Windows does not permit safe replacement of an open memory-mapped LMDB file, so `overwrite=True` is not supported there. A failed conversion or validation leaves the previous destination unchanged. All dpdata readers of the destination must be closed before replacement; readers owned by other libraries must likewise be closed because their environments cannot be coordinated through dpdata's process-local cache.

## Preserving the system partition

The `frame_system_ids` entry of the metadata records, for every frame, the index of the source system it belongs to. DeePMD-kit uses this partition for system-wise sampling, for example through `prob_sys_size`.

When a {class}`dpdata.MultiSystems` is dumped through `to("lmdb", ...)`, its frames are first grouped by chemical formula, and systems that share a formula are merged into a single entry. The resulting `frame_system_ids` therefore reflect the formula grouping rather than the original sources, and the number of systems may be smaller than the number of inputs.

When the original partition must be retained, the function {func}`dpdata.formats.lmdb.dump_systems` writes an ordered sequence of systems without formula merging; each input becomes one system, numbered in iteration order.

```python
import dpdata
from dpdata.formats.lmdb import dump_systems

systems = [dpdata.LabeledSystem(d, fmt="deepmd/npy") for d in directories]
dump_systems(systems, "data.lmdb", type_map=["H", "C", "N", "O"])
```

If `type_map` is supplied, the systems may be provided as a generator, so that the sequence need not be held in memory.

```python
def load():
    for d in directories:
        yield dpdata.LabeledSystem(d, fmt="deepmd/npy")


dump_systems(load(), "data.lmdb", type_map=["H", "C", "N", "O"])
```

Reading an LMDB through dpdata subsequently groups frames by composition and does not reconstruct the original `frame_system_ids` partition. The partition remains available to the DeePMD-kit LMDB data loader.

## Load data

```python
import dpdata

ms = dpdata.MultiSystems.from_file("data.lmdb", fmt="lmdb")
```

Frames are grouped by composition, and each composition is returned as one system. Atom order is canonicalized by a stable sort on the global atom type; coordinates and all registered atomic fields are permuted consistently. Frames of one composition must have identical field sets and a consistent periodic-boundary condition. By default (`mixed_type=False`) the element table of each resulting system is restricted to the elements that the system contains. When `mixed_type=True`, every system retains the complete element set from the database. Loading through {class}`dpdata.MultiSystems` may normalize the order of `atom_names`, so callers should use element names rather than assume that stored numerical indices are retained.

The general {class}`dpdata.System` constructor applies its `type_map` argument after format parsing. Consequently, a direct single-system load with an explicit `type_map` retains that complete requested table even when `mixed_type=False`; this is standard `System` behavior. The compact/full distinction above describes the dictionaries yielded to {class}`dpdata.MultiSystems`.

```python
ms = dpdata.MultiSystems.from_file("data.lmdb", fmt="lmdb", mixed_type=True)
```

A database that holds a single composition may also be read into a {class}`dpdata.LabeledSystem`. If the database contains several compositions, only the first can be represented by a single system and a warning is issued.

```python
ls = dpdata.LabeledSystem("data.lmdb", fmt="lmdb")
```

The reader loads all selected frames into memory. The default `max_frames=100000` guard rejects larger data sets before decoding; set it to `None` only when sufficient memory is available. Large training data sets should normally be consumed directly by the DeePMD-kit data loader.
