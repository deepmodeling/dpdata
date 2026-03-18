---
name: dpdata-driver
description: Use dpdata Python Driver plugins to label systems (energies/forces/virials) via System.predict(), list available drivers, and build Driver objects (ase/deepmd/gaussian/sqm/hybrid). Use when working with dpdata Python API (not CLI) and you need driver-based energy/force prediction, plugin registration keys, or examples of using dpdata with ASE calculators or DeePMD models.
---

# dpdata-driver

Use dpdata “driver plugins” to **label** a `dpdata.System` (predict energies/forces/virials) and obtain a `dpdata.LabeledSystem`.

## Key idea

- A **Driver** converts an unlabeled `System` into a `LabeledSystem` by computing:
  - `energies` (required)
  - `forces` (optional but common)
  - `virials` (optional)

In dpdata, this is exposed as:

- `System.predict(*args, driver="dp", **kwargs) -> LabeledSystem`

`driver` can be:

- a **string key** (plugin name), e.g. `"ase"`, `"dp"`, `"gaussian"`
- a **Driver object**, e.g. `Driver.get_driver("ase")(...)`

## List supported driver keys (runtime)

When unsure what drivers exist in *this* dpdata version/env, query them at runtime:

```python
from dpdata.driver import Driver

print(sorted(Driver.get_drivers().keys()))
```

In the current repo state, keys include:

- `ase`
- `dp` / `deepmd` / `deepmd-kit`
- `gaussian`
- `sqm`
- `hybrid`

(Exact set depends on dpdata version and installed extras.)

## Minimal workflow

```python
import dpdata
from dpdata.system import System

sys = System("input.xyz", fmt="xyz")
ls = sys.predict(driver="ase", calculator=...)  # returns dpdata.LabeledSystem
```

### Verify you got a labeled system

```python
assert "energies" in ls.data
# optional:
# assert "forces" in ls.data
# assert "virials" in ls.data
```

## Example: use the ASE driver with an ASE calculator (runnable)

This is the easiest *fully runnable* example because it doesn’t require external QM software.

Dependencies (recommended): declare script dependencies with uv inline metadata, then run with `uv run`.

```python
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "dpdata",
#   "numpy",
#   "ase",
# ]
# ///
```

Script:

```python
import numpy as np
from ase.calculators.emt import EMT
from dpdata.system import System

# write a tiny molecule
open("tmp.xyz", "w").write("""2\n\nH 0 0 0\nH 0 0 0.74\n""")

sys = System("tmp.xyz", fmt="xyz")
ls = sys.predict(driver="ase", calculator=EMT())

print("energies", np.array(ls.data["energies"]))
print("forces shape", np.array(ls.data["forces"]).shape)
if "virials" in ls.data:
    print("virials shape", np.array(ls.data["virials"]).shape)
else:
    print("virials: <not provided by this driver/calculator>")
```

## Example: pass a Driver object instead of a string

```python
from ase.calculators.emt import EMT
from dpdata.driver import Driver
from dpdata.system import System

sys = System("tmp.xyz", fmt="xyz")
ase_driver = Driver.get_driver("ase")(calculator=EMT())
ls = sys.predict(driver=ase_driver)
```

## Hybrid driver

Use `driver="hybrid"` to sum energies/forces/virials from multiple drivers.

The `HybridDriver` accepts `drivers=[ ... ]` where each item is either:

- a `Driver` instance
- a dict like `{"type": "sqm", ...}` (type is the driver key)

Example (structure only; may require external executables):

```python
from dpdata.driver import Driver

hyb = Driver.get_driver("hybrid")(
    drivers=[
        {"type": "sqm", "qm_theory": "DFTB3"},
        {"type": "dp", "dp": "frozen_model.pb"},
    ]
)
# ls = sys.predict(driver=hyb)
```

## Notes / gotchas

- Many drivers require extra dependencies or external programs:
  - `dp` requires `deepmd-kit` + a model file
  - `gaussian` requires Gaussian and a valid executable (default `g16`)
  - `sqm` requires AmberTools `sqm`
- If you just need file format conversion, use the existing **dpdata CLI** skill instead.
