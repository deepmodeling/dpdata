---
name: dpdata-minimizer
description: Minimize geometries with dpdata minimizer plugins via System.minimize(), including how minimizers relate to drivers (ASEMinimizer needs a dpdata Driver) and how to list supported minimizers (ase/sqm). Use when doing geometry optimization/minimization through dpdata Python API.
---

# dpdata-minimizer

Use dpdata “minimizer plugins” to **optimize/minimize geometry** and return a `dpdata.LabeledSystem`.

## Key idea

- A **Minimizer** performs geometry optimization (updates coordinates) and returns labeled results.
- dpdata exposes this as:

```python
System.minimize(*args, minimizer: str|Minimizer, **kwargs) -> LabeledSystem
```

## List supported minimizer keys (runtime)

```python
from dpdata.driver import Minimizer

print(sorted(Minimizer.get_minimizers().keys()))
```

In the current dpdata repo, minimizer keys include:

- `ase`
- `sqm`

## Relationship to drivers (important)

Some minimizers require a **dpdata Driver object**.

Example: `ASEMinimizer` takes a dpdata `Driver` in its constructor:

- `minimizer="ase"` requires `driver=<dpdata Driver>` (e.g. the ASE driver wrapping an ASE calculator).

So you generally do:

1. Construct a driver
1. Construct a minimizer (or let dpdata do it by passing the right kwargs)
1. Call `System.minimize(...)`

## Runnable example: ASE minimizer with an ASE calculator

Use uv inline script metadata so the example runs reproducibly with `uv run`.

```python
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "dpdata",
#   "numpy",
#   "ase",
# ]
# ///

import numpy as np
from ase.calculators.emt import EMT

from dpdata.driver import Driver
from dpdata.system import System

open("tmp.xyz", "w").write("""2\n\nH 0 0 0\nH 0 0 0.74\n""")

sys = System("tmp.xyz", fmt="xyz")

# Build a dpdata driver that can provide energies/forces to ASE optimizers.
ase_driver = Driver.get_driver("ase")(calculator=EMT())

# Minimize using the ASE minimizer plugin.
# NOTE: ASEMinimizer expects `driver` (not `calculator`) as input.
ls = sys.minimize(minimizer="ase", driver=ase_driver, fmax=0.05, max_steps=5)

print("coords", np.array(ls.data["coords"]).shape)
print("energies", np.array(ls.data["energies"]))
print("forces", np.array(ls.data["forces"]).shape)
```

## Notes / gotchas

- `System.minimize(...)` accepts either a minimizer key string or a Minimizer object.
- If you previously used `System.predict(driver="ase", calculator=...)`, be aware that minimization is different: you need to pass a **driver** into the minimizer (ASEMinimizer does not accept `calculator=`).
- `sqm` minimizer requires AmberTools `sqm` executable and typically won’t be runnable in CI.
