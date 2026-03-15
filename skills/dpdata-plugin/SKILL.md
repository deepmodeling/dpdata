---
name: dpdata-plugin
description: Create and install dpdata plugins (especially custom Format readers/writers) using Format.register(...) and pyproject.toml entry_points under 'dpdata.plugins'. Use when extending dpdata with new formats or distributing plugins as separate Python packages.
---

# dpdata-plugin

dpdata loads plugins in two ways:

1. **Built-in plugins** in `dpdata.plugins.*` (imported automatically)
1. **External plugins** exposed via Python package entry points: `dpdata.plugins`

This skill focuses on **external plugin packages**, the recommended way to add new formats without modifying dpdata itself.

## What can be extended?

Most commonly: add a new **Format** (file reader/writer) via:

```python
from dpdata.format import Format


@Format.register("myfmt")
class MyFormat(Format): ...
```

## How dpdata discovers plugins

dpdata imports `dpdata.plugins` during normal use (e.g. `dpdata.system` imports it). That module:

- imports every built-in module in `dpdata/plugins/*.py`
- then loads all **entry points** in group `dpdata.plugins`

So an external plugin package only needs to ensure that importing the entry-point target triggers the `@Format.register(...)` side effects.

## Minimal external plugin package (based on plugin_example/)

### 1) Create a new Python package

Example layout:

```text
dpdata_random/
  pyproject.toml
  dpdata_random/
    __init__.py
```

### 2) Implement and register your Format

In `dpdata_random/__init__.py` (shortened example):

```python
from __future__ import annotations

import numpy as np
from dpdata.format import Format


@Format.register("random")
class RandomFormat(Format):
    def from_system(self, N, **kwargs):
        return {
            "atom_numbs": [20],
            "atom_names": ["X"],
            "atom_types": [0] * 20,
            "cells": np.repeat(np.eye(3)[None, ...], N, axis=0) * 100.0,
            "coords": np.random.rand(N, 20, 3) * 100.0,
            "orig": np.zeros(3),
            "nopbc": False,
        }
```

Return dicts must match dpdata’s expected schema (cells/coords/atom_names/atom_types/...).

### 3) Expose an entry point

In `pyproject.toml`:

```toml
[project]
name = "dpdata_random"
version = "0.0.0"
dependencies = ["numpy", "dpdata"]

[project.entry-points.'dpdata.plugins']
random = "dpdata_random:RandomFormat"
```

Any importable target works; this pattern points directly at the class.

### 4) Install and test

In a clean env (recommended via `uv`):

```bash
uv run --with dpdata --with numpy python3 - <<'PY'
import dpdata
from dpdata.format import Format

# importing dpdata will load entry points (dpdata.plugins)
print('random' in Format.get_formats())
PY
```

If it prints `True`, your plugin was discovered.

## Debug checklist

- Did you install the plugin package into the same environment where you run dpdata?
- Does `pyproject.toml` contain `[project.entry-points.'dpdata.plugins']`?
- Does importing the entry point module/class execute the `@Format.register(...)` decorator?
- If using `uv run`, remember each command runs in its own environment unless you’re in a `uv` project (or you rely on `uv run --with ...`).
