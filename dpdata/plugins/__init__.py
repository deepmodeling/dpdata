from __future__ import annotations

import importlib
from importlib import metadata
from pathlib import Path

PACKAGE_BASE = "dpdata.plugins"
NOT_LOADABLE = ("__init__.py",)

for module_file in Path(__file__).parent.glob("*.py"):
    if module_file.name not in NOT_LOADABLE:
        module_name = f".{module_file.stem}"
        importlib.import_module(module_name, PACKAGE_BASE)

# https://setuptools.readthedocs.io/en/latest/userguide/entry_point.html
try:
    eps = metadata.entry_points(group="dpdata.plugins")
except TypeError:
    eps = metadata.entry_points().get("dpdata.plugins", [])
for ep in eps:
    plugin = ep.load()
