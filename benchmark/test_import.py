from __future__ import annotations

import importlib
import subprocess
import sys
from unittest import mock

import pytest


@pytest.mark.parametrize("mod_name", ["dpdata", "dpdata.cli"])
def test_bench_module_import(benchmark, mod_name):
    """Benchmark the import time."""

    @benchmark
    def _():
        with mock.patch("sys.modules", {}):
            importlib.import_module(mod_name, "test_bench_imports")
