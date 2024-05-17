from __future__ import annotations

import subprocess
import sys

import pytest


@pytest.mark.benchmark
def test_import():
    """Test import dpdata."""
    subprocess.check_output(
        [sys.executable, "-c", "'from dpdata import LabeledSystem'"]
    ).decode("ascii")


@pytest.mark.benchmark
def test_cli():
    """Test dpdata command."""
    subprocess.check_output([sys.executable, "-m", "dpdata", "-h"]).decode("ascii")
