from __future__ import annotations

import io
import unittest
from pathlib import Path

from dpdata.utils import read_text_file


class TestReadFile(unittest.TestCase):
    def test_read_text_file_from_string_io(self):
        string_io = io.StringIO("Hello, world!")
        with read_text_file(string_io) as file:
            self.assertEqual(file.read(), "Hello, world!")

    def test_read_text_file_from_file_str(self):
        with read_text_file("/proc/cpuinfo") as file:
            self.assertEqual(file.read(), Path("/proc/cpuinfo").read_text())
        
    def test_read_text_file_from_file_path(self):
        with read_text_file(Path("/proc/cpuinfo")) as file:
            self.assertEqual(file.read(), Path("/proc/cpuinfo").read_text())
