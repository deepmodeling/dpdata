from __future__ import annotations

import os
import shutil
import tempfile
import unittest

from context import dpdata

FIXTURE = os.path.join("cp2k", "aimd")


class TestCp2kAimdMissingInputs(unittest.TestCase):
    """``cp2k/aimd_output`` locates its inputs by globbing a directory.

    A directory missing one of them used to fail with a bare ``IndexError``
    from subscripting the empty match list.
    """

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.work = os.path.join(self.tmp_dir, "aimd")
        shutil.copytree(FIXTURE, self.work)

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def _remove(self, name):
        os.remove(os.path.join(self.work, name))

    def test_complete_directory_still_loads(self):
        system = dpdata.LabeledSystem(self.work, fmt="cp2k/aimd_output")
        self.assertGreater(system.get_nframes(), 0)

    def test_missing_trajectory_names_the_pattern(self):
        self._remove("DPGEN-pos-1.xyz")
        with self.assertRaises(FileNotFoundError) as caught:
            dpdata.LabeledSystem(self.work, fmt="cp2k/aimd_output")
        message = str(caught.exception)
        self.assertIn("*pos*.xyz", message)
        self.assertIn("trajectory file", message)
        # The listing tells the user what the parser did see.
        self.assertIn("cp2k.log", message)

    def test_missing_log_names_the_pattern(self):
        self._remove("cp2k.log")
        with self.assertRaises(FileNotFoundError) as caught:
            dpdata.LabeledSystem(self.work, fmt="cp2k/aimd_output")
        message = str(caught.exception)
        self.assertIn("*.log", message)
        self.assertIn("output log", message)

    def test_pointing_at_a_file_says_so(self):
        # A frequent mistake: passing the log file rather than its directory.
        log = os.path.join(self.work, "cp2k.log")
        with self.assertRaises(FileNotFoundError) as caught:
            dpdata.LabeledSystem(log, fmt="cp2k/aimd_output")
        self.assertIn("is not a directory", str(caught.exception))

    def test_missing_directory_says_so(self):
        absent = os.path.join(self.tmp_dir, "absent")
        with self.assertRaises(FileNotFoundError) as caught:
            dpdata.LabeledSystem(absent, fmt="cp2k/aimd_output")
        self.assertIn("is not a directory", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
