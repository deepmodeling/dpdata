"""Discoverable workflow guidance belongs to help, not command output."""

from __future__ import annotations

import contextlib
import io
import unittest
from unittest.mock import patch

from dpdata.cli import dpdata_parser
from dpdata.skills import workflow_help


class TestSkillHelp(unittest.TestCase):
    def test_contextual_help(self) -> None:
        for command, skill in {"": "dpdata-cli"}.items():
            with self.subTest(command=command):
                output = io.StringIO()
                with (
                    contextlib.redirect_stdout(output),
                    self.assertRaises(SystemExit) as result,
                ):
                    dpdata_parser().parse_args(
                        ([command] if command else []) + ["--help"]
                    )
                self.assertEqual(result.exception.code, 0)
                self.assertIn(f"dpdata skills read {skill}", output.getvalue())

    def test_missing_guide_has_no_pointer(self) -> None:
        self.assertEqual(workflow_help("not-a-bundled-skill"), "")
        with patch("dpdata.skills._root", side_effect=ValueError("missing resources")):
            self.assertEqual(workflow_help("dpdata-cli"), "")
