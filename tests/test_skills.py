from __future__ import annotations

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from dpdata.cli import dpdata_cli, dpdata_parser
from dpdata.skills import _read, _root, skills_cli


class TestSkills(unittest.TestCase):
    def test_list_and_read(self):
        root = _root()
        expected = sorted(
            p.name for p in root.iterdir() if p.joinpath("SKILL.md").is_file()
        )
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            skills_cli(["list"])
        self.assertEqual(output.getvalue().splitlines(), expected)
        self.assertTrue(expected)
        for name in expected:
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                skills_cli(["read", name])
            self.assertEqual(
                output.getvalue(),
                root.joinpath(name, "SKILL.md").read_bytes().decode("utf-8"),
            )

    def test_reference(self):
        self.assertEqual(
            _read(_root(), "dpdata-cli", "README.md"),
            _root().joinpath("dpdata-cli", "README.md").read_bytes().decode("utf-8"),
        )

    def test_bad_paths(self):
        for path in [
            "/etc/passwd",
            "../SKILL.md",
            "references/../../SKILL.md",
            "a\\b",
            "C:/a",
            "a\x00b",
            "./SKILL.md",
            "missing.md",
        ]:
            with self.subTest(path=path), self.assertRaises(ValueError):
                _read(_root(), "dpdata-cli", path)
        for name in ["../dpdata-cli", "/dpdata-cli", "missing"]:
            with self.subTest(name=name), self.assertRaises(ValueError):
                _read(_root(), name, "SKILL.md")

    def test_symlink_escape(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "skills"
            skill = root / "example"
            skill.mkdir(parents=True)
            (skill / "SKILL.md").write_text("example", encoding="utf-8")
            secret = Path(directory) / "secret"
            secret.write_text("outside", encoding="utf-8")
            (skill / "escape.md").symlink_to(secret)
            with self.assertRaises(ValueError):
                _read(root, "example", "escape.md")

    def test_help_and_legacy_parser(self):
        self.assertIn("skills list", dpdata_parser().format_help())
        for argv in [["POSCAR", "-n"], ["./skills"], ["--", "skills"]]:
            self.assertEqual(
                dpdata_parser().parse_args(argv).from_file,
                argv[0] if argv[0] != "--" else "skills",
            )

    def test_no_conversion(self):
        with (
            patch("sys.argv", ["dpdata", "skills", "list"]),
            patch("dpdata.cli.convert") as convert,
            contextlib.redirect_stdout(io.StringIO()),
        ):
            dpdata_cli()
        convert.assert_not_called()

    def test_unknown_is_usage_error(self):
        with (
            contextlib.redirect_stderr(io.StringIO()),
            self.assertRaises(SystemExit) as error,
        ):
            skills_cli(["read", "missing"])
        self.assertEqual(error.exception.code, 2)
