"""Read the public workflow skills shipped with this dpdata version."""

from __future__ import annotations

import argparse
import re
import sys
from importlib.resources import files
from pathlib import Path


def _root():
    bundled = files("dpdata").joinpath("_skills")
    if bundled.is_dir():
        return bundled
    # Editable/source checkouts use the same canonical directory as the build.
    source = Path(__file__).resolve().parent.parent
    if (source / "pyproject.toml").is_file() and (source / "skills").is_dir():
        return source / "skills"
    raise ValueError("Bundled skills are missing; reinstall dpdata.")


def _parts(value):
    if "\\" in value or "\x00" in value or ":" in value:
        raise ValueError("Use a relative skill path without backslashes or colons.")
    parts = value.split("/")
    if any(part in ("", ".", "..") for part in parts):
        raise ValueError("Skill paths cannot be absolute or contain '.' or '..'.")
    return parts


def _read(root, name, path):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", name):
        raise ValueError(f"Invalid skill name: {name!r}")
    skill = root.joinpath(name)
    if not skill.joinpath("SKILL.md").is_file():
        raise ValueError(f"Unknown skill: {name}")
    target = skill.joinpath(*_parts(path))
    if isinstance(target, Path) and (
        not skill.resolve().is_relative_to(root.resolve())
        or not target.resolve().is_relative_to(skill.resolve())
    ):
        raise ValueError("Skill paths must stay within the named bundled skill.")
    if not target.is_file():
        raise ValueError(f"Unknown skill file: {name}/{path}")
    try:
        return target.read_bytes().decode("utf-8")
    except UnicodeError as exc:
        raise ValueError("Only UTF-8 text skill files can be read.") from exc


def skills_cli(argv):
    """List and read bundled instructions without running a workflow."""
    parser = argparse.ArgumentParser(
        prog="dpdata skills", description="Read offline, bundled workflow instructions."
    )
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("list", help="list available workflows")
    reading = commands.add_parser("read", help="read a workflow or its reference file")
    reading.add_argument("name", help="skill name from 'dpdata skills list'")
    reading.add_argument(
        "path", nargs="?", default="SKILL.md", help="relative text file"
    )
    args = parser.parse_args(argv)
    try:
        root = _root()
        if args.command == "list":
            names = sorted(
                item.name
                for item in root.iterdir()
                if item.joinpath("SKILL.md").is_file()
            )
            print("\n".join(names))
        else:
            content = _read(root, args.name, args.path)
            sys.stdout.write(content)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


def workflow_help(name: str) -> str:
    """Point help readers to an available, version-matched workflow guide."""
    try:
        available = _root().joinpath(name, "SKILL.md").is_file()
    except (ValueError, OSError):
        return ""
    if not available:
        return ""
    return f"Workflow guide (read first): dpdata skills read {name}"
