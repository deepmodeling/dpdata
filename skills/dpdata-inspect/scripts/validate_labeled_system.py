#!/usr/bin/env python3
"""Validate dpdata systems against the reusable labeled-data contract.

The validator intentionally checks structure and finite values only. Unit and
stress conventions are producer contracts, so callers must state them with the
metadata options when running in strict mode.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

import numpy as np


@dataclass
class ValidationReport:
    """Machine-readable validation result for one or more systems."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    summary: dict[str, Any] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        """Whether no contract violations were found."""
        return not self.errors

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation."""
        return {"ok": self.ok, **asdict(self)}


def _as_array(
    data: Mapping[str, Any], key: str, errors: list[str]
) -> np.ndarray | None:
    """Get *key* as an array, recording malformed values instead of raising."""
    if key not in data:
        return None
    try:
        return np.asarray(data[key])
    except (TypeError, ValueError) as exc:
        errors.append(f"{key}: cannot convert to an array ({exc})")
        return None


def _check_finite(name: str, value: np.ndarray, errors: list[str]) -> None:
    """Reject NaN and infinity in a numeric field."""
    if not np.issubdtype(value.dtype, np.number):
        errors.append(f"{name}: expected numeric values, got {value.dtype}")
        return
    if not np.all(np.isfinite(value)):
        errors.append(f"{name}: contains NaN or infinity")


def _check_shape(
    name: str,
    value: np.ndarray,
    shape: tuple[int, ...],
    errors: list[str],
) -> None:
    """Check an exact array shape."""
    if value.shape != shape:
        errors.append(f"{name}: shape {value.shape}, expected {shape}")


def validate_data(
    data: Mapping[str, Any],
    *,
    require_labels: bool = True,
    require_pbc: bool = False,
    expected_type_map: Sequence[str] | None = None,
    strict_metadata: bool = False,
    energy_unit: str | None = None,
    force_unit: str | None = None,
    virial_unit: str | None = None,
    stress_sign: str | None = None,
    provenance: str | None = None,
    producer: str | None = None,
    producer_version: str | None = None,
    parser_version: str | None = None,
) -> ValidationReport:
    """Validate one dpdata ``System.data`` mapping.

    Parameters are deliberately explicit because dpdata cannot infer units,
    stress/virial sign conventions, or the producer/version contract from the
    arrays alone.  Set ``strict_metadata`` to require those declarations.
    """
    report = ValidationReport()
    errors = report.errors
    warnings = report.warnings

    missing = [
        key
        for key in ("coords", "atom_types", "atom_names", "atom_numbs")
        if key not in data
    ]
    if missing:
        errors.append("missing required fields: " + ", ".join(missing))
        return report

    coords = _as_array(data, "coords", errors)
    atom_types = _as_array(data, "atom_types", errors)
    cells = _as_array(data, "cells", errors)
    if coords is None or atom_types is None:
        return report

    if coords.ndim != 3 or coords.shape[-1] != 3:
        errors.append(f"coords: shape {coords.shape}, expected (nframes, natoms, 3)")
        return report
    nframes, natoms = coords.shape[:2]
    _check_finite("coords", coords, errors)

    if atom_types.shape != (natoms,):
        errors.append(f"atom_types: shape {atom_types.shape}, expected ({natoms},)")
    elif not np.issubdtype(atom_types.dtype, np.integer):
        errors.append(f"atom_types: expected integer values, got {atom_types.dtype}")
    elif np.any(atom_types < 0):
        errors.append("atom_types: contains a negative type index")

    atom_names = data["atom_names"]
    atom_numbs = data["atom_numbs"]
    try:
        ntypes = len(atom_names)
        counts = np.asarray(atom_numbs)
    except (TypeError, ValueError) as exc:
        errors.append(f"atom_names/atom_numbs: malformed ({exc})")
        return report
    if ntypes != len(counts):
        errors.append(
            f"atom_names/atom_numbs: lengths differ ({ntypes} versus {len(counts)})"
        )
    if counts.ndim != 1 or not np.issubdtype(counts.dtype, np.integer):
        errors.append("atom_numbs: expected a one-dimensional integer sequence")
    elif np.any(counts < 0):
        errors.append("atom_numbs: contains a negative atom count")
    elif int(np.sum(counts)) != natoms:
        errors.append(
            f"atom_numbs: sum {int(np.sum(counts))}, expected natoms {natoms}"
        )
    if atom_types.size and not ntypes:
        errors.append("atom_names: empty type_map for non-empty atom_types")
    elif (
        atom_types.size
        and np.issubdtype(atom_types.dtype, np.integer)
        and np.max(atom_types) >= ntypes
    ):
        errors.append("atom_types: index is outside atom_names/type_map")
    if expected_type_map is not None and list(atom_names) != list(expected_type_map):
        errors.append(
            "atom_names: does not match type_map "
            f"({list(atom_names)!r} versus {list(expected_type_map)!r})"
        )

    nopbc = bool(data.get("nopbc", False))
    if require_pbc and nopbc:
        errors.append("nopbc: system is non-periodic while periodic cells are required")
    if cells is None:
        if not nopbc or require_pbc:
            errors.append("cells: missing for a periodic system")
    else:
        _check_shape("cells", cells, (nframes, 3, 3), errors)
        if cells.shape == (nframes, 3, 3):
            _check_finite("cells", cells, errors)
            if require_pbc:
                try:
                    determinants = np.linalg.det(cells)
                except (np.linalg.LinAlgError, TypeError, ValueError) as exc:
                    errors.append(f"cells: cannot compute determinants ({exc})")
                else:
                    if np.any(np.isclose(determinants, 0.0)):
                        errors.append(
                            "cells: singular cell while periodic cells are required"
                        )

    label_shapes = {
        "energies": (nframes,),
        "forces": (nframes, natoms, 3),
        "virials": (nframes, 3, 3),
        "atom_pref": (nframes, natoms),
    }
    if require_labels and "energies" not in data:
        errors.append("energies: missing required energy labels")
    for name, expected_shape in label_shapes.items():
        if name not in data:
            continue
        value = _as_array(data, name, errors)
        if value is None:
            continue
        _check_shape(name, value, expected_shape, errors)
        _check_finite(name, value, errors)

    if (require_labels or "energies" in data) and energy_unit is None:
        warnings.append("energy unit was not declared")
    if "forces" in data and force_unit is None:
        warnings.append("force unit was not declared")
    if "virials" in data:
        if virial_unit is None:
            warnings.append("virial/stress unit was not declared")
        if stress_sign is None:
            warnings.append("virial/stress sign convention was not declared")
    if provenance is None:
        warnings.append("raw-data provenance was not supplied")
    if producer is None:
        warnings.append("producer contract was not declared")
    if producer_version is None:
        warnings.append("producer version was not declared")
    if parser_version is None:
        warnings.append("parser version was not declared")
    if provenance is not None:
        provenance_path = Path(provenance)
        if not provenance_path.is_file():
            errors.append(f"provenance: file does not exist ({provenance})")
        elif provenance_path.stat().st_size == 0:
            errors.append(f"provenance: file is empty ({provenance})")

    if strict_metadata:
        metadata_missing = []
        if (require_labels or "energies" in data) and energy_unit is None:
            metadata_missing.append("energy unit")
        if "forces" in data and force_unit is None:
            metadata_missing.append("force unit")
        if "virials" in data and virial_unit is None:
            metadata_missing.append("virial/stress unit")
        if "virials" in data and stress_sign is None:
            metadata_missing.append("virial/stress sign")
        if provenance is None:
            metadata_missing.append("provenance")
        if producer is None:
            metadata_missing.append("producer")
        if producer_version is None:
            metadata_missing.append("producer version")
        if parser_version is None:
            metadata_missing.append("parser version")
        if metadata_missing:
            errors.append(
                "strict metadata is enabled; missing " + ", ".join(metadata_missing)
            )

    report.summary = {
        "frames": nframes,
        "atoms": natoms,
        "types": ntypes,
        "labeled": "energies" in data,
        "forces": "forces" in data,
        "virials": "virials" in data,
        "periodic": not nopbc,
        "atom_names": list(atom_names),
        "producer": producer,
        "producer_version": producer_version,
        "parser_version": parser_version,
    }
    return report


def validate_systems(
    systems: Sequence[Any],
    **kwargs: Any,
) -> list[ValidationReport]:
    """Validate systems and enforce one shared type map across a collection."""
    reports = [validate_data(system.data, **kwargs) for system in systems]
    names: list[str] | None = None
    for index, (system, report) in enumerate(zip(systems, reports, strict=True)):
        current = list(system.data.get("atom_names", []))
        if names is None:
            names = current
        elif current != names:
            report.errors.append(
                "atom_names: differs from system 0 in a multi-system dataset "
                f"(system {index}: {current!r}, system 0: {names!r})"
            )
    return reports


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate dpdata frame, label, atom-order, and metadata contracts."
    )
    parser.add_argument("input", help="Input path accepted by dpdata")
    parser.add_argument(
        "-f",
        "--format",
        default="auto",
        help="dpdata format (default: auto)",
    )
    parser.add_argument(
        "--multi",
        action="store_true",
        help="Load a directory containing multiple systems",
    )
    parser.add_argument(
        "--allow-unlabeled",
        action="store_true",
        help="Allow a dataset without energy labels",
    )
    parser.add_argument(
        "--require-pbc",
        action="store_true",
        help="Require nonsingular periodic cells for every frame",
    )
    parser.add_argument(
        "--type-map",
        nargs="+",
        help="Expected atom names in type-map order (also passed to the loader)",
    )
    parser.add_argument("--energy-unit", help="Declared energy unit, e.g. eV")
    parser.add_argument("--force-unit", help="Declared force unit, e.g. eV/A")
    parser.add_argument("--virial-unit", help="Declared virial/stress unit")
    parser.add_argument("--stress-sign", help="Declared virial/stress sign convention")
    parser.add_argument(
        "--provenance",
        help="Path to a non-empty raw input or provenance record",
    )
    parser.add_argument("--producer", help="Producer name, e.g. CP2K or VASP")
    parser.add_argument("--producer-version", help="Producer version or release")
    parser.add_argument("--parser-version", help="dpdata/parser version used")
    parser.add_argument(
        "--strict-metadata",
        action="store_true",
        help="Turn undeclared units, stress sign, and provenance into errors",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of text")
    return parser.parse_args(argv)


def _load_systems(args: argparse.Namespace) -> list[Any]:
    # Import lazily so validate_data can be used in a minimal Python environment.
    from dpdata.system import LabeledSystem, MultiSystems, System

    loader_kwargs = {"type_map": args.type_map}
    if args.multi:
        collection = MultiSystems.from_file(
            args.input,
            fmt=args.format,
            labeled=not args.allow_unlabeled,
            **loader_kwargs,
        )
        return list(collection.systems.values())
    cls = System if args.allow_unlabeled else LabeledSystem
    return [cls(args.input, fmt=args.format, **loader_kwargs)]


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point."""
    args = _parse_args(argv)
    try:
        systems = _load_systems(args)
        if not systems:
            raise ValueError("no systems found")
        reports = validate_systems(
            systems,
            require_labels=not args.allow_unlabeled,
            require_pbc=args.require_pbc,
            expected_type_map=args.type_map,
            strict_metadata=args.strict_metadata,
            energy_unit=args.energy_unit,
            force_unit=args.force_unit,
            virial_unit=args.virial_unit,
            stress_sign=args.stress_sign,
            provenance=args.provenance,
            producer=args.producer,
            producer_version=args.producer_version,
            parser_version=args.parser_version,
        )
    except Exception as exc:  # noqa: BLE001 - expose parser failures as CLI errors
        report = ValidationReport(errors=[f"load: {exc}"])
        reports = [report]

    result = {
        "ok": all(report.ok for report in reports),
        "systems": [report.to_dict() for report in reports],
    }
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        status = "PASS" if result["ok"] else "FAIL"
        print(f"{status}: validated {len(reports)} system(s)")
        for index, report in enumerate(reports):
            prefix = f"system {index}: "
            for message in report.errors:
                print(f"ERROR {prefix}{message}")
            for message in report.warnings:
                print(f"WARNING {prefix}{message}")
            if report.summary:
                summary = report.summary
                print(
                    f"  frames={summary['frames']} atoms={summary['atoms']} "
                    f"labeled={summary['labeled']} periodic={summary['periodic']}"
                )
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
