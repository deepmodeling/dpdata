from __future__ import annotations

import bz2
import datetime
import gzip
import importlib
import json
from enum import Enum
from pathlib import Path
from typing import Any, BinaryIO, TextIO, cast
from uuid import UUID

import numpy as np


def _detect_format(filename: str | Path, fmt: str | None = None) -> str:
    if fmt is not None:
        return fmt
    suffixes = [suffix.lower() for suffix in Path(filename).suffixes]
    if suffixes and suffixes[-1] in {".gz", ".z", ".bz2"}:
        suffixes.pop()
    suffix = suffixes[-1] if suffixes else ""
    if suffix == ".mpk":
        return "mpk"
    if suffix in {".yaml", ".yml"}:
        return "yaml"
    return "json"


def _open_text(filename: str | Path, mode: str) -> TextIO:
    path = str(filename)
    lower_path = path.lower()
    if lower_path.endswith((".gz", ".z")):
        return cast("TextIO", gzip.open(path, mode, encoding="utf-8"))
    if lower_path.endswith(".bz2"):
        return cast("TextIO", bz2.open(path, mode, encoding="utf-8"))
    return cast("TextIO", open(path, mode, encoding="utf-8"))


def _open_binary(filename: str | Path, mode: str) -> BinaryIO:
    path = str(filename)
    lower_path = path.lower()
    if lower_path.endswith((".gz", ".z")):
        return cast("BinaryIO", gzip.open(path, mode))
    if lower_path.endswith(".bz2"):
        return cast("BinaryIO", bz2.open(path, mode))
    return cast("BinaryIO", open(path, mode))


def _yaml_dump(obj: Any, fp, *args: Any, **kwargs: Any) -> None:
    try:
        yaml = importlib.import_module("yaml")

        if "sort_keys" not in kwargs:
            kwargs["sort_keys"] = False
        getattr(yaml, "safe_dump")(obj, fp, *args, **kwargs)
    except ModuleNotFoundError:
        try:
            ruamel_yaml = importlib.import_module("ruamel.yaml")
        except ModuleNotFoundError as e:
            raise RuntimeError(
                "Dumping YAML files requires PyYAML or ruamel.yaml."
            ) from e
        yaml = getattr(ruamel_yaml, "YAML")()
        if "indent" in kwargs:
            indent = kwargs.pop("indent")
            yaml.indent(mapping=indent, sequence=indent, offset=2)
        yaml.dump(obj, fp, *args, **kwargs)


def _yaml_load(fp, *args: Any, **kwargs: Any) -> Any:
    try:
        yaml = importlib.import_module("yaml")

        return getattr(yaml, "safe_load")(fp, *args, **kwargs)
    except ModuleNotFoundError:
        try:
            ruamel_yaml = importlib.import_module("ruamel.yaml")
        except ModuleNotFoundError as e:
            raise RuntimeError(
                "Loading YAML files requires PyYAML or ruamel.yaml."
            ) from e
        yaml = getattr(ruamel_yaml, "YAML")(typ="safe")
        return yaml.load(fp, *args, **kwargs)


def _encode_ndarray(obj: np.ndarray) -> dict[str, Any]:
    if str(obj.dtype).startswith("complex"):
        data = [obj.real.tolist(), obj.imag.tolist()]
    else:
        data = obj.tolist()
    return {
        "@module": "numpy",
        "@class": "array",
        "dtype": str(obj.dtype),
        "data": data,
    }


def to_serializable(obj: Any) -> Any:
    """Convert common dpdata objects to monty-compatible plain data."""
    if isinstance(obj, dict):
        return {to_serializable(k): to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_serializable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _encode_ndarray(obj)
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, datetime.datetime):
        return {
            "@module": "datetime",
            "@class": "datetime",
            "string": str(obj),
        }
    if isinstance(obj, UUID):
        return {"@module": "uuid", "@class": "UUID", "string": str(obj)}
    if isinstance(obj, Path):
        return {"@module": "pathlib", "@class": "Path", "string": str(obj)}
    if isinstance(obj, Enum):
        return {
            "@module": obj.__class__.__module__,
            "@class": obj.__class__.__name__,
            "value": to_serializable(obj.value),
        }
    if hasattr(obj, "as_dict"):
        data = obj.as_dict()
        if "@module" not in data:
            data["@module"] = obj.__class__.__module__
        if "@class" not in data:
            data["@class"] = obj.__class__.__name__
        return to_serializable(data)
    return obj


def _decode_ndarray(data: dict[str, Any]) -> np.ndarray:
    dtype = data["dtype"]
    if dtype.startswith("complex"):
        real, imag = data["data"]
        return np.array(real, dtype=dtype) + np.array(imag, dtype=dtype) * 1j
    return np.array(data["data"], dtype=dtype)


def process_decoded(obj: Any) -> Any:
    """Decode monty-style dictionaries used by existing dpdata JSON files."""
    if isinstance(obj, dict):
        if "@module" in obj and "@class" in obj:
            module_name = obj["@module"]
            class_name = obj["@class"]
            if module_name == "numpy" and class_name == "array":
                return _decode_ndarray(obj)
            if module_name == "datetime" and class_name == "datetime":
                try:
                    return datetime.datetime.fromisoformat(obj["string"])
                except ValueError:
                    value = obj["string"].split("+")[0]
                    try:
                        return datetime.datetime.strptime(value, "%Y-%m-%d %H:%M:%S.%f")
                    except ValueError:
                        return datetime.datetime.strptime(value, "%Y-%m-%d %H:%M:%S")
            if module_name == "uuid" and class_name == "UUID":
                return UUID(obj["string"])
            if module_name == "pathlib" and class_name == "Path":
                return Path(obj["string"])
            try:
                module = importlib.import_module(module_name)
                cls = getattr(module, class_name)
            except (AttributeError, ImportError, ModuleNotFoundError):
                cls = None
            if cls is not None:
                # Decode fields before construction so custom objects receive the
                # same nested Python values as top-level serialized objects.
                data = {
                    k: process_decoded(v)
                    for k, v in obj.items()
                    if not k.startswith("@")
                }
                if hasattr(cls, "from_dict"):
                    return cls.from_dict(data)
                if isinstance(cls, type) and issubclass(cls, Enum):
                    return cls(data["value"])
        return {process_decoded(k): process_decoded(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [process_decoded(v) for v in obj]
    return obj


def dumpfn(
    obj: Any,
    filename: str | Path,
    *args: Any,
    fmt: str | None = None,
    **kwargs: Any,
) -> None:
    """Dump an object to JSON, YAML, or msgpack without requiring monty."""
    fmt = _detect_format(filename, fmt)
    obj = to_serializable(obj)
    if fmt == "json":
        with _open_text(filename, "wt") as fp:
            json.dump(obj, fp, *args, **kwargs)
        return
    if fmt == "yaml":
        with _open_text(filename, "wt") as fp:
            _yaml_dump(obj, fp, *args, **kwargs)
        return
    if fmt == "mpk":
        try:
            import msgpack
        except ModuleNotFoundError as e:
            raise RuntimeError("Dumping msgpack files requires msgpack.") from e
        with _open_binary(filename, "wb") as fp:
            msgpack.dump(obj, fp, *args, **kwargs)
        return
    raise TypeError(f"Invalid format: {fmt}")


def loadfn(
    filename: str | Path,
    *args: Any,
    fmt: str | None = None,
    **kwargs: Any,
) -> Any:
    """Load JSON, YAML, or msgpack data and decode monty-style objects."""
    fmt = _detect_format(filename, fmt)
    if fmt == "json":
        with _open_text(filename, "rt") as fp:
            return process_decoded(json.load(fp, *args, **kwargs))
    if fmt == "yaml":
        with _open_text(filename, "rt") as fp:
            return process_decoded(_yaml_load(fp, *args, **kwargs))
    if fmt == "mpk":
        try:
            import msgpack
        except ModuleNotFoundError as e:
            raise RuntimeError("Loading msgpack files requires msgpack.") from e
        with _open_binary(filename, "rb") as fp:
            return process_decoded(msgpack.load(fp, *args, **kwargs))
    raise TypeError(f"Invalid format: {fmt}")
