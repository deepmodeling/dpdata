"""LMDB format that is fully interoperable with DeePMD-kit.

The on-disk layout is a *flat* sequence of frames, identical to what
DeePMD-kit's ``LmdbDataReader`` consumes and what the community
``npy_to_lmdb`` / ``json_to_lmdb`` converters produce:

- ``__metadata__`` (msgpack dict, string keys)::

      {
          "nframes": int,
          "frame_idx_fmt": "012d",
          "type_map": [str, ...],          # global element table
          "frame_nlocs": [int, ...],       # atoms per frame
          "frame_system_ids": [int, ...],  # source-system index per frame
      }

- one entry per frame, keyed by the zero-padded global frame index
  (``"000000000000"`` by default), whose value is a msgpack dict with
  string keys. Array values use the manual encoding
  ``{"type": str(dtype), "shape": [...], "data": <bytes>}`` and
  ``atom_numbs`` is a plain list of per-type counts over ``type_map``.

  Core per-frame keys use the plural names consumed by DeePMD-kit
  (``coords``, ``cells``, ``energies``, ``forces``, ``virials``,
  ``atom_types``). Registered additional fields use
  :attr:`dpdata.data_type.DataType.deepmd_name`.

``atom_types`` is stored as ``int32`` global indices into ``type_map``.

Notes
-----
This reader loads frames into memory, mirroring ``dpdata``'s in-memory
``System`` model. A configurable frame-count guard prevents accidental
decoding of very large databases; training datasets should normally be
consumed directly by DeePMD-kit's streaming dataloader.
"""

from __future__ import annotations

import itertools
import numbers
import os
import shutil
import tempfile
import threading
import warnings
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import lmdb
import msgpack
import numpy as np

from dpdata.data_type import Axis, DataType
from dpdata.format import Format

DEFAULT_MAP_SIZE = 1024**4  # 1 TiB sparse, matches the reference converters
DEFAULT_FRAME_FMT = "012d"
DEFAULT_MAX_FRAMES = 100_000
DEFAULT_WRITE_BATCH_SIZE = 1_000
METADATA_KEY = b"__metadata__"

_RESERVED_KEYS = frozenset(
    {
        "atom_numbs",
        "atom_names",
        "atom_types",
        "orig",
        "real_atom_types",
        "real_atom_names",
        "nopbc",
    }
)
_CORE_DISK_NAMES = {
    "cells": "cells",
    "coords": "coords",
    "energies": "energies",
    "forces": "forces",
    "virials": "virials",
}
_CORE_FRAME_KEYS = frozenset(_CORE_DISK_NAMES.values())
_DEEPMD_CORE_REMAP = {
    "coords": "coord",
    "cells": "box",
    "energies": "energy",
    "forces": "force",
    "atom_types": "atype",
    "virials": "virial",
}
_FORBIDDEN_CUSTOM_DISK_NAMES = frozenset(
    {
        METADATA_KEY.decode(),
        "atom_types",
        "atom_numbs",
        "atom_names",
        "orig",
        "type",
        "natoms",
        "real_natoms_vec",
        "fid",
        "sid",
        *_CORE_FRAME_KEYS,
        *_DEEPMD_CORE_REMAP.values(),
    }
)

_TOKEN_TO_AXIS = {
    "nframes": Axis.NFRAMES,
    "natoms": Axis.NATOMS,
    "ntypes": Axis.NTYPES,
    "nbonds": Axis.NBONDS,
}


@dataclass(frozen=True)
class _FieldSpec:
    """Describe one dpdata field and its LMDB representation."""

    dp_name: str
    disk_name: str
    shape: tuple[int | Axis, ...]
    deepmd_name: str

    @property
    def frame_axis(self) -> int | None:
        """Return the frame axis in the in-memory array."""
        return (
            self.shape.index(Axis.NFRAMES)
            if Axis.NFRAMES in self.shape
            else None
        )

    @property
    def atom_axes(self) -> tuple[int, ...]:
        """Return atom axes in the in-memory array."""
        return tuple(i for i, dim in enumerate(self.shape) if dim is Axis.NATOMS)

    @property
    def payload_atom_axes(self) -> tuple[int, ...]:
        """Return atom axes after the frame axis has been removed."""
        frame_axis = self.frame_axis
        if frame_axis is None:
            return self.atom_axes
        return tuple(axis - (axis > frame_axis) for axis in self.atom_axes)


@dataclass(frozen=True)
class _AtomSource:
    """Validated atom-type data for one source system."""

    names: tuple[str, ...]
    types: np.ndarray
    masks: np.ndarray | None


_PROTOCOL_DTYPES = (
    DataType(
        "spins",
        np.ndarray,
        (Axis.NFRAMES, Axis.NATOMS, 3),
        required=False,
        deepmd_name="spin",
    ),
    DataType(
        "force_mags",
        np.ndarray,
        (Axis.NFRAMES, Axis.NATOMS, 3),
        required=False,
        deepmd_name="force_mag",
    ),
    DataType(
        "charges",
        np.ndarray,
        (Axis.NFRAMES, Axis.NATOMS),
        required=False,
        deepmd_name="charge",
    ),
    DataType(
        "move",
        np.ndarray,
        (Axis.NFRAMES, Axis.NATOMS, 3),
        required=False,
        deepmd_name="move",
    ),
    DataType(
        "hessian",
        np.ndarray,
        (Axis.NFRAMES, Axis.NATOMS, 3, Axis.NATOMS, 3),
        required=False,
        deepmd_name="hessian",
    ),
    DataType(
        "fparam",
        np.ndarray,
        (Axis.NFRAMES, -1),
        required=False,
        deepmd_name="fparam",
    ),
    DataType(
        "aparam",
        np.ndarray,
        (Axis.NFRAMES, Axis.NATOMS, -1),
        required=False,
        deepmd_name="aparam",
    ),
    DataType(
        "atom_ener",
        np.ndarray,
        (Axis.NFRAMES, Axis.NATOMS),
        required=False,
        deepmd_name="atom_ener",
    ),
    DataType(
        "charge_spin",
        np.ndarray,
        (Axis.NFRAMES, -1),
        required=False,
        deepmd_name="charge_spin",
    ),
)
_PROTOCOL_FIELD_NAMES = frozenset(
    {
        *_CORE_DISK_NAMES,
        "atom_pref",
        *(dtype.name for dtype in _PROTOCOL_DTYPES),
    }
)
_PROTOCOL_DISK_NAMES = {
    dtype.deepmd_name: dtype.name for dtype in _PROTOCOL_DTYPES
}

_READ_ENV_CACHE: dict[str, tuple[lmdb.Environment, int]] = {}
_READ_ENV_LOCK = threading.Lock()
_PUBLISHING_PATHS: set[str] = set()
_READ_ENV_PID = os.getpid()


def _shape_tokens(shape: tuple[int | Axis, ...]) -> list[int | str]:
    """Serialize a (possibly symbolic) DataType shape to msgpack-able tokens."""
    tokens: list[int | str] = []
    for s in shape:
        if isinstance(s, Axis):
            tokens.append(s.value)
        else:
            tokens.append(int(s))
    return tokens


def _tokens_to_shape(tokens: list[int | str | bytes]) -> tuple[int | Axis, ...]:
    """Inverse of :func:`_shape_tokens`."""
    out: list[int | Axis] = []
    for tok in tokens:
        t = tok.decode() if isinstance(tok, bytes) else tok
        if isinstance(t, str) and t in _TOKEN_TO_AXIS:
            out.append(_TOKEN_TO_AXIS[t])
        else:
            out.append(int(t))
    return tuple(out)


class LMDBError(Exception):
    """Base class for LMDB errors."""


class LMDBMetadataError(LMDBError):
    """Metadata not found or malformed in LMDB."""


class LMDBFrameError(LMDBError):
    """Frame data not found in LMDB."""


def _encode_array(arr: np.ndarray) -> dict[str, Any]:
    """Encode a numpy array as ``{"type", "shape", "data"}`` (string keys).

    ``shape`` is taken from the original array (so 0-d scalars stay 0-d);
    ``data`` is the C-ordered byte buffer.
    """
    a = np.asarray(arr)
    if a.dtype.hasobject or a.dtype.fields is not None:
        raise LMDBError(
            f"LMDB array dtype {a.dtype} is not a portable raw-byte dtype."
        )
    try:
        restored_dtype = np.dtype(str(a.dtype))
    except TypeError as exc:
        raise LMDBError(
            f"LMDB array dtype {a.dtype} cannot be reconstructed from its "
            "serialized name."
        ) from exc
    if restored_dtype != a.dtype:
        raise LMDBError(
            f"LMDB array dtype {a.dtype} is not preserved by its serialized "
            f"name {str(a.dtype)!r}."
        )
    return {
        "type": str(a.dtype),
        "shape": list(a.shape),
        "data": np.ascontiguousarray(a).tobytes(),
    }


def _is_encoded_array(val) -> bool:
    """Whether ``val`` is a manually encoded array dict (str or byte keys)."""
    if not isinstance(val, dict):
        return False
    return ("type" in val and "data" in val) or (b"type" in val and b"data" in val)


def _decode_array(d: dict[str | bytes, Any]) -> np.ndarray:
    """Reconstruct a numpy array from the manual encoding.

    Tolerates both string and byte keys so files written by any compliant
    producer can be read back.
    """
    type_key = "type" if "type" in d else b"type"
    shape_key = "shape" if "shape" in d else b"shape"
    data_key = "data" if "data" in d else b"data"
    raw_type = d[type_key]
    dtype = np.dtype(raw_type.decode() if isinstance(raw_type, bytes) else raw_type)
    if dtype.hasobject or dtype.fields is not None:
        raise ValueError(f"dtype {dtype} is not a portable raw-byte dtype")
    buf = d[data_key]
    if shape_key in d:
        shape = tuple(d[shape_key])
    else:
        shape = (len(buf) // dtype.itemsize,)
    # frombuffer returns a read-only view; copy so downstream code can mutate.
    return np.frombuffer(buf, dtype=dtype).reshape(shape).copy()


def _decode_frame(raw: bytes) -> dict[str, Any]:
    """Decode a msgpack frame into a dict of numpy arrays / python scalars."""
    try:
        frame = msgpack.unpackb(raw, raw=False)
    except (
        msgpack.exceptions.ExtraData,
        msgpack.exceptions.FormatError,
        msgpack.exceptions.StackError,
        ValueError,
    ) as exc:
        raise LMDBFrameError(f"Cannot decode frame: {exc}") from exc
    if not isinstance(frame, dict):
        raise LMDBFrameError("Each frame must contain a mapping.")
    out = {}
    for key, val in frame.items():
        if isinstance(key, bytes):
            name = key.decode()
        elif isinstance(key, str):
            name = key
        else:
            raise LMDBFrameError("Frame keys must be strings or bytes.")
        try:
            out[name] = _decode_array(val) if _is_encoded_array(val) else val
        except (TypeError, ValueError) as exc:
            raise LMDBFrameError(
                f"Cannot decode array field '{name}': {exc}"
            ) from exc
    return out


def _meta_get(meta: dict[str | bytes, Any], key: str) -> Any:
    """Read a metadata value, tolerating string or byte keys."""
    if key in meta:
        return meta[key]
    bkey = key.encode()
    if bkey in meta:
        return meta[bkey]
    return None


def _read_metadata(txn: lmdb.Transaction) -> dict[str | bytes, Any]:
    raw = txn.get(METADATA_KEY)
    if raw is None:
        raise LMDBMetadataError("LMDB database does not contain __metadata__.")
    try:
        metadata = msgpack.unpackb(raw, raw=False)
    except (
        msgpack.exceptions.ExtraData,
        msgpack.exceptions.FormatError,
        msgpack.exceptions.StackError,
        ValueError,
    ) as exc:
        raise LMDBMetadataError(f"Cannot decode __metadata__: {exc}") from exc
    if not isinstance(metadata, dict):
        raise LMDBMetadataError("__metadata__ must contain a mapping.")
    return metadata


def _packb(value: Any) -> bytes:
    """Pack a value and enforce the byte-oriented LMDB contract."""
    packed = msgpack.packb(value, use_bin_type=True)
    if not isinstance(packed, bytes):
        raise TypeError("msgpack.packb did not return bytes.")
    return packed


def _require_integer(
    value: Any, label: str, *, minimum: int
) -> int:
    """Return a strictly validated integer metadata value."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, numbers.Integral
    ):
        raise LMDBMetadataError(f"{label} must be an integer.")
    result = int(value)
    if result < minimum:
        raise LMDBMetadataError(f"{label} must be at least {minimum}.")
    return result


def _all_data_types() -> dict[str, DataType]:
    """Return registered and protocol-level data types by dpdata name."""
    from dpdata.system import LabeledSystem, System

    dtypes = {dt.name: dt for dt in _PROTOCOL_DTYPES}
    for dt in (*System.DTYPES, *LabeledSystem.DTYPES):
        dtypes[dt.name] = dt
    return dtypes


def _disk_name(dtype: DataType) -> str:
    """Return the LMDB key associated with a dpdata data type."""
    return _CORE_DISK_NAMES.get(dtype.name, dtype.deepmd_name)


def _validate_field_namespace(dp_name: str, disk_name: str) -> None:
    """Prevent custom fields from shadowing DeePMD structural data."""
    if dp_name in _CORE_DISK_NAMES:
        expected = _CORE_DISK_NAMES[dp_name]
        if disk_name != expected:
            raise LMDBError(
                f"Core field '{dp_name}' must use LMDB key '{expected}', not "
                f"'{disk_name}'."
            )
        return
    protocol_owner = _PROTOCOL_DISK_NAMES.get(disk_name)
    if protocol_owner is not None and protocol_owner != dp_name:
        raise LMDBError(
            f"Custom field '{dp_name}' cannot use protocol LMDB key "
            f"'{disk_name}', which belongs to '{protocol_owner}'."
        )
    if disk_name in _FORBIDDEN_CUSTOM_DISK_NAMES or disk_name.startswith("find_"):
        raise LMDBError(
            f"Custom field '{dp_name}' cannot use reserved LMDB key "
            f"'{disk_name}'."
        )


def _field_spec(
    dtype: DataType, *, validate_namespace: bool = True
) -> _FieldSpec:
    """Build a field specification from a registered data type."""
    if dtype.shape is None:
        raise LMDBError(
            f"Data type '{dtype.name}' has no declared shape and cannot be "
            "represented unambiguously in per-frame LMDB records."
        )
    disk_name = _disk_name(dtype)
    if validate_namespace:
        _validate_field_namespace(dtype.name, disk_name)
    return _FieldSpec(
        dp_name=dtype.name,
        disk_name=disk_name,
        shape=tuple(dtype.shape),
        deepmd_name=dtype.deepmd_name,
    )


def _field_registries(
    *, prefer_protocol: bool = False
) -> tuple[
    dict[str, _FieldSpec], dict[str, _FieldSpec]
]:
    """Return field specifications indexed by dpdata and disk names."""
    data_types = _all_data_types()
    if prefer_protocol:
        for dtype in _PROTOCOL_DTYPES:
            data_types[dtype.name] = dtype
    by_dp: dict[str, _FieldSpec] = {}
    by_disk: dict[str, _FieldSpec] = {}
    ambiguous_disk_names: set[str] = set()
    for dtype in data_types.values():
        if dtype.name in _RESERVED_KEYS or dtype.shape is None:
            continue
        spec = _field_spec(dtype, validate_namespace=False)
        by_dp[spec.dp_name] = spec
        if spec.disk_name in ambiguous_disk_names:
            continue
        previous = by_disk.get(spec.disk_name)
        if previous is not None and previous.dp_name != spec.dp_name:
            del by_disk[spec.disk_name]
            ambiguous_disk_names.add(spec.disk_name)
            continue
        by_disk[spec.disk_name] = spec
    return by_dp, by_disk


def _validate_names(names: list[str] | tuple[str, ...], label: str) -> tuple[str, ...]:
    """Validate and normalize an element table."""
    if isinstance(names, (str, bytes)):
        raise LMDBError(f"{label} must be a sequence of element names.")
    normalized = tuple(str(name) for name in names)
    if not normalized:
        raise LMDBError(f"{label} must contain at least one element.")
    if len(set(normalized)) != len(normalized):
        raise LMDBError(f"{label} contains duplicate element names: {normalized}.")
    return normalized


def _validate_integer_types(
    values: np.ndarray,
    *,
    name: str,
    upper_bound: int,
    allow_virtual: bool,
) -> np.ndarray:
    """Validate an atom-type array without silently truncating values."""
    array = np.asarray(values)
    if (
        np.issubdtype(array.dtype, np.bool_)
        or not np.issubdtype(array.dtype, np.integer)
    ):
        raise LMDBError(f"{name} must use an integer dtype, got {array.dtype}.")
    lower_bound = -1 if allow_virtual else 0
    if np.issubdtype(array.dtype, np.unsignedinteger):
        out_of_range = np.any(array >= upper_bound)
    else:
        out_of_range = np.any(array < lower_bound) or np.any(array >= upper_bound)
    if out_of_range:
        constraint = (
            f"-1 or indices in [0, {upper_bound})"
            if allow_virtual
            else f"indices in [0, {upper_bound})"
        )
        raise LMDBError(f"{name} must contain only {constraint}.")
    return array.astype(np.int64, copy=False)


def _prepare_atom_source(
    data: dict[str, Any], nframes: int, natoms: int
) -> _AtomSource:
    """Validate standard or mixed atom types for one source system."""
    if "real_atom_types" in data or "real_atom_names" in data:
        if "real_atom_types" not in data or "real_atom_names" not in data:
            raise LMDBError(
                "Mixed-type data requires both real_atom_types and real_atom_names."
            )
        names = _validate_names(data["real_atom_names"], "real_atom_names")
        real_types = _validate_integer_types(
            np.asarray(data["real_atom_types"]),
            name="real_atom_types",
            upper_bound=len(names),
            allow_virtual=True,
        )
        if real_types.shape != (nframes, natoms):
            raise LMDBError(
                "real_atom_types must have shape "
                f"({nframes}, {natoms}), got {real_types.shape}."
            )
        masks = real_types >= 0
        if np.any(np.sum(masks, axis=1) == 0):
            raise LMDBError("Every frame must contain at least one real atom.")
        return _AtomSource(names, real_types, masks)

    names = _validate_names(data["atom_names"], "atom_names")
    atom_types = _validate_integer_types(
        np.asarray(data["atom_types"]),
        name="atom_types",
        upper_bound=len(names),
        allow_virtual=False,
    )
    if atom_types.shape != (natoms,):
        raise LMDBError(
            f"atom_types must have shape ({natoms},), got {atom_types.shape}."
        )
    expected_counts = np.bincount(atom_types, minlength=len(names)).tolist()
    if list(data["atom_numbs"]) != expected_counts:
        raise LMDBError(
            "atom_numbs is inconsistent with atom_types: "
            f"expected {expected_counts}, got {data['atom_numbs']}."
        )
    return _AtomSource(names, atom_types, None)


def _source_data(system: Any) -> dict[str, Any]:
    """Return a raw data dictionary from a System-like object."""
    if isinstance(system, dict):
        return system
    if hasattr(system, "data") and isinstance(system.data, dict):
        return system.data
    raise TypeError("Each input must be a dpdata System or a data dictionary.")


def _source_active_names(data: dict[str, Any]) -> tuple[str, ...]:
    """Return active elements in first-appearance order for type-map discovery."""
    coords = np.asarray(data["coords"])
    if coords.ndim != 3 or coords.shape[2] != 3:
        raise LMDBError(
            f"coords must have shape (nframes, natoms, 3), got {coords.shape}."
        )
    source = _prepare_atom_source(data, coords.shape[0], coords.shape[1])
    types = source.types[source.types >= 0]
    active_indices = set(int(value) for value in np.unique(types))
    return tuple(
        name for index, name in enumerate(source.names) if index in active_indices
    )


def _remove_path(path: Path) -> None:
    """Remove a file or directory if it exists."""
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    elif path.exists() or path.is_symlink():
        path.unlink()


def _normalized_path(path: str | Path) -> str:
    """Return one cache key for relative, absolute, and case aliases."""
    return os.path.normcase(os.path.realpath(os.fspath(path)))


def _reset_read_cache_after_fork() -> None:
    """Discard inherited LMDB environments and synchronization state."""
    global _READ_ENV_CACHE, _READ_ENV_LOCK, _PUBLISHING_PATHS, _READ_ENV_PID

    for env, _ in _READ_ENV_CACHE.values():
        try:
            env.close()
        except lmdb.Error:
            pass
    _READ_ENV_CACHE = {}
    _PUBLISHING_PATHS = set()
    _READ_ENV_LOCK = threading.Lock()
    _READ_ENV_PID = os.getpid()


if hasattr(os, "register_at_fork"):
    os.register_at_fork(after_in_child=_reset_read_cache_after_fork)


def _open_read_env(path: str) -> tuple[str, lmdb.Environment]:
    """Open or share a read-only LMDB environment within dpdata."""
    global _READ_ENV_PID

    if _READ_ENV_PID != os.getpid():
        _reset_read_cache_after_fork()
    resolved = _normalized_path(path)
    with _READ_ENV_LOCK:
        if resolved in _PUBLISHING_PATHS:
            raise LMDBError(
                f"LMDB '{path}' is being replaced and cannot be opened."
            )
        cached = _READ_ENV_CACHE.get(resolved)
        if cached is not None:
            env, count = cached
            _READ_ENV_CACHE[resolved] = (env, count + 1)
            return resolved, env
        try:
            env = lmdb.open(
                path,
                readonly=True,
                lock=False,
                subdir=True,
                readahead=False,
                meminit=False,
            )
        except lmdb.Error as exc:
            raise LMDBError(
                f"Cannot open LMDB '{path}'. If another library has already opened "
                "this path in the current process, close that reader before using "
                "dpdata."
            ) from exc
        _READ_ENV_CACHE[resolved] = (env, 1)
        return resolved, env


def _close_read_env(resolved: str) -> None:
    """Release a cached read-only LMDB environment."""
    with _READ_ENV_LOCK:
        env, count = _READ_ENV_CACHE[resolved]
        if count == 1:
            del _READ_ENV_CACHE[resolved]
            env.close()
        else:
            _READ_ENV_CACHE[resolved] = (env, count - 1)


def _begin_publish(path: Path) -> str:
    """Reserve a destination for publication after checking active readers."""
    resolved = _normalized_path(path)
    with _READ_ENV_LOCK:
        if resolved in _READ_ENV_CACHE:
            raise LMDBError(
                f"Cannot replace LMDB '{path}' while a dpdata reader is active."
            )
        if resolved in _PUBLISHING_PATHS:
            raise LMDBError(f"LMDB '{path}' is already being replaced.")
        _PUBLISHING_PATHS.add(resolved)
    return resolved


def _end_publish(resolved: str) -> None:
    """Release a publication reservation."""
    with _READ_ENV_LOCK:
        _PUBLISHING_PATHS.discard(resolved)


def _open_publish_guard(path: Path) -> lmdb.Environment:
    """Hold the process-wide LMDB slot until atomic replacement completes."""
    try:
        return lmdb.open(
            str(path),
            readonly=True,
            lock=False,
            subdir=True,
            readahead=False,
            meminit=False,
        )
    except lmdb.Error as exc:
        raise LMDBError(
            f"Cannot replace LMDB '{path}' while another LMDB environment is "
            "open in this process. Close DeePMD-kit and other readers first."
        ) from exc


class _LMDBWriter:
    """Write a complete LMDB to a temporary sibling and publish it safely."""

    def __init__(
        self,
        directory: str,
        *,
        map_size: int,
        frame_idx_fmt: str,
        type_map: list[str] | None,
        write_batch_size: int,
        overwrite: bool,
    ) -> None:
        if write_batch_size <= 0:
            raise ValueError("write_batch_size must be positive.")
        normalized_type_map = (
            list(_validate_names(type_map, "type_map"))
            if type_map is not None
            else None
        )
        self.directory = Path(directory).resolve()
        if self.directory.exists() and not overwrite:
            raise FileExistsError(
                f"LMDB destination '{directory}' already exists; pass "
                "overwrite=True to replace it."
            )
        self.directory.parent.mkdir(parents=True, exist_ok=True)
        self.temp_directory = Path(
            tempfile.mkdtemp(
                prefix=f".{self.directory.name}.tmp-",
                dir=self.directory.parent,
            )
        )
        self.overwrite = overwrite
        self.frame_idx_fmt = frame_idx_fmt
        self.type_map = normalized_type_map
        self.write_batch_size = write_batch_size
        self.frame_index = 0
        self.system_index = 0
        self.frame_nlocs: list[int] = []
        self.frame_system_ids: list[int] = []
        self.data_shapes: dict[str, list[int | str]] = {}
        self.data_names: dict[str, str] = {}
        self._closed = False
        self._published = False
        try:
            self.env = lmdb.open(
                str(self.temp_directory), map_size=map_size, subdir=True
            )
        except BaseException:
            _remove_path(self.temp_directory)
            raise

    def _build_specs(
        self,
        data: dict[str, Any],
        *,
        nframes: int,
        natoms: int,
        has_virtual_atoms: bool,
    ) -> list[_FieldSpec]:
        """Build and validate field specifications for one system."""
        by_dp, _ = _field_registries()
        all_dtypes = _all_data_types()
        specs: list[_FieldSpec] = []
        used_disk_names: dict[str, str] = {}
        for name, value in data.items():
            if name in _RESERVED_KEYS or not isinstance(value, np.ndarray):
                continue
            array = np.asarray(value)
            registered_dtype = all_dtypes.get(name)
            if registered_dtype is not None and registered_dtype.shape is None:
                raise LMDBError(
                    f"Data type '{name}' has no declared shape and cannot be "
                    "represented unambiguously in per-frame LMDB records."
                )
            registered = name in by_dp
            if registered:
                spec = by_dp[name]
            else:
                if array.ndim > 0 and array.shape[0] == nframes:
                    shape: tuple[int | Axis, ...] = (
                        Axis.NFRAMES,
                        *array.shape[1:],
                    )
                else:
                    shape = tuple(array.shape)
                spec = _FieldSpec(name, name, shape, name)
                payload_shape = list(shape)
                if Axis.NFRAMES in payload_shape:
                    payload_shape.remove(Axis.NFRAMES)
                if has_virtual_atoms and natoms in payload_shape:
                    raise LMDBError(
                        f"Unregistered field '{name}' may contain an atom axis, "
                        "but its shape is ambiguous for mixed data with virtual "
                        "atoms. Register a DataType with Axis.NATOMS before writing."
                    )
            _validate_field_namespace(spec.dp_name, spec.disk_name)
            previous = used_disk_names.get(spec.disk_name)
            if previous is not None and previous != spec.dp_name:
                raise LMDBError(
                    f"Fields '{previous}' and '{spec.dp_name}' both map to LMDB "
                    f"key '{spec.disk_name}'."
                )
            used_disk_names[spec.disk_name] = spec.dp_name
            self._validate_source_shape(spec, array, nframes, natoms)
            specs.append(spec)
            if spec.disk_name not in _CORE_FRAME_KEYS:
                existing_shape = self.data_shapes.get(spec.disk_name)
                tokens = _shape_tokens(spec.shape)
                if existing_shape is not None and existing_shape != tokens:
                    raise LMDBError(
                        f"Field '{spec.dp_name}' has incompatible shapes across "
                        f"systems: {existing_shape} and {tokens}."
                    )
                self.data_shapes[spec.disk_name] = tokens
            if spec.disk_name != spec.dp_name:
                existing_name = self.data_names.get(spec.disk_name)
                if existing_name is not None and existing_name != spec.dp_name:
                    raise LMDBError(
                        f"LMDB key '{spec.disk_name}' maps to conflicting dpdata "
                        f"names '{existing_name}' and '{spec.dp_name}'."
                    )
                self.data_names[spec.disk_name] = spec.dp_name
        if "coords" not in used_disk_names:
            raise LMDBError("coords is required for every system.")
        return specs

    @staticmethod
    def _validate_source_shape(
        spec: _FieldSpec, array: np.ndarray, nframes: int, natoms: int
    ) -> None:
        """Validate an in-memory field against its symbolic shape."""
        if array.ndim != len(spec.shape):
            raise LMDBError(
                f"Field '{spec.dp_name}' has rank {array.ndim}, expected "
                f"{len(spec.shape)} from {spec.shape}."
            )
        for axis, (actual, expected) in enumerate(zip(array.shape, spec.shape)):
            if expected is Axis.NFRAMES and actual != nframes:
                raise LMDBError(
                    f"Field '{spec.dp_name}' axis {axis} must have {nframes} "
                    f"frames, got {actual}."
                )
            if expected is Axis.NATOMS and actual != natoms:
                raise LMDBError(
                    f"Field '{spec.dp_name}' axis {axis} must have {natoms} "
                    f"atoms, got {actual}."
                )
            if isinstance(expected, int) and expected != -1 and actual != expected:
                raise LMDBError(
                    f"Field '{spec.dp_name}' axis {axis} must have length "
                    f"{expected}, got {actual}."
                )

    @staticmethod
    def _frame_payload(
        array: np.ndarray,
        spec: _FieldSpec,
        frame_index: int,
        atom_mask: np.ndarray | None,
    ) -> np.ndarray:
        """Extract one frame and remove virtual atoms on every atom axis."""
        frame_axis = spec.frame_axis
        payload = (
            np.take(array, frame_index, axis=frame_axis)
            if frame_axis is not None
            else array
        )
        if atom_mask is not None:
            for axis in spec.payload_atom_axes:
                payload = np.compress(atom_mask, payload, axis=axis)
        return payload

    def _put_records(self, records: list[tuple[bytes, bytes]]) -> None:
        """Commit one batch, growing the LMDB map and retrying when necessary."""
        while True:
            try:
                with self.env.begin(write=True) as txn:
                    for key, value in records:
                        txn.put(key, value)
                return
            except lmdb.MapFullError:
                current_size = int(self.env.info()["map_size"])
                self.env.set_mapsize(
                    max(current_size * 2, current_size + (1 << 20))
                )

    def write_system(self, data: dict[str, Any]) -> None:
        """Append one validated source system."""
        coords = np.asarray(data["coords"])
        if coords.ndim != 3 or coords.shape[2] != 3:
            raise LMDBError(
                f"coords must have shape (nframes, natoms, 3), got {coords.shape}."
            )
        nframes, natoms = coords.shape[:2]
        if nframes == 0:
            raise LMDBError("A source system must contain at least one frame.")
        atom_source = _prepare_atom_source(data, nframes, natoms)

        if self.type_map is None:
            self.type_map = list(atom_source.names)
        type_map = _validate_names(self.type_map, "type_map")
        name_to_global = {name: index for index, name in enumerate(type_map)}
        active_indices = np.unique(atom_source.types[atom_source.types >= 0])
        active_names = {atom_source.names[int(index)] for index in active_indices}
        missing_active = sorted(active_names - set(type_map))
        if missing_active:
            raise LMDBError(
                f"Elements {missing_active} are not present in type_map {type_map}."
            )
        local_to_global = np.array(
            [name_to_global.get(name, -1) for name in atom_source.names],
            dtype=np.int64,
        )
        has_virtual_atoms = (
            atom_source.masks is not None and not np.all(atom_source.masks)
        )
        specs = self._build_specs(
            data,
            nframes=nframes,
            natoms=natoms,
            has_virtual_atoms=has_virtual_atoms,
        )
        nopbc = bool(data.get("nopbc", False))
        fields = [(spec, np.asarray(data[spec.dp_name])) for spec in specs]

        for batch_start in range(0, nframes, self.write_batch_size):
            batch_end = min(batch_start + self.write_batch_size, nframes)
            records: list[tuple[bytes, bytes]] = []
            batch_nlocs: list[int] = []
            for offset, frame_i in enumerate(range(batch_start, batch_end)):
                if atom_source.masks is None:
                    local_types = atom_source.types
                    atom_mask = None
                else:
                    atom_mask = atom_source.masks[frame_i]
                    local_types = atom_source.types[frame_i, atom_mask]
                global_types = local_to_global[local_types]
                if np.any(global_types < 0):
                    missing_indices = np.unique(local_types[global_types < 0])
                    missing_names = [
                        atom_source.names[int(index)] for index in missing_indices
                    ]
                    raise LMDBError(
                        f"Elements {missing_names} are not present in type_map "
                        f"{list(type_map)}."
                    )
                global_types = global_types.astype(np.int32)
                frame: dict[str, Any] = {
                    "atom_types": _encode_array(global_types)
                }
                for spec, array in fields:
                    if spec.dp_name == "cells" and nopbc:
                        continue
                    payload = self._frame_payload(
                        array, spec, frame_i, atom_mask
                    )
                    frame[spec.disk_name] = _encode_array(payload)
                counts = np.bincount(
                    global_types, minlength=len(type_map)
                )[: len(type_map)]
                frame["atom_numbs"] = [int(count) for count in counts]
                key = format(
                    self.frame_index + offset, self.frame_idx_fmt
                ).encode("ascii")
                records.append((key, _packb(frame)))
                batch_nlocs.append(int(global_types.size))
            self._put_records(records)
            self.frame_nlocs.extend(batch_nlocs)
            self.frame_system_ids.extend(
                [self.system_index] * len(records)
            )
            self.frame_index += len(records)
        self.system_index += 1

    def _metadata(self) -> dict[str, Any]:
        """Return validated metadata for the completed database."""
        if self.system_index == 0 or self.frame_index == 0:
            raise LMDBError("Cannot write an empty LMDB dataset.")
        if self.type_map is None:
            raise LMDBError("Cannot write LMDB metadata without a type_map.")
        metadata: dict[str, Any] = {
            "nframes": self.frame_index,
            "frame_idx_fmt": self.frame_idx_fmt,
            "type_map": self.type_map,
            "frame_nlocs": self.frame_nlocs,
            "frame_system_ids": self.frame_system_ids,
        }
        if self.data_shapes:
            metadata["dp_data_shapes"] = self.data_shapes
        if self.data_names:
            metadata["dp_data_names"] = self.data_names
        return metadata

    def finish(self) -> None:
        """Commit metadata, close the temporary database, and publish it."""
        self._put_records([(METADATA_KEY, _packb(self._metadata()))])
        self.env.sync()
        self.env.close()
        self._closed = True
        self._validate_staged_database()
        self._publish()
        self._published = True

    def _validate_staged_database(self) -> None:
        """Reopen the staged database and verify its complete key sequence."""
        env = lmdb.open(
            str(self.temp_directory),
            readonly=True,
            lock=False,
            subdir=True,
            readahead=False,
            meminit=False,
        )
        try:
            with env.begin() as txn:
                entries = int(txn.stat()["entries"])
                if entries != self.frame_index + 1:
                    raise LMDBError(
                        f"Staged LMDB has {entries} entries; expected "
                        f"{self.frame_index + 1}."
                    )
                metadata = _read_metadata(txn)
                if _meta_get(metadata, "nframes") != self.frame_index:
                    raise LMDBError("Staged LMDB nframes does not match written frames.")
                if len(_meta_get(metadata, "frame_nlocs") or []) != self.frame_index:
                    raise LMDBError(
                        "Staged LMDB frame_nlocs length is inconsistent."
                    )
                if (
                    len(_meta_get(metadata, "frame_system_ids") or [])
                    != self.frame_index
                ):
                    raise LMDBError(
                        "Staged LMDB frame_system_ids length is inconsistent."
                    )
                for index in range(self.frame_index):
                    key = format(index, self.frame_idx_fmt).encode("ascii")
                    if txn.get(key) is None:
                        raise LMDBError(
                            f"Staged LMDB is missing frame key {key!r}."
                        )
        finally:
            env.close()

    def _publish(self) -> None:
        """Publish by one atomic rename of a directory or its data file."""
        reservation = _begin_publish(self.directory)
        try:
            if not self.directory.exists():
                os.replace(self.temp_directory, self.directory)
                return
            if not self.overwrite:
                raise FileExistsError(
                    f"LMDB destination '{self.directory}' already exists."
                )
            target_data = self.directory / "data.mdb"
            staged_data = self.temp_directory / "data.mdb"
            if not self.directory.is_dir() or not target_data.is_file():
                raise LMDBError(
                    f"Existing destination '{self.directory}' is not an LMDB "
                    "directory containing data.mdb."
                )
            guard = _open_publish_guard(self.directory)
            try:
                os.replace(staged_data, target_data)
            finally:
                guard.close()
            _remove_path(self.temp_directory)
        finally:
            _end_publish(reservation)

    def abort(self) -> None:
        """Close and discard an incomplete temporary database."""
        if not self._closed:
            self.env.close()
            self._closed = True
        if not self._published:
            _remove_path(self.temp_directory)


class LMDBFormat(Format):
    """DeePMD-kit compatible LMDB format.

    A single flat LMDB stores all frames from one or many systems. The
    same on-disk format is produced regardless of whether the source is a
    standard or a mixed-type system, so the output is always readable by
    DeePMD-kit's ``LmdbDataReader``.

    The ``mixed_type`` keyword controls only how frames are mapped back to
    ``dpdata`` objects on read (see :meth:`from_multi_systems`).

    Examples
    --------
    Write a single labeled system::

        >>> import dpdata
        >>> ls = dpdata.LabeledSystem("OUTCAR", fmt="vasp/outcar")
        >>> ls.to("lmdb", "data.lmdb")

    Write many systems into one LMDB, forcing a global type map::

        >>> ms = dpdata.MultiSystems(s1, s2, s3)
        >>> ms.to("lmdb", "data.lmdb", type_map=["H", "C", "N", "O"])

    Read back as standard (per-composition) systems::

        >>> ms = dpdata.MultiSystems.from_file("data.lmdb", fmt="lmdb")

    Read back keeping the full global type map on every system::

        >>> ms = dpdata.MultiSystems.from_file(
        ...     "data.lmdb", fmt="lmdb", mixed_type=True
        ... )

    Note that loading through :class:`dpdata.MultiSystems` normalises the
    ``atom_names`` order (the element set is kept, but reordered); a direct
    single-system load preserves the stored order.
    """

    # ------------------------------------------------------------------ #
    # Writing
    # ------------------------------------------------------------------ #
    def to_multi_systems(
        self,
        formulas,
        directory,
        map_size: int = DEFAULT_MAP_SIZE,
        frame_idx_fmt: str = DEFAULT_FRAME_FMT,
        type_map: list[str] | None = None,
        write_batch_size: int = DEFAULT_WRITE_BATCH_SIZE,
        overwrite: bool = False,
        **kwargs,
    ):
        """Write multiple dpdata systems to one LMDB.

        Parameters
        ----------
        formulas : list[str]
            One handle per system (the value is not used on disk).
        directory : str
            Output LMDB directory.
        map_size : int, optional
            Maximum LMDB size in bytes. Default is 1 TiB (sparse).
        frame_idx_fmt : str, optional
            Format used for the per-frame integer key. Default ``"012d"``.
        type_map : list[str], optional
            Global element table. If ``None``, the element list of the
            first system written is used (for a ``MultiSystems`` this is
            the union of all systems' elements).
        write_batch_size : int, optional
            Number of frames committed per LMDB write transaction.
        overwrite : bool, optional
            Whether to replace an existing destination after the new database
            has been written and validated. The default is ``False``.
        **kwargs : dict
            other parameters

        Yields
        ------
        tuple
            ``(self, formula)`` handle consumed by :meth:`to_system`.
        """
        writer = _LMDBWriter(
            directory,
            map_size=map_size,
            frame_idx_fmt=frame_idx_fmt,
            type_map=type_map,
            write_batch_size=write_batch_size,
            overwrite=overwrite,
        )
        self._writer = writer
        completed = False
        try:
            for ff in formulas:
                yield (self, ff)
            writer.finish()
            completed = True
        finally:
            if not completed:
                writer.abort()
            self._writer = None

    def dump_systems(
        self,
        systems,
        directory,
        map_size: int = DEFAULT_MAP_SIZE,
        frame_idx_fmt: str = DEFAULT_FRAME_FMT,
        type_map: list[str] | None = None,
        write_batch_size: int = DEFAULT_WRITE_BATCH_SIZE,
        overwrite: bool = False,
    ):
        """Write an ordered sequence of systems, one ``frame_system_id`` each.

        Unlike :meth:`to_multi_systems` (the path used by
        ``MultiSystems.to('lmdb', ...)``), the systems are **not** merged by
        formula: every element of ``systems`` becomes exactly one source
        system in the database, numbered ``0, 1, 2, ...`` in iteration order.
        This preserves the system partition recorded in ``frame_system_ids``,
        which DeePMD-kit uses for ``prob_sys_size`` based sampling.

        Parameters
        ----------
        systems : Iterable[System]
            An ordered iterable of :class:`dpdata.System` /
            :class:`dpdata.LabeledSystem` objects (or raw data dicts). A
            :class:`dpdata.MultiSystems` must not be used here, because its
            frames are already merged by formula.
        directory : str
            Output LMDB directory.
        map_size : int, optional
            Maximum LMDB size in bytes. Default is 1 TiB (sparse).
        frame_idx_fmt : str, optional
            Format used for the per-frame integer key. Default ``"012d"``.
        type_map : list[str], optional
            Global element table. If ``None``, the union of the elements of
            all systems is used, in first-appearance order. When the systems
            are produced lazily (a generator) a ``type_map`` should be given
            so that the whole sequence need not be held in memory.
        write_batch_size : int, optional
            Number of frames committed per LMDB write transaction.
        overwrite : bool, optional
            Whether to replace an existing destination after the new database
            has been written and validated. The default is ``False``.
        """
        if hasattr(systems, "systems"):
            raise TypeError(
                "dump_systems expects the original ordered systems, not a "
                "MultiSystems object whose systems are already merged by formula."
            )
        if type_map is None:
            systems = list(systems)
            if not systems:
                raise LMDBError("Cannot write an empty LMDB dataset.")
            resolved: list[str] = []
            seen: set[str] = set()
            for ss in systems:
                data = _source_data(ss)
                for name in _source_active_names(data):
                    if name not in seen:
                        seen.add(name)
                        resolved.append(name)
        else:
            resolved = [str(n) for n in type_map]
            iterator = iter(systems)
            try:
                first = next(iterator)
            except StopIteration:
                raise LMDBError("Cannot write an empty LMDB dataset.") from None
            systems = itertools.chain((first,), iterator)

        writer = _LMDBWriter(
            directory,
            map_size=map_size,
            frame_idx_fmt=frame_idx_fmt,
            type_map=resolved,
            write_batch_size=write_batch_size,
            overwrite=overwrite,
        )
        try:
            for ss in systems:
                writer.write_system(_source_data(ss))
            writer.finish()
        except BaseException:
            writer.abort()
            raise

    def to_system(self, data, file_name, **kwargs):
        """Save a single (unlabeled) System to an LMDB database."""
        self._to_any(data, file_name, **kwargs)

    def to_labeled_system(self, data, file_name, **kwargs):
        """Save a single LabeledSystem to an LMDB database."""
        self._to_any(data, file_name, **kwargs)

    def _to_any(self, data, file_name, **kwargs):
        if isinstance(file_name, tuple) and file_name[0] is self:
            if self._writer is None:
                raise LMDBError("LMDB writer is not active.")
            self._writer.write_system(data)
            return
        gen = self.to_multi_systems(["system"], file_name, **kwargs)
        try:
            next(gen)
            writer = self._writer
            if writer is None:
                raise LMDBError("LMDB writer is not active.")
            writer.write_system(data)
            next(gen)
        except StopIteration:
            return
        except BaseException:
            gen.close()
            raise

    # ------------------------------------------------------------------ #
    # Reading
    # ------------------------------------------------------------------ #
    def from_multi_systems(
        self,
        directory,
        mixed_type: bool = False,
        type_map: list[str] | None = None,
        max_frames: int | None = DEFAULT_MAX_FRAMES,
        **kwargs,
    ):
        """Load systems from a flat LMDB.

        Frames are grouped by atom-count composition. Atom order is
        canonicalized by a stable sort on the global atom type, and every
        registered atomic field follows the same permutation. Each
        composition becomes one ``dpdata`` system.

        Parameters
        ----------
        directory : str
            Path to the LMDB directory.
        mixed_type : bool, optional
            If ``False`` (default) each system's ``atom_names`` is the
            compact set of elements it actually contains. If ``True`` every
            system keeps the full global ``type_map`` as ``atom_names``
            (with zero counts for absent elements).
        type_map : list[str], optional
            Requested element table for the returned systems. When the file
            stores a ``type_map``, the stored global indices are remapped to
            this table **by element name** (consistent with DeePMD-kit); every
            element in the file must be present in ``type_map``. When the file
            has no ``type_map``, the indices are named positionally from this
            argument. Defaults to the ``type_map`` stored in the file.
        max_frames : int or None, optional
            Maximum number of frames loaded into memory. The default is
            100,000. Set to ``None`` only when sufficient memory is available.
        **kwargs : dict
            other parameters

        Yields
        ------
        dict
            system data dictionary for each composition group.
        """
        frames, meta = self._read_all_frames(directory, max_frames=max_frames)
        file_type_map = _meta_get(meta, "type_map")
        if file_type_map:
            names = list(
                _validate_names(
                    [
                        n.decode() if isinstance(n, bytes) else n
                        for n in file_type_map
                    ],
                    "file type_map",
                )
            )
            self._validate_frame_types(frames, names, meta)
            if type_map is not None:
                requested = list(_validate_names(type_map, "type_map"))
                missing = [name for name in names if name not in requested]
                if missing:
                    raise LMDBError(
                        f"Elements {missing} from the file type_map are missing "
                        f"from the requested type_map {requested}."
                    )
                remap = np.array([requested.index(n) for n in names], dtype=int)
                for fr in frames:
                    fr["atom_types"] = remap[np.asarray(fr["atom_types"]).astype(int)]
                    fr["atom_numbs"] = np.bincount(
                        fr["atom_types"], minlength=len(requested)
                    ).tolist()
                names = requested
        elif type_map is not None:
            names = list(_validate_names(type_map, "type_map"))
            self._validate_frame_types(frames, names, meta)
        else:
            names = self._infer_type_names(frames)
            self._validate_frame_types(frames, names, meta)
        self._normalize_reference_fields(frames)
        specs = self._resolve_read_specs(frames, meta)
        frames = self._rename_frame_fields(frames, specs)
        datasets = self._group_frames(frames, names, mixed_type, specs)
        dp_specs = {spec.dp_name: spec for spec in specs.values()}
        from dpdata.system import LabeledSystem, System

        original_system_dtypes = System.DTYPES
        original_labeled_dtypes = LabeledSystem.DTYPES
        self._register_specs(datasets, dp_specs)
        completed = False
        try:
            yield from datasets
            completed = True
        finally:
            if not completed:
                System.DTYPES = original_system_dtypes
                LabeledSystem.DTYPES = original_labeled_dtypes

    @staticmethod
    def _normalize_reference_fields(frames: list[dict[str, Any]]) -> None:
        """Normalize flattened atomic parameters emitted by reference converters."""
        for index, frame in enumerate(frames):
            if "aparam" not in frame:
                continue
            value = np.asarray(frame["aparam"])
            natoms = len(frame["atom_types"])
            if value.ndim == 1:
                if value.size == 0 or value.size % natoms:
                    raise LMDBFrameError(
                        f"Frame {index} aparam length {value.size} is not a "
                        f"positive multiple of natoms={natoms}."
                    )
                frame["aparam"] = value.reshape(natoms, value.size // natoms)
            elif value.ndim != 2 or value.shape[0] != natoms:
                raise LMDBFrameError(
                    f"Frame {index} aparam must have shape (natoms, ndof) or "
                    f"a flattened length divisible by natoms={natoms}, got "
                    f"{value.shape}."
                )

    def from_system(self, file_name, **kwargs):
        """Load data for a single System from an LMDB database."""
        if isinstance(file_name, dict):
            return file_name
        return self._first_system(file_name, require_labeled=False, **kwargs)

    def from_labeled_system(self, file_name, **kwargs):
        """Load data for a single LabeledSystem from an LMDB database."""
        if isinstance(file_name, dict):
            return file_name
        return self._first_system(file_name, require_labeled=True, **kwargs)

    def _first_system(self, file_name, *, require_labeled: bool, **kwargs):
        """Return the first composition group, warning if more than one exists.

        A single ``System`` can only hold one composition; an LMDB usually
        stores several. Callers that need every system should use
        ``MultiSystems.from_file``.
        """
        from dpdata.system import LabeledSystem, System

        original_system_dtypes = System.DTYPES
        original_labeled_dtypes = LabeledSystem.DTYPES
        systems = list(self.from_multi_systems(file_name, **kwargs))
        first = systems[0]
        if require_labeled and "energies" not in first:
            System.DTYPES = original_system_dtypes
            LabeledSystem.DTYPES = original_labeled_dtypes
            raise LMDBFrameError(
                "LMDB does not contain energies required by LabeledSystem. "
                "Load it as System or pass labeled=False to MultiSystems."
            )
        if len(systems) > 1:
            warnings.warn(
                f"LMDB '{file_name}' contains more than one composition; only the "
                "first is loaded into a single System. Use "
                "dpdata.MultiSystems.from_file(..., fmt='lmdb') to load all of them.",
                stacklevel=2,
            )
        return first

    def _read_all_frames(
        self, directory: str, *, max_frames: int | None
    ) -> tuple[list[dict[str, Any]], dict[str | bytes, Any]]:
        resolved, env = _open_read_env(directory)
        try:
            with env.begin() as txn:
                meta = _read_metadata(txn)
                frame_fmt = _meta_get(meta, "frame_idx_fmt") or DEFAULT_FRAME_FMT
                if not isinstance(frame_fmt, str):
                    raise LMDBMetadataError("frame_idx_fmt must be a string.")
                try:
                    format(0, frame_fmt)
                except (TypeError, ValueError) as exc:
                    raise LMDBMetadataError(
                        f"Invalid frame_idx_fmt {frame_fmt!r}: {exc}"
                    ) from exc
                nframes = _meta_get(meta, "nframes")
                if nframes is None:
                    raise LMDBMetadataError("Metadata does not contain 'nframes'.")
                nframes = _require_integer(nframes, "nframes", minimum=1)
                if max_frames is not None:
                    if max_frames <= 0:
                        raise ValueError("max_frames must be positive or None.")
                    if nframes > max_frames:
                        raise LMDBError(
                            f"LMDB '{directory}' contains {nframes:,} frames, "
                            f"exceeding max_frames={max_frames:,}. Use DeePMD-kit's "
                            "streaming reader, or set max_frames=None when sufficient "
                            "memory is available."
                        )
                frame_nlocs = _meta_get(meta, "frame_nlocs")
                if frame_nlocs is not None:
                    if not isinstance(frame_nlocs, (list, tuple)):
                        raise LMDBMetadataError(
                            "frame_nlocs must be a sequence of integers."
                        )
                    if len(frame_nlocs) != nframes:
                        raise LMDBMetadataError(
                            "frame_nlocs length does not match nframes."
                        )
                    frame_nlocs = [
                        _require_integer(
                            value,
                            f"frame_nlocs[{index}]",
                            minimum=1,
                        )
                        for index, value in enumerate(frame_nlocs)
                    ]
                    meta["frame_nlocs"] = frame_nlocs
                frame_system_ids = _meta_get(meta, "frame_system_ids")
                if frame_system_ids is not None:
                    if not isinstance(frame_system_ids, (list, tuple)):
                        raise LMDBMetadataError(
                            "frame_system_ids must be a sequence of integers."
                        )
                    if len(frame_system_ids) != nframes:
                        raise LMDBMetadataError(
                            "frame_system_ids length does not match nframes."
                        )
                    frame_system_ids = [
                        _require_integer(
                            value,
                            f"frame_system_ids[{index}]",
                            minimum=0,
                        )
                        for index, value in enumerate(frame_system_ids)
                    ]
                    meta["frame_system_ids"] = frame_system_ids
                frames: list[dict[str, Any]] = []
                for i in range(int(nframes)):
                    key = format(i, frame_fmt).encode("ascii")
                    raw = txn.get(key)
                    if raw is None:
                        raise LMDBFrameError(f"Frame data not found for key: {key!r}")
                    frames.append(_decode_frame(bytes(raw)))
        finally:
            _close_read_env(resolved)
        return frames, meta

    @staticmethod
    def _infer_type_names(frames: list[dict[str, Any]]) -> list[str]:
        """Create positional names for a file without type_map metadata."""
        maximum = -1
        for frame in frames:
            raw = np.asarray(frame["atom_types"])
            if not np.issubdtype(raw.dtype, np.integer) or raw.ndim != 1:
                raise LMDBFrameError("atom_types must be a one-dimensional integer array.")
            if raw.size and int(raw.min()) < 0:
                raise LMDBFrameError("atom_types must not contain negative indices.")
            if raw.size:
                maximum = max(maximum, int(raw.max()))
        if maximum < 0:
            raise LMDBFrameError("Every frame must contain at least one atom.")
        return [f"Type_{index}" for index in range(maximum + 1)]

    @staticmethod
    def _validate_frame_types(
        frames: list[dict[str, Any]],
        names: list[str],
        metadata: dict[str | bytes, Any],
    ) -> None:
        """Validate per-frame atom types and count metadata."""
        frame_nlocs = _meta_get(metadata, "frame_nlocs")
        for index, frame in enumerate(frames):
            try:
                atom_types = _validate_integer_types(
                    np.asarray(frame["atom_types"]),
                    name=f"frame {index} atom_types",
                    upper_bound=len(names),
                    allow_virtual=False,
                )
            except (KeyError, LMDBError) as exc:
                raise LMDBFrameError(str(exc)) from exc
            if atom_types.ndim != 1 or atom_types.size == 0:
                raise LMDBFrameError(
                    f"Frame {index} atom_types must be a non-empty 1-D array."
                )
            if frame_nlocs is not None and int(frame_nlocs[index]) != atom_types.size:
                raise LMDBMetadataError(
                    f"frame_nlocs[{index}] does not match atom_types length."
                )
            expected = np.bincount(atom_types, minlength=len(names)).tolist()
            stored = frame.get("atom_numbs")
            if stored is not None and list(stored) != expected:
                raise LMDBFrameError(
                    f"Frame {index} atom_numbs is inconsistent with atom_types."
                )
            frame["atom_types"] = atom_types
            frame["atom_numbs"] = expected

    @staticmethod
    def _metadata_mapping(
        metadata: dict[str | bytes, Any], key: str
    ) -> dict[str, Any]:
        """Return a metadata mapping with decoded string keys."""
        raw = _meta_get(metadata, key) or {}
        if not isinstance(raw, dict):
            raise LMDBMetadataError(f"{key} must be a mapping.")
        return {
            k.decode() if isinstance(k, bytes) else str(k): v
            for k, v in raw.items()
        }

    def _resolve_read_specs(
        self,
        frames: list[dict[str, Any]],
        metadata: dict[str | bytes, Any],
    ) -> dict[str, _FieldSpec]:
        """Resolve every disk field to a dpdata name and symbolic shape."""
        by_dp, by_disk = _field_registries(prefer_protocol=True)
        shape_hints = self._metadata_mapping(metadata, "dp_data_shapes")
        name_hints = self._metadata_mapping(metadata, "dp_data_names")
        disk_names = {
            key
            for frame in frames
            for key in frame
            if key not in {"atom_types", "atom_numbs"}
        }
        specs: dict[str, _FieldSpec] = {}
        used_dp_names: dict[str, str] = {}
        for disk_name in sorted(disk_names):
            hinted_name = name_hints.get(disk_name)
            if isinstance(hinted_name, bytes):
                hinted_name = hinted_name.decode()
            known = by_disk.get(disk_name)
            dp_name = str(hinted_name) if hinted_name is not None else (
                known.dp_name if known is not None else disk_name
            )
            if known is not None and hinted_name is not None and dp_name != known.dp_name:
                raise LMDBMetadataError(
                    f"dp_data_names maps '{disk_name}' to '{dp_name}', but the "
                    f"registered protocol maps it to '{known.dp_name}'."
                )
            expected_spec = by_dp.get(dp_name)
            if (
                expected_spec is not None
                and disk_name != expected_spec.disk_name
            ):
                raise LMDBMetadataError(
                    f"LMDB key '{disk_name}' cannot map to registered field "
                    f"'{dp_name}', whose protocol key is "
                    f"'{expected_spec.disk_name}'."
                )
            try:
                _validate_field_namespace(dp_name, disk_name)
            except LMDBError as exc:
                raise LMDBFrameError(str(exc)) from exc
            previous_disk = used_dp_names.get(dp_name)
            if previous_disk is not None and previous_disk != disk_name:
                raise LMDBMetadataError(
                    f"LMDB keys '{previous_disk}' and '{disk_name}' both map to "
                    f"dpdata field '{dp_name}'."
                )
            used_dp_names[dp_name] = disk_name

            if disk_name in shape_hints:
                shape = _tokens_to_shape(shape_hints[disk_name])
                if (
                    known is not None
                    and not self._shape_hint_compatible(known.shape, shape)
                ):
                    raise LMDBMetadataError(
                        f"dp_data_shapes for known field '{disk_name}' changes "
                        f"its protocol shape from {known.shape} to {shape}."
                    )
            elif (
                known is not None
                and known.dp_name in _PROTOCOL_FIELD_NAMES
            ):
                shape = known.shape
            else:
                shape = self._infer_field_shape(disk_name, frames)
            if shape is None:
                raise LMDBMetadataError(
                    f"No shape is available for LMDB field '{disk_name}'."
                )
            deepmd_name = (
                known.deepmd_name
                if known is not None
                else by_dp[dp_name].deepmd_name
                if dp_name in by_dp
                else disk_name
            )
            specs[disk_name] = _FieldSpec(
                dp_name=dp_name,
                disk_name=disk_name,
                shape=tuple(shape),
                deepmd_name=deepmd_name,
            )
        return specs

    @staticmethod
    def _shape_hint_compatible(
        protocol: tuple[int | Axis, ...],
        hinted: tuple[int | Axis, ...],
    ) -> bool:
        """Allow metadata to refine wildcard dimensions but not protocol axes."""
        if len(protocol) != len(hinted):
            return False
        for expected, actual in zip(protocol, hinted):
            if expected == -1:
                if isinstance(actual, int):
                    continue
                return False
            if expected != actual:
                return False
        return True

    @staticmethod
    def _infer_field_shape(
        disk_name: str, frames: list[dict[str, Any]]
    ) -> tuple[int | Axis, ...]:
        """Infer unknown shapes across all atom counts without local coincidences."""
        observed = [
            (np.asarray(frame[disk_name]).shape, len(frame["atom_types"]))
            for frame in frames
            if disk_name in frame
        ]
        ranks = {len(shape) for shape, _ in observed}
        if len(ranks) != 1:
            raise LMDBFrameError(
                f"Field '{disk_name}' has inconsistent ranks across frames."
            )
        rank = ranks.pop()
        inferred: list[int | Axis] = [Axis.NFRAMES]
        for axis in range(rank):
            dimensions = [shape[axis] for shape, _ in observed]
            follows_natoms = all(
                dim == natoms
                for dim, (_, natoms) in zip(dimensions, observed)
            )
            if follows_natoms:
                if len({natoms for _, natoms in observed}) == 1:
                    raise LMDBFrameError(
                        f"Field '{disk_name}' axis {axis} has the same length as "
                        "the atom count in every observed frame; it is ambiguous "
                        "without dp_data_shapes or a registered DataType."
                    )
                inferred.append(Axis.NATOMS)
            elif len(set(dimensions)) == 1:
                inferred.append(dimensions[0])
            else:
                raise LMDBFrameError(
                    f"Field '{disk_name}' axis {axis} varies independently of "
                    "the atom count and cannot be represented by one DataType."
                )
        return tuple(inferred)

    @staticmethod
    def _rename_frame_fields(
        frames: list[dict[str, Any]],
        specs: dict[str, _FieldSpec],
    ) -> list[dict[str, Any]]:
        """Map disk field names to dpdata names and validate payload shapes."""
        renamed_frames: list[dict[str, Any]] = []
        for frame_index, frame in enumerate(frames):
            renamed: dict[str, Any] = {
                "atom_types": frame["atom_types"],
                "atom_numbs": frame["atom_numbs"],
            }
            natoms = len(frame["atom_types"])
            for disk_name, spec in specs.items():
                if disk_name not in frame:
                    continue
                if spec.dp_name in renamed:
                    raise LMDBFrameError(
                        f"Frame {frame_index} contains duplicate field "
                        f"'{spec.dp_name}'."
                    )
                value = np.asarray(frame[disk_name])
                LMDBFormat._validate_payload_shape(
                    spec, value, natoms, frame_index
                )
                renamed[spec.dp_name] = value
            renamed_frames.append(renamed)
        return renamed_frames

    @staticmethod
    def _validate_payload_shape(
        spec: _FieldSpec,
        value: np.ndarray,
        natoms: int,
        frame_index: int,
    ) -> None:
        """Validate a decoded per-frame payload against its field specification."""
        expected = list(spec.shape)
        if spec.frame_axis is not None:
            del expected[spec.frame_axis]
        if value.ndim != len(expected):
            raise LMDBFrameError(
                f"Frame {frame_index} field '{spec.disk_name}' has shape "
                f"{value.shape}, incompatible with {spec.shape}."
            )
        for axis, (actual, dimension) in enumerate(zip(value.shape, expected)):
            if dimension is Axis.NATOMS and actual != natoms:
                raise LMDBFrameError(
                    f"Frame {frame_index} field '{spec.disk_name}' axis {axis} "
                    f"must have {natoms} atoms, got {actual}."
                )
            if isinstance(dimension, int) and dimension != -1 and actual != dimension:
                raise LMDBFrameError(
                    f"Frame {frame_index} field '{spec.disk_name}' axis {axis} "
                    f"must have length {dimension}, got {actual}."
                )

    def _group_frames(
        self,
        frames: list[dict[str, Any]],
        names: list[str],
        mixed_type: bool,
        specs_by_disk: dict[str, _FieldSpec],
    ) -> list[dict[str, Any]]:
        """Canonicalize atom order and group frames by composition."""
        specs = {spec.dp_name: spec for spec in specs_by_disk.values()}
        canonical_frames: list[dict[str, Any]] = []
        groups: OrderedDict[tuple[int, ...], list[int]] = OrderedDict()
        for index, frame in enumerate(frames):
            canonical = self._canonicalize_frame(frame, specs)
            canonical_frames.append(canonical)
            composition = tuple(
                np.bincount(
                    canonical["atom_types"], minlength=len(names)
                ).tolist()
            )
            groups.setdefault(composition, []).append(index)

        results: list[dict[str, Any]] = []
        for composition, indices in groups.items():
            self._validate_group_schema(canonical_frames, indices)
            global_atom_types = canonical_frames[indices[0]]["atom_types"]
            if mixed_type:
                atom_names = list(names)
                atom_types = global_atom_types.copy()
            else:
                present = [index for index, count in enumerate(composition) if count]
                atom_names = [names[index] for index in present]
                remap = np.full(len(names), -1, dtype=int)
                remap[present] = np.arange(len(present))
                atom_types = remap[global_atom_types]
            data: dict[str, Any] = {
                "atom_names": atom_names,
                "atom_numbs": np.bincount(
                    atom_types, minlength=len(atom_names)
                ).tolist(),
                "atom_types": atom_types,
                "orig": np.zeros(3),
            }
            self._aggregate_group(data, canonical_frames, indices, specs)
            results.append(data)
        return results

    @staticmethod
    def _canonicalize_frame(
        frame: dict[str, Any], specs: dict[str, _FieldSpec]
    ) -> dict[str, Any]:
        """Sort atoms stably by global type and reorder every atomic field."""
        atom_types = np.asarray(frame["atom_types"], dtype=np.int64)
        permutation = np.argsort(atom_types, kind="stable")
        canonical: dict[str, Any] = {
            "atom_types": atom_types[permutation],
            "atom_numbs": frame["atom_numbs"],
        }
        for name, value in frame.items():
            if name in {"atom_types", "atom_numbs"}:
                continue
            spec = specs[name]
            reordered = np.asarray(value)
            for axis in spec.payload_atom_axes:
                reordered = np.take(reordered, permutation, axis=axis)
            canonical[name] = reordered
        return canonical

    @staticmethod
    def _validate_group_schema(
        frames: list[dict[str, Any]], indices: list[int]
    ) -> None:
        """Reject label or PBC inconsistencies that a System cannot represent."""
        expected = set(frames[indices[0]]) - {"atom_types", "atom_numbs"}
        for index in indices[1:]:
            actual = set(frames[index]) - {"atom_types", "atom_numbs"}
            if actual != expected:
                missing = sorted(expected - actual)
                extra = sorted(actual - expected)
                raise LMDBFrameError(
                    "Frames with the same composition must contain identical "
                    f"fields; frame {index} is missing {missing} and adds {extra}."
                )

    @staticmethod
    def _aggregate_group(
        data: dict[str, Any],
        frames: list[dict[str, Any]],
        indices: list[int],
        specs: dict[str, _FieldSpec],
    ) -> None:
        """Aggregate canonical per-frame payloads into one System data dictionary."""
        fields = set(frames[indices[0]]) - {"atom_types", "atom_numbs"}
        for name in fields:
            values = [np.asarray(frames[index][name]) for index in indices]
            data[name] = LMDBFormat._aggregate_values(values, specs[name])
        if "cells" not in fields:
            data["nopbc"] = True
            data["cells"] = np.zeros((len(indices), 3, 3))

    @staticmethod
    def _aggregate_values(
        values: list[np.ndarray], spec: _FieldSpec
    ) -> np.ndarray:
        """Aggregate values while preserving byte order and static semantics."""
        first = values[0]
        for value in values[1:]:
            if value.dtype != first.dtype or value.shape != first.shape:
                raise LMDBFrameError(
                    f"Field '{spec.disk_name}' has inconsistent dtype or shape "
                    "within one composition."
                )
        frame_axis = spec.frame_axis
        if frame_axis is None:
            if any(not np.array_equal(first, value, equal_nan=True) for value in values[1:]):
                raise LMDBFrameError(
                    f"Static field '{spec.disk_name}' differs between frames."
                )
            return first.copy()
        output_shape = list(first.shape)
        output_shape.insert(frame_axis, len(values))
        output = np.empty(output_shape, dtype=first.dtype)
        for index, value in enumerate(values):
            selection: list[slice | int] = [slice(None)] * output.ndim
            selection[frame_axis] = index
            output[tuple(selection)] = value
        return output

    @staticmethod
    def _register_specs(
        datasets: list[dict[str, Any]], specs: dict[str, _FieldSpec]
    ) -> None:
        """Validate all schemas, then commit process-global registrations."""
        from dpdata.system import LabeledSystem, System

        existing = {dt.name: dt for dt in (*System.DTYPES, *LabeledSystem.DTYPES)}
        protocol_dtypes = {dtype.name: dtype for dtype in _PROTOCOL_DTYPES}
        used_names = {
            name
            for data in datasets
            for name in data
            if name in specs
        }
        pending: list[DataType] = []
        for name in sorted(used_names):
            if name not in specs:
                continue
            spec = specs[name]
            current = existing.get(name)
            if current is not None:
                merged_shape = LMDBFormat._merge_shapes(
                    current.shape, spec.shape
                )
                protocol_dtype = protocol_dtypes.get(name)
                if (
                    merged_shape is None
                    and protocol_dtype is not None
                    and LMDBFormat._shapes_compatible(
                        current.shape, protocol_dtype.shape
                    )
                    and LMDBFormat._shapes_compatible(
                        spec.shape, protocol_dtype.shape
                    )
                ):
                    merged_shape = protocol_dtype.shape
                if (
                    merged_shape is None
                    or current.deepmd_name != spec.deepmd_name
                ):
                    raise LMDBError(
                        f"Data type '{name}' conflicts with the process-global "
                        f"definition: existing shape/name "
                        f"{current.shape}/{current.deepmd_name}, LMDB shape/name "
                        f"{spec.shape}/{spec.deepmd_name}."
                    )
                if merged_shape != current.shape:
                    dtype = DataType(
                        name,
                        np.ndarray,
                        shape=merged_shape,
                        required=False,
                        deepmd_name=spec.deepmd_name,
                    )
                    pending.append(dtype)
                    existing[name] = dtype
                continue
            dtype = DataType(
                name,
                np.ndarray,
                shape=spec.shape,
                required=False,
                deepmd_name=spec.deepmd_name,
            )
            pending.append(dtype)
            existing[name] = dtype
        if pending:
            replacements = {dtype.name: dtype for dtype in pending}
            for cls in (System, LabeledSystem):
                current = {dtype.name: dtype for dtype in cls.DTYPES}
                current.update(replacements)
                cls.DTYPES = tuple(current.values())

    @staticmethod
    def _shapes_compatible(
        first: tuple[int | Axis, ...] | None,
        second: tuple[int | Axis, ...] | None,
    ) -> bool:
        """Return whether two shapes differ only through wildcard dimensions."""
        if first is None or second is None or len(first) != len(second):
            return first == second
        return all(
            left == right
            or left == -1
            or right == -1
            for left, right in zip(first, second)
        )

    @staticmethod
    def _merge_shapes(
        first: tuple[int | Axis, ...] | None,
        second: tuple[int | Axis, ...] | None,
    ) -> tuple[int | Axis, ...] | None:
        """Return the least restrictive compatible shape."""
        if not LMDBFormat._shapes_compatible(first, second):
            return None
        if first is None or second is None:
            return first
        return tuple(
            -1 if left == -1 or right == -1 else left
            for left, right in zip(first, second)
        )


def dump_systems(
    systems,
    file_name,
    type_map: list[str] | None = None,
    map_size: int = DEFAULT_MAP_SIZE,
    frame_idx_fmt: str = DEFAULT_FRAME_FMT,
    write_batch_size: int = DEFAULT_WRITE_BATCH_SIZE,
    overwrite: bool = False,
):
    """Write an ordered sequence of systems to one LMDB, preserving identity.

    Each element of ``systems`` is stored as a distinct source system,
    numbered ``0, 1, 2, ...`` in iteration order, and recorded in the
    ``frame_system_ids`` metadata. In contrast to
    ``MultiSystems.to('lmdb', ...)``, systems are not merged by formula, so
    the system partition used by DeePMD-kit's ``prob_sys_size`` is kept.

    Parameters
    ----------
    systems : Iterable[System]
        An ordered iterable of :class:`dpdata.System` /
        :class:`dpdata.LabeledSystem` objects or raw system data dictionaries.
        Do not pass a :class:`dpdata.MultiSystems`, whose frames are already
        merged by formula.
    file_name : str
        Output LMDB directory.
    type_map : list[str], optional
        Global element table. If ``None``, the union of the elements of all
        systems is used. Provide it explicitly to stream a generator without
        materialising the whole sequence.
    map_size : int, optional
        Maximum LMDB size in bytes. Default is 1 TiB (sparse).
    frame_idx_fmt : str, optional
        Format used for the per-frame integer key. Default ``"012d"``.
    write_batch_size : int, optional
        Number of frames committed per LMDB write transaction.
    overwrite : bool, optional
        Whether to replace an existing destination after the new database has
        been written and validated. The default is ``False``.

    Examples
    --------
    >>> import dpdata
    >>> from dpdata.formats.lmdb import dump_systems
    >>> systems = [
    ...     dpdata.LabeledSystem(d, fmt="deepmd/npy") for d in directories
    ... ]
    >>> dump_systems(systems, "data.lmdb", type_map=["H", "C", "N", "O"])
    """
    LMDBFormat().dump_systems(
        systems,
        file_name,
        map_size=map_size,
        frame_idx_fmt=frame_idx_fmt,
        type_map=type_map,
        write_batch_size=write_batch_size,
        overwrite=overwrite,
    )
