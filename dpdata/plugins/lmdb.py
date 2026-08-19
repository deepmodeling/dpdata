from __future__ import annotations

from dpdata.format import Format

try:
    from dpdata.formats.deepmd.lmdb.format import LMDBFormat
except ModuleNotFoundError as exc:
    if exc.name not in {"lmdb", "msgpack"}:
        raise

    _lmdb_import_error = exc

    class LMDBFormat(Format):
        """Placeholder used when the LMDB backend dependencies are unavailable."""

        def __init__(self, *args, **kwargs) -> None:
            raise ModuleNotFoundError(
                "The deepmd/lmdb format requires the 'lmdb' and 'msgpack' packages."
            ) from _lmdb_import_error


# Canonical name; ``lmdb`` is kept as a backward-compatible alias.
Format.register("deepmd/lmdb")(LMDBFormat)
Format.register("lmdb")(LMDBFormat)
