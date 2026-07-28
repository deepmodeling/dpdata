from __future__ import annotations

from dpdata.format import Format
from dpdata.formats.deepmd.lmdb.format import LMDBFormat

# Canonical name; ``lmdb`` is kept as a backward-compatible alias.
Format.register("deepmd/lmdb")(LMDBFormat)
Format.register("lmdb")(LMDBFormat)
