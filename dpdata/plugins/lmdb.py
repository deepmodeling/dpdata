from __future__ import annotations

from dpdata.format import Format
from dpdata.formats.lmdb.format import LMDBFormat

Format.register("lmdb")(LMDBFormat)
