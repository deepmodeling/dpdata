from __future__ import annotations

from dpdata.format import Format
from dpdata.lmdb.format import LMDBFormat

Format.register("lmdb")(LMDBFormat)
