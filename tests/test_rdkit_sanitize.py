from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path
from unittest.mock import patch

MODULE_PATH = Path(__file__).resolve().parents[1] / "dpdata" / "rdkit" / "sanitize.py"
SPEC = importlib.util.spec_from_file_location("dpdata_rdkit_sanitize", MODULE_PATH)
sanitize = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(sanitize)


class _Bond:
    def __init__(self, order):
        self.order = order

    def GetBondTypeAsDouble(self):
        return self.order


class _AtomNewAPI:
    def __init__(self, explicit_valence=4, bond_orders=None):
        self.explicit_valence = explicit_valence
        self.bond_orders = bond_orders or [1, 1, 1, 1]

    def GetBonds(self):
        return [_Bond(order) for order in self.bond_orders]

    def GetValence(self, valence_type):
        self.valence_type = valence_type
        return self.explicit_valence

    def GetExplicitValence(self):
        raise AssertionError("legacy API should not be used when new API is available")

    def GetSymbol(self):
        return "C"

    def GetIdx(self):
        return 0


class _AtomOldAPI:
    def __init__(self, explicit_valence=4, bond_orders=None):
        self.explicit_valence = explicit_valence
        self.bond_orders = bond_orders or [1, 1, 1, 1]

    def GetBonds(self):
        return [_Bond(order) for order in self.bond_orders]

    def GetExplicitValence(self):
        return self.explicit_valence

    def GetSymbol(self):
        return "C"

    def GetIdx(self):
        return 0


class _ChemNewAPI:
    class ValenceType:
        EXPLICIT = object()


class TestGetExplicitValence(unittest.TestCase):
    def test_prefers_new_rdkit_valence_api(self):
        atom = _AtomNewAPI(explicit_valence=4)
        with patch.dict(
            "sys.modules", {"rdkit": type("_Rdkit", (), {"Chem": _ChemNewAPI})}
        ):
            self.assertEqual(sanitize.get_explicit_valence(atom), 4)
            self.assertIs(atom.valence_type, _ChemNewAPI.ValenceType.EXPLICIT)

    def test_falls_back_to_legacy_api_when_new_api_missing(self):
        atom = _AtomOldAPI(explicit_valence=4)
        with patch.dict(
            "sys.modules", {"rdkit": type("_Rdkit", (), {"Chem": object()})}
        ):
            self.assertEqual(sanitize.get_explicit_valence(atom), 4)


if __name__ == "__main__":
    unittest.main()
