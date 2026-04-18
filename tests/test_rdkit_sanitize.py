from __future__ import annotations

import importlib.util
import os
import sys
import unittest
from unittest.mock import patch


_SKIP_REASON = None
if importlib.util.find_spec("rdkit") is None:
    _SKIP_REASON = "requires rdkit"
    sanitize = None
else:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    try:
        from dpdata.rdkit import sanitize
    except ModuleNotFoundError as exc:
        _SKIP_REASON = f"missing test dependency: {exc.name}"
        sanitize = None


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


class _ChemOldAPI:
    pass


class _RdkitModule:
    def __init__(self, chem):
        self.Chem = chem


@unittest.skipIf(_SKIP_REASON is not None, _SKIP_REASON or "skip")
class TestGetExplicitValence(unittest.TestCase):
    def test_prefers_new_rdkit_valence_api(self):
        atom = _AtomNewAPI(explicit_valence=4)
        fake_rdkit = _RdkitModule(_ChemNewAPI)
        with patch.dict(
            "sys.modules", {"rdkit": fake_rdkit, "rdkit.Chem": _ChemNewAPI}
        ):
            self.assertEqual(sanitize.get_explicit_valence(atom), 4)
            self.assertIs(atom.valence_type, _ChemNewAPI.ValenceType.EXPLICIT)

    def test_falls_back_to_legacy_api_when_new_api_missing(self):
        atom = _AtomOldAPI(explicit_valence=4)
        fake_rdkit = _RdkitModule(_ChemOldAPI)
        with patch.object(
            atom, "GetExplicitValence", wraps=atom.GetExplicitValence
        ) as mock_get_explicit_valence:
            with patch.dict(
                "sys.modules", {"rdkit": fake_rdkit, "rdkit.Chem": _ChemOldAPI}
            ):
                self.assertEqual(sanitize.get_explicit_valence(atom), 4)
            mock_get_explicit_valence.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()
