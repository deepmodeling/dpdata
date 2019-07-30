from pathlib import Path
from monty.serialization import loadfn,dumpfn

fpdt=str(Path(__file__).absolute().parent / "periodic_table.json")
_pdt=loadfn(fpdt)

class Element:

    def __init__(self, symbol: str):
        assert symbol in _pdt.keys()
        self.symbol = "%s" % symbol
        d = _pdt[symbol]
        self._Z = d['atomic_no']
        self._name = d['name']
        self._X = d['X']
        self._mass = d['atomic_mass']
        self._r = d['radius']
        self._cr = d["calculated_radius"]

    def __str__(self):
        return self.symbol

    def __repr__(self):
        return "Element : %s"%self.symbol

    @property
    def Z(self):
        return self._Z

    @property
    def name(self):
        return self._name

    @property
    def X(self):
        return self._X

    @property
    def mass(self):
        return self._mass

    @property
    def radius(self):
        return self._r

    @property
    def calculated_radius(self):
        return self._cr
