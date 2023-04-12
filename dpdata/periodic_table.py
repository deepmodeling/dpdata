from pathlib import Path

from monty.serialization import loadfn

fpdt = str(Path(__file__).absolute().parent / "periodic_table.json")
_pdt = loadfn(fpdt)
ELEMENTS = [
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
]


class Element:
    def __init__(self, symbol: str):
        assert symbol in ELEMENTS
        self.symbol = "%s" % symbol
        d = _pdt[symbol]
        self._Z = d["atomic_no"]
        self._name = d["name"]
        self._X = d["X"]
        self._mass = d["atomic_mass"]
        self._r = d["radius"]
        self._cr = d["calculated_radius"]

    def __str__(self):
        return self.symbol

    def __repr__(self):
        return "Element : %s" % self.symbol

    @classmethod
    def from_Z(cls, Z):
        assert Z > 0
        assert Z < len(ELEMENTS)
        return cls(ELEMENTS[Z - 1])

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
