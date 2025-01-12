from __future__ import annotations

from abc import ABC

# scipy physical constant version 2018
physical_constants = {
    "Avogadro constant": 6.02214076e23,  # mol^-1
    "elementary charge": 1.602176634e-19,  # C
    "atomic unit of length": 5.29177210903e-11,  # m
    "atomic unit of energy": 4.3597447222071e-18,  # J
    "Rydberg constant": 10973731.568160,  # m^-1
    "Planck constant": 6.62607015e-34,  # J·s
    "speed of light in vacuum": 299792458,  # m·s^-1
}


def scipy_constant_value(key: str) -> float:
    return physical_constants[key]


AVOGADRO = scipy_constant_value("Avogadro constant")  # Avagadro constant
ELE_CHG = scipy_constant_value("elementary charge")  # Elementary Charge, in C
BOHR = scipy_constant_value("atomic unit of length")  # Bohr, in m
HARTREE = scipy_constant_value("atomic unit of energy")  # Hartree, in Jole
RYDBERG = (
    scipy_constant_value("Rydberg constant")
    * scipy_constant_value("Planck constant")
    * scipy_constant_value("speed of light in vacuum")
)  # Rydberg, in Jole

# energy conversions
econvs = {
    "eV": 1.0,
    "hartree": HARTREE / ELE_CHG,
    "kJ_mol": 1 / (ELE_CHG * AVOGADRO / 1000),
    "kcal_mol": 1 / (ELE_CHG * AVOGADRO / 1000 / 4.184),
    "rydberg": RYDBERG / ELE_CHG,
    "J": 1 / ELE_CHG,
    "kJ": 1000 / ELE_CHG,
}

# length conversions
lconvs = {
    "angstrom": 1.0,
    "bohr": BOHR * 1e10,
    "nm": 10.0,
    "m": 1e10,
}


def check_unit(unit):
    if unit not in econvs.keys() and unit not in lconvs.keys():
        try:
            eunit = unit.split("/")[0]
            lunit = unit.split("/")[1]
            if eunit not in econvs.keys():
                raise RuntimeError(f"Invaild unit: {unit}")
            if lunit not in lconvs.keys():
                raise RuntimeError(f"Invalid unit: {unit}")
        except Exception:
            raise RuntimeError(f"Invalid unit: {unit}")


class Conversion(ABC):
    def __init__(self, unitA, unitB, check=True):
        """Parent class for unit conversion.

        Parameters
        ----------
        unitA : str
            unit to be converted
        unitB : str
            unit which unitA is converted to, i.e. `1 unitA = self._value unitB`
        check : bool
            whether to check unit validity

        Examples
        --------
        >>> conv = Conversion("foo", "bar", check=False)
        >>> conv.set_value("10.0")
        >>> print(conv)
        1 foo = 10.0 bar
        >>> conv.value()
        10.0
        """
        if check:
            check_unit(unitA)
            check_unit(unitB)
        self.unitA = unitA
        self.unitB = unitB
        self._value = 0.0

    def value(self):
        return self._value

    def set_value(self, value):
        self._value = value

    def __repr__(self):
        return f"1 {self.unitA} = {self._value} {self.unitB}"

    def __str__(self):
        return self.__repr__()


class EnergyConversion(Conversion):
    def __init__(self, unitA, unitB):
        """Class for energy conversion.

        Examples
        --------
        >>> conv = EnergyConversion("eV", "kcal_mol")
        >>> conv.value()
        23.06054783061903
        """
        super().__init__(unitA, unitB)
        self.set_value(econvs[unitA] / econvs[unitB])


class LengthConversion(Conversion):
    def __init__(self, unitA, unitB):
        """Class for length conversion.

        Examples
        --------
        >>> conv = LengthConversion("angstrom", "nm")
        >>> conv.value()
        0.1
        """
        super().__init__(unitA, unitB)
        self.set_value(lconvs[unitA] / lconvs[unitB])


class ForceConversion(Conversion):
    def __init__(self, unitA, unitB):
        """Class for force conversion.

        Parameters
        ----------
        unitA, unitB : str
            in format of "energy_unit/length_unit"

        Examples
        --------
        >>> conv = ForceConversion("kJ_mol/nm", "eV/angstrom")
        >>> conv.value()
        0.0010364269656262175
        """
        super().__init__(unitA, unitB)
        econv = EnergyConversion(unitA.split("/")[0], unitB.split("/")[0]).value()
        lconv = LengthConversion(unitA.split("/")[1], unitB.split("/")[1]).value()
        self.set_value(econv / lconv)


class PressureConversion(Conversion):
    def __init__(self, unitA, unitB):
        """Class for pressure conversion.

        Parameters
        ----------
        unitA, unitB : str
            in format of "energy_unit/length_unit^3", or in `["Pa", "pa", "kPa", "kpa", "bar", "kbar"]`

        Examples
        --------
        >>> conv = PressureConversion("kbar", "eV/angstrom^3")
        >>> conv.value()
        0.0006241509074460763
        """
        super().__init__(unitA, unitB, check=False)
        unitA, factorA = self._convert_unit(unitA)
        unitB, factorB = self._convert_unit(unitB)
        eunitA, lunitA = self._split_unit(unitA)
        eunitB, lunitB = self._split_unit(unitB)
        econv = EnergyConversion(eunitA, eunitB).value() * factorA / factorB
        lconv = LengthConversion(lunitA, lunitB).value()
        self.set_value(econv / lconv**3)

    def _convert_unit(self, unit):
        if unit == "Pa" or unit == "pa":
            return "J/m^3", 1
        elif unit == "kPa" or unit == "kpa":
            return "kJ/m^3", 1
        elif unit == "GPa" or unit == "Gpa":
            return "kJ/m^3", 1e6
        elif unit == "bar":
            return "J/m^3", 1e5
        elif unit == "kbar":
            return "kJ/m^3", 1e5
        else:
            return unit, 1

    def _split_unit(self, unit):
        eunit = unit.split("/")[0]
        lunit = unit.split("/")[1][:-2]
        return eunit, lunit
