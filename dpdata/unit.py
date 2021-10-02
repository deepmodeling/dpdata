from abc import ABC

AVOGADRO = 6.02214179E23    # Avagadro constant
ELE_CHG  = 1.602176487E-19  # Elementary Charge, in C

# energy conversions
econvs = {
    "eV": 1.0,
    "hartree": 4.35974394E-18 / ELE_CHG,
    "kJ_mol": 1 / (ELE_CHG * AVOGADRO / 1000),
    "kcal_mol": 1 / (ELE_CHG * AVOGADRO / 1000 / 4.184),
    "rydberg": 2.17987197E-18 / ELE_CHG,
    "J": 1 / ELE_CHG,
    "kJ": 1000 / ELE_CHG
}

# length conversions
lconvs = {
    "angstrom": 1.0,
    "bohr": 0.52917720859,
    "nm": 10.0,
    "m": 1E10,
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
        except:
            raise RuntimeError(f"Invalid unit: {unit}")

class Conversion(ABC):
    def __init__(self, unitA, unitB, check=True):
        """
        1 unitA = self.value unitB
        """
        if check:
            check_unit(unitA)
            check_unit(unitB)
        self.unitA = unitA
        self.unitB = unitB
        self._value = 0.0
    
    def value(self):
        return self._value
    
    def __repr__(self):
        return f"1 {self.unitA} = {self._value} {self.unitB}"
    
    def __str__(self):
        return self.__repr__()

class EnergyConversion(Conversion):
    def __init__(self, unitA, unitB):
        super().__init__(unitA, unitB)
        self._value = econvs[unitA] / econvs[unitB]

class LengthConversion(Conversion):
    def __init__(self, unitA, unitB):
        super().__init__(unitA, unitB)
        self._value = lconvs[unitA] / lconvs[unitB]

class ForceConversion(Conversion):
    def __init__(self, unitA, unitB):
        super().__init__(unitA, unitB)
        econv = EnergyConversion(unitA.split("/")[0], unitB.split("/")[0]).value()
        lconv = LengthConversion(unitA.split("/")[1], unitB.split("/")[1]).value()
        self._value = econv / lconv

class VirialConversion(Conversion):
    def __init__(self, unitA, unitB):
        super().__init__(unitA, unitB, check=False)
        unitA, factorA = self._convert_unit(unitA)
        unitB, factorB = self._convert_unit(unitB)
        eunitA, lunitA = self._split_unit(unitA)
        eunitB, lunitB = self._split_unit(unitB)
        econv = EnergyConversion(eunitA, eunitB).value() * factorA / factorB
        lconv = LengthConversion(lunitA, lunitB).value()
        self._value = econv / lconv**3
    
    def _convert_unit(self, unit):
        if unit == "Pa" or unit == "pa":
            return "J/m^3", 1
        elif unit == "kPa" or unit == "kpa":
            return "kJ/m^3", 1
        elif unit == "bar":
            return "J/m^3", 1E5
        elif unit == "kbar":
            return "kJ/m^3", 1E5
        else:
            return unit, 1

    def _split_unit(self, unit):
        eunit = unit.split("/")[0]
        lunit = unit.split("/")[1][:-2]
        return eunit, lunit



