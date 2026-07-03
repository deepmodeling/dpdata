"""Unit conversion helpers for extended XYZ (extxyz) format.

This module provides a table-driven approach to convert energy, force, and
stress/pressure values from various unit systems commonly found in extxyz
files into dpdata's internal units:

- Energy: eV
- Force: eV/angstrom
- Stress/Pressure: eV/angstrom^3 (before volume multiplication to get virial)
"""

from __future__ import annotations

from dpdata.unit import EnergyConversion, ForceConversion, PressureConversion

# ---------------------------------------------------------------------------
# Unit alias mapping tables
# Keys are LOWERCASE; values are canonical names recognized by dpdata.unit
# ---------------------------------------------------------------------------

_ENERGY_UNIT_MAP: dict[str, str] = {
    "ev": "eV",
    "hartree": "hartree",
    "ha": "hartree",
    "ry": "rydberg",
    "rydberg": "rydberg",
    "kcal/mol": "kcal_mol",
    "kcal_mol": "kcal_mol",
    "kj/mol": "kJ_mol",
    "kj_mol": "kJ_mol",
}

_LENGTH_UNIT_MAP: dict[str, str] = {
    "angstrom": "angstrom",
    "ang": "angstrom",
    "ang.": "angstrom",
    "a": "angstrom",
    "bohr": "bohr",
    "nm": "nm",
}

_PRESSURE_UNIT_MAP: dict[str, str] = {
    "gpa": "GPa",
    "kbar": "kbar",
    "bar": "bar",
    "ev/angstrom^3": "eV/angstrom^3",
    "ev/ang^3": "eV/angstrom^3",
    "ev/a^3": "eV/angstrom^3",
    "ha/bohr^3": "hartree/bohr^3",
    "hartree/bohr^3": "hartree/bohr^3",
}

# dpdata internal unit strings
_INTERNAL_ENERGY = "eV"
_INTERNAL_FORCE = "eV/angstrom"
_INTERNAL_PRESSURE = "eV/angstrom^3"


def _parse_force_unit(raw: str) -> tuple[str, str]:
    """Split a composite force unit string into (energy_part, length_part).

    Examples
    --------
    >>> _parse_force_unit("kcal/mol/angstrom")
    ('kcal/mol', 'angstrom')
    >>> _parse_force_unit("hartree/bohr")
    ('hartree', 'bohr')
    >>> _parse_force_unit("ev/ang")
    ('ev', 'ang')
    """
    # Try matching known energy prefixes (longest first) to handle
    # composite names like "kcal/mol" that themselves contain "/".
    for e_key in sorted(_ENERGY_UNIT_MAP.keys(), key=len, reverse=True):
        prefix = e_key + "/"
        if raw.startswith(prefix):
            l_part = raw[len(prefix) :]
            if l_part:
                return e_key, l_part
    # Fallback: split on last "/"
    parts = raw.rsplit("/", 1)
    if len(parts) == 2 and parts[0] and parts[1]:
        return parts[0], parts[1]
    raise ValueError(f"Cannot parse force unit string: '{raw}'")


def _get_unit_factor(unit_str: str | None, quantity: str) -> float:
    """Return the multiplicative factor to convert from the given unit to dpdata internals.

    Parameters
    ----------
    unit_str : str or None
        The unit string read from the extxyz header (e.g. "hartree", "kcal/mol/angstrom").
        If None, returns 1.0 (assumes data is already in internal units).
    quantity : str
        One of "energy", "force", or "stress".

    Returns
    -------
    float
        Conversion factor such that ``value_internal = value_file * factor``.

    Raises
    ------
    ValueError
        If the unit string is not recognized or the quantity type is invalid.
    """
    if unit_str is None:
        return 1.0

    key = unit_str.lower().strip()

    if quantity == "energy":
        canonical = _ENERGY_UNIT_MAP.get(key)
        if canonical is None:
            raise ValueError(
                f"Unsupported energy unit: '{unit_str}'. "
                f"Supported: {list(_ENERGY_UNIT_MAP.keys())}"
            )
        return EnergyConversion(canonical, _INTERNAL_ENERGY).value()

    elif quantity == "force":
        e_part, l_part = _parse_force_unit(key)
        e_canonical = _ENERGY_UNIT_MAP.get(e_part)
        l_canonical = _LENGTH_UNIT_MAP.get(l_part)
        if e_canonical is None:
            raise ValueError(
                f"Unsupported energy part in force unit: '{e_part}' "
                f"(from '{unit_str}'). Supported: {list(_ENERGY_UNIT_MAP.keys())}"
            )
        if l_canonical is None:
            raise ValueError(
                f"Unsupported length part in force unit: '{l_part}' "
                f"(from '{unit_str}'). Supported: {list(_LENGTH_UNIT_MAP.keys())}"
            )
        src_unit = f"{e_canonical}/{l_canonical}"
        return ForceConversion(src_unit, _INTERNAL_FORCE).value()

    elif quantity == "stress":
        canonical = _PRESSURE_UNIT_MAP.get(key)
        if canonical is None:
            raise ValueError(
                f"Unsupported stress/pressure unit: '{unit_str}'. "
                f"Supported: {list(_PRESSURE_UNIT_MAP.keys())}"
            )
        return PressureConversion(canonical, _INTERNAL_PRESSURE).value()

    else:
        raise ValueError(
            f"Unknown quantity type: '{quantity}'. Must be 'energy', 'force', or 'stress'."
        )
