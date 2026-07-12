from __future__ import annotations

import tempfile
import unittest
import warnings

import numpy as np
from comp_sys import CompLabeledSys, IsPBC
from context import dpdata


class TestQuipGapxyz1(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.multi_systems = dpdata.MultiSystems.from_file(
            "xyz/xyz_unittest.xyz", "quip/gap/xyz"
        )
        self.system_1 = self.multi_systems.systems["B1C9"]
        self.system_2 = dpdata.LabeledSystem("xyz/B1C9", fmt="deepmd")
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


class TestQuipGapxyz2(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_temp0 = dpdata.MultiSystems.from_file(
            file_name="xyz/xyz_unittest.xyz", fmt="quip/gap/xyz"
        )
        self.system_1 = self.system_temp0.systems["B5C7"]  # .sort_atom_types()
        self.system_temp1 = dpdata.LabeledSystem("xyz/B1C9", fmt="deepmd")
        self.system_temp2 = dpdata.LabeledSystem("xyz/B5C7", fmt="deepmd")
        self.system_temp3 = dpdata.MultiSystems(self.system_temp2, self.system_temp1)
        self.system_2 = self.system_temp3.systems["B5C7"]
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


class TestQuipGapxyzsort1(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.multi_systems_1 = dpdata.MultiSystems.from_file(
            "xyz/xyz_unittest.sort.xyz", "quip/gap/xyz"
        )
        self.system_1 = self.multi_systems_1.systems["B5C7"]
        self.system_1.sort_atom_types()
        self.multi_systems_2 = dpdata.MultiSystems.from_file(
            "xyz/xyz_unittest.xyz", "quip/gap/xyz"
        )
        self.system_2 = self.multi_systems_2.systems["B5C7"]
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


class TestQuipGapxyzsort2(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.multi_systems_1 = dpdata.MultiSystems.from_file(
            "xyz/xyz_unittest.sort.xyz", "quip/gap/xyz"
        )
        self.system_1 = self.multi_systems_1.systems["B1C9"]
        self.system_1.sort_atom_types()
        self.multi_systems_2 = dpdata.MultiSystems.from_file(
            "xyz/xyz_unittest.xyz", "quip/gap/xyz"
        )
        self.system_2 = self.multi_systems_2.systems["B1C9"]
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


class TestQuipGapxyzfield(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.multi_systems_1 = dpdata.MultiSystems.from_file(
            "xyz/xyz_unittest.field.xyz", "quip/gap/xyz"
        )
        self.system_1 = self.multi_systems_1.systems["B1C9"]
        self.system_1.sort_atom_types()
        self.multi_systems_2 = dpdata.MultiSystems.from_file(
            "xyz/xyz_unittest.xyz", "quip/gap/xyz"
        )
        self.system_2 = self.multi_systems_2.systems["B1C9"]
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


class TestQuipGapxyzfield2(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.multi_systems_1 = dpdata.MultiSystems.from_file(
            "xyz/xyz_unittest.field.xyz", "quip/gap/xyz"
        )
        self.system_1 = self.multi_systems_1.systems["B5C7"]
        self.system_1.sort_atom_types()
        self.multi_systems_2 = dpdata.MultiSystems.from_file(
            "xyz/xyz_unittest.xyz", "quip/gap/xyz"
        )
        self.system_2 = self.multi_systems_2.systems["B5C7"]
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


class TestQuipGapxyzNoVirials(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.multi_systems_1 = dpdata.MultiSystems.from_file(
            "xyz/xyz_B5C7_novirials.xyz", fmt="quip/gap/xyz"
        )
        self.system_1 = self.multi_systems_1.systems["B5C7"]
        self.system_1.sort_atom_types()
        self.system_2 = dpdata.LabeledSystem("xyz/B5C7_novirials", fmt="deepmd/raw")
        self.places = 6
        self.e_places = 6
        self.f_places = 6


# ---------- stress / virial conversion (fixes #973) ----------


class TestStressToVirial(unittest.TestCase):
    """Read extxyz with stress= header → virials via -V*stress."""

    def setUp(self):
        self.ms = dpdata.MultiSystems.from_file("xyz/stress_only.xyz", fmt="extxyz")
        self.system = list(self.ms.systems.values())[0]

    def test_has_virials(self):
        self.assertIn("virials", self.system.data)

    def test_virial_values(self):
        """Virial = -V * stress, V=27, stress=diag(0.01,0.02,0.03)."""
        expected = np.array([[[-0.27, 0, 0], [0, -0.54, 0], [0, 0, -0.81]]])
        np.testing.assert_allclose(self.system.data["virials"], expected, atol=1e-10)

    def test_energy(self):
        np.testing.assert_allclose(self.system.data["energies"], [-1.5])

    def test_forces(self):
        expected = np.array([[[0.1, 0, 0], [-0.1, 0, 0], [0, 0, 0]]])
        np.testing.assert_allclose(self.system.data["forces"], expected)


class TestStressVoigt(unittest.TestCase):
    """Read extxyz with 6-component Voigt stress."""

    def setUp(self):
        self.ms = dpdata.MultiSystems.from_file("xyz/stress_voigt.xyz", fmt="extxyz")
        self.system = list(self.ms.systems.values())[0]

    def test_has_virials(self):
        self.assertIn("virials", self.system.data)

    def test_virial_values(self):
        """Voigt stress=[0.01,0.02,0.03,0.004,0.005,0.006] → 3×3 → virial=-V*stress."""
        stress = np.array(
            [[0.01, 0.006, 0.005], [0.006, 0.02, 0.004], [0.005, 0.004, 0.03]]
        )
        expected = np.array([-27.0 * stress])
        np.testing.assert_allclose(self.system.data["virials"], expected, atol=1e-10)


class TestStressSignPositive(unittest.TestCase):
    """Test stress_sign=1 (opposite convention)."""

    def setUp(self):
        self.ms = dpdata.MultiSystems.from_file(
            "xyz/stress_only.xyz", fmt="extxyz", stress_sign=1
        )
        self.system = list(self.ms.systems.values())[0]

    def test_virial_values_positive(self):
        """Virial = +V * stress with stress_sign=1."""
        expected = np.array([[[0.27, 0, 0], [0, 0.54, 0], [0, 0, 0.81]]])
        np.testing.assert_allclose(self.system.data["virials"], expected, atol=1e-10)


class TestVirialsKey(unittest.TestCase):
    """Read extxyz with 'virials' (plural) instead of 'virial'."""

    def setUp(self):
        self.ms = dpdata.MultiSystems.from_file("xyz/virials_key.xyz", fmt="extxyz")
        self.system = list(self.ms.systems.values())[0]

    def test_has_virials(self):
        self.assertIn("virials", self.system.data)

    def test_virial_values(self):
        expected = np.array([[[0.27, 0, 0], [0, 0.54, 0], [0, 0, 0.81]]])
        np.testing.assert_allclose(self.system.data["virials"], expected, atol=1e-10)


class TestStressesKey(unittest.TestCase):
    """Read extxyz with 'stresses' (plural) instead of 'stress'."""

    def setUp(self):
        self.ms = dpdata.MultiSystems.from_file("xyz/stresses_key.xyz", fmt="extxyz")
        self.system = list(self.ms.systems.values())[0]

    def test_has_virials(self):
        self.assertIn("virials", self.system.data)

    def test_virial_values(self):
        """Same as stress_only: virial = -27 * diag(0.01,0.02,0.03)."""
        expected = np.array([[[-0.27, 0, 0], [0, -0.54, 0], [0, 0, -0.81]]])
        np.testing.assert_allclose(self.system.data["virials"], expected, atol=1e-10)


# ---------- robustness ----------


class TestUnknownProperties(unittest.TestCase):
    """Read extxyz with extra per-atom props (magmom) — should warn, not crash."""

    def test_parses_without_crash(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ms = dpdata.MultiSystems.from_file("xyz/unknown_props.xyz", fmt="extxyz")
            system = list(ms.systems.values())[0]
            self.assertEqual(system.get_nframes(), 1)
            warn_msgs = [str(x.message) for x in w]
            self.assertTrue(
                any("magmom" in msg for msg in warn_msgs),
                f"Expected warning about 'magmom', got: {warn_msgs}",
            )

    def test_forces_correct(self):
        ms = dpdata.MultiSystems.from_file("xyz/unknown_props.xyz", fmt="extxyz")
        system = list(ms.systems.values())[0]
        expected = np.array([[[0.1, 0, 0], [-0.1, 0, 0], [0, 0, 0]]])
        np.testing.assert_allclose(system.data["forces"], expected)


class TestForcesKey(unittest.TestCase):
    """Read extxyz using 'forces' (ASE style) instead of 'force'."""

    def setUp(self):
        self.ms = dpdata.MultiSystems.from_file("xyz/forces_key.xyz", fmt="extxyz")
        self.system = list(self.ms.systems.values())[0]

    def test_forces_parsed(self):
        expected = np.array([[[0.1, 0, 0], [-0.1, 0, 0], [0, 0, 0]]])
        np.testing.assert_allclose(self.system.data["forces"], expected)


class TestNoPBC(unittest.TestCase):
    """Read extxyz without Lattice → nopbc=True, dummy 100Å box."""

    def setUp(self):
        self.ms = dpdata.MultiSystems.from_file("xyz/nopbc.xyz", fmt="extxyz")
        self.system = list(self.ms.systems.values())[0]

    def test_nopbc(self):
        self.assertTrue(self.system.data.get("nopbc", False))

    def test_dummy_cell(self):
        expected = np.diag([100.0, 100.0, 100.0])
        np.testing.assert_allclose(self.system.data["cells"][0], expected)

    def test_energy(self):
        np.testing.assert_allclose(self.system.data["energies"], [-1.5])


class TestEnergyKeyVariants(unittest.TestCase):
    """Read extxyz with 'Energy' instead of 'energy'."""

    def setUp(self):
        self.ms = dpdata.MultiSystems.from_file(
            "xyz/energy_key_variants.xyz", fmt="extxyz"
        )
        self.system = list(self.ms.systems.values())[0]

    def test_energy_parsed(self):
        np.testing.assert_allclose(self.system.data["energies"], [-1.5])


class TestPBCFalseHeader(unittest.TestCase):
    """Read extxyz with pbc='F F F' → nopbc=True even with Lattice present."""

    def setUp(self):
        self.ms = dpdata.MultiSystems.from_file("xyz/pbc_false.xyz", fmt="extxyz")
        self.system = list(self.ms.systems.values())[0]

    def test_nopbc(self):
        self.assertTrue(self.system.data.get("nopbc", False))

    def test_cell_still_parsed(self):
        """Lattice should still be read even if pbc=F."""
        expected = np.diag([3.0, 3.0, 3.0])
        np.testing.assert_allclose(self.system.data["cells"][0], expected)


# ---------- writer ----------


class TestWriteForcesKey(unittest.TestCase):
    """Writer should output 'forces' (not 'force') for ASE compat."""

    def test_output_has_forces_key(self):
        ms = dpdata.MultiSystems.from_file("xyz/stress_only.xyz", fmt="extxyz")
        system = list(ms.systems.values())[0]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
            system.to("extxyz", f.name)
            f.flush()
            with open(f.name) as fread:
                content = fread.read()
        self.assertIn("forces:R:3", content)
        self.assertNotIn("force:R:3", content)


class TestWriteStressField(unittest.TestCase):
    """Writer should output stress= alongside virial= when virials present."""

    def test_output_has_stress(self):
        ms = dpdata.MultiSystems.from_file("xyz/xyz_unittest.xyz", fmt="quip/gap/xyz")
        system = list(ms.systems.values())[0]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
            system.to("extxyz", f.name)
            f.flush()
            with open(f.name) as fread:
                content = fread.read()
        self.assertIn("virial=", content)
        self.assertIn("stress=", content)


class TestWritePBC(unittest.TestCase):
    """Writer should output pbc field."""

    def test_pbc_true(self):
        ms = dpdata.MultiSystems.from_file("xyz/xyz_unittest.xyz", fmt="quip/gap/xyz")
        system = list(ms.systems.values())[0]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
            system.to("extxyz", f.name)
            f.flush()
            with open(f.name) as fread:
                content = fread.read()
        self.assertIn('pbc="T T T"', content)

    def test_pbc_false(self):
        ms = dpdata.MultiSystems.from_file("xyz/nopbc.xyz", fmt="extxyz")
        system = list(ms.systems.values())[0]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
            system.to("extxyz", f.name)
            f.flush()
            with open(f.name) as fread:
                content = fread.read()
        self.assertIn('pbc="F F F"', content)


# ---------- roundtrip ----------


class TestRoundtripStress(unittest.TestCase):
    """Write extxyz with stress → read back → virials preserved."""

    def test_roundtrip(self):
        ms1 = dpdata.MultiSystems.from_file("xyz/stress_only.xyz", fmt="extxyz")
        sys1 = list(ms1.systems.values())[0]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
            sys1.to("extxyz", f.name)
            f.flush()
            ms2 = dpdata.MultiSystems.from_file(f.name, fmt="extxyz")
        sys2 = list(ms2.systems.values())[0]
        np.testing.assert_allclose(
            sys1.data["virials"], sys2.data["virials"], atol=1e-6
        )
        np.testing.assert_allclose(
            sys1.data["energies"], sys2.data["energies"], atol=1e-10
        )
        np.testing.assert_allclose(sys1.data["forces"], sys2.data["forces"], atol=1e-10)


class TestFromLabeledSystemDirect(unittest.TestCase):
    """LabeledSystem('file.xyz', fmt='extxyz') should work directly."""

    def test_direct_read(self):
        system = dpdata.LabeledSystem("xyz/stress_only.xyz", fmt="extxyz")
        self.assertEqual(system.get_nframes(), 1)
        np.testing.assert_allclose(system.data["energies"], [-1.5])
        self.assertIn("virials", system.data)


# ============================================================================
# Unit conversion tests (Phase 5 of PR #678 redesign)
# ============================================================================


class TestUnitConvertModule(unittest.TestCase):
    """Direct tests for dpdata.formats.xyz._unit_convert helpers."""

    def test_get_unit_factor_none_returns_1(self):
        from dpdata.formats.xyz._unit_convert import _get_unit_factor

        self.assertEqual(_get_unit_factor(None, "energy"), 1.0)
        self.assertEqual(_get_unit_factor(None, "force"), 1.0)
        self.assertEqual(_get_unit_factor(None, "stress"), 1.0)

    def test_get_unit_factor_energy_ev(self):
        from dpdata.formats.xyz._unit_convert import _get_unit_factor

        self.assertAlmostEqual(_get_unit_factor("eV", "energy"), 1.0)

    def test_get_unit_factor_energy_hartree(self):
        from dpdata.formats.xyz._unit_convert import _get_unit_factor
        from dpdata.unit import EnergyConversion

        expected = EnergyConversion("hartree", "eV").value()
        self.assertAlmostEqual(_get_unit_factor("hartree", "energy"), expected)
        self.assertAlmostEqual(_get_unit_factor("Ha", "energy"), expected)

    def test_get_unit_factor_energy_kcal_mol(self):
        from dpdata.formats.xyz._unit_convert import _get_unit_factor
        from dpdata.unit import EnergyConversion

        expected = EnergyConversion("kcal_mol", "eV").value()
        self.assertAlmostEqual(_get_unit_factor("kcal/mol", "energy"), expected)

    def test_get_unit_factor_force_hartree_bohr(self):
        from dpdata.formats.xyz._unit_convert import _get_unit_factor
        from dpdata.unit import ForceConversion

        expected = ForceConversion("hartree/bohr", "eV/angstrom").value()
        self.assertAlmostEqual(_get_unit_factor("hartree/bohr", "force"), expected)

    def test_get_unit_factor_force_kcal_mol_ang(self):
        from dpdata.formats.xyz._unit_convert import _get_unit_factor
        from dpdata.unit import ForceConversion

        expected = ForceConversion("kcal_mol/angstrom", "eV/angstrom").value()
        self.assertAlmostEqual(_get_unit_factor("kcal/mol/angstrom", "force"), expected)

    def test_get_unit_factor_stress_gpa(self):
        from dpdata.formats.xyz._unit_convert import _get_unit_factor
        from dpdata.unit import PressureConversion

        expected = PressureConversion("GPa", "eV/angstrom^3").value()
        self.assertAlmostEqual(_get_unit_factor("GPa", "stress"), expected)

    def test_get_unit_factor_stress_kbar(self):
        from dpdata.formats.xyz._unit_convert import _get_unit_factor
        from dpdata.unit import PressureConversion

        expected = PressureConversion("kbar", "eV/angstrom^3").value()
        self.assertAlmostEqual(_get_unit_factor("kbar", "stress"), expected)

    def test_unsupported_unit_raises(self):
        from dpdata.formats.xyz._unit_convert import _get_unit_factor

        with self.assertRaises(ValueError):
            _get_unit_factor("nonsense", "energy")
        with self.assertRaises(ValueError):
            _get_unit_factor("bad/unit", "force")
        with self.assertRaises(ValueError):
            _get_unit_factor("megapascal", "stress")

    def test_parse_force_unit(self):
        from dpdata.formats.xyz._unit_convert import _parse_force_unit

        self.assertEqual(
            _parse_force_unit("kcal/mol/angstrom"), ("kcal/mol", "angstrom")
        )
        self.assertEqual(_parse_force_unit("hartree/bohr"), ("hartree", "bohr"))
        self.assertEqual(_parse_force_unit("ev/ang"), ("ev", "ang"))

    def test_parse_force_unit_invalid(self):
        from dpdata.formats.xyz._unit_convert import _parse_force_unit

        with self.assertRaises(ValueError):
            _parse_force_unit("noslash")


class TestEnergyUnitHartree(unittest.TestCase):
    """Test reading extxyz with energy-unit=hartree and force-unit=hartree/bohr."""

    def test_conversion(self):
        from dpdata.unit import EnergyConversion, ForceConversion

        system = dpdata.LabeledSystem("xyz/energy_hartree.xyz", fmt="extxyz")
        e_factor = EnergyConversion("hartree", "eV").value()
        f_factor = ForceConversion("hartree/bohr", "eV/angstrom").value()

        # energy: -0.5 hartree → eV
        np.testing.assert_allclose(
            system.data["energies"], [-0.5 * e_factor], rtol=1e-10
        )
        # forces: first atom [0.01, 0.02, 0.03] hartree/bohr → eV/angstrom
        np.testing.assert_allclose(
            system.data["forces"][0, 0],
            np.array([0.01, 0.02, 0.03]) * f_factor,
            rtol=1e-10,
        )


class TestEnergyUnitKcalMol(unittest.TestCase):
    """Test reading extxyz with energy-unit=kcal/mol and force-unit=kcal/mol/angstrom."""

    def test_conversion(self):
        from dpdata.unit import EnergyConversion, ForceConversion

        system = dpdata.LabeledSystem("xyz/energy_kcal_mol.xyz", fmt="extxyz")
        e_factor = EnergyConversion("kcal_mol", "eV").value()
        f_factor = ForceConversion("kcal_mol/angstrom", "eV/angstrom").value()

        np.testing.assert_allclose(
            system.data["energies"], [-10.0 * e_factor], rtol=1e-10
        )
        np.testing.assert_allclose(
            system.data["forces"][0, 0],
            np.array([1.0, 2.0, 3.0]) * f_factor,
            rtol=1e-10,
        )


class TestStressUnitGPa(unittest.TestCase):
    """Test reading extxyz with stress-unit=GPa."""

    def test_conversion(self):
        from dpdata.unit import PressureConversion

        system = dpdata.LabeledSystem("xyz/stress_gpa.xyz", fmt="extxyz")
        s_factor = PressureConversion("GPa", "eV/angstrom^3").value()
        volume = 5.0**3  # cubic cell 5x5x5

        # stress is identity * 1.0 GPa → virial = -V * stress_internal
        expected_virial_diag = -1.0 * volume * 1.0 * s_factor
        np.testing.assert_allclose(
            np.diag(system.data["virials"][0]),
            [expected_virial_diag] * 3,
            rtol=1e-10,
        )


class TestStressUnitKbar(unittest.TestCase):
    """Test reading extxyz with stress-unit=kbar."""

    def test_conversion(self):
        from dpdata.unit import PressureConversion

        system = dpdata.LabeledSystem("xyz/stress_kbar.xyz", fmt="extxyz")
        s_factor = PressureConversion("kbar", "eV/angstrom^3").value()
        volume = 5.0**3

        # stress is identity * 10.0 kbar → virial = -V * stress_internal
        expected_virial_diag = -1.0 * volume * 10.0 * s_factor
        np.testing.assert_allclose(
            np.diag(system.data["virials"][0]),
            [expected_virial_diag] * 3,
            rtol=1e-10,
        )


class TestNoUnitHeaderDefaults(unittest.TestCase):
    """Without unit headers, values should be unchanged (factor=1.0)."""

    def test_default_no_conversion(self):
        # stress_only.xyz has no unit headers
        system = dpdata.LabeledSystem("xyz/stress_only.xyz", fmt="extxyz")
        # Energy should be exactly as in file
        np.testing.assert_allclose(system.data["energies"], [-1.5])


class TestUnsupportedUnitRaises(unittest.TestCase):
    """Unsupported unit string should raise ValueError, not silently pass."""

    def test_raises_on_bad_energy_unit(self):
        with self.assertRaises(ValueError):
            dpdata.LabeledSystem("xyz/unsupported_unit.xyz", fmt="extxyz")


class TestVirialUnitConversion(unittest.TestCase):
    """Virial has energy dimensions; must be converted when energy-unit is set."""

    def test_virial_hartree(self):
        from dpdata.unit import EnergyConversion

        system = dpdata.LabeledSystem("xyz/virial_hartree.xyz", fmt="extxyz")
        e_factor = EnergyConversion("hartree", "eV").value()

        # virial in file is identity matrix in hartree → should be identity * e_factor in eV
        expected_diag = 1.0 * e_factor
        np.testing.assert_allclose(
            np.diag(system.data["virials"][0]),
            [expected_diag] * 3,
            rtol=1e-10,
        )


class TestAtomicUnitAliases(unittest.TestCase):
    """energy-unit=au and force-unit=a.u. should map to hartree and hartree/bohr."""

    def test_au_alias(self):
        from dpdata.unit import EnergyConversion, ForceConversion

        system = dpdata.LabeledSystem("xyz/energy_au.xyz", fmt="extxyz")
        e_factor = EnergyConversion("hartree", "eV").value()
        f_factor = ForceConversion("hartree/bohr", "eV/angstrom").value()

        np.testing.assert_allclose(
            system.data["energies"], [-0.5 * e_factor], rtol=1e-10
        )
        np.testing.assert_allclose(
            system.data["forces"][0, 0],
            np.array([0.01, 0.02, 0.03]) * f_factor,
            rtol=1e-10,
        )


class TestStressUnitAtomicUnitAlias(unittest.TestCase):
    """stress-unit=au/bohr^3 should map to hartree/bohr^3."""

    def test_stress_au_bohr3(self):
        from dpdata.unit import PressureConversion

        system = dpdata.LabeledSystem("xyz/stress_au_bohr3.xyz", fmt="extxyz")
        s_factor = PressureConversion("hartree/bohr^3", "eV/angstrom^3").value()
        volume = 5.0**3

        # stress is identity * 0.001 au/bohr^3 → virial = -V * stress_internal
        expected_virial_diag = -1.0 * volume * 0.001 * s_factor
        np.testing.assert_allclose(
            np.diag(system.data["virials"][0]),
            [expected_virial_diag] * 3,
            rtol=1e-10,
        )


if __name__ == "__main__":
    unittest.main()
