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


if __name__ == "__main__":
    unittest.main()
