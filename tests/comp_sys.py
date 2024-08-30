from __future__ import annotations

import numpy as np


class CompSys:
    def test_len_func(self):
        self.assertEqual(len(self.system_1), len(self.system_2))

    def test_add_func(self):
        self.assertEqual(
            len(self.system_1 + self.system_1), len(self.system_2 + self.system_2)
        )

    def test_atom_numbs(self):
        self.assertEqual(
            self.system_1.data["atom_numbs"], self.system_2.data["atom_numbs"]
        )

    def test_atom_names(self):
        self.assertEqual(
            self.system_1.data["atom_names"], self.system_2.data["atom_names"]
        )

    def test_atom_types(self):
        np.testing.assert_array_equal(
            self.system_1.data["atom_types"], self.system_2.data["atom_types"]
        )

    def test_orig(self):
        for d0 in range(3):
            self.assertEqual(
                self.system_1.data["orig"][d0], self.system_2.data["orig"][d0]
            )

    def test_nframs(self):
        self.assertEqual(self.system_1.get_nframes(), self.system_2.get_nframes())

    def test_cell(self):
        self.assertEqual(self.system_1.get_nframes(), self.system_2.get_nframes())
        if not self.system_1.nopbc and not self.system_2.nopbc:
            np.testing.assert_almost_equal(
                self.system_1.data["cells"],
                self.system_2.data["cells"],
                decimal=self.places,
                err_msg="cell failed",
            )

    def test_coord(self):
        self.assertEqual(self.system_1.get_nframes(), self.system_2.get_nframes())
        # think about direct coord
        if self.system_1.nopbc:
            # nopbc doesn't need to test cells
            return
        tmp_cell = self.system_1.data["cells"]
        tmp_cell = np.reshape(tmp_cell, [-1, 3])
        tmp_cell_norm = np.reshape(np.linalg.norm(tmp_cell, axis=1), [-1, 1, 3])
        np.testing.assert_almost_equal(
            self.system_1.data["coords"] / tmp_cell_norm,
            self.system_2.data["coords"] / tmp_cell_norm,
            decimal=self.places,
            err_msg="coord failed",
        )

    def test_nopbc(self):
        self.assertEqual(self.system_1.nopbc, self.system_2.nopbc)

    def test_data_check(self):
        self.system_1.check_data()
        self.system_2.check_data()


class CompLabeledSys(CompSys):
    def test_energy(self):
        self.assertEqual(self.system_1.get_nframes(), self.system_2.get_nframes())
        np.testing.assert_almost_equal(
            self.system_1.data["energies"],
            self.system_2.data["energies"],
            decimal=self.e_places,
            err_msg="energies failed",
        )

    def test_force(self):
        self.assertEqual(self.system_1.get_nframes(), self.system_2.get_nframes())
        np.testing.assert_almost_equal(
            self.system_1.data["forces"],
            self.system_2.data["forces"],
            decimal=self.f_places,
            err_msg="forces failed",
        )

    def test_virial(self):
        self.assertEqual(self.system_1.get_nframes(), self.system_2.get_nframes())
        # if len(self.system_1['virials']) == 0:
        #     self.assertEqual(len(self.system_1['virials']), 0)
        #     return
        if not self.system_1.has_virial():
            self.assertFalse(self.system_2.has_virial())
            return
        np.testing.assert_almost_equal(
            self.system_1["virials"],
            self.system_2["virials"],
            decimal=self.v_places,
            err_msg="virials failed",
        )


def _make_comp_ms_test_func(comp_sys_test_func):
    """Functional making test function for multi-sys from the that for
    single sys.
    """

    def comp_ms_test_func(iobj):
        iobj.assertEqual(len(iobj.ms_1), len(iobj.ms_2))
        keys = [ii.formula for ii in iobj.ms_1]
        for kk in keys:
            iobj.system_1 = iobj.ms_1[kk]
            iobj.system_2 = iobj.ms_2[kk]
            comp_sys_test_func(iobj)
        del iobj.system_1
        del iobj.system_2

    return comp_ms_test_func


def _make_comp_ms_class(comp_class):
    """Class functinal making multi-sys comparison case from single-sys
    comparison case.
    """

    class CompMS:
        pass

    test_methods = [
        func
        for func in dir(comp_class)
        if callable(getattr(comp_class, func)) and func.startswith("test_")
    ]

    for func in test_methods:
        setattr(CompMS, func, _make_comp_ms_test_func(getattr(comp_class, func)))

    return CompMS


# MultiSystems comparison from single System comparison
CompMultiSys = _make_comp_ms_class(CompSys)

# LabeledMultiSystems comparison from single LabeledSystem comparison
CompLabeledMultiSys = _make_comp_ms_class(CompLabeledSys)


class MultiSystems:
    def test_systems_name(self):
        self.assertEqual(set(self.systems.systems), set(self.system_names))

    def test_systems_size(self):
        for name, size in self.system_sizes.items():
            self.assertEqual(self.systems[name].get_nframes(), size)

    def test_atom_names(self):
        self.assertEqual(self.atom_names, self.systems.atom_names)


class IsPBC:
    def test_is_pbc(self):
        self.assertFalse(self.system_1.nopbc)
        self.assertFalse(self.system_2.nopbc)


class IsNoPBC:
    def test_is_nopbc(self):
        self.assertTrue(self.system_1.nopbc)
        self.assertTrue(self.system_2.nopbc)


class MSAllIsPBC:
    def test_is_pbc(self):
        self.assertTrue(all([not ss.nopbc for ss in self.ms_1]))
        self.assertTrue(all([not ss.nopbc for ss in self.ms_2]))


class MSAllIsNoPBC:
    def test_is_nopbc(self):
        self.assertTrue(all([ss.nopbc for ss in self.ms_1]))
        self.assertTrue(all([ss.nopbc for ss in self.ms_2]))
