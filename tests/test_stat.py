from context import dpdata

import unittest


class TestStat(unittest.TestCase):
    def test_errors(self):
        system1 = dpdata.LabeledSystem("gaussian/methane.gaussianlog", fmt="gaussian/log")
        system2 = dpdata.LabeledSystem("amber/sqm_opt.out", fmt="sqm/out")

        e = dpdata.stat.Errors(system1, system2)
        self.assertAlmostEqual(e.e_mae, 1014.7946598792427, 6)
        self.assertAlmostEqual(e.e_rmse, 1014.7946598792427, 6)
        self.assertAlmostEqual(e.f_mae, 0.004113640526088011, 6)
        self.assertAlmostEqual(e.f_rmse, 0.005714011247538185, 6)

    def test_multi_errors(self):
        system1 = dpdata.MultiSystems(dpdata.LabeledSystem("gaussian/methane.gaussianlog", fmt="gaussian/log"))
        system2 = dpdata.MultiSystems(dpdata.LabeledSystem("amber/sqm_opt.out", fmt="sqm/out"))

        e = dpdata.stat.MultiErrors(system1, system2)
        self.assertAlmostEqual(e.e_mae, 1014.7946598792427, 6)
        self.assertAlmostEqual(e.e_rmse, 1014.7946598792427, 6)
        self.assertAlmostEqual(e.f_mae, 0.004113640526088011, 6)
        self.assertAlmostEqual(e.f_rmse, 0.005714011247538185, 6)
