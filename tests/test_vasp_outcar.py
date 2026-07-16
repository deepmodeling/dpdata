from __future__ import annotations

import io
import unittest
import warnings

import numpy as np
from comp_sys import CompLabeledSys, IsPBC
from context import dpdata

from dpdata.formats.vasp.outcar import _get_frames_lower
from dpdata.utils import uniq_atom_names


class TestVaspOUTCAR(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem()
        self.system_1.from_vasp_xml("poscars/vasprun.h2o.md.xml")
        self.system_2 = dpdata.LabeledSystem()
        self.system_2.from_vasp_outcar("poscars/OUTCAR.h2o.md")
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 4


class TestVaspOUTCARIncompleteForceTable(unittest.TestCase):
    """Regression tests for force tables cut short before all atom rows."""

    @staticmethod
    def _block(force_rows, *, include_header=False, energy=-1.0):
        lines = []
        if include_header:
            lines.extend(
                [
                    " TITEL  = PAW_PBE H 15Jun2001",
                    " NELM = 60; maximum number of electronic SC steps",
                    " ions per type = 2",
                ]
            )
        lines.extend(
            [
                " VOLUME and BASIS-vectors are now :",
                " filler",
                " filler",
                " filler",
                " filler",
                " 1.0 0.0 0.0",
                " 0.0 1.0 0.0",
                " 0.0 0.0 1.0",
                " POSITION                                       TOTAL-FORCE (eV/Angst)",
                " -----------------------------------------------------------------------------------",
                *force_rows,
                f" free  energy   TOTEN  =       {energy:.6f} eV",
            ]
        )
        return "\n".join(lines) + "\n"

    def test_separator_before_all_atoms_skips_only_incomplete_frame(self):
        # The first frame is valid, but the second table reaches VASP's
        # separator after only one of the two expected atoms.  Older tests
        # exclusively used complete tables, so the unconditional float
        # conversion of the separator was never exercised.
        valid_rows = ["0 0 0 1 2 3", "1 1 1 4 5 6"]
        incomplete_rows = ["2 2 2 7 8 9"]
        contents = self._block(valid_rows, include_header=True) + self._block(
            incomplete_rows, energy=-2.0
        )

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            frames = _get_frames_lower(io.StringIO(contents), "OUTCAR")

        self.assertEqual(frames[3].shape, (1, 3, 3))
        self.assertEqual(frames[4].shape, (1, 2, 3))
        self.assertEqual(frames[6].shape, (1, 2, 3))
        self.assertTrue(
            any("incomplete labels in frame 2" in str(item.message) for item in caught)
        )

    def test_truncated_table_does_not_index_past_block(self):
        # A file truncated immediately after its first atom previously raised
        # IndexError while the parser blindly indexed all ``ntot`` rows.
        valid_rows = ["0 0 0 1 2 3", "1 1 1 4 5 6"]
        truncated_block = "\n".join(
            [
                " POSITION                                       TOTAL-FORCE (eV/Angst)",
                " -----------------------------------------------------------------------------------",
                "0 0 0 1 2 3",
            ]
        )
        contents = self._block(valid_rows, include_header=True) + truncated_block

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            frames = _get_frames_lower(io.StringIO(contents), "OUTCAR")

        self.assertEqual(frames[4].shape, (1, 2, 3))
        self.assertTrue(
            any("expected 2 atom rows, found 1" in str(item.message) for item in caught)
        )


class TestVaspOUTCARTypeMap(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        sys0 = dpdata.LabeledSystem("poscars/OUTCAR.ch4.unconverged", fmt="vasp/outcar")
        sys0.data["atom_names"] = ["A", "C", "B", "H", "D"]
        sys0.data["atom_numbs"] = [0, 1, 0, 4, 0]
        sys0.data["atom_types"] = np.array([3, 3, 3, 3, 1], dtype=int)
        sys1 = dpdata.LabeledSystem(
            "poscars/OUTCAR.ch4.unconverged",
            fmt="vasp/outcar",
            type_map=["A", "C", "B", "H", "D"],
        )
        self.system_1 = sys0
        self.system_2 = sys1
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


class TestVaspOUTCARSkip(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        begin = 1
        step = 3
        end = 10
        self.system_1 = dpdata.LabeledSystem(
            "poscars/OUTCAR.h2o.md.10", fmt="vasp/outcar", begin=begin, step=step
        )
        self.system_2 = dpdata.LabeledSystem(
            "poscars/OUTCAR.h2o.md.10", fmt="vasp/outcar"
        ).sub_system(np.arange(begin, end, step))
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


class TestVaspOUTCARVdw(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem("poscars/OUTCAR.Ge.vdw", fmt="vasp/outcar")
        self.system_2 = dpdata.LabeledSystem()
        self.system_2.from_vasp_xml("poscars/vasprun.Ge.vdw.xml")
        self.places = 5
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6


class TestDuplicatedAtomNames(unittest.TestCase):
    def test(self):
        system = dpdata.LabeledSystem("poscars/6362_OUTCAR", fmt="vasp/outcar")
        expected_types = [0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1]
        self.assertEqual(list(system["atom_types"]), expected_types)
        self.assertEqual(system["atom_names"], ["B", "O"])
        self.assertEqual(system["atom_numbs"], [8, 6])

    def test_type_map(self):
        system = dpdata.LabeledSystem(
            "poscars/6362_OUTCAR", fmt="vasp/outcar", type_map=["O", "B"]
        )
        expected_types = [1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0]
        self.assertEqual(list(system["atom_types"]), expected_types)
        self.assertEqual(system["atom_names"], ["O", "B"])
        self.assertEqual(system["atom_numbs"], [6, 8])


class TestUniqAtomNames(unittest.TestCase):
    def test(self):
        data = {}
        data["atom_names"] = ["O", "H", "O", "H"]
        data["atom_types"] = np.array([0, 1, 2, 3, 3, 2, 1], dtype=int)

        data = uniq_atom_names(data)
        self.assertEqual(list(data["atom_types"]), [0, 1, 0, 1, 1, 0, 1])
        self.assertEqual(list(data["atom_names"]), ["O", "H"])
        self.assertEqual(list(data["atom_numbs"]), [3, 4])


class TestVaspOUTCARML(unittest.TestCase):
    def test(self):
        system1 = dpdata.LabeledSystem(
            "poscars/OUTCAR.ch4.ml", fmt="vasp/outcar", ml=True
        )
        system2 = dpdata.LabeledSystem(
            "poscars/OUTCAR.ch4.ml", fmt="vasp/outcar", ml=False
        )
        expected_types = [0, 0, 0, 0, 1]
        self.assertEqual(list(system1["atom_types"]), expected_types)
        self.assertEqual(system1["atom_names"], ["H", "C"])
        self.assertEqual(len(system1["energies"]), 10)
        self.assertEqual(list(system2["atom_types"]), expected_types)
        self.assertEqual(system2["atom_names"], ["H", "C"])
        self.assertEqual(len(system2["energies"]), 4)


class TestVaspOUTCARNWRITE0(unittest.TestCase):
    def test(self):
        # only the first and last frames that have forces are read
        ss = dpdata.LabeledSystem("poscars/Ti-aimd-nwrite0/OUTCAR")
        self.assertEqual(ss.get_nframes(), 2)


class TestVaspAtomNamesV6(unittest.TestCase):
    def test(self):
        # in vasp v6, the key TITEL is removed. check if the atom names
        # are correctly parsed.
        ss = dpdata.LabeledSystem("poscars/Ti-O-Ti-v6/OUTCAR")
        self.assertEqual(ss.get_atom_names(), ["Ti", "O"])
        self.assertEqual(ss.get_atom_numbs(), [6, 2])
        np.testing.assert_equal(ss.get_atom_types(), [0, 0, 0, 1, 1, 0, 0, 0])


class TestVaspOUTCARLongIonTypes(unittest.TestCase):
    def test(self):
        # vasp<=6.3 only print ions per type for the first 10 types of atoms
        # raise exception when the bug is triggered.
        with self.assertRaises(RuntimeError) as c:
            ss = dpdata.LabeledSystem("poscars/outcar.longit/OUTCAR")
        self.assertTrue(
            "The number of the atom numbers per each type" in str(c.exception)
        )
        self.assertTrue("does not match that of the atom types" in str(c.exception))


if __name__ == "__main__":
    unittest.main()
