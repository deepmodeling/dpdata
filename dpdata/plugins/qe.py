from __future__ import annotations

import dpdata.formats.md.pbc
import dpdata.formats.qe.scf
import dpdata.formats.qe.traj
from dpdata.format import Format


@Format.register("qe/cp/traj")
class QECPTrajFormat(Format):
    @Format.post("rot_lower_triangular")
    def from_system(self, file_name, begin=0, step=1, **kwargs):
        data, _ = dpdata.formats.qe.traj.to_system_data(
            file_name + ".in", file_name, begin=begin, step=step
        )
        data["coords"] = dpdata.formats.md.pbc.apply_pbc(
            data["coords"],
            data["cells"],
        )
        return data

    @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, begin=0, step=1, **kwargs):
        data, cs = dpdata.formats.qe.traj.to_system_data(
            file_name + ".in", file_name, begin=begin, step=step
        )
        data["coords"] = dpdata.formats.md.pbc.apply_pbc(
            data["coords"],
            data["cells"],
        )
        data["energies"], data["forces"], es = dpdata.formats.qe.traj.to_system_label(
            file_name + ".in", file_name, begin=begin, step=step
        )
        assert cs == es, "the step key between files are not consistent"
        return data


@Format.register("qe/pw/scf")
class QECPPWSCFFormat(Format):
    @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, **kwargs):
        data = {}
        (
            data["atom_names"],
            data["atom_numbs"],
            data["atom_types"],
            data["cells"],
            data["coords"],
            data["energies"],
            data["forces"],
            tmp_virial,
        ) = dpdata.formats.qe.scf.get_frame(file_name)
        if tmp_virial is not None:
            data["virials"] = tmp_virial
        return data
