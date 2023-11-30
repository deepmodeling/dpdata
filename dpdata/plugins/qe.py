import dpdata.md.pbc
import dpdata.qe.scf
import dpdata.qe.traj
from dpdata.format import Format
import numpy as np


@Format.register("qe/cp/traj")
class QECPTrajFormat(Format):
    @Format.post("rot_lower_triangular")
    def from_system(self, file_name, begin=0, step=1, **kwargs):
        data, _ = dpdata.qe.traj.to_system_data(
            file_name + ".in", file_name, begin=begin, step=step
        )
        data["coords"] = dpdata.md.pbc.apply_pbc(
            data["coords"],
            data["cells"],
        )
        return data

    @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, begin=0, step=1, **kwargs):
        data, cs = dpdata.qe.traj.to_system_data(
            file_name + ".in", file_name, begin=begin, step=step
        )
        data["coords"] = dpdata.md.pbc.apply_pbc(
            data["coords"],
            data["cells"],
        )
        data["energies"], data["forces"], stress, es = dpdata.qe.traj.to_system_label(
            file_name + ".in", file_name, begin=begin, step=step
        )
        if stress is not None:
            # Compute virials of all frames by: virial = stress * volume
            virial = stress * np.linalg.det(data["cells"])[:, np.newaxis, np.newaxis]
            data["virials"] = virial
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
        ) = dpdata.qe.scf.get_frame(file_name)
        if tmp_virial is not None:
            data["virials"] = tmp_virial
        return data
