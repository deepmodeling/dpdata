import dpdata.fhi_aims.output
from dpdata.format import Format


@Format.register("fhi_aims/md")
@Format.register("fhi_aims/output")
class FhiMDFormat(Format):
    def from_labeled_system(
        self, file_name, md=True, begin=0, step=1, convergence_check=True, **kwargs
    ):
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
        ) = dpdata.fhi_aims.output.get_frames(
            file_name,
            md=md,
            begin=begin,
            step=step,
            convergence_check=convergence_check,
        )
        if tmp_virial is not None:
            data["virials"] = tmp_virial
        return data


@Format.register("fhi_aims/scf")
class FhiSCFFormat(Format):
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
        ) = dpdata.fhi_aims.output.get_frames(file_name, md=False, begin=0, step=1)
        if tmp_virial is not None:
            data["virials"] = tmp_virial
        return data
