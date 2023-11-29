import dpdata.md.pbc
import dpdata.openmx.omx
from dpdata.format import Format


@Format.register("openmx")
class QECPTrajFormat(Format):
    @Format.post("rot_lower_triangular")
    def from_system(self, file_name, begin=0, step=1, **kwargs):
        fname=f"{file_name}.dat"
        mdname=f"{file_name}.md"

        data, _ = dpdata.openmx.omx.to_system_data(
            fname, mdname, begin=begin, step=step
        )
        data["coords"] = dpdata.md.pbc.apply_pbc(
            data["coords"],
            data["cells"],
        )
        return data

    @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, begin=0, step=1, **kwargs):
        fname=f"{file_name}.dat"
        mdname=f"{file_name}.md"

        data, cs = dpdata.openmx.omx.to_system_data(
            fname, mdname, begin=begin, step=step
        )
        data["coords"] = dpdata.md.pbc.apply_pbc(
            data["coords"],
            data["cells"],
        )
        data["energies"], data["forces"] = dpdata.openmx.omx.to_system_label(
            fname, mdname, data, begin=begin, step=step
        )
        return data
