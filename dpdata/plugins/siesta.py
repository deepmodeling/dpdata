import dpdata.siesta.aiMD_output
import dpdata.siesta.output
from dpdata.format import Format


@Format.register("siesta/output")
class SiestaOutputFormat(Format):
    def from_system(self, file_name, **kwargs):
        data = {}
        (
            data["atom_names"],
            data["atom_numbs"],
            data["atom_types"],
            data["cells"],
            data["coords"],
            _e,
            _f,
            _v,
        ) = dpdata.siesta.output.obtain_frame(file_name)
        return data

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
            data["virials"],
        ) = dpdata.siesta.output.obtain_frame(file_name)
        return data


@Format.register("siesta/aimd_output")
@Format.register_from("from_siesta_aiMD_output")
class SiestaAIMDOutputFormat(Format):
    def from_system(self, file_name, **kwargs):
        data = {}
        (
            data["atom_names"],
            data["atom_numbs"],
            data["atom_types"],
            data["cells"],
            data["coords"],
            _e,
            _f,
            _v,
        ) = dpdata.siesta.aiMD_output.get_aiMD_frame(file_name)
        return data

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
            data["virials"],
        ) = dpdata.siesta.aiMD_output.get_aiMD_frame(file_name)
        return data
