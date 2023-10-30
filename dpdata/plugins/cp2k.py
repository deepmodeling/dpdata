import glob

import dpdata.cp2k.output
from dpdata.cp2k.output import Cp2kSystems
from dpdata.format import Format

string_warning = """
Hi, you got an error from dpdata,
please check if your cp2k files include full information,
otherwise its version is not supported by dpdata.
Try use dpdata plugin from cp2kdata package,
for details, please refer to
https://robinzyb.github.io/cp2kdata/
"""


@Format.register("cp2k/aimd_output")
class CP2KAIMDOutputFormat(Format):
    def from_labeled_system(self, file_name, restart=False, **kwargs):
        xyz_file = sorted(glob.glob(f"{file_name}/*pos*.xyz"))[0]
        log_file = sorted(glob.glob(f"{file_name}/*.log"))[0]
        try:
            return tuple(Cp2kSystems(log_file, xyz_file, restart))
        except (StopIteration, RuntimeError) as e:
            # StopIteration is raised when pattern match is failed
            raise PendingDeprecationWarning(string_warning) from e


@Format.register("cp2k/output")
class CP2KOutputFormat(Format):
    def from_labeled_system(self, file_name, restart=False, **kwargs):
        try:
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
            ) = dpdata.cp2k.output.get_frames(file_name)
            if tmp_virial is not None:
                data["virials"] = tmp_virial
            return data
        # TODO: in the future, we should add exact error type here
        # TODO: when pattern match is failed
        # TODO: For now just use RuntimeError as a placeholder.
        except RuntimeError as e:
            raise PendingDeprecationWarning(string_warning) from e
