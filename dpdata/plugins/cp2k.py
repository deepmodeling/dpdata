import dpdata.cp2k.output
import glob
from dpdata.cp2k.output import Cp2kSystems
from dpdata.format import Format


@Format.register("cp2k/aimd_output")
class CP2KAIMDOutputFormat(Format):
    def from_labeled_system(self, file_name, restart=False, **kwargs):
        xyz_file = sorted(glob.glob("{}/*pos*.xyz".format(file_name)))[0]
        log_file = sorted(glob.glob("{}/*.log".format(file_name)))[0]
        return tuple(Cp2kSystems(log_file, xyz_file, restart))


@Format.register("cp2k/output")
class CP2KOutputFormat(Format):
    def from_labeled_system(self, file_name, restart=False, **kwargs):
        data = {}
        data['atom_names'], \
            data['atom_numbs'], \
            data['atom_types'], \
            data['cells'], \
            data['coords'], \
            data['energies'], \
            data['forces'], \
            tmp_virial \
            = dpdata.cp2k.output.get_frames(file_name)
        if tmp_virial is not None:
            data['virials'] = tmp_virial
        return data
