from dpdata.xyz.quip_gap_xyz import QuipGapxyzSystems
from dpdata.format import Format


@Format.register("quip/gap/xyz")
@Format.register("quip/gap/xyz_file")
class QuipGapXYZFormat(Format):
    def from_labeled_system(self, data, **kwargs):
        return data

    def from_multi_systems(self, file_name, **kwargs):
        # here directory is the file_name
        return QuipGapxyzSystems(file_name)