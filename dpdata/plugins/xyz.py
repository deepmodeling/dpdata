import numpy as np

from dpdata.xyz.quip_gap_xyz import QuipGapxyzSystems
from dpdata.xyz.xyz import coord_to_xyz
from dpdata.format import Format

@Format.register("xyz")
class XYZFormat(Format):
    """XYZ foramt.

    Examples
    --------
    >>> s.to("xyz", "a.xyz")
    """
    def to_system(self, data, file_name, **kwargs):
        buff = []
        types = np.array(data['atom_names'])[data['atom_types']]
        for cc in data['coords']:
            buff.append(coord_to_xyz(cc, types))
        with open(file_name, 'w') as fp:
            fp.write("\n".join(buff))


@Format.register("quip/gap/xyz")
@Format.register("quip/gap/xyz_file")
class QuipGapXYZFormat(Format):
    def from_labeled_system(self, data, **kwargs):
        return data

    def from_multi_systems(self, file_name, **kwargs):
        # here directory is the file_name
        return QuipGapxyzSystems(file_name)