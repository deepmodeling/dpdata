import numpy as np

from dpdata.format import Format
from dpdata.xyz.xyz import coord_to_xyz


@Format.register("3dmol")
class AmberMDFormat(Format):
    """See 3D structure of a frame in jupyter.
    
    py3Dmol should be installed in advance.
    """
    def to_system(self, data: dict, f_idx: int = 0, **kwargs):
        import py3Dmol
        types = np.array(data['atom_names'])[data['atom_types']]
        xyz = coord_to_xyz(data['coords'][f_idx], types)
        viewer = py3Dmol.view(width=300, height=300)
        viewer.addModel(xyz, 'xyz')
        viewer.zoomTo()
        return viewer
