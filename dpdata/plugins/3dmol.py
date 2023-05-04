from typing import Tuple

import numpy as np

from dpdata.format import Format
from dpdata.xyz.xyz import coord_to_xyz


@Format.register("3dmol")
class Py3DMolFormat(Format):
    """3DMol format.

    To use this format,  py3Dmol should be installed in advance.
    """

    def to_system(
        self,
        data: dict,
        f_idx: int = 0,
        size: Tuple[int] = (300, 300),
        style: dict = {"stick": {}, "sphere": {"radius": 0.4}},
        **kwargs,
    ):
        """Show 3D structure of a frame in jupyter.

        Parameters
        ----------
        data : dict
            system data
        f_idx : int
            frame index to show
        size : tuple[int]
            (width, height) of the widget
        style : dict
            style of 3DMol. Read 3DMol documentation for details.
        **kwargs : dict
            other parameters

        Examples
        --------
        >>> system.to_3dmol()
        """
        import py3Dmol

        types = np.array(data["atom_names"])[data["atom_types"]]
        xyz = coord_to_xyz(data["coords"][f_idx], types)
        viewer = py3Dmol.view(width=size[0], height=size[1])
        viewer.addModel(xyz, "xyz")
        viewer.setStyle(style.copy())
        viewer.zoomTo()
        return viewer
