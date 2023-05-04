import numpy as np

import dpdata.pwmat.atomconfig
import dpdata.pwmat.movement
from dpdata.format import Format


@Format.register("movement")
@Format.register("mlmd")
@Format.register("pwmat/movement")
@Format.register("pwmat/mlmd")
@Format.register("pwmat/output")
class PwmatOutputFormat(Format):
    @Format.post("rot_lower_triangular")
    def from_labeled_system(
        self, file_name, begin=0, step=1, convergence_check=True, **kwargs
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
        ) = dpdata.pwmat.movement.get_frames(
            file_name, begin=begin, step=step, convergence_check=convergence_check
        )
        if tmp_virial is not None:
            data["virials"] = tmp_virial
        # scale virial to the unit of eV
        if "virials" in data:
            v_pref = 1 * 1e3 / 1.602176621e6
            for ii in range(data["coords"].shape[0]):
                vol = np.linalg.det(np.reshape(data["cells"][ii], [3, 3]))
                data["virials"][ii] *= v_pref * vol
        return data


@Format.register("atom.config")
@Format.register("final.config")
@Format.register("pwmat/atom.config")
@Format.register("pwmat/final.config")
class PwmatAtomconfigFormat(Format):
    @Format.post("rot_lower_triangular")
    def from_system(self, file_name, **kwargs):
        with open(file_name) as fp:
            lines = [line.rstrip("\n") for line in fp]
        return dpdata.pwmat.atomconfig.to_system_data(lines)

    def to_system(self, data, file_name, frame_idx=0, *args, **kwargs):
        """Dump the system in pwmat atom.config format.

        Parameters
        ----------
        data : dict
            The system data
        file_name : str
            The output file name
        frame_idx : int
            The index of the frame to dump
        *args : list
            other parameters
        **kwargs : dict
            other parameters
        """
        assert frame_idx < len(data["coords"])
        w_str = dpdata.pwmat.atomconfig.from_system_data(data, frame_idx)
        with open(file_name, "w") as fp:
            fp.write(w_str)
