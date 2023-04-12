import numpy as np

import dpdata.vasp.outcar
import dpdata.vasp.poscar
import dpdata.vasp.xml
from dpdata.format import Format
from dpdata.utils import uniq_atom_names


@Format.register("poscar")
@Format.register("contcar")
@Format.register("vasp/poscar")
@Format.register("vasp/contcar")
class VASPPoscarFormat(Format):
    @Format.post("rot_lower_triangular")
    def from_system(self, file_name, **kwargs):
        with open(file_name) as fp:
            lines = [line.rstrip("\n") for line in fp]
        data = dpdata.vasp.poscar.to_system_data(lines)
        data = uniq_atom_names(data)
        return data

    def to_system(self, data, file_name, frame_idx=0, **kwargs):
        """Dump the system in vasp POSCAR format.

        Parameters
        ----------
        data : dict
            The system data
        file_name : str
            The output file name
        frame_idx : int
            The index of the frame to dump
        **kwargs : dict
            other parameters
        """
        w_str = VASPStringFormat().to_system(data, frame_idx=frame_idx)
        with open(file_name, "w") as fp:
            fp.write(w_str)


@Format.register("vasp/string")
class VASPStringFormat(Format):
    def to_system(self, data, frame_idx=0, **kwargs):
        """Dump the system in vasp POSCAR format string.

        Parameters
        ----------
        data : dict
            The system data
        frame_idx : int
            The index of the frame to dump
        **kwargs : dict
            other parameters
        """
        assert frame_idx < len(data["coords"])
        return dpdata.vasp.poscar.from_system_data(data, frame_idx)


# rotate the system to lammps convention
@Format.register("outcar")
@Format.register("vasp/outcar")
class VASPOutcarFormat(Format):
    @Format.post("rot_lower_triangular")
    def from_labeled_system(
        self, file_name, begin=0, step=1, convergence_check=True, **kwargs
    ):
        data = {}
        ml = kwargs.get("ml", False)
        (
            data["atom_names"],
            data["atom_numbs"],
            data["atom_types"],
            data["cells"],
            data["coords"],
            data["energies"],
            data["forces"],
            tmp_virial,
        ) = dpdata.vasp.outcar.get_frames(
            file_name,
            begin=begin,
            step=step,
            ml=ml,
            convergence_check=convergence_check,
        )
        if tmp_virial is not None:
            data["virials"] = tmp_virial
        # scale virial to the unit of eV
        if "virials" in data:
            v_pref = 1 * 1e3 / 1.602176621e6
            for ii in range(data["cells"].shape[0]):
                vol = np.linalg.det(np.reshape(data["cells"][ii], [3, 3]))
                data["virials"][ii] *= v_pref * vol
        data = uniq_atom_names(data)
        return data


# rotate the system to lammps convention
@Format.register("xml")
@Format.register("vasp/xml")
class VASPXMLFormat(Format):
    @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, begin=0, step=1, **kwargs):
        data = {}
        (
            data["atom_names"],
            data["atom_types"],
            data["cells"],
            data["coords"],
            data["energies"],
            data["forces"],
            data["virials"],
        ) = dpdata.vasp.xml.analyze(
            file_name, type_idx_zero=True, begin=begin, step=step
        )
        data["atom_numbs"] = []
        for ii in range(len(data["atom_names"])):
            data["atom_numbs"].append(sum(data["atom_types"] == ii))
        # the vasp xml assumes the direct coordinates
        # apply the transform to the cartesan coordinates
        for ii in range(data["cells"].shape[0]):
            data["coords"][ii] = np.matmul(data["coords"][ii], data["cells"][ii])
        # scale virial to the unit of eV
        v_pref = 1 * 1e3 / 1.602176621e6
        for ii in range(data["cells"].shape[0]):
            vol = np.linalg.det(np.reshape(data["cells"][ii], [3, 3]))
            data["virials"][ii] *= v_pref * vol
        data = uniq_atom_names(data)
        return data
