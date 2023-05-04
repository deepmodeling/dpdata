import dpdata.gromacs.gro
from dpdata.format import Format


@Format.register("gro")
@Format.register("gromacs/gro")
class GromacsGroFormat(Format):
    def from_system(self, file_name, format_atom_name=True, **kwargs):
        """Load gromacs .gro file.

        Parameters
        ----------
        file_name : str
            The input file name
        format_atom_name : bool
            Whether to format the atom name
        **kwargs : dict
            other parameters
        """
        return dpdata.gromacs.gro.file_to_system_data(
            file_name, format_atom_name=format_atom_name, **kwargs
        )

    def to_system(self, data, file_name=None, frame_idx=-1, **kwargs):
        """Dump the system in gromacs .gro format.

        Parameters
        ----------
        data : dict
            System data
        file_name : str or None
            The output file name. If None, return the file content as a string
        frame_idx : int
            The index of the frame to dump
        **kwargs : dict
            other parameters
        """
        assert frame_idx < len(data["coords"])
        if frame_idx == -1:
            strs = []
            for idx in range(data["coords"].shape[0]):
                gro_str = dpdata.gromacs.gro.from_system_data(data, f_idx=idx, **kwargs)
                strs.append(gro_str)
            gro_str = "\n".join(strs)
        else:
            gro_str = dpdata.gromacs.gro.from_system_data(
                data, f_idx=frame_idx, **kwargs
            )

        if file_name is None:
            return gro_str
        else:
            with open(file_name, "w+") as fp:
                fp.write(gro_str)
