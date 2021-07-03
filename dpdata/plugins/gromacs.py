import dpdata.gromacs.gro
from dpdata.format import Format


@Format.register("gro")
@Format.register("gromacs/gro")
@Format.register_from("from_gromacs_gro")
@Format.register_to("to_gromacs_gro")
class PwmatOutputFormat(Format):
    def from_system(self, file_name, format_atom_name=True, **kwargs):
        """
        Load gromacs .gro file

        Parameters
        ----------
        file_name : str
            The input file name
        """
        return dpdata.gromacs.gro.file_to_system_data(file_name, format_atom_name=format_atom_name)

    def to_system(self, data, file_name=None, frame_idx=-1, **kwargs):
        """
        Dump the system in gromacs .gro format

        Parameters
        ----------
        file_name : str or None
            The output file name. If None, return the file content as a string
        frame_idx : int
            The index of the frame to dump
        """
        assert(frame_idx < len(data['coords']))
        if frame_idx == -1:
            strs = []
            for idx in range(data['coords'].shape[0]):
                gro_str = dpdata.gromacs.gro.from_system_data(data, f_idx=idx)
                strs.append(gro_str)
            gro_str = "\n".join(strs)
        else:
            gro_str = dpdata.gromacs.gro.from_system_data(
                data, f_idx=frame_idx)

        if file_name is None:
            return gro_str
        else:
            with open(file_name, 'w+') as fp:
                fp.write(gro_str)
