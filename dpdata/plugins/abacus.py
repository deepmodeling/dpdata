import dpdata.abacus.md
import dpdata.abacus.relax
import dpdata.abacus.scf
from dpdata.format import Format


@Format.register("abacus/stru")
@Format.register("stru")
class AbacusSTRUFormat(Format):
    def from_system(self, file_name, **kwargs):
        return dpdata.abacus.scf.get_frame_from_stru(file_name)

    def to_system(self, data, file_name, frame_idx=0, **kwargs):
        """Dump the system into ABACUS STRU format file.

        Parameters
        ----------
        data : dict
            System data
        file_name : str
            The output file name
        frame_idx : int
            The index of the frame to dump
        pp_file : list of string, optional
            List of pseudo potential files
        numerical_orbital : list of string, optional
            List of orbital files
        mass : list of float, optional
            List of atomic masses
        numerical_descriptor : str, optional
            numerical descriptor file
        **kwargs : dict
            other parameters
        """
        pp_file = kwargs.get("pp_file")
        numerical_orbital = kwargs.get("numerical_orbital")
        mass = kwargs.get("mass")
        numerical_descriptor = kwargs.get("numerical_descriptor")
        stru_string = dpdata.abacus.scf.make_unlabeled_stru(
            data=data,
            frame_idx=frame_idx,
            pp_file=pp_file,
            numerical_orbital=numerical_orbital,
            numerical_descriptor=numerical_descriptor,
            mass=mass,
        )
        with open(file_name, "w") as fp:
            fp.write(stru_string)


@Format.register("abacus/scf")
@Format.register("abacus/pw/scf")
@Format.register("abacus/lcao/scf")
class AbacusSCFFormat(Format):
    # @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, **kwargs):
        return dpdata.abacus.scf.get_frame(file_name)


@Format.register("abacus/md")
@Format.register("abacus/pw/md")
@Format.register("abacus/lcao/md")
class AbacusMDFormat(Format):
    # @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, **kwargs):
        return dpdata.abacus.md.get_frame(file_name)


@Format.register("abacus/relax")
@Format.register("abacus/pw/relax")
@Format.register("abacus/lcao/relax")
class AbacusRelaxFormat(Format):
    # @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, **kwargs):
        return dpdata.abacus.relax.get_frame(file_name)
