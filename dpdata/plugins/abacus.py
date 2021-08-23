import dpdata.abacus.scf
from dpdata.format import Format


@Format.register("abacus/scf")
@Format.register("abacus/pw/scf")
class AbacusSCFFormat(Format):
    @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, **kwargs):
        return dpdata.abacus.scf.get_frame(file_name)
