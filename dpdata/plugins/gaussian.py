import dpdata.gaussian.log
from dpdata.format import Format


@Format.register("gaussian/log")
class GaussianLogFormat(Format):
    def from_labeled_system(self, file_name, md=False, **kwargs):
        try:
            return dpdata.gaussian.log.to_system_data(file_name, md=md)
        except AssertionError:
            return {
                'energies': [],
                'forces': [],
                'nopbc': True
            }


@Format.register("gaussian/md")
class GaussianMDFormat(Format):
    def from_labeled_system(self, file_name, **kwargs):
        return GaussianLogFormat().from_labeled_system(file_name, md=True)
