import dpdata.deepmd.raw
import dpdata.deepmd.comp
import numpy as np
from dpdata.format import Format


@Format.register("deepmd")
@Format.register("deepmd/raw")
class DeePMDRawFormat(Format):
    def from_system(self, file_name, type_map=None, **kwargs):
        return dpdata.deepmd.raw.to_system_data(file_name, type_map=type_map, labels=False)

    def to_system(self, data, file_name, **kwargs):
        """Dump the system in deepmd raw format to directory `file_name`
        """
        dpdata.deepmd.raw.dump(file_name, data)

    def from_labeled_system(self, file_name, type_map, **kwargs):
        return dpdata.deepmd.raw.to_system_data(file_name, type_map=type_map, labels=True)

    MultiMode = Format.MultiModes.Directory


@Format.register("deepmd/npy")
@Format.register("deepmd/comp")
class DeePMDCompFormat(Format):
    def from_system(self, file_name, type_map=None, **kwargs):
        return dpdata.deepmd.comp.to_system_data(file_name, type_map=type_map, labels=False)

    def to_system(self, data, file_name, set_size=5000, prec=np.float32, **kwargs):
        """
        Dump the system in deepmd compressed format (numpy binary) to `folder`.

        The frames are firstly split to sets, then dumped to seperated subfolders named as `folder/set.000`, `folder/set.001`, ....

        Each set contains `set_size` frames.
        The last set may have less frames than `set_size`.

        Parameters
        ----------
        data: dict
            System data
        file_name : str
            The output folder
        set_size : int
            The size of each set.
        prec: {numpy.float32, numpy.float64}
            The floating point precision of the compressed data
        """
        dpdata.deepmd.comp.dump(
            file_name, data, set_size=set_size, comp_prec=prec)

    def from_labeled_system(self, file_name, type_map, **kwargs):
        return dpdata.deepmd.comp.to_system_data(file_name, type_map=type_map, labels=True)

    MultiMode = Format.MultiModes.Directory
