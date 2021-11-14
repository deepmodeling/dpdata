import dpdata.deepmd.raw
import dpdata.deepmd.comp
import dpdata.deepmd.hdf5
import numpy as np
import h5py
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

    def from_labeled_system(self, file_name, type_map=None, **kwargs):
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

    def from_labeled_system(self, file_name, type_map=None, **kwargs):
        return dpdata.deepmd.comp.to_system_data(file_name, type_map=type_map, labels=True)

    MultiMode = Format.MultiModes.Directory

@Format.register("deepmd/hdf5")
class DeePMDCompFormat(Format):
    """HDF5 format for DeePMD-kit.
    
    Examples
    --------
    Dump a MultiSystems to a HDF5 file:
    >>> import dpdata
    >>> dpdata.MultiSystems().from_deepmd_npy("data").to_deepmd_hdf5("data.hdf5")
    """
    def from_system(self, file_name, type_map=None, **kwargs):
        s = file_name.split("#")
        name = s[1] if len(s) > 1 else ""
        with h5py.File(s[0], 'r') as f:
            return dpdata.deepmd.hdf5.to_system_data(f, name, type_map=type_map, labels=False)

    def from_labeled_system(self, file_name, type_map=None, **kwargs):
        s = file_name.split("#")
        name = s[1] if len(s) > 1 else ""
        with h5py.File(s[0], 'r') as f:
            return dpdata.deepmd.hdf5.to_system_data(f, name, type_map=type_map, labels=True)
    
    def to_system(self,
                  data : dict,
                  file_name : str,
                  set_size : int = 5000, 
                  comp_prec : np.dtype = np.float32,
                  **kwargs):
        s = file_name.split("#")
        name = s[1] if len(s) > 1 else ""
        mode = 'a' if name else 'w'
        with h5py.File(s[0], mode) as f:
            dpdata.deepmd.hdf5.dump(f, name, data, set_size = set_size, comp_prec = comp_prec)

    def from_multi_systems(self,
                  directory,
                  **kwargs):
        with h5py.File(directory, 'r') as f:
            return ["%s#%s" % (directory, ff) for ff in f.keys()]

    def to_multi_systems(self,
                  formulas,
                  directory,
                  **kwargs):
        return ["%s#%s" % (directory, ff) for ff in formulas]
