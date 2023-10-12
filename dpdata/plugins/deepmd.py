from __future__ import annotations

import os

try:
    import h5py
except ImportError:
    pass
import numpy as np

import dpdata
import dpdata.deepmd.comp
import dpdata.deepmd.hdf5
import dpdata.deepmd.mixed
import dpdata.deepmd.raw
from dpdata.driver import Driver
from dpdata.format import Format


@Format.register("deepmd")
@Format.register("deepmd/raw")
class DeePMDRawFormat(Format):
    def from_system(self, file_name, type_map=None, **kwargs):
        return dpdata.deepmd.raw.to_system_data(
            file_name, type_map=type_map, labels=False
        )

    def to_system(self, data, file_name, **kwargs):
        """Dump the system in deepmd raw format to directory `file_name`."""
        dpdata.deepmd.raw.dump(file_name, data)

    def from_labeled_system(self, file_name, type_map=None, **kwargs):
        return dpdata.deepmd.raw.to_system_data(
            file_name, type_map=type_map, labels=True
        )

    MultiMode = Format.MultiModes.Directory


@Format.register("deepmd/npy")
@Format.register("deepmd/comp")
class DeePMDCompFormat(Format):
    def from_system(self, file_name, type_map=None, **kwargs):
        return dpdata.deepmd.comp.to_system_data(
            file_name, type_map=type_map, labels=False
        )

    def to_system(self, data, file_name, set_size=5000, prec=np.float64, **kwargs):
        """Dump the system in deepmd compressed format (numpy binary) to `folder`.

        The frames are firstly split to sets, then dumped to seperated subfolders named as `folder/set.000`, `folder/set.001`, ....

        Each set contains `set_size` frames.
        The last set may have less frames than `set_size`.

        Parameters
        ----------
        data : dict
            System data
        file_name : str
            The output folder
        set_size : int
            The size of each set.
        prec : {numpy.float32, numpy.float64}
            The floating point precision of the compressed data
        **kwargs : dict
            other parameters
        """
        dpdata.deepmd.comp.dump(file_name, data, set_size=set_size, comp_prec=prec)

    def from_labeled_system(self, file_name, type_map=None, **kwargs):
        return dpdata.deepmd.comp.to_system_data(
            file_name, type_map=type_map, labels=True
        )

    MultiMode = Format.MultiModes.Directory


@Format.register("deepmd/npy/mixed")
class DeePMDMixedFormat(Format):
    """Mixed type numpy format for DeePMD-kit.
    Under this format, systems with the same number of atoms but different formula can be put together
    for a larger system, especially when the frame numbers in systems are sparse.
    This also helps to mixture the type information together for model training with type embedding network.

    Examples
    --------
    Dump a MultiSystems into a mixed type numpy directory:

    >>> import dpdata
    >>> dpdata.MultiSystems(*systems).to_deepmd_npy_mixed("mixed_dir")

    Load a mixed type data into a MultiSystems:

    >>> import dpdata
    >>> dpdata.MultiSystems().load_systems_from_file("mixed_dir", fmt="deepmd/npy/mixed")
    """

    def from_system_mix(self, file_name, type_map=None, **kwargs):
        return dpdata.deepmd.mixed.to_system_data(
            file_name, type_map=type_map, labels=False
        )

    def to_system(
        self, data, file_name, set_size: int = 2000, prec=np.float64, **kwargs
    ):
        """Dump the system in deepmd mixed type format (numpy binary) to `folder`.

        The frames were already split to different systems, so these frames can be dumped to one single subfolders
            named as `folder/set.000`, containing less than `set_size` frames.

        Parameters
        ----------
        data : dict
            System data
        file_name : str
            The output folder
        set_size : int, default=2000
            set size
        prec : {numpy.float32, numpy.float64}
            The floating point precision of the compressed data
        **kwargs : dict
            other parameters
        """
        dpdata.deepmd.mixed.dump(file_name, data, set_size=set_size, comp_prec=prec)

    def from_labeled_system_mix(self, file_name, type_map=None, **kwargs):
        return dpdata.deepmd.mixed.to_system_data(
            file_name, type_map=type_map, labels=True
        )

    def mix_system(self, *system, type_map, **kwargs):
        """Mix the systems into mixed_type ones according to the unified given type_map.

        Parameters
        ----------
        *system : System
            The systems to mix
        type_map : list of str
            Maps atom type to name
        **kwargs : dict
            other parameters

        Returns
        -------
        mixed_systems: dict
            dict of mixed system with key 'atom_numbs'
        """
        return dpdata.deepmd.mixed.mix_system(*system, type_map=type_map, **kwargs)

    def from_multi_systems(self, directory, **kwargs):
        sys_dir = []
        for root, dirs, files in os.walk(directory):
            if (
                "type_map.raw" in files
            ):  # mixed_type format systems must have type_map.raw
                sys_dir.append(root)
        return sys_dir

    MultiMode = Format.MultiModes.Directory


@Format.register("deepmd/hdf5")
class DeePMDHDF5Format(Format):
    """HDF5 format for DeePMD-kit.

    Examples
    --------
    Dump a MultiSystems to a HDF5 file:

    >>> import dpdata
    >>> dpdata.MultiSystems().from_deepmd_npy("data").to_deepmd_hdf5("data.hdf5")
    """

    def _from_system(
        self,
        file_name: str | (h5py.Group | h5py.File),
        type_map: list[str],
        labels: bool,
    ):
        """Convert HDF5 file to System or LabeledSystem data.

        This method is used to switch from labeled or non-labeled options.

        Parameters
        ----------
        file_name : str or h5py.Group or h5py.File
            file name of the HDF5 file or HDF5 object. If it is a string,
            hashtag is used to split path to the HDF5 file and the HDF5 group
        type_map : dict[str]
            type map
        labels : bool
            if Labeled

        Returns
        -------
        dict
            System or LabeledSystem data

        Raises
        ------
        TypeError
            file_name is not str or h5py.Group or h5py.File
        """
        if isinstance(file_name, (h5py.Group, h5py.File)):
            return dpdata.deepmd.hdf5.to_system_data(
                file_name, "", type_map=type_map, labels=labels
            )
        elif isinstance(file_name, str):
            s = file_name.split("#")
            name = s[1] if len(s) > 1 else ""
            with h5py.File(s[0], "r") as f:
                return dpdata.deepmd.hdf5.to_system_data(
                    f, name, type_map=type_map, labels=labels
                )
        else:
            raise TypeError("Unsupported file_name")

    def from_system(
        self,
        file_name: str | (h5py.Group | h5py.File),
        type_map: list[str] | None = None,
        **kwargs,
    ) -> dict:
        """Convert HDF5 file to System data.

        Parameters
        ----------
        file_name : str or h5py.Group or h5py.File
            file name of the HDF5 file or HDF5 object. If it is a string,
            hashtag is used to split path to the HDF5 file and the HDF5 group
        type_map : dict[str]
            type map
        **kwargs : dict
            other parameters

        Returns
        -------
        dict
            System data

        Raises
        ------
        TypeError
            file_name is not str or h5py.Group or h5py.File
        """
        return self._from_system(file_name, type_map=type_map, labels=False)

    def from_labeled_system(
        self,
        file_name: str | (h5py.Group | h5py.File),
        type_map: list[str] | None = None,
        **kwargs,
    ) -> dict:
        """Convert HDF5 file to LabeledSystem data.

        Parameters
        ----------
        file_name : str or h5py.Group or h5py.File
            file name of the HDF5 file or HDF5 object. If it is a string,
            hashtag is used to split path to the HDF5 file and the HDF5 group
        type_map : dict[str]
            type map
        **kwargs : dict
            other parameters

        Returns
        -------
        dict
            LabeledSystem data

        Raises
        ------
        TypeError
            file_name is not str or h5py.Group or h5py.File
        """
        return self._from_system(file_name, type_map=type_map, labels=True)

    def to_system(
        self,
        data: dict,
        file_name: str | (h5py.Group | h5py.File),
        set_size: int = 5000,
        comp_prec: np.dtype = np.float64,
        **kwargs,
    ):
        """Convert System data to HDF5 file.

        Parameters
        ----------
        data : dict
            data dict
        file_name : str or h5py.Group or h5py.File
            file name of the HDF5 file or HDF5 object. If it is a string,
            hashtag is used to split path to the HDF5 file and the HDF5 group
        set_size : int, default=5000
            set size
        comp_prec : np.dtype
            data precision
        **kwargs : dict
            other parameters
        """
        if isinstance(file_name, (h5py.Group, h5py.File)):
            dpdata.deepmd.hdf5.dump(
                file_name, "", data, set_size=set_size, comp_prec=comp_prec
            )
        elif isinstance(file_name, str):
            s = file_name.split("#")
            name = s[1] if len(s) > 1 else ""
            with h5py.File(s[0], "w") as f:
                dpdata.deepmd.hdf5.dump(
                    f, name, data, set_size=set_size, comp_prec=comp_prec
                )
        else:
            raise TypeError("Unsupported file_name")

    def from_multi_systems(self, directory: str, **kwargs) -> h5py.Group:
        """Generate HDF5 groups from a HDF5 file, which will be
        passed to `from_system`.

        Parameters
        ----------
        directory : str
            HDF5 file name
        **kwargs : dict
            other parameters

        Yields
        ------
        h5py.Group
            a HDF5 group in the HDF5 file
        """
        with h5py.File(directory, "r") as f:
            for ff in f.keys():
                yield f[ff]

    def to_multi_systems(
        self, formulas: list[str], directory: str, **kwargs
    ) -> h5py.Group:
        """Generate HDF5 groups, which will be passed to `to_system`.

        Parameters
        ----------
        formulas : list[str]
            formulas of MultiSystems
        directory : str
            HDF5 file name
        **kwargs : dict
            other parameters

        Yields
        ------
        h5py.Group
            a HDF5 group with the name of formula
        """
        with h5py.File(directory, "w") as f:
            for ff in formulas:
                yield f.create_group(ff)


@Driver.register("dp")
@Driver.register("deepmd")
@Driver.register("deepmd-kit")
class DPDriver(Driver):
    """DeePMD-kit driver.

    Parameters
    ----------
    dp : deepmd.DeepPot or str
        The deepmd-kit potential class or the filename of the model.

    Examples
    --------
    >>> DPDriver("frozen_model.pb")
    """

    def __init__(self, dp: str) -> None:
        try:
            # DP 1.x
            import deepmd.DeepPot as DeepPot
        except ModuleNotFoundError:
            # DP 2.x
            from deepmd.infer import DeepPot
        if not isinstance(dp, DeepPot):
            self.dp = DeepPot(dp)
        else:
            self.dp = dp
        self.enable_auto_batch_size = (
            "auto_batch_size" in DeepPot.__init__.__code__.co_varnames
        )

    def label(self, data: dict) -> dict:
        """Label a system data by deepmd-kit. Returns new data with energy, forces, and virials.

        Parameters
        ----------
        data : dict
            data with coordinates and atom types

        Returns
        -------
        dict
            labeled data with energies and forces
        """
        type_map = self.dp.get_type_map()

        ori_sys = dpdata.System.from_dict({"data": data})
        ori_sys_copy = ori_sys.copy()
        ori_sys.sort_atom_names(type_map=type_map)
        atype = ori_sys["atom_types"]
        ori_sys = ori_sys_copy

        if not self.enable_auto_batch_size:
            labeled_sys = dpdata.LabeledSystem()
            for ss in ori_sys:
                coord = ss["coords"].reshape((1, ss.get_natoms() * 3))
                if not ss.nopbc:
                    cell = ss["cells"].reshape((1, 9))
                else:
                    cell = None
                e, f, v = self.dp.eval(coord, cell, atype)
                data = ss.data
                data["energies"] = e.reshape((1,))
                data["forces"] = f.reshape((1, ss.get_natoms(), 3))
                data["virials"] = v.reshape((1, 3, 3))
                this_sys = dpdata.LabeledSystem.from_dict({"data": data})
                labeled_sys.append(this_sys)
            data = labeled_sys.data
        else:
            # since v2.0.2, auto batch size is supported
            coord = ori_sys.data["coords"].reshape(
                (ori_sys.get_nframes(), ori_sys.get_natoms() * 3)
            )
            if not ori_sys.nopbc:
                cell = ori_sys.data["cells"].reshape((ori_sys.get_nframes(), 9))
            else:
                cell = None
            e, f, v = self.dp.eval(coord, cell, atype)
            data = ori_sys.data.copy()
            data["energies"] = e.reshape((ori_sys.get_nframes(),))
            data["forces"] = f.reshape((ori_sys.get_nframes(), ori_sys.get_natoms(), 3))
            data["virials"] = v.reshape((ori_sys.get_nframes(), 3, 3))
        return data
