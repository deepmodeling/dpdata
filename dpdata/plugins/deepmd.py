from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np

import dpdata
import dpdata.formats.deepmd.comp
import dpdata.formats.deepmd.hdf5
import dpdata.formats.deepmd.mixed
import dpdata.formats.deepmd.raw
from dpdata.data_type import Axis, DataType
from dpdata.driver import Driver
from dpdata.format import Format

if TYPE_CHECKING:
    import h5py


def register_spin():
    dt = DataType(
        "spins",
        np.ndarray,
        (Axis.NFRAMES, Axis.NATOMS, 3),
        required=False,
        deepmd_name="spin",
    )
    dpdata.System.register_data_type(dt)
    dpdata.LabeledSystem.register_data_type(dt)

    dt = DataType(
        "force_mags",
        np.ndarray,
        (Axis.NFRAMES, Axis.NATOMS, 3),
        required=False,
        deepmd_name="force_mag",
    )
    dpdata.System.register_data_type(dt)
    dpdata.LabeledSystem.register_data_type(dt)


@Format.register("deepmd")
@Format.register("deepmd/raw")
class DeePMDRawFormat(Format):
    def from_system(self, file_name, type_map=None, **kwargs):
        register_spin()
        return dpdata.formats.deepmd.raw.to_system_data(
            file_name, type_map=type_map, labels=False
        )

    def to_system(self, data, file_name, **kwargs):
        """Dump the system in deepmd raw format to directory `file_name`."""
        dpdata.formats.deepmd.raw.dump(file_name, data)

    def from_labeled_system(self, file_name, type_map=None, **kwargs):
        register_spin()
        return dpdata.formats.deepmd.raw.to_system_data(
            file_name, type_map=type_map, labels=True
        )

    MultiMode = Format.MultiModes.Directory


@Format.register("deepmd/npy")
@Format.register("deepmd/comp")
class DeePMDCompFormat(Format):
    def from_system(self, file_name, type_map=None, **kwargs):
        register_spin()
        return dpdata.formats.deepmd.comp.to_system_data(
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
        dpdata.formats.deepmd.comp.dump(
            file_name, data, set_size=set_size, comp_prec=prec
        )

    def from_labeled_system(self, file_name, type_map=None, **kwargs):
        register_spin()
        return dpdata.formats.deepmd.comp.to_system_data(
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

    Dump with ``atom_numb_pad`` to reduce the number of subdirectories.
    Systems are padded with virtual atoms (type -1) so that atom counts are
    rounded up to the nearest multiple of the given number:

    >>> dpdata.MultiSystems(*systems).to_deepmd_npy_mixed("mixed_dir", atom_numb_pad=8)

    Load a mixed type data into a MultiSystems:

    >>> import dpdata
    >>> dpdata.MultiSystems().load_systems_from_file("mixed_dir", fmt="deepmd/npy/mixed")
    """

    def from_system_mix(self, file_name, type_map=None, **kwargs):
        return dpdata.formats.deepmd.mixed.to_system_data(
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
        dpdata.formats.deepmd.mixed.dump(
            file_name, data, set_size=set_size, comp_prec=prec
        )

    def from_labeled_system_mix(self, file_name, type_map=None, **kwargs):
        return dpdata.formats.deepmd.mixed.to_system_data(
            file_name, type_map=type_map, labels=True
        )

    def mix_system(self, *system, type_map, atom_numb_pad=None, **kwargs):
        """Mix the systems into mixed_type ones according to the unified given type_map.

        Parameters
        ----------
        *system : System
            The systems to mix
        type_map : list of str
            Maps atom type to name
        atom_numb_pad : int, optional
            If provided, pad atom counts to the next multiple of this number
            using virtual atoms (type -1 in real_atom_types). This reduces the
            number of subdirectories when systems have many different atom counts.
            For example, ``atom_numb_pad=8`` groups systems into multiples of 8:
            a 5-atom system is padded to 8, a 9-atom system is padded to 16, etc.
            Virtual atoms are transparently removed when loading the data back.
        **kwargs : dict
            other parameters

        Returns
        -------
        mixed_systems: dict
            dict of mixed system with key 'atom_numbs'

        Examples
        --------
        Dump with padding so that atom counts are rounded up to multiples of 8:

        >>> import dpdata
        >>> dpdata.MultiSystems(*systems).to_deepmd_npy_mixed("mixed_dir", atom_numb_pad=8)
        """
        return dpdata.formats.deepmd.mixed.mix_system(
            *system, type_map=type_map, atom_numb_pad=atom_numb_pad, **kwargs
        )

    def from_multi_systems(self, directory, **kwargs):
        register_spin()
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
        import h5py

        register_spin()

        if isinstance(file_name, (h5py.Group, h5py.File)):
            return dpdata.formats.deepmd.hdf5.to_system_data(
                file_name, "", type_map=type_map, labels=labels
            )
        elif isinstance(file_name, str):
            s = file_name.split("#")
            name = s[1] if len(s) > 1 else ""
            with h5py.File(s[0], "r") as f:
                return dpdata.formats.deepmd.hdf5.to_system_data(
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
        import h5py

        if isinstance(file_name, (h5py.Group, h5py.File)):
            dpdata.formats.deepmd.hdf5.dump(
                file_name, "", data, set_size=set_size, comp_prec=comp_prec
            )
        elif isinstance(file_name, str):
            s = file_name.split("#")
            name = s[1] if len(s) > 1 else ""
            with h5py.File(s[0], "w") as f:
                dpdata.formats.deepmd.hdf5.dump(
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
        import h5py

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
        import h5py

        with h5py.File(directory, "w") as f:
            for ff in formulas:
                yield f.create_group(ff)


@Format.register("deepmd/hdf5/mixed")
class DeePMDHDF5MixedFormat(DeePMDMixedFormat):
    """Mixed type HDF5 format for DeePMD-kit.

    Mixed type data stores frames with the same atom count in one dataset even
    when their formulas differ. The placeholder ``type.raw`` contains only the
    mixed token type, while ``set.*/real_atom_types.npy`` stores the real atom
    type layout for each frame. Loading reconstructs regular Systems by
    splitting frames with different ``real_atom_types`` rows.

    The HDF5 layout mirrors ``deepmd/npy/mixed`` inside HDF5 groups. For
    :class:`dpdata.MultiSystems`, each top-level mixed group is keyed by the
    number of atoms after optional padding, such as ``"4"`` or ``"8"``. A
    string path may include ``"#group/path"`` to read or write mixed data under
    a nested HDF5 group.

    Examples
    --------
    Dump a :class:`dpdata.MultiSystems` object to a mixed HDF5 file:

    >>> systems.to_deepmd_hdf5_mixed("mixed.hdf5")

    Dump with atom-count padding:

    >>> systems.to_deepmd_hdf5_mixed("mixed.hdf5", atom_numb_pad=8)

    Load a mixed HDF5 file into :class:`dpdata.MultiSystems`:

    >>> dpdata.MultiSystems().from_deepmd_hdf5_mixed("mixed.hdf5")
    """

    @staticmethod
    def _load_hdf5_mixed_data(group, type_map=None, labels=True):
        """Load one mixed HDF5 group as a backend data dict.

        Parameters
        ----------
        group : h5py.Group or h5py.File
            HDF5 object containing one mixed DeePMD system group. The group must
            contain ``type.raw``, ``type_map.raw`` and ``set.*`` children.
        type_map : list[str], optional
            Type map used by the generic HDF5 loader.
        labels : bool, default=True
            Whether labeled data such as energies and forces should be loaded.

        Returns
        -------
        dict
            Mixed-type data dict consumed by
            :func:`dpdata.formats.deepmd.mixed.to_system_data`.
        """
        return dpdata.formats.deepmd.hdf5.to_system_data(
            group, "", type_map=type_map, labels=labels
        )

    @staticmethod
    def _dump_hdf5_mixed_data(group, data, set_size, comp_prec, remove_sets=True):
        """Dump one mixed data dict to an HDF5 group.

        Parameters
        ----------
        group : h5py.Group or h5py.File
            Destination HDF5 object.
        data : dict
            Mixed-type data dict prepared by
            :func:`dpdata.formats.deepmd.mixed.dump`.
        set_size : int
            Maximum number of frames per ``set.*`` group.
        comp_prec : numpy.dtype
            Floating point precision for dumped frame data.
        remove_sets : bool, default=True
            Accepted for backend compatibility. HDF5 groups are recreated by the
            caller, so this argument is not used.
        """
        dpdata.formats.deepmd.hdf5.dump(
            group, "", data, set_size=set_size, comp_prec=comp_prec
        )

    @staticmethod
    def _iter_mixed_groups(group):
        """Yield mixed DeePMD HDF5 groups under ``group``.

        A group is considered a mixed system group when it contains
        ``type.raw``, ``type_map.raw`` and at least one ``set.*`` group with a
        ``real_atom_types.npy`` dataset. If the current group is not a system
        group, nested HDF5 groups are searched recursively. This supports files
        written either as a single mixed system at the file root or as
        MultiSystems groups such as ``/4`` and ``/8``.

        Parameters
        ----------
        group : h5py.Group or h5py.File
            HDF5 group or file to scan.

        Yields
        ------
        h5py.Group or h5py.File
            Mixed system groups to pass to ``from_system_mix``.
        """
        import h5py

        set_groups = [
            item
            for key, item in group.items()
            if key.startswith("set.") and isinstance(item, h5py.Group)
        ]
        is_mixed_group = (
            "type.raw" in group
            and "type_map.raw" in group
            and any("real_atom_types.npy" in set_group for set_group in set_groups)
        )
        if is_mixed_group:
            yield group
            return
        for item in group.values():
            if isinstance(item, h5py.Group):
                yield from DeePMDHDF5MixedFormat._iter_mixed_groups(item)

    @staticmethod
    def _get_group(file, name):
        """Return ``file`` or a named child group.

        Parameters
        ----------
        file : h5py.File or h5py.Group
            Root HDF5 object.
        name : str
            Child group path. An empty string selects ``file`` itself.

        Returns
        -------
        h5py.File or h5py.Group
            Selected HDF5 object.
        """
        if not name:
            return file
        return file[name]

    @staticmethod
    def _create_group(file, name):
        """Create a named child group.

        Parameters
        ----------
        file : h5py.File or h5py.Group
            Root HDF5 object.
        name : str
            Child group path. An empty string selects ``file`` itself.

        Returns
        -------
        h5py.File or h5py.Group
            Created group, or ``file`` when ``name`` is empty.
        """
        if not name:
            return file
        return file.create_group(name)

    def from_system_mix(self, file_name, type_map=None, **kwargs):
        """Load unlabeled mixed HDF5 data and split it into Systems.

        Parameters
        ----------
        file_name : str or h5py.Group or h5py.File
            HDF5 file, HDF5 group, or string in ``"file.hdf5#group"`` form.
        type_map : list[str], optional
            Type map used to remap real atom types while loading.
        **kwargs : dict
            Additional keyword arguments accepted for format API compatibility.

        Returns
        -------
        list[dict]
            Unlabeled System data dicts reconstructed from the mixed data.
        """
        return self._from_system_mix(file_name, type_map=type_map, labels=False)

    def from_labeled_system_mix(self, file_name, type_map=None, **kwargs):
        """Load labeled mixed HDF5 data and split it into LabeledSystems.

        Parameters
        ----------
        file_name : str or h5py.Group or h5py.File
            HDF5 file, HDF5 group, or string in ``"file.hdf5#group"`` form.
        type_map : list[str], optional
            Type map used to remap real atom types while loading.
        **kwargs : dict
            Additional keyword arguments accepted for format API compatibility.

        Returns
        -------
        list[dict]
            LabeledSystem data dicts reconstructed from the mixed data.
        """
        return self._from_system_mix(file_name, type_map=type_map, labels=True)

    def _from_system_mix(self, file_name, type_map=None, labels=True):
        """Load mixed HDF5 data through the shared mixed backend.

        Parameters
        ----------
        file_name : str or h5py.Group or h5py.File
            HDF5 file, HDF5 group, or string in ``"file.hdf5#group"`` form.
            When a file object is given, the object itself is interpreted as the
            mixed system group.
        type_map : list[str], optional
            Type map used to remap real atom types while loading.
        labels : bool, default=True
            Whether labeled data such as energies and forces should be loaded.

        Returns
        -------
        list[dict]
            System or LabeledSystem data dicts split out of the mixed HDF5 data.

        Raises
        ------
        TypeError
            If ``file_name`` is not a string, HDF5 group, or HDF5 file.
        """
        import h5py

        register_spin()

        if isinstance(file_name, (h5py.Group, h5py.File)):
            return dpdata.formats.deepmd.mixed.to_system_data(
                file_name,
                type_map=type_map,
                labels=labels,
                load_func=self._load_hdf5_mixed_data,
            )
        elif isinstance(file_name, str):
            s = file_name.split("#")
            name = s[1] if len(s) > 1 else ""
            with h5py.File(s[0], "r") as f:
                return dpdata.formats.deepmd.mixed.to_system_data(
                    self._get_group(f, name),
                    type_map=type_map,
                    labels=labels,
                    load_func=self._load_hdf5_mixed_data,
                )
        else:
            raise TypeError("Unsupported file_name")

    def to_system(
        self,
        data,
        file_name,
        set_size: int = 2000,
        prec=np.float64,
        comp_prec=None,
        **kwargs,
    ):
        """Dump a System data dict in mixed HDF5 format.

        Parameters
        ----------
        data : dict
            System or LabeledSystem data dict. If it is not already in mixed
            type form, it is copied and converted before dumping.
        file_name : str or h5py.Group or h5py.File
            HDF5 file, HDF5 group, or string in ``"file.hdf5#group"`` form.
            Strings open the target file in write mode. HDF5 objects are written
            in place.
        set_size : int, default=2000
            Maximum number of frames per ``set.*`` group.
        prec : numpy.dtype, default=numpy.float64
            Floating point precision for dumped frame data. Kept for
            consistency with ``deepmd/npy/mixed``.
        comp_prec : numpy.dtype, optional
            Explicit floating point precision. When provided, this overrides
            ``prec``.
        **kwargs : dict
            Additional keyword arguments accepted for format API compatibility.

        Raises
        ------
        TypeError
            If ``file_name`` is not a string, HDF5 group, or HDF5 file.
        """
        import h5py

        if comp_prec is None:
            comp_prec = prec

        if isinstance(file_name, (h5py.Group, h5py.File)):
            dpdata.formats.deepmd.mixed.dump(
                file_name,
                data,
                set_size=set_size,
                comp_prec=comp_prec,
                dump_func=self._dump_hdf5_mixed_data,
            )
        elif isinstance(file_name, str):
            s = file_name.split("#")
            name = s[1] if len(s) > 1 else ""
            with h5py.File(s[0], "w") as f:
                dpdata.formats.deepmd.mixed.dump(
                    self._create_group(f, name),
                    data,
                    set_size=set_size,
                    comp_prec=comp_prec,
                    dump_func=self._dump_hdf5_mixed_data,
                )
        else:
            raise TypeError("Unsupported file_name")

    def from_multi_systems(self, directory, **kwargs):
        """Generate mixed HDF5 groups for MultiSystems loading.

        Parameters
        ----------
        directory : str or h5py.Group or h5py.File
            HDF5 file, HDF5 group, or string in ``"file.hdf5#group"`` form. The
            selected object may be either one mixed system group or a container
            of mixed groups.
        **kwargs : dict
            Additional keyword arguments accepted for format API compatibility.

        Yields
        ------
        h5py.Group or h5py.File
            Mixed HDF5 groups that will be passed to ``from_system_mix``.

        Raises
        ------
        TypeError
            If ``directory`` is not a string, HDF5 group, or HDF5 file.
        """
        import h5py

        register_spin()

        if isinstance(directory, (h5py.Group, h5py.File)):
            yield from self._iter_mixed_groups(directory)
        elif isinstance(directory, str):
            s = directory.split("#")
            name = s[1] if len(s) > 1 else ""
            with h5py.File(s[0], "r") as f:
                yield from self._iter_mixed_groups(self._get_group(f, name))
        else:
            raise TypeError("Unsupported directory")

    def to_multi_systems(self, formulas, directory, **kwargs):
        """Generate HDF5 groups for MultiSystems mixed dumping.

        Parameters
        ----------
        formulas : list[str]
            Mixed group names produced by ``mix_system``. For mixed HDF5 these
            names are atom counts after optional padding.
        directory : str or h5py.Group or h5py.File
            HDF5 file, HDF5 group, or string in ``"file.hdf5#group"`` form.
            Strings open the target file in write mode.
        **kwargs : dict
            Additional keyword arguments accepted for format API compatibility.

        Yields
        ------
        h5py.Group
            Destination groups that will be passed to ``to_system``.

        Raises
        ------
        TypeError
            If ``directory`` is not a string, HDF5 group, or HDF5 file.
        """
        import h5py

        if isinstance(directory, (h5py.Group, h5py.File)):
            for ff in formulas:
                if ff in directory:
                    del directory[ff]
                yield directory.create_group(ff)
        elif isinstance(directory, str):
            s = directory.split("#")
            name = s[1] if len(s) > 1 else ""
            with h5py.File(s[0], "w") as f:
                root = self._create_group(f, name)
                for ff in formulas:
                    yield root.create_group(ff)
        else:
            raise TypeError("Unsupported directory")


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
