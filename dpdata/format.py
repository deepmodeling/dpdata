"""Implement the format plugin system."""
import os
from abc import ABC

from .plugin import Plugin


class Format(ABC):
    """The abstract base class for all formats.

    To add a new format, one should create a new class inherited from this class, and then

    - implement several methods, such as :meth:`from_system`;
    - register the format with a key;
    - add documentation in the class docstring;

    The new format can be either insider or outside the package.
    """

    __FormatPlugin = Plugin()
    __FromPlugin = Plugin()
    __ToPlugin = Plugin()

    @staticmethod
    def register(key):
        """Register a format plugin.

        By default, after a format plugin is registered, the following methods
        will be registered as well for :meth:`System`, :meth:`LabeledSystem`, :meth:`MultiSystems`, and
        :meth:`BondOrderSystem`:

        - from_{key.replace('/', '_')}
        - to_{key.replace('/', '_')}
        - from({key}, ...)
        - to({key}, ...)

        The decorator should be explicitly executed before :mod:`dpdata.system`
        is imported. A module will be imported automatically if it

        - is a submodule of :mod:`dpdata.plugins`;
        - is registered at the `dpdata.plugins` entry point

        Parameters
        ----------
        key : str
            The key to register the plugin.

        Returns
        -------
        function
            The decorator function.

        Examples
        --------
        Register a format plugin:

        >>> @Format.register('test')
        ... @Format.register('test2')
        ... class TestFormat(Format):
        ...     pass
        """
        return Format.__FormatPlugin.register(key)

    @staticmethod
    def register_from(key):
        """Register a from method if the target method name is not default.

        Parameters
        ----------
        key : str
            The key to register the plugin.

        Returns
        -------
        function
            The decorator function.

        Examples
        --------
        Register a from method:

        >>> @Format.register_from('from_test_haha')
        ... @Format.register('test)
        ... class TestFormat(Format):
        ...     pass

        This will register a from method named from_test_haha, although the
        format name is test.
        """
        return Format.__FromPlugin.register(key)

    @staticmethod
    def register_to(key):
        """Register a to method if the target method name is not default.

        Parameters
        ----------
        key : str
            The key to register the plugin.

        Returns
        -------
        function
            The decorator function.

        Examples
        --------
        Register a to method:

        >>> @Format.register_to('to_test_haha')
        ... @Format.register('test')
        ... class TestFormat(Format):
        ...     pass

        This will register a to method named to_test_haha, although the
        format name is test.
        """
        return Format.__ToPlugin.register(key)

    @staticmethod
    def get_formats():
        """Get all registered formats."""
        return Format.__FormatPlugin.plugins

    @staticmethod
    def get_from_methods():
        """Get all registered from methods."""
        return Format.__FromPlugin.plugins

    @staticmethod
    def get_to_methods():
        """Get all registered to methods."""
        return Format.__ToPlugin.plugins

    @staticmethod
    def post(func_name):
        """Register a post function for from method.

        Such function will be called after the "from" method is called.

        Parameters
        ----------
        func_name : str or list of str
            The name of the post function.

        Returns
        -------
        function
            The decorator function.

        Examples
        --------
        Register a post function:

        >>> @Format.post('remove_pbc')
        ... @Format.register('test')
        ... class TestFormat(Format):
        ...     pass
        """

        def decorator(object):
            if not isinstance(func_name, (list, tuple, set)):
                object.post_func = (func_name,)
            else:
                object.post_func = func_name
            return object

        return decorator

    def from_system(self, file_name, **kwargs):
        """Implement System.from that converts from this format to System.

        Parameters
        ----------
        file_name : str
            file name, i.e. the first argument
        **kwargs : dict
            keyword arguments that will be passed from the method

        Returns
        -------
        data : dict
            system data, whose keys are defined in System.DTYPES
        """
        raise NotImplementedError(
            "%s doesn't support System.from" % (self.__class__.__name__)
        )

    def to_system(self, data, *args, **kwargs):
        """Implement System.to that converts from System to this format.

        Parameters
        ----------
        data : dict
            system data, whose keys are defined in System.DTYPES
        *args : list
            arguments that will be passed from the method
        **kwargs : dict
            keyword arguments that will be passed from the method
        """
        raise NotImplementedError(
            "%s doesn't support System.to" % (self.__class__.__name__)
        )

    def from_labeled_system(self, file_name, **kwargs):
        """Implement LabeledSystem.from that converts from this format to LabeledSystem.

        Parameters
        ----------
        file_name : str
            file name, i.e. the first argument
        **kwargs : dict
            keyword arguments that will be passed from the method

        Returns
        -------
        data : dict
            system data, whose keys are defined in LabeledSystem.DTYPES
        """
        raise NotImplementedError(
            "%s doesn't support LabeledSystem.from" % (self.__class__.__name__)
        )

    def to_labeled_system(self, data, *args, **kwargs):
        """Implement LabeledSystem.to that converts from LabeledSystem to this format.

        By default, LabeledSystem.to will fallback to System.to.

        Parameters
        ----------
        data : dict
            system data, whose keys are defined in LabeledSystem.DTYPES
        *args : list
            arguments that will be passed from the method
        **kwargs : dict
            keyword arguments that will be passed from the method
        """
        return self.to_system(data, *args, **kwargs)

    def from_bond_order_system(self, file_name, **kwargs):
        """Implement BondOrderSystem.from that converts from this format to BondOrderSystem.

        Parameters
        ----------
        file_name : str
            file name, i.e. the first argument
        **kwargs : dict
            keyword arguments that will be passed from the method

        Returns
        -------
        data : dict
            system data
        """
        raise NotImplementedError(
            "%s doesn't support BondOrderSystem.from" % (self.__class__.__name__)
        )

    def to_bond_order_system(self, data, rdkit_mol, *args, **kwargs):
        """Implement BondOrderSystem.to that converts from BondOrderSystem to this format.

        By default, BondOrderSystem.to will fallback to LabeledSystem.to.

        Parameters
        ----------
        data : dict
            system data
        rdkit_mol : rdkit.Chem.rdchem.Mol
            rdkit mol object
        *args : list
            arguments that will be passed from the method
        **kwargs : dict
            keyword arguments that will be passed from the method
        """
        return self.to_system(data, *args, **kwargs)

    class MultiModes:
        """File mode for MultiSystems.

        The current implemented modes are:

            - 0 (default): not implemented
            - 1: every directory under the top-level directory is a system.
        """

        NotImplemented = 0
        Directory = 1

    MultiMode = MultiModes.NotImplemented

    def from_multi_systems(self, directory, **kwargs):
        """Implement MultiSystems.from that converts from this format to MultiSystems.

        By default, this method follows MultiMode to implement the conversion.

        Parameters
        ----------
        directory : str
            directory of system
        **kwargs : dict
            keyword arguments that will be passed from the method

        Returns
        -------
        filenames: list[str]
            list of filenames
        """
        if self.MultiMode == self.MultiModes.Directory:
            return [
                os.path.join(directory, name)
                for name in os.listdir(directory)
                if os.path.isdir(os.path.join(directory, name))
            ]
        raise NotImplementedError(
            "%s doesn't support MultiSystems.from" % (self.__class__.__name__)
        )

    def to_multi_systems(self, formulas, directory, **kwargs):
        """Implement MultiSystems.to that converts from MultiSystems to this format.

        By default, this method follows MultiMode to implement the conversion.

        Parameters
        ----------
        formulas : list[str]
            list of formulas
        directory : str
            directory of system
        **kwargs : dict
            keyword arguments that will be passed from the method
        """
        if self.MultiMode == self.MultiModes.Directory:
            return [os.path.join(directory, ff) for ff in formulas]
        raise NotImplementedError(
            "%s doesn't support MultiSystems.to" % (self.__class__.__name__)
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
            keyword arguments that will be passed from the method

        Returns
        -------
        mixed_systems: dict
            dict of mixed system with key 'atom_numbs'
        """
        raise NotImplementedError(
            "%s doesn't support System.from" % (self.__class__.__name__)
        )
