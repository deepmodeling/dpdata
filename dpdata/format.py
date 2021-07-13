"""Implement the format plugin system."""
import os
from collections import abc
from abc import ABC

from .plugin import Plugin


class Format(ABC):
    __FormatPlugin = Plugin()
    __FromPlugin = Plugin()
    __ToPlugin = Plugin()

    @staticmethod
    def register(key):
        return Format.__FormatPlugin.register(key)

    @staticmethod
    def register_from(key):
        return Format.__FromPlugin.register(key)

    @staticmethod
    def register_to(key):
        return Format.__ToPlugin.register(key)
    
    @staticmethod
    def get_formats():
        return Format.__FormatPlugin.plugins

    @staticmethod
    def get_from_methods():
        return Format.__FromPlugin.plugins

    @staticmethod
    def get_to_methods():
        return Format.__ToPlugin.plugins
    
    @staticmethod
    def post(func_name):
        def decorator(object):
            if not isinstance(func_name, (list, tuple, set)):
                object.post_func = (func_name,)
            else:
                object.post_func = func_name
            return object
        return decorator

    def from_system(self, file_name, **kwargs):
        """System.from

        Parameters
        ----------
        file_name: str
            file name

        Returns
        -------
        data: dict
            system data
        """
        raise NotImplementedError("%s doesn't support System.from" %(self.__class__.__name__))

    def to_system(self, data, *args, **kwargs):
        """System.to

        Parameters
        ----------
        data: dict
            system data
        """
        raise NotImplementedError("%s doesn't support System.to" %(self.__class__.__name__))

    def from_labeled_system(self, file_name, **kwargs):
        raise NotImplementedError("%s doesn't support LabeledSystem.from" %(self.__class__.__name__))

    def to_labeled_system(self, data, *args, **kwargs):
        return self.to_system(data, *args, **kwargs)

    def from_bond_order_system(self, file_name, **kwargs):
        raise NotImplementedError("%s doesn't support BondOrderSystem.from" %(self.__class__.__name__))

    def to_bond_order_system(self, data, rdkit_mol, *args, **kwargs):
        return self.to_system(data, *args, **kwargs)

    class MultiModes:
        """File mode for MultiSystems
        0 (default): not implemented
        1: every directory under the top-level directory is a system
        """
        NotImplemented = 0
        Directory = 1

    MultiMode = MultiModes.NotImplemented

    def from_multi_systems(self, directory, **kwargs):
        """MultiSystems.from
        
        Parameters
        ----------
        directory: str
            directory of system
        
        Returns
        -------
        filenames: list[str]
            list of filenames
        """
        if self.MultiMode == self.MultiModes.Directory:
            return [os.path.join(directory, name) for name in os.listdir(directory) if os.path.isdir(os.path.join(directory, name))]
        raise NotImplementedError("%s doesn't support MultiSystems.from" %(self.__class__.__name__))

    def to_multi_systems(self, formulas, directory, **kwargs):
        if self.MultiMode == self.MultiModes.Directory:
            return [os.path.join(directory, ff) for ff in formulas]
        raise NotImplementedError("%s doesn't support MultiSystems.to" %(self.__class__.__name__))
        
