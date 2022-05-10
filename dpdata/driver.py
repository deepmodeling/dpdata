"""Driver plugin system."""
from typing import Callable
from .plugin import Plugin
from abc import ABC, abstractmethod


class Driver(ABC):
    """The base class for a driver plugin. A driver can
    label a pure System to generate the LabeledSystem.

    See Also
    --------
    dpdata.plugins.deepmd.DPDriver : an example of Driver
    """
    __DriverPlugin = Plugin()

    @staticmethod
    def register(key: str) -> Callable:
        """Register a driver plugin. Used as decorators.
        
        Parameter
        ---------
        key: str
            key of the plugin.
        
        Returns
        -------
        Callable
            decorator of a class

        Examples
        --------
        >>> @Driver.register("some_driver")
        ... class SomeDriver(Driver):
        ...     pass
        """
        return Driver.__DriverPlugin.register(key)

    @staticmethod
    def get_driver(key: str) -> "Driver":
        """Get a driver plugin.
        
        Parameter
        ---------
        key: str
            key of the plugin.
        
        Returns
        -------
        Driver
            the specific driver class
        
        Raises
        ------
        RuntimeError
            if the requested driver is not implemented
        """
        try:
            return Driver.__DriverPlugin.plugins[key]
        except KeyError as e:
            raise RuntimeError('Unknown driver: ' + key) from e
    
    def __init__(self, *args, **kwargs) -> None:
        """Setup the driver."""

    @abstractmethod
    def label(self, data: dict) -> dict:
        """Label a system data. Returns new data with energy, forces, and virials.
        
        Parameters
        ----------
        data : dict
            data with coordinates and atom types
        
        Returns
        -------
        dict
            labeled data with energies and forces
        """
        return NotImplemented
