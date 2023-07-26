from enum import Enum, unique
from typing import TYPE_CHECKING, Tuple

import numpy as np

from dpdata.plugin import Plugin

if TYPE_CHECKING:
    from dpdata.system import System


@unique
class Axis(Enum):
    """Data axis."""

    NFRAMES = "nframes"
    NATOMS = "natoms"
    NTYPES = "ntypes"
    NBONDS = "nbonds"


class AnyInt(int):
    """AnyInt equals to any other integer."""

    def __eq__(self, other):
        return True


class DataError(Exception):
    """Data is not correct."""


class DataType:
    """DataType represents a type of data, like coordinates, energies, etc.

    Parameters
    ----------
    name : str
        name of data
    dtype : type or tuple[type]
        data type, e.g. np.ndarray
    shape : tuple[int], optional
        shape of data. Used when data is list or np.ndarray. Use Axis to
        represents numbers
    required : bool, default=True
        whether this data is required
    """

    def __init__(
        self,
        name: str,
        dtype: type,
        shape: Tuple[int, Axis] = None,
        required: bool = True,
    ) -> None:
        self.name = name
        self.dtype = dtype
        self.shape = shape
        self.required = required

    def real_shape(self, system: "System") -> Tuple[int]:
        """Returns expected real shape of a system."""
        shape = []
        for ii in self.shape:
            if ii is Axis.NFRAMES:
                shape.append(system.get_nframes())
            elif ii is Axis.NTYPES:
                shape.append(system.get_ntypes())
            elif ii is Axis.NATOMS:
                shape.append(system.get_natoms())
            elif ii is Axis.NBONDS:
                # BondOrderSystem
                shape.append(system.get_nbonds())
            elif ii == -1:
                shape.append(AnyInt(-1))
            elif isinstance(ii, int):
                shape.append(ii)
            else:
                raise RuntimeError("Shape is not an int!")
        return tuple(shape)

    def check(self, system: "System"):
        """Check if a system has correct data of this type.

        Parameters
        ----------
        system : System
            checked system

        Raises
        ------
        DataError
            type or shape of data is not correct
        """
        # check if exists
        if self.name in system.data:
            data = system.data[self.name]
            # check dtype
            # allow list for empty np.ndarray
            if isinstance(data, list) and not len(data):
                pass
            elif not isinstance(data, self.dtype):
                raise DataError(
                    f"Type of {self.name} is {type(data).__name__}, but expected {self.dtype.__name__}"
                )
            # check shape
            if self.shape is not None:
                shape = self.real_shape(system)
                # skip checking empty list of np.ndarray
                if isinstance(data, np.ndarray):
                    if data.size and shape != data.shape:
                        raise DataError(
                            f"Shape of {self.name} is {data.shape}, but expected {shape}"
                        )
                elif isinstance(data, list):
                    if len(shape) and shape[0] != len(data):
                        raise DataError(
                            "Length of %s is %d, but expected %d"
                            % (self.name, len(data), shape[0])
                        )
                else:
                    raise RuntimeError("Unsupported type to check shape")
        elif self.required:
            raise DataError("%s not found in data" % self.name)


__system_data_type_plugin = Plugin()
__labeled_system_data_type_plugin = Plugin()


def register_data_type(data_type: DataType, labeled: bool):
    """Register a data type.

    Parameters
    ----------
    data_type : DataType
        data type to be registered
    labeled : bool
        whether this data type is for LabeledSystem
    """
    plugin = __labeled_system_data_type_plugin if labeled else __system_data_type_plugin
    plugin.register(data_type.name)(data_type)


def get_data_types(labeled: bool):
    """Get all registered data types.

    Parameters
    ----------
    labeled : bool
        whether this data type is for LabeledSystem
    """
    plugin = __labeled_system_data_type_plugin if labeled else __system_data_type_plugin
    return tuple(plugin.plugins.values())
