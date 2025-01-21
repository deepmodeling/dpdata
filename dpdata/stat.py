from __future__ import annotations

from abc import ABCMeta, abstractmethod
from functools import lru_cache
from typing import Any

import numpy as np

from dpdata.system import LabeledSystem, MultiSystems


def mae(errors: np.ndarray) -> np.floating[Any]:
    """Compute the mean absolute error (MAE).

    Parameters
    ----------
    errors : np.ndarray
        errors between two values

    Returns
    -------
    floating[Any]
        mean absolute error (MAE)
    """
    return np.mean(np.abs(errors))


def rmse(errors: np.ndarray) -> np.floating[Any]:
    """Compute the root mean squared error (RMSE).

    Parameters
    ----------
    errors : np.ndarray
        errors between two values

    Returns
    -------
    floating[Any]
        root mean squared error (RMSE)
    """
    return np.sqrt(np.mean(np.square(errors)))


class ErrorsBase(metaclass=ABCMeta):
    """Compute errors (deviations) between two systems. The type of system is assigned by SYSTEM_TYPE.

    Parameters
    ----------
    system_1 : object
        system 1
    system_2 : object
        system 2
    """

    SYSTEM_TYPE = object

    def __init__(self, system_1: SYSTEM_TYPE, system_2: SYSTEM_TYPE) -> None:
        assert isinstance(system_1, self.SYSTEM_TYPE), (
            f"system_1 should be {self.SYSTEM_TYPE.__name__}"
        )
        assert isinstance(system_2, self.SYSTEM_TYPE), (
            f"system_2 should be {self.SYSTEM_TYPE.__name__}"
        )
        self.system_1 = system_1
        self.system_2 = system_2

    @property
    @abstractmethod
    def e_errors(self) -> np.ndarray:
        """Energy errors."""

    @property
    @abstractmethod
    def f_errors(self) -> np.ndarray:
        """Force errors."""

    @property
    def e_mae(self) -> np.floating[Any]:
        """Energy MAE."""
        return mae(self.e_errors)

    @property
    def e_rmse(self) -> np.floating[Any]:
        """Energy RMSE."""
        return rmse(self.e_errors)

    @property
    def f_mae(self) -> np.floating[Any]:
        """Force MAE."""
        return mae(self.f_errors)

    @property
    def f_rmse(self) -> np.floating[Any]:
        """Force RMSE."""
        return rmse(self.f_errors)


class Errors(ErrorsBase):
    """Compute errors (deviations) between two LabeledSystems.

    Parameters
    ----------
    system_1 : object
        system 1
    system_2 : object
        system 2

    Examples
    --------
    Get errors between referenced system and predicted system:

    >>> e = dpdata.stat.Errors(system_1, system_2)
    >>> print("%.4f %.4f %.4f %.4f" % (e.e_mae, e.e_rmse, e.f_mae, e.f_rmse))
    """

    SYSTEM_TYPE = LabeledSystem

    @property
    @lru_cache
    def e_errors(self) -> np.ndarray:
        """Energy errors."""
        assert isinstance(self.system_1, self.SYSTEM_TYPE)
        assert isinstance(self.system_2, self.SYSTEM_TYPE)
        return self.system_1["energies"] - self.system_2["energies"]

    @property
    @lru_cache
    def f_errors(self) -> np.ndarray:
        """Force errors."""
        assert isinstance(self.system_1, self.SYSTEM_TYPE)
        assert isinstance(self.system_2, self.SYSTEM_TYPE)
        return (self.system_1["forces"] - self.system_2["forces"]).ravel()


class MultiErrors(ErrorsBase):
    """Compute errors (deviations) between two MultiSystems.

    Parameters
    ----------
    system_1 : object
        system 1
    system_2 : object
        system 2

    Examples
    --------
    Get errors between referenced system and predicted system:

    >>> e = dpdata.stat.MultiErrors(system_1, system_2)
    >>> print("%.4f %.4f %.4f %.4f" % (e.e_mae, e.e_rmse, e.f_mae, e.f_rmse))
    """

    SYSTEM_TYPE = MultiSystems

    @property
    @lru_cache
    def e_errors(self) -> np.ndarray:
        """Energy errors."""
        assert isinstance(self.system_1, self.SYSTEM_TYPE)
        assert isinstance(self.system_2, self.SYSTEM_TYPE)
        errors = []
        for nn in self.system_1.systems.keys():
            ss1 = self.system_1[nn]
            ss2 = self.system_2[nn]
            errors.append(Errors(ss1, ss2).e_errors.ravel())
        return np.concatenate(errors)

    @property
    @lru_cache
    def f_errors(self) -> np.ndarray:
        """Force errors."""
        assert isinstance(self.system_1, self.SYSTEM_TYPE)
        assert isinstance(self.system_2, self.SYSTEM_TYPE)
        errors = []
        for nn in self.system_1.systems.keys():
            ss1 = self.system_1[nn]
            ss2 = self.system_2[nn]
            errors.append(Errors(ss1, ss2).f_errors.ravel())
        return np.concatenate(errors)
