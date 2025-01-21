from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.calculator import (  # noqa: TID253
    Calculator,
    PropertyNotImplementedError,
    all_changes,
)

import dpdata

from .driver import Driver

if TYPE_CHECKING:
    from ase import Atoms


class DPDataCalculator(Calculator):
    """Implementation of ASE deepmd calculator based on a driver.

    Parameters
    ----------
    driver : Driver
        dpdata driver
    """

    @property
    def name(self) -> str:
        return "dpdata"

    implemented_properties = ["energy", "free_energy", "forces", "virial", "stress"]

    def __init__(self, driver: Driver, **kwargs) -> None:
        Calculator.__init__(self, label=Driver.__name__, **kwargs)
        self.driver = driver

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: list[str] = ["energy", "forces"],
        system_changes: list[str] = all_changes,
    ):
        """Run calculation with a driver.

        Parameters
        ----------
        atoms : Optional[Atoms], optional
            atoms object to run the calculation on, by default None
        properties : List[str], optional
            unused, only for function signature compatibility,
            by default ["energy", "forces"]
        system_changes : List[str], optional
            unused, only for function signature compatibility, by default all_changes
        """
        assert atoms is not None
        atoms = atoms.copy()

        system = dpdata.System(atoms, fmt="ase/structure")
        data = system.predict(driver=self.driver).data

        self.results["energy"] = data["energies"][0]
        # see https://gitlab.com/ase/ase/-/merge_requests/2485
        self.results["free_energy"] = data["energies"][0]
        if "forces" in data:
            self.results["forces"] = data["forces"][0]
        if "virials" in data:
            self.results["virial"] = data["virials"][0].reshape(3, 3)

        # convert virial into stress for lattice relaxation
        if "stress" in properties:
            if sum(atoms.get_pbc()) > 0:
                # the usual convention (tensile stress is positive)
                # stress = -virial / volume
                stress = (
                    -0.5
                    * (data["virials"][0].copy() + data["virials"][0].copy().T)
                    / atoms.get_volume()
                )
                # Voigt notation
                self.results["stress"] = stress.flat[[0, 4, 8, 5, 2, 1]]
            else:
                raise PropertyNotImplementedError
