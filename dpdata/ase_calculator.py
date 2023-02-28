from typing import TYPE_CHECKING, List, Optional

from ase.calculators.calculator import (
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

    name = "dpdata"
    implemented_properties = ["energy", "free_energy", "forces", "virial", "stress"]

    def __init__(self, driver: Driver, **kwargs) -> None:
        Calculator.__init__(self, label=Driver.__name__, **kwargs)
        self.driver = driver

    def calculate(
        self,
        atoms: Optional["Atoms"] = None,
        properties: List[str] = ["energy", "forces"],
        system_changes: List[str] = all_changes,
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
        if atoms is not None:
            self.atoms = atoms.copy()

        system = dpdata.System(self.atoms, fmt="ase/structure")
        data = system.predict(driver=self.driver).data

        self.results["energy"] = data["energies"][0]
        # see https://gitlab.com/ase/ase/-/merge_requests/2485
        self.results["free_energy"] = data["energies"][0]
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
