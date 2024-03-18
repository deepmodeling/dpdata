from typing import TYPE_CHECKING, Optional, Type

import numpy as np

import dpdata
from dpdata.driver import Driver, Minimizer
from dpdata.format import Format

try:
    import ase.io
    from ase.calculators.calculator import PropertyNotImplementedError
    from ase.io import Trajectory

    if TYPE_CHECKING:
        from ase.optimize.optimize import Optimizer
except ImportError:
    pass


@Format.register("ase/structure")
class ASEStructureFormat(Format):
    """Format for the `Atomic Simulation Environment <https://wiki.fysik.dtu.dk/ase/>`_ (ase).

    ASE supports parsing a few dozen of data formats. As described in i
    `the documentation <ihttps://wiki.fysik.dtu.dk/ase/ase/io/io.html>`_,
    many of these formats can be determined automatically.
    Use the `ase_fmt` keyword argument to supply the format if
    automatic detection fails.
    """

    def from_system(self, atoms: "ase.Atoms", **kwargs) -> dict:
        """Convert ase.Atoms to a System.

        Parameters
        ----------
        atoms : ase.Atoms
            an ASE Atoms, containing a structure
        **kwargs : dict
            other parameters

        Returns
        -------
        dict
            data dict
        """
        symbols = atoms.get_chemical_symbols()
        atom_names = list(dict.fromkeys(symbols))
        atom_numbs = [symbols.count(symbol) for symbol in atom_names]
        atom_types = np.array([atom_names.index(symbol) for symbol in symbols]).astype(
            int
        )
        cells = atoms.cell[:]
        coords = atoms.get_positions()
        info_dict = {
            "atom_names": atom_names,
            "atom_numbs": atom_numbs,
            "atom_types": atom_types,
            "cells": np.array([cells]),
            "coords": np.array([coords]),
            "orig": np.zeros(3),
            "nopbc": not np.any(atoms.get_pbc()),
        }
        return info_dict

    def from_labeled_system(self, atoms: "ase.Atoms", **kwargs) -> dict:
        """Convert ase.Atoms to a LabeledSystem. Energies and forces
        are calculated by the calculator.

        Parameters
        ----------
        atoms : ase.Atoms
            an ASE Atoms, containing a structure
        **kwargs : dict
            other parameters

        Returns
        -------
        dict
            data dict

        Raises
        ------
        RuntimeError
            ASE will raise RuntimeError if the atoms does not
            have a calculator
        """
        info_dict = self.from_system(atoms)
        try:
            energies = atoms.get_potential_energy(force_consistent=True)
        except PropertyNotImplementedError:
            energies = atoms.get_potential_energy()
        forces = atoms.get_forces()
        info_dict = {
            **info_dict,
            "energies": np.array([energies]),
            "forces": np.array([forces]),
        }
        try:
            stress = atoms.get_stress(False)
        except PropertyNotImplementedError:
            pass
        else:
            virials = np.array([-atoms.get_volume() * stress])
            info_dict["virials"] = virials
        return info_dict

    def from_multi_systems(
        self,
        file_name: str,
        begin: Optional[int] = None,
        end: Optional[int] = None,
        step: Optional[int] = None,
        ase_fmt: Optional[str] = None,
        **kwargs,
    ) -> "ase.Atoms":
        """Convert a ASE supported file to ASE Atoms.

        It will finally be converted to MultiSystems.

        Parameters
        ----------
        file_name : str
            path to file
        begin : int, optional
            begin frame index
        end : int, optional
            end frame index
        step : int, optional
            frame index step
        ase_fmt : str, optional
            ASE format. See the ASE documentation about supported formats
        **kwargs : dict
            other parameters

        Yields
        ------
        ase.Atoms
            ASE atoms in the file
        """
        frames = ase.io.read(file_name, format=ase_fmt, index=slice(begin, end, step))
        yield from frames

    def to_system(self, data, **kwargs):
        """Convert System to ASE Atom obj."""
        from ase import Atoms

        structures = []
        species = [data["atom_names"][tt] for tt in data["atom_types"]]

        for ii in range(data["coords"].shape[0]):
            structure = Atoms(
                symbols=species,
                positions=data["coords"][ii],
                pbc=not data.get("nopbc", False),
                cell=data["cells"][ii],
            )
            structures.append(structure)

        return structures

    def to_labeled_system(self, data, *args, **kwargs):
        """Convert System to ASE Atoms object."""
        from ase import Atoms
        from ase.calculators.singlepoint import SinglePointCalculator

        structures = []
        species = [data["atom_names"][tt] for tt in data["atom_types"]]

        for ii in range(data["coords"].shape[0]):
            structure = Atoms(
                symbols=species,
                positions=data["coords"][ii],
                pbc=not data.get("nopbc", False),
                cell=data["cells"][ii],
            )

            results = {"energy": data["energies"][ii], "forces": data["forces"][ii]}
            if "virials" in data:
                # convert to GPa as this is ase convention
                # v_pref = 1 * 1e4 / 1.602176621e6
                vol = structure.get_volume()
                # results['stress'] = data["virials"][ii] / (v_pref * vol)
                results["stress"] = -data["virials"][ii] / vol

            structure.calc = SinglePointCalculator(structure, **results)
            structures.append(structure)

        return structures


@Format.register("ase/traj")
class ASETrajFormat(Format):
    """Format for the ASE's trajectory format <https://wiki.fysik.dtu.dk/ase/ase/io/trajectory.html#module-ase.io.trajectory>`_ (ase).'
    a `traj' contains a sequence of frames, each of which is an `Atoms' object.
    """

    def from_system(
        self,
        file_name: str,
        begin: Optional[int] = 0,
        end: Optional[int] = None,
        step: Optional[int] = 1,
        **kwargs,
    ) -> dict:
        """Read ASE's trajectory file to `System` of multiple frames.

        Parameters
        ----------
        file_name : str
            ASE's trajectory file
        begin : int, optional
            begin frame index
        end : int, optional
            end frame index
        step : int, optional
            frame index step
        **kwargs : dict
            other parameters

        Returns
        -------
        dict_frames: dict
            a dictionary containing data of multiple frames
        """
        traj = Trajectory(file_name)
        sub_traj = traj[begin:end:step]
        dict_frames = ASEStructureFormat().from_system(sub_traj[0])
        for atoms in sub_traj[1:]:
            tmp = ASEStructureFormat().from_system(atoms)
            dict_frames["cells"] = np.append(dict_frames["cells"], tmp["cells"][0])
            dict_frames["coords"] = np.append(dict_frames["coords"], tmp["coords"][0])

        ## Correct the shape of numpy arrays
        dict_frames["cells"] = dict_frames["cells"].reshape(-1, 3, 3)
        dict_frames["coords"] = dict_frames["coords"].reshape(len(sub_traj), -1, 3)

        return dict_frames

    def from_labeled_system(
        self,
        file_name: str,
        begin: Optional[int] = 0,
        end: Optional[int] = None,
        step: Optional[int] = 1,
        **kwargs,
    ) -> dict:
        """Read ASE's trajectory file to `System` of multiple frames.

        Parameters
        ----------
        file_name : str
            ASE's trajectory file
        begin : int, optional
            begin frame index
        end : int, optional
            end frame index
        step : int, optional
            frame index step
        **kwargs : dict
            other parameters

        Returns
        -------
        dict_frames: dict
            a dictionary containing data of multiple frames
        """
        traj = Trajectory(file_name)
        sub_traj = traj[begin:end:step]

        ## check if the first frame has a calculator
        if sub_traj[0].calc is None:
            raise ValueError(
                "The input trajectory does not contain energies and forces, may not be a labeled system."
            )

        dict_frames = ASEStructureFormat().from_labeled_system(sub_traj[0])
        for atoms in sub_traj[1:]:
            tmp = ASEStructureFormat().from_labeled_system(atoms)
            dict_frames["cells"] = np.append(dict_frames["cells"], tmp["cells"][0])
            dict_frames["coords"] = np.append(dict_frames["coords"], tmp["coords"][0])
            dict_frames["energies"] = np.append(
                dict_frames["energies"], tmp["energies"][0]
            )
            dict_frames["forces"] = np.append(dict_frames["forces"], tmp["forces"][0])
            if "virials" in tmp.keys() and "virials" in dict_frames.keys():
                dict_frames["virials"] = np.append(
                    dict_frames["virials"], tmp["virials"][0]
                )

        ## Correct the shape of numpy arrays
        dict_frames["cells"] = dict_frames["cells"].reshape(-1, 3, 3)
        dict_frames["coords"] = dict_frames["coords"].reshape(len(sub_traj), -1, 3)
        dict_frames["forces"] = dict_frames["forces"].reshape(len(sub_traj), -1, 3)
        if "virials" in dict_frames.keys():
            dict_frames["virials"] = dict_frames["virials"].reshape(-1, 3, 3)

        return dict_frames


@Driver.register("ase")
class ASEDriver(Driver):
    """ASE Driver.

    Parameters
    ----------
    calculator : ase.calculators.calculator.Calculato
        ASE calculator
    """

    def __init__(self, calculator: "ase.calculators.calculator.Calculator") -> None:
        """Setup the driver."""
        self.calculator = calculator

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
        # convert data to ase data
        system = dpdata.System(data=data)
        # list[Atoms]
        structures = system.to_ase_structure()
        labeled_system = dpdata.LabeledSystem()
        for atoms in structures:
            atoms.calc = self.calculator
            ls = dpdata.LabeledSystem(
                atoms, fmt="ase/structure", type_map=data["atom_names"]
            )
            labeled_system.append(ls)
        return labeled_system.data


@Minimizer.register("ase")
class ASEMinimizer(Minimizer):
    """ASE minimizer.

    Parameters
    ----------
    driver : Driver
        dpdata driver
    optimizer : type, optional
        ase optimizer class
    fmax : float, optional, default=5e-3
        force convergence criterion
    max_steps : int, optional
        max steps to optimize
    optimizer_kwargs : dict, optional
        other parameters for optimizer
    """

    def __init__(
        self,
        driver: Driver,
        optimizer: Optional[Type["Optimizer"]] = None,
        fmax: float = 5e-3,
        max_steps: Optional[int] = None,
        optimizer_kwargs: dict = {},
    ) -> None:
        self.calculator = driver.ase_calculator
        if optimizer is None:
            from ase.optimize import LBFGS

            self.optimizer = LBFGS
        else:
            self.optimizer = optimizer
        self.optimizer_kwargs = {
            "logfile": None,
            **optimizer_kwargs.copy(),
        }
        self.fmax = fmax
        self.max_steps = max_steps

    def minimize(self, data: dict) -> dict:
        """Minimize the geometry.

        Parameters
        ----------
        data : dict
            data with coordinates and atom types

        Returns
        -------
        dict
            labeled data with minimized coordinates, energies, and forces
        """
        system = dpdata.System(data=data)
        # list[Atoms]
        structures = system.to_ase_structure()
        labeled_system = dpdata.LabeledSystem()
        for atoms in structures:
            atoms.calc = self.calculator
            dyn = self.optimizer(atoms, **self.optimizer_kwargs)
            dyn.run(fmax=self.fmax, steps=self.max_steps)
            ls = dpdata.LabeledSystem(
                atoms, fmt="ase/structure", type_map=data["atom_names"]
            )
            labeled_system.append(ls)
        return labeled_system.data
