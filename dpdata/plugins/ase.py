from dpdata.driver import Driver
from dpdata.format import Format
import numpy as np
import dpdata
try:
    import ase.io
    from ase.calculators.calculator import PropertyNotImplementedError
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

        Returns
        -------
        dict
            data dict
        """
        symbols = atoms.get_chemical_symbols()
        atom_names = list(set(symbols))
        atom_numbs = [symbols.count(symbol) for symbol in atom_names]
        atom_types = np.array([atom_names.index(symbol) for symbol in symbols]).astype(int)
        cells = atoms.cell[:]
        coords = atoms.get_positions()
        info_dict = {
            'atom_names': atom_names,
            'atom_numbs': atom_numbs,
            'atom_types': atom_types,
            'cells': np.array([cells]).astype('float32'),
            'coords': np.array([coords]).astype('float32'),
            'orig': np.zeros(3),
        }
        return info_dict

    def from_labeled_system(self, atoms: "ase.Atoms", **kwargs) -> dict:
        """Convert ase.Atoms to a LabeledSystem. Energies and forces
        are calculated by the calculator.

        Parameters
        ----------
        atoms : ase.Atoms
            an ASE Atoms, containing a structure

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
            ** info_dict,
            'energies': np.array([energies]).astype('float32'),
            'forces': np.array([forces]).astype('float32'),
        }
        try:
            stress = atoms.get_stress(False)
        except PropertyNotImplementedError:
            pass
        else:
            virials = np.array([-atoms.get_volume() * stress]).astype('float32')
            info_dict['virials'] = virials
        return info_dict

    def from_multi_systems(self, file_name: str, begin: int=None, end: int=None, step: int=None, ase_fmt: str=None, **kwargs) -> "ase.Atoms":
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

        Yields
        ------
        ase.Atoms
            ASE atoms in the file
        """
        frames = ase.io.read(file_name, format=ase_fmt, index=slice(begin, end, step))
        for atoms in frames:
            yield atoms

    def to_system(self, data, **kwargs):
        '''
        convert System to ASE Atom obj

        '''
        from ase import Atoms

        structures = []
        species = [data['atom_names'][tt] for tt in data['atom_types']]

        for ii in range(data['coords'].shape[0]):
            structure = Atoms(
                symbols=species, positions=data['coords'][ii], pbc=not data.get('nopbc', False), cell=data['cells'][ii])
            structures.append(structure)

        return structures

    def to_labeled_system(self, data, *args, **kwargs):
        '''Convert System to ASE Atoms object.'''
        from ase import Atoms
        from ase.calculators.singlepoint import SinglePointCalculator

        structures = []
        species = [data['atom_names'][tt] for tt in data['atom_types']]

        for ii in range(data['coords'].shape[0]):
            structure = Atoms(
                symbols=species,
                positions=data['coords'][ii],
                pbc=not data.get('nopbc', False),
                cell=data['cells'][ii]
            )

            results = {
                'energy': data["energies"][ii],
                'forces': data["forces"][ii]
            }
            if "virials" in data:
                # convert to GPa as this is ase convention
                # v_pref = 1 * 1e4 / 1.602176621e6
                vol = structure.get_volume()
                # results['stress'] = data["virials"][ii] / (v_pref * vol)
                results['stress'] = -data["virials"][ii] / vol

            structure.calc = SinglePointCalculator(structure, **results)
            structures.append(structure)

        return structures


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
            ls = dpdata.LabeledSystem(atoms, fmt="ase/structure", type_map=data['atom_names'])
            labeled_system.append(ls)
        return labeled_system.data
