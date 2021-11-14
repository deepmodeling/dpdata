from dpdata.format import Format
import numpy as np
try:
    import ase.io
    from ase.calculators.calculator import PropertyNotImplementedError
except ImportError:
    pass


@Format.register("ase/structure")
class ASEStructureFormat(Format):
    def from_labeled_system(self, data, **kwargs):
        return data

    def from_multi_systems(self, file_name, begin=None, end=None, step=None, fmt='traj', **kwargs):
        frames = ase.io.read(file_name, format=fmt, index=slice(begin, end, step))
        for atoms in frames:
            symbols = atoms.get_chemical_symbols()
            atom_names = list(set(symbols))
            atom_numbs = [symbols.count(symbol) for symbol in atom_names]
            atom_types = np.array([atom_names.index(symbol) for symbol in symbols]).astype(int)
    
            cells = atoms.cell[:]
            coords = atoms.get_positions()
            try:
                energies = atoms.get_potential_energy(force_consistent=True)
            except PropertyNotImplementedError:
                energies = atoms.get_potential_energy()
            forces = atoms.get_forces()
            info_dict = {
                'atom_names': atom_names,
                'atom_numbs': atom_numbs,
                'atom_types': atom_types,
                'cells': np.array([cells]).astype('float32'),
                'coords': np.array([coords]).astype('float32'),
                'energies': np.array([energies]).astype('float32'),
                'forces': np.array([forces]).astype('float32'),
                'orig': [0,0,0],
            }
            try:
                stress = atoms.get_stress(False)
                virials = np.array([-atoms.get_volume() * stress]).astype('float32')
                info_dict['virials'] = virials
            except PropertyNotImplementedError:
                pass
            yield info_dict

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
