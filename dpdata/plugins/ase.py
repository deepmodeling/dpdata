from dpdata.format import Format


@Format.register("ase/structure")
@Format.register_to("to_ase_structure")
class ASEStructureFormat(Format):
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
                v_pref = 1 * 1e4 / 1.602176621e6
                vol = structure.get_volume()
                results['stress'] = data["virials"][ii] / (v_pref * vol)

            structure.calc = SinglePointCalculator(structure, **results)
            structures.append(structure)

        return structures
