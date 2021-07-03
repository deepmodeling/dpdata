from dpdata.format import Format


@Format.register("pymatgen/structure")
@Format.register_to("to_pymatgen_structure")
class PyMatgenStructureFormat(Format):
    def to_system(self, data, **kwargs):
        """convert System to Pymatgen Structure obj
        """
        structures = []
        try:
            from pymatgen import Structure
        except ModuleNotFoundError as e:
            raise ImportError('No module pymatgen.Structure') from e

        species = []
        for name, numb in zip(data['atom_names'], data['atom_numbs']):
            species.extend([name]*numb)
        for ii in range(data['coords'].shape[0]):
            structure = Structure(
                data['cells'][ii], species, data['coords'][ii], coords_are_cartesian=True)
            structures.append(structure)
        return structures


@Format.register("pymatgen/computedstructureentry")
@Format.register_to("to_pymatgen_ComputedStructureEntry")
class PyMatgenCSEFormat(Format):
    def to_labeled_system(self, data, *args, **kwargs):
        """convert System to Pymagen ComputedStructureEntry obj
        """
        try:
            from pymatgen.entries.computed_entries import ComputedStructureEntry
        except ModuleNotFoundError as e:
            raise ImportError(
                'No module ComputedStructureEntry in pymatgen.entries.computed_entries') from e

        entries = []

        for ii, structure in enumerate(PyMatgenStructureFormat().to_system(data)):
            energy = data['energies'][ii]
            csedata = {'forces': data['forces'][ii],
                    'virials': data['virials'][ii]}

            entry = ComputedStructureEntry(structure, energy, data=csedata)
            entries.append(entry)
        return entries
