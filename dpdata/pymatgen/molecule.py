from pymatgen.core import Molecule

def to_system_data(file_name) :
    mol = Molecule.from_file(file_name)
    elem_mol = list(str(site.species.elements[0]) for site in mol.sites)
    elem_counter = Counter(elem_mol)
    atom_names = elem_counter.keys()
    atom_numbs = elem_counter.values()

    atom_types = []
    for idx,ii in enumerate(atom_numbs) :
        for jj in range(ii) :
            atom_types.append(idx)
    
    system = {}
    system['atom_names'] = atom_names
    system['atom_numbs'] = atom_numbs
    system['atom_types'] = np.array(atom_types, dtype = int)
    system['orig'] = np.array([0, 0, 0])
    system['cells'] = []

    mol.get_centered_molecule()
    system['coords'] = [np.copy(mol.cart_coords)]

    return system
