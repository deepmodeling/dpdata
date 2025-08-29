from __future__ import annotations

import pickle

from dpdata.format import Format
from dpdata.unit import EnergyConversion


@Format.register("fennol")
class FeNNolFormat(Format):
    """The FeNNol format plugin for dpdata.
    
    FeNNol (https://github.com/thomasple/FeNNol/) uses a pickle format
    for training machine learning models. This plugin supports exporting
    dpdata LabeledSystem to FeNNol format.
    
    The format consists of a dictionary with 'training' and 'validation' keys,
    where each contains a list of structures with:
    - 'species': atomic species/elements
    - 'coordinates': atomic positions in Angstroms
    - 'formation_energy': energy in kcal/mol
    - 'shifted_energy': energy in kcal/mol (same as formation_energy in this implementation)
    - 'forces': atomic forces in kcal/mol/Angstrom
    
    Examples
    --------
    Export a LabeledSystem to FeNNol format:
    
    >>> import dpdata
    >>> ls = dpdata.LabeledSystem("OUTCAR", fmt="vasp/outcar")
    >>> ls.to("fennol", "data.pkl")
    """

    def to_labeled_system(self, data, file_name, train_size=0.8, **kwargs):
        """Convert dpdata LabeledSystem to FeNNol format.
        
        Parameters
        ----------
        data : dict
            LabeledSystem data
        file_name : str
            Output pickle file name
        train_size : float, optional
            Fraction of data to use for training (default: 0.8)
        **kwargs : dict
            Other parameters
        """
        # Unit conversions
        energy_conv = EnergyConversion("eV", "kcal_mol").value()
        force_conv = EnergyConversion("eV", "kcal_mol").value()  # eV/Angstrom to kcal/mol/Angstrom
        
        # Extract data
        atom_names = data["atom_names"]
        atom_types = data["atom_types"]
        coords = data["coords"]  # shape: (nframes, natoms, 3)
        energies = data["energies"]  # shape: (nframes,)
        forces = data["forces"]  # shape: (nframes, natoms, 3)
        
        nframes = coords.shape[0]
        natoms = coords.shape[1]
        
        # Create species array from atom_types and atom_names
        species = [atom_names[atom_types[i]] for i in range(natoms)]
        
        # Prepare data structures
        structures = []
        
        for i in range(nframes):
            structure = {
                "species": species,
                "coordinates": coords[i].copy(),  # Already in Angstroms
                "formation_energy": energies[i] * energy_conv,  # Convert eV to kcal/mol
                "shifted_energy": energies[i] * energy_conv,  # Same as formation_energy
                "forces": forces[i] * force_conv,  # Convert eV/Angstrom to kcal/mol/Angstrom
            }
            structures.append(structure)
        
        # Split into training and validation sets
        n_train = int(nframes * train_size)
        training_data = structures[:n_train]
        validation_data = structures[n_train:]
        
        # Create FeNNol format dictionary
        fennol_data = {
            "training": training_data,
            "validation": validation_data,
            "description": f"Generated from dpdata with {nframes} frames, {n_train} training, {nframes - n_train} validation"
        }
        
        # Save to pickle file
        with open(file_name, 'wb') as f:
            pickle.dump(fennol_data, f)