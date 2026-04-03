from __future__ import annotations

import pickle
from typing import Any

from dpdata.format import Format
from dpdata.unit import EnergyConversion


class _MultiSystemCollector:
    """Helper class to collect data from multiple systems."""

    def __init__(self, filename: str, **kwargs):
        self.filename = filename
        self.kwargs = kwargs
        self.systems_data = []
        self.expected_count = 0
        self.processed_count = 0
    
    def add_system_data(self, system_name: str, data: dict):
        """Add data from a single system."""
        self.systems_data.append((system_name, data))
        self.processed_count += 1
        
        # Write combined data when all systems are processed
        if self.processed_count == self.expected_count:
            self._write_combined_data()
    
    def _write_combined_data(self):
        """Combine data from all systems and write to FeNNol format."""
        if not self.systems_data:
            return
            
        # Unit conversions
        energy_conv = EnergyConversion("eV", "kcal_mol").value()
        force_conv = EnergyConversion("eV", "kcal_mol").value()
        
        all_structures = []
        total_frames = 0
        
        # Process each system
        for system_name, data in self.systems_data:
            atom_names = data["atom_names"]
            atom_types = data["atom_types"] 
            coords = data["coords"]
            energies = data["energies"]
            forces = data["forces"]
            
            nframes = coords.shape[0]
            natoms = coords.shape[1]
            total_frames += nframes
            
            # Create species array from atom_types and atom_names
            species = [atom_names[atom_types[i]] for i in range(natoms)]
            
            # Process each frame
            for i in range(nframes):
                structure = {
                    "species": species,
                    "coordinates": coords[i].copy(),
                    "formation_energy": energies[i] * energy_conv,
                    "shifted_energy": energies[i] * energy_conv,
                    "forces": forces[i] * force_conv,
                    "system_name": system_name,  # Track which system this came from
                }
                all_structures.append(structure)
        
        # Split into training and validation sets
        train_size = self.kwargs.get('train_size', 0.8)
        n_train = int(total_frames * train_size)
        training_data = all_structures[:n_train]
        validation_data = all_structures[n_train:]
        
        # Create FeNNol format dictionary
        fennol_data = {
            "training": training_data,
            "validation": validation_data,
            "description": f"Generated from dpdata MultiSystems with {len(self.systems_data)} systems, "
                          f"{total_frames} frames, {n_train} training, {total_frames - n_train} validation"
        }
        
        # Save to pickle file
        with open(self.filename, 'wb') as f:
            pickle.dump(fennol_data, f)


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
    
    Export multiple systems to a single FeNNol file:
    
    >>> ms = dpdata.MultiSystems(ls1, ls2)
    >>> ms.to("fennol", "combined_data.pkl")
    """
    
    def __init__(self):
        super().__init__()
        self._multi_collector = None

    def to_multi_systems(self, formulas: list[str], directory: str, **kwargs: Any):
        """Generate collectors for writing multiple systems to the same FeNNol file.
        
        Parameters
        ----------
        formulas : list[str]
            formulas/names of systems
        directory : str
            FeNNol pickle file name
        **kwargs : dict
            other parameters (e.g., train_size)
            
        Yields
        ------
        _MultiSystemCollector
            collector object that systems will write their data to
        """
        # Create shared collector for all systems
        self._multi_collector = _MultiSystemCollector(directory, **kwargs)
        self._multi_collector.expected_count = len(formulas)
        
        # Yield the same collector for each system
        for formula in formulas:
            yield self._multi_collector

    def to_labeled_system(self, data, file_name, train_size=0.8, **kwargs):
        """Convert dpdata LabeledSystem to FeNNol format.
        
        Parameters
        ----------
        data : dict
            LabeledSystem data
        file_name : str or _MultiSystemCollector
            Output pickle file name or multi-system collector
        train_size : float, optional
            Fraction of data to use for training (default: 0.8)
        **kwargs : dict
            Other parameters
        """
        # Check if this is being called from MultiSystems
        if isinstance(file_name, _MultiSystemCollector):
            # Add data to the collector instead of writing directly
            system_name = kwargs.get('system_name', 'unnamed_system')
            file_name.add_system_data(system_name, data)
            return
        
        # Original single-system implementation
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