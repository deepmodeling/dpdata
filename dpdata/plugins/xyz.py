from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from dpdata.format import Format
from dpdata.periodic_table import Element
from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType
from dpdata.xyz.quip_gap_xyz import QuipGapxyzSystems
from dpdata.xyz.xyz import coord_to_xyz, xyz_to_coord


@Format.register("xyz")
class XYZFormat(Format):
    """XYZ foramt.

    Examples
    --------
    >>> s.to("xyz", "a.xyz")
    """

    def to_system(self, data, file_name: FileType, **kwargs):
        buff = []
        types = np.array(data["atom_names"])[data["atom_types"]]
        for cc in data["coords"]:
            buff.append(coord_to_xyz(cc, types))
        with open_file(file_name, "w") as fp:
            fp.write("\n".join(buff))

    def from_system(self, file_name: FileType, **kwargs):
        with open_file(file_name) as fp:
            coords, types = xyz_to_coord(fp.read())
        atom_names, atom_types, atom_numbs = np.unique(
            types, return_inverse=True, return_counts=True
        )
        return {
            "atom_names": list(atom_names),
            "atom_numbs": list(atom_numbs),
            "atom_types": atom_types,
            "coords": coords.reshape((1, *coords.shape)),
            "cells": np.eye(3).reshape((1, 3, 3)) * 100,
            "nopbc": True,
            "orig": np.zeros(3),
        }


@Format.register("quip/gap/xyz")
@Format.register("quip/gap/xyz_file")
class QuipGapXYZFormat(Format):
    # Class variable to track which files have been written to
    _written_files = set()
    
    def from_labeled_system(self, data, **kwargs):
        return data

    def from_multi_systems(self, file_name, **kwargs):
        # here directory is the file_name
        return QuipGapxyzSystems(file_name)

    def to_labeled_system(self, data, file_name: FileType, **kwargs):
        """Write LabeledSystem data to QUIP/GAP XYZ format file.
        
        Parameters
        ----------
        data : dict
            system data
        file_name : FileType
            output file name
        **kwargs : dict
            additional arguments
        """
        frames = []
        nframes = len(data["energies"])
        
        for frame_idx in range(nframes):
            frame_lines = self._format_single_frame(data, frame_idx)
            frames.append("\n".join(frame_lines))
        
        # Determine if we should append based on whether file was written before
        file_path = str(file_name)
        if file_path in self._written_files:
            mode = "a"
        else:
            mode = "w"
            self._written_files.add(file_path)
        
        with open_file(file_name, mode) as fp:
            if mode == "a":
                # Add newline separator if appending
                fp.write("\n")
            fp.write("\n".join(frames))

    def to_multi_systems(self, formulas, directory, **kwargs):
        """Return single filename for all systems in QUIP/GAP XYZ format.
        
        For QUIP/GAP XYZ format, all systems are written to a single file.
        We use a class variable to track which files have been written to.
        
        Parameters
        ----------
        formulas : list[str]
            list of system names/formulas  
        directory : str
            output filename 
        **kwargs : dict
            additional arguments
            
        Returns
        -------
        list[str]
            list with same filename for all systems
        """
        # Clear the written files set for this new operation
        file_path = str(directory)
        if file_path in self._written_files:
            self._written_files.remove(file_path)
        
        # Return the same filename for all systems
        return [directory] * len(formulas)

    def _format_single_frame(self, data, frame_idx):
        """Format a single frame of system data into QUIP/GAP XYZ format lines.
        
        Parameters
        ----------
        data : dict
            system data
        frame_idx : int
            frame index
            
        Returns
        -------
        list[str]
            lines for the frame
        """
        # Number of atoms
        natoms = len(data["atom_types"])
        
        # Build header line with metadata
        header_parts = []
        
        # Energy
        energy = data["energies"][frame_idx]
        header_parts.append(f"energy={energy:.12e}")
        
        # Virial (if present)
        if "virials" in data:
            virial = data["virials"][frame_idx]
            virial_str = "    ".join(f"{v:.12e}" for v in virial.flatten())
            header_parts.append(f'virial="{virial_str}"')
        
        # Lattice
        cell = data["cells"][frame_idx]
        lattice_str = "   ".join(f"{c:.12e}" for c in cell.flatten())
        header_parts.append(f'Lattice="{lattice_str}"')
        
        # Properties
        header_parts.append("Properties=species:S:1:pos:R:3:Z:I:1:force:R:3")
        
        header_line = "    ".join(header_parts)
        
        # Format atom lines
        atom_lines = []
        coords = data["coords"][frame_idx]
        forces = data["forces"][frame_idx]
        atom_names = np.array(data["atom_names"])
        atom_types = data["atom_types"]
        
        for i in range(natoms):
            atom_type_idx = atom_types[i]
            species = atom_names[atom_type_idx]
            x, y, z = coords[i]
            fx, fy, fz = forces[i]
            atomic_number = Element(species).Z
            
            atom_line = f"{species}    {x:.11e}   {y:.11e}   {z:.11e}   {atomic_number}    {fx:.11e}  {fy:.11e}   {fz:.11e}"
            atom_lines.append(atom_line)
        
        # Combine all lines for this frame
        frame_lines = [str(natoms), header_line] + atom_lines
        return frame_lines
