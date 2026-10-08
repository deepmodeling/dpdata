from __future__ import annotations

import numpy as np

from dpdata.format import Format
from dpdata.formats.dftbplus.output import parse_dftb_plus
from dpdata.unit import EnergyConversion, ForceConversion, PressureConversion

energy_convert = EnergyConversion("hartree", "eV").value()
force_convert = ForceConversion("hartree/bohr", "eV/angstrom").value()
pressure_convert = PressureConversion("hartree/bohr^3", "eV/angstrom^3").value()


@Format.register("dftbplus")
class DFTBplusFormat(Format):
    """DFTB+ input/output pair for one labeled configuration.

    Pass a tuple containing the DFTB+ input geometry file and output result
    file. The reader combines symbols and coordinates from the input with the
    energy and forces from the output. A periodic GenFormat geometry (``S`` or
    ``F``) keeps its cell, and the ``Total stress tensor`` in the output, when
    present, is converted to the virial. See the `DFTB+ documentation
    <https://dftbplus.org/documentation>`_ for the underlying file formats.
    """

    def from_labeled_system(self, file_paths, **kwargs):
        """Reads system information from the given DFTB+ file paths.

        Parameters
        ----------
        file_paths : tuple
            A tuple containing the input and output file paths.
            - Input file (file_in): Contains information about symbols, coord and cell.
            - Output file (file_out): Contains information about energy, force and stress.
        **kwargs : dict
            other parameters

        """
        file_in, file_out = file_paths
        frame = parse_dftb_plus(file_in, file_out)
        symbols = frame["symbols"]
        coord = frame["coords"]
        last_occurrence = {v: i for i, v in enumerate(symbols)}
        atom_names = np.array(sorted(np.unique(symbols), key=last_occurrence.get))
        atom_types = np.array([np.where(atom_names == s)[0][0] for s in symbols])
        atom_numbs = np.array([np.sum(atom_types == i) for i in range(len(atom_names))])
        natoms = coord.shape[0]

        data = {
            "atom_types": atom_types,
            "atom_names": list(atom_names),
            "atom_numbs": list(atom_numbs),
            "coords": coord.reshape((1, natoms, 3)),
            "energies": np.array([frame["energy"] * energy_convert]),
            "forces": (frame["forces"] * force_convert).reshape((1, natoms, 3)),
            "orig": np.zeros(3),
        }
        cell = frame["cell"]
        if cell is None:
            data["cells"] = np.zeros((1, 3, 3))
            data["nopbc"] = True
        else:
            data["cells"] = cell.reshape((1, 3, 3))
            if frame["stress"] is not None:
                # DFTB+ reports the stress with the opposite sign to ASE's
                # (see ase.calculators.dftb), so the virial is +volume * stress.
                volume = abs(np.linalg.det(cell))
                virial = frame["stress"] * volume * pressure_convert
                data["virials"] = virial.reshape((1, 3, 3))
        return data
