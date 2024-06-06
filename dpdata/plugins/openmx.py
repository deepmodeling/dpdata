from __future__ import annotations

import dpdata.md.pbc
import dpdata.openmx.omx
from dpdata.format import Format


@Format.register("openmx/md")
class OPENMXFormat(Format):
    """Format for the `OpenMX <https://www.openmx-square.org/>`.

    OpenMX (Open source package for Material eXplorer) is a nano-scale material simulation package based on DFT, norm-conserving pseudopotentials, and pseudo-atomic localized basis functions.

    Note that two output files, System.Name.dat and System.Name.md, are required.

    Use the `openmx/md` keyword argument to supply this format.
    """

    @Format.post("rot_lower_triangular")
    def from_system(self, file_name: str, **kwargs) -> dict:
        """Read from OpenMX output.

        Parameters
        ----------
        file_name : str
            file name, which is specified by a input file, i.e. System.Name.dat
        **kwargs : dict
            other parameters

        Returns
        -------
        dict
            data dict
        """
        fname = f"{file_name}.dat"
        mdname = f"{file_name}.md"

        data, _ = dpdata.openmx.omx.to_system_data(fname, mdname)
        data["coords"] = dpdata.md.pbc.apply_pbc(
            data["coords"],
            data["cells"],
        )
        return data

    @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name: str, **kwargs) -> dict:
        """Read from OpenMX output.

        Parameters
        ----------
        file_name : str
            file name, which is specified by a input file, i.e. System.Name.dat
        **kwargs : dict
            other parameters

        Returns
        -------
        dict
            data dict
        """
        fname = f"{file_name}.dat"
        mdname = f"{file_name}.md"

        data, cs = dpdata.openmx.omx.to_system_data(fname, mdname)
        data["coords"] = dpdata.md.pbc.apply_pbc(
            data["coords"],
            data["cells"],
        )
        data["energies"], data["forces"] = dpdata.openmx.omx.to_system_label(
            fname, mdname
        )
        return data
