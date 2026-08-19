from __future__ import annotations

import dpdata.formats.openmx.omx
import dpdata.md.pbc
from dpdata.format import Format


@Format.register("openmx/md")
class OPENMXFormat(Format):
    """Output pair from `OpenMX <https://www.openmx-square.org/>`_.

    OpenMX (Open source package for Material eXplorer) is a nano-scale material simulation package based on DFT, norm-conserving pseudopotentials, and pseudo-atomic localized basis functions.

    Note that two output files, System.Name.dat and System.Name.md, are required.

    Use the ``openmx/md`` alias and pass the shared ``System.Name`` prefix;
    dpdata appends ``.dat`` and ``.md`` automatically.
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

        data, _ = dpdata.formats.openmx.omx.to_system_data(fname, mdname)
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

        data, cs = dpdata.formats.openmx.omx.to_system_data(fname, mdname)
        data["coords"] = dpdata.md.pbc.apply_pbc(
            data["coords"],
            data["cells"],
        )
        data["energies"], data["forces"] = dpdata.formats.openmx.omx.to_system_label(
            fname, mdname
        )
        return data
