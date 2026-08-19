from __future__ import annotations

import dpdata.formats.siesta.aiMD_output
import dpdata.formats.siesta.output
from dpdata.format import Format


@Format.register("siesta/output")
class SiestaOutputFormat(Format):
    """SIESTA single-step text output.

    `SIESTA <https://siesta-project.org/>`_ (Spanish Initiative for
    Electronic Simulations with Thousands of Atoms) is an open-source DFT
    package based on linear-scaling methods and LCAO basis sets.

    The format can be loaded as an unlabeled structure or as a labeled system
    containing the energy, forces, and virial reported by SIESTA.
    """

    def from_system(self, file_name, **kwargs):
        """Load geometry from a SIESTA output file.

        Parameters
        ----------
        file_name : str or os.PathLike
            SIESTA text output.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Unlabeled system data.
        """
        data = {}
        (
            data["atom_names"],
            data["atom_numbs"],
            data["atom_types"],
            data["cells"],
            data["coords"],
            _e,
            _f,
            _v,
        ) = dpdata.formats.siesta.output.obtain_frame(file_name)
        return data

    def from_labeled_system(self, file_name, **kwargs):
        """Load geometry and labels from a SIESTA output file.

        Parameters
        ----------
        file_name : str or os.PathLike
            SIESTA text output.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled system data with energy, forces, and virial.
        """
        data = {}
        (
            data["atom_names"],
            data["atom_numbs"],
            data["atom_types"],
            data["cells"],
            data["coords"],
            data["energies"],
            data["forces"],
            data["virials"],
        ) = dpdata.formats.siesta.output.obtain_frame(file_name)
        return data


@Format.register("siesta/aimd_output")
@Format.register_from("from_siesta_aiMD_output")
class SiestaAIMDOutputFormat(Format):
    """SIESTA ab initio molecular-dynamics output.

    `SIESTA <https://siesta-project.org/>`_ (Spanish Initiative for
    Electronic Simulations with Thousands of Atoms) is an open-source DFT
    package based on LCAO basis sets.

    This reader handles the multi-frame layout emitted by SIESTA AIMD runs and
    can return either the trajectory geometry alone or all available labels.
    """

    def from_system(self, file_name, **kwargs):
        """Load geometry from a SIESTA AIMD output.

        Parameters
        ----------
        file_name : str or os.PathLike
            SIESTA AIMD output file.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Unlabeled trajectory data.
        """
        data = {}
        (
            data["atom_names"],
            data["atom_numbs"],
            data["atom_types"],
            data["cells"],
            data["coords"],
            _e,
            _f,
            _v,
        ) = dpdata.formats.siesta.aiMD_output.get_aiMD_frame(file_name)
        return data

    def from_labeled_system(self, file_name, **kwargs):
        """Load geometry and labels from a SIESTA AIMD output.

        Parameters
        ----------
        file_name : str or os.PathLike
            SIESTA AIMD output file.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled trajectory data with energies, forces, and virials.
        """
        data = {}
        (
            data["atom_names"],
            data["atom_numbs"],
            data["atom_types"],
            data["cells"],
            data["coords"],
            data["energies"],
            data["forces"],
            data["virials"],
        ) = dpdata.formats.siesta.aiMD_output.get_aiMD_frame(file_name)
        return data
