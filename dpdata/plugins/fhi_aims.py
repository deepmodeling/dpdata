from __future__ import annotations

import dpdata.formats.fhi_aims.output
from dpdata.format import Format


@Format.register("fhi_aims/md")
@Format.register("fhi_aims/output")
class FhiMDFormat(Format):
    """FHI-aims molecular-dynamics or multi-step output.

    `FHI-aims <https://fhi-aims.org/>`_ is an all-electron electronic
    structure code based on numeric atom-centered orbitals.

    The text-output reader extracts geometries, energies, forces, and optional
    virials from converged FHI-aims calculation steps.
    """

    def from_labeled_system(
        self, file_name, md=True, begin=0, step=1, convergence_check=True, **kwargs
    ):
        """Load labeled frames from FHI-aims output.

        Parameters
        ----------
        file_name : str or os.PathLike
            FHI-aims output file.
        md : bool, default=True
            Parse the output as a multi-step molecular-dynamics calculation.
        begin : int, default=0
            Index of the first frame to load.
        step : int, default=1
            Load every ``step``-th frame.
        convergence_check : bool, default=True
            Exclude unconverged calculation steps when enabled.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled trajectory data.
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
            tmp_virial,
        ) = dpdata.formats.fhi_aims.output.get_frames(
            file_name,
            md=md,
            begin=begin,
            step=step,
            convergence_check=convergence_check,
        )
        if tmp_virial is not None:
            data["virials"] = tmp_virial
        return data


@Format.register("fhi_aims/scf")
class FhiSCFFormat(Format):
    """FHI-aims single-point self-consistent-field output.

    `FHI-aims <https://fhi-aims.org/>`_ is an all-electron electronic
    structure code. The text reader loads the calculation geometry, total
    energy, forces, and an optional virial into one
    :class:`dpdata.LabeledSystem` frame. Use ``fhi_aims/output`` for
    molecular-dynamics or other multi-step output.
    """

    def from_labeled_system(self, file_name, **kwargs):
        """Load the first labeled frame from FHI-aims SCF output.

        Parameters
        ----------
        file_name : str or os.PathLike
            FHI-aims output file.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled system data with energy, forces, and optional virial.
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
            tmp_virial,
        ) = dpdata.formats.fhi_aims.output.get_frames(
            file_name, md=False, begin=0, step=1
        )
        if tmp_virial is not None:
            data["virials"] = tmp_virial
        return data
