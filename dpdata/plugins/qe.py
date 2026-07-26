from __future__ import annotations

import dpdata.formats.qe.scf
import dpdata.formats.qe.traj
import dpdata.md.pbc
from dpdata.format import Format


@Format.register("qe/cp/traj")
class QECPTrajFormat(Format):
    """Quantum ESPRESSO CP trajectory files sharing a common prefix.

    `Quantum ESPRESSO <https://www.quantum-espresso.org/>`_ is an integrated
    suite of open-source codes for electronic-structure calculations based on
    DFT, plane waves, and pseudopotentials.

    Given ``file_name='run'``, dpdata reads ``run.in`` together with the CP
    trajectory files rooted at ``run``. Loading as a labeled system also reads
    the matching energy and force records.
    """

    @Format.post("rot_lower_triangular")
    def from_system(self, file_name, begin=0, step=1, **kwargs):
        """Load coordinates and cells from a Quantum ESPRESSO CP trajectory.

        Parameters
        ----------
        file_name : str
            Common prefix of the CP input and trajectory files.
        begin : int, default=0
            Index of the first frame to load.
        step : int, default=1
            Load every ``step``-th frame.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Unlabeled trajectory data.
        """
        data, _ = dpdata.formats.qe.traj.to_system_data(
            file_name + ".in", file_name, begin=begin, step=step
        )
        data["coords"] = dpdata.md.pbc.apply_pbc(
            data["coords"],
            data["cells"],
        )
        return data

    @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, begin=0, step=1, **kwargs):
        """Load a labeled Quantum ESPRESSO CP trajectory.

        Parameters
        ----------
        file_name : str
            Common prefix of the CP input, trajectory, energy, and force files.
        begin : int, default=0
            Index of the first frame to load.
        step : int, default=1
            Load every ``step``-th frame.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Coordinates, cells, energies, and forces for the selected frames.
        """
        data, cs = dpdata.formats.qe.traj.to_system_data(
            file_name + ".in", file_name, begin=begin, step=step
        )
        data["coords"] = dpdata.md.pbc.apply_pbc(
            data["coords"],
            data["cells"],
        )
        data["energies"], data["forces"], es = dpdata.formats.qe.traj.to_system_label(
            file_name + ".in", file_name, begin=begin, step=step
        )
        assert cs == es, "the step key between files are not consistent"
        return data


@Format.register("qe/pw/scf")
class QECPPWSCFFormat(Format):
    """Quantum ESPRESSO PWscf self-consistent-field output.

    `Quantum ESPRESSO <https://www.quantum-espresso.org/>`_ is an integrated
    suite of open-source codes for DFT calculations.

    The reader extracts the final cell, coordinates, total energy, forces, and
    optional stress/virial from a ``pw.x`` text output.
    """

    @Format.post("rot_lower_triangular")
    def from_labeled_system(self, file_name, **kwargs):
        """Load a labeled Quantum ESPRESSO PWscf calculation.

        Parameters
        ----------
        file_name : str or list[str]
            Quantum ESPRESSO ``pw.x`` output file. The matching input file is
            inferred by replacing ``out`` with ``in`` in the base name; pass
            ``[input_file, output_file]`` to give both paths explicitly.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled system data for the calculation.
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
        ) = dpdata.formats.qe.scf.get_frame(file_name)
        if tmp_virial is not None:
            data["virials"] = tmp_virial
        return data
