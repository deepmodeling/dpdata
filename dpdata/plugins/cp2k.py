from __future__ import annotations

import glob

import dpdata.formats.cp2k.output
from dpdata.format import Format
from dpdata.formats.cp2k.output import Cp2kSystems

string_warning = """
Hi, you got an error from dpdata,
please check if your cp2k files include full information,
otherwise its version is not supported by dpdata.
Try use dpdata plugin from cp2kdata package,
for details, please refer to
https://robinzyb.github.io/cp2kdata/
"""


@Format.register("cp2k/aimd_output")
class CP2KAIMDOutputFormat(Format):
    """CP2K AIMD calculation directory.

    `CP2K <https://www.cp2k.org/>`_ is a quantum chemistry and solid state
    physics software package that can perform atomistic simulations.

    The reader pairs the first ``*pos*.xyz`` trajectory with the first CP2K
    ``.log`` file in the directory and extracts coordinates, cells, energies,
    forces, and virials where available.
    """

    def from_labeled_system(self, file_name, restart=False, **kwargs):
        """Load a labeled CP2K AIMD trajectory.

        Parameters
        ----------
        file_name : str or os.PathLike
            Directory containing CP2K position and log files.
        restart : bool, default=False
            Whether the trajectory is from a restarted CP2K calculation.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        tuple[dict, ...]
            One or more labeled system-data dictionaries parsed from the run.
        """
        xyz_file = sorted(glob.glob(f"{file_name}/*pos*.xyz"))[0]
        log_file = sorted(glob.glob(f"{file_name}/*.log"))[0]
        try:
            return tuple(Cp2kSystems(log_file, xyz_file, restart))
        except (StopIteration, RuntimeError) as e:
            # StopIteration is raised when pattern match is failed
            raise PendingDeprecationWarning(string_warning) from e


@Format.register("cp2k/output")
class CP2KOutputFormat(Format):
    """Single CP2K output file containing coordinates and calculation labels.

    `CP2K <https://www.cp2k.org/>`_ is a quantum chemistry and solid state
    physics software package. This legacy reader targets standard CP2K text
    output. For newer or unsupported CP2K layouts, use the separately
    maintained ``cp2kdata`` plugin referenced by the warning raised on parse
    failure.
    """

    def from_labeled_system(self, file_name, restart=False, **kwargs):
        """Load frames from a CP2K text output.

        Parameters
        ----------
        file_name : str or os.PathLike
            CP2K output file.
        restart : bool, default=False
            Reserved for compatibility with the AIMD reader.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled system data.
        """
        try:
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
            ) = dpdata.formats.cp2k.output.get_frames(file_name)
            if tmp_virial is not None:
                data["virials"] = tmp_virial
            return data
        # TODO: in the future, we should add exact error type here
        # TODO: when pattern match is failed
        # TODO: For now just use RuntimeError as a placeholder.
        except RuntimeError as e:
            raise PendingDeprecationWarning(string_warning) from e
