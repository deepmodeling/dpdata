from __future__ import annotations

import os
import subprocess as sp
import tempfile
from typing import TYPE_CHECKING

import numpy as np

import dpdata.formats.gaussian.fchk
import dpdata.formats.gaussian.gjf
import dpdata.formats.gaussian.log
from dpdata.data_type import Axis, DataType
from dpdata.driver import Driver
from dpdata.format import Format
from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType


def register_hessian_data(data):
    if "hessian" in data:
        dt = DataType(
            "hessian",
            np.ndarray,
            (Axis.NFRAMES, Axis.NATOMS, 3, Axis.NATOMS, 3),
            required=False,
            deepmd_name="hessian",
        )
        dpdata.LabeledSystem.register_data_type(dt)


@Format.register("gaussian/log")
class GaussianLogFormat(Format):
    """Gaussian text output containing energies, coordinates, and forces.

    `Gaussian <https://gaussian.com/>`_ is a general-purpose electronic
    structure package for molecules.

    Standard single-point or optimization output is read by default. Set
    ``md=True`` when the file contains a Gaussian molecular-dynamics run.
    """

    def from_labeled_system(self, file_name: FileType, md=False, **kwargs):
        """Load labeled frames from a Gaussian log file.

        Parameters
        ----------
        file_name : str or os.PathLike or file-like object
            Gaussian output file.
        md : bool, default=False
            Parse multiple molecular-dynamics frames instead of the standard
            calculation layout.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled molecular data, or empty label arrays when parsing fails.
        """
        try:
            return dpdata.formats.gaussian.log.to_system_data(file_name, md=md)
        except AssertionError:
            return {"energies": [], "forces": [], "nopbc": True}


@Format.register("gaussian/fchk")
class GaussianFChkFormat(Format):
    """Gaussian formatted checkpoint (``.fchk``) file.

    `Gaussian <https://gaussian.com/>`_ is a general-purpose electronic
    structure package. Formatted checkpoint files provide molecular geometry
    and energy and may also contain gradients and Cartesian force constants
    (the Hessian).
    """

    def from_labeled_system(
        self, file_name: FileType, has_forces=True, has_hessian=True, **kwargs
    ):
        """Load a Gaussian formatted checkpoint file.

        Parameters
        ----------
        file_name : str or os.PathLike or file-like object
            Gaussian ``.fchk`` file.
        has_forces : bool, default=True
            Expect and parse Cartesian gradients as forces.
        has_hessian : bool, default=True
            Expect and parse Cartesian force constants as a Hessian.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled molecular data and, when requested, Hessian data.
        """
        try:
            data = dpdata.formats.gaussian.fchk.to_system_data(
                file_name, has_forces=has_forces, has_hessian=has_hessian
            )
            register_hessian_data(data)
            return data
        except AssertionError:
            return {"energies": [], "forces": [], "hessian": [], "nopbc": True}


@Format.register("gaussian/md")
class GaussianMDFormat(Format):
    """Gaussian molecular-dynamics text output.

    `Gaussian <https://gaussian.com/>`_ is a general-purpose electronic
    structure package. This alias uses the Gaussian log reader with
    multi-frame MD parsing enabled.
    """

    def from_labeled_system(self, file_name: FileType, **kwargs):
        """Load a Gaussian molecular-dynamics trajectory.

        Parameters
        ----------
        file_name : str or os.PathLike or file-like object
            Gaussian MD log file.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled molecular-dynamics frames.
        """
        return GaussianLogFormat().from_labeled_system(file_name, md=True)


@Format.register("gaussian/gjf")
class GaussiaGJFFormat(Format):
    """Gaussian input (``.gjf``/``.com``) file.

    `Gaussian <https://gaussian.com/>`_ is a general-purpose electronic
    structure package. The reader extracts the molecular geometry from an
    input deck. The writer creates a Gaussian job for the supplied frames
    using keyword arguments documented by
    :func:`dpdata.formats.gaussian.gjf.make_gaussian_input`.
    """

    def from_system(self, file_name: FileType, **kwargs):
        """Read Gaussian input file.

        Parameters
        ----------
        file_name : str
            file name
        **kwargs : dict
            keyword arguments
        """
        with open_file(file_name) as fp:
            text = fp.read()
        return dpdata.formats.gaussian.gjf.read_gaussian_input(text)

    def to_system(self, data: dict, file_name: FileType, **kwargs):
        """Generate Gaussian input file.

        Parameters
        ----------
        data : dict
            system data
        file_name : str
            file name
        **kwargs : dict
            Other parameters to make input files. See :meth:`dpdata.formats.gaussian.gjf.make_gaussian_input`
        """
        text = dpdata.formats.gaussian.gjf.make_gaussian_input(data, **kwargs)
        with open_file(file_name, "w") as fp:
            fp.write(text)


@Driver.register("gaussian")
class GaussianDriver(Driver):
    """Gaussian driver.

    Note that "force" keyword must be added. If the number of atoms is large,
    "Geom=PrintInputOrient" should be added.

    Parameters
    ----------
    gaussian_exec : str, default=g16
        path to gaussian program
    **kwargs : dict
        other arguments to make input files. See :meth:`dpdata.formats.gaussian.gjf.make_gaussian_input`

    Examples
    --------
    Use B3LYP method to calculate potential energy of a methane molecule:

    >>> labeled_system = system.predict(keywords="force b3lyp/6-31g**", driver="gaussian")
    >>> labeled_system['energies'][0]
    -1102.714590995794
    """

    def __init__(self, gaussian_exec: str = "g16", **kwargs) -> None:
        self.gaussian_exec = gaussian_exec
        self.kwargs = kwargs

    def label(self, data: dict) -> dict:
        """Label a system data. Returns new data with energy, forces, and virials.

        Parameters
        ----------
        data : dict
            data with coordinates and atom types

        Returns
        -------
        dict
            labeled data with energies and forces
        """
        ori_system = dpdata.System(data=data)
        labeled_system = dpdata.LabeledSystem()
        with tempfile.TemporaryDirectory() as d:
            for ii, ss in enumerate(ori_system):
                inp_fn = os.path.join(d, "%d.gjf" % ii)  # noqa: UP031
                out_fn = os.path.join(d, "%d.log" % ii)  # noqa: UP031
                ss.to("gaussian/gjf", inp_fn, **self.kwargs)
                try:
                    sp.check_output([*self.gaussian_exec.split(), inp_fn])
                except sp.CalledProcessError as e:
                    with open_file(out_fn) as f:
                        out = f.read()
                    raise RuntimeError("Run gaussian failed! Output:\n" + out) from e
                labeled_system.append(dpdata.LabeledSystem(out_fn, fmt="gaussian/log"))
        return labeled_system.data
