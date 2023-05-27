import os
import subprocess as sp
import tempfile

import dpdata.gaussian.gjf
import dpdata.gaussian.log
from dpdata.driver import Driver
from dpdata.format import Format


@Format.register("gaussian/log")
class GaussianLogFormat(Format):
    def from_labeled_system(self, file_name, md=False, **kwargs):
        try:
            return dpdata.gaussian.log.to_system_data(file_name, md=md)
        except AssertionError:
            return {"energies": [], "forces": [], "nopbc": True}


@Format.register("gaussian/md")
class GaussianMDFormat(Format):
    def from_labeled_system(self, file_name, **kwargs):
        return GaussianLogFormat().from_labeled_system(file_name, md=True)


@Format.register("gaussian/gjf")
class GaussiaGJFFormat(Format):
    """Gaussian input file."""

    def from_system(self, file_name: str, **kwargs):
        """Read Gaussian input file.

        Parameters
        ----------
        file_name : str
            file name
        **kwargs : dict
            keyword arguments
        """
        with open(file_name) as fp:
            text = fp.read()
        return dpdata.gaussian.gjf.read_gaussian_input(text)

    def to_system(self, data: dict, file_name: str, **kwargs):
        """Generate Gaussian input file.

        Parameters
        ----------
        data : dict
            system data
        file_name : str
            file name
        **kwargs : dict
            Other parameters to make input files. See :meth:`dpdata.gaussian.gjf.make_gaussian_input`
        """
        text = dpdata.gaussian.gjf.make_gaussian_input(data, **kwargs)
        with open(file_name, "w") as fp:
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
        other arguments to make input files. See :meth:`dpdata.gaussian.gjf.make_gaussian_input`

    Examples
    --------
    Use B3LYP method to calculate potential energy of a methane molecule:

    >>> labeled_system = system.predict(keywords="force b3lyp/6-31g**", driver="gaussian")
    >>> labeled_system['energies'][0]
    -1102.714590995794
    """

    def __init__(self, gaussian_exec: str = "g16", **kwargs: dict) -> None:
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
                inp_fn = os.path.join(d, "%d.gjf" % ii)
                out_fn = os.path.join(d, "%d.log" % ii)
                ss.to("gaussian/gjf", inp_fn, **self.kwargs)
                try:
                    sp.check_output([*self.gaussian_exec.split(), inp_fn])
                except sp.CalledProcessError as e:
                    with open(out_fn) as f:
                        out = f.read()
                    raise RuntimeError("Run gaussian failed! Output:\n" + out) from e
                labeled_system.append(dpdata.LabeledSystem(out_fn, fmt="gaussian/log"))
        return labeled_system.data
