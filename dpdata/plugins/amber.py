import os
import subprocess as sp
import tempfile

import dpdata.amber.md
import dpdata.amber.sqm
from dpdata.driver import Driver, Minimizer
from dpdata.format import Format


@Format.register("amber/md")
class AmberMDFormat(Format):
    def from_system(
        self,
        file_name=None,
        parm7_file=None,
        nc_file=None,
        use_element_symbols=None,
        **kwargs,
    ):
        # assume the prefix is the same if the spefic name is not given
        if parm7_file is None:
            parm7_file = file_name + ".parm7"
        if nc_file is None:
            nc_file = file_name + ".nc"
        return dpdata.amber.md.read_amber_traj(
            parm7_file=parm7_file,
            nc_file=nc_file,
            use_element_symbols=use_element_symbols,
            labeled=False,
        )

    def from_labeled_system(
        self,
        file_name=None,
        parm7_file=None,
        nc_file=None,
        mdfrc_file=None,
        mden_file=None,
        mdout_file=None,
        use_element_symbols=None,
        **kwargs,
    ):
        # assume the prefix is the same if the spefic name is not given
        if parm7_file is None:
            parm7_file = file_name + ".parm7"
        if nc_file is None:
            nc_file = file_name + ".nc"
        if mdfrc_file is None:
            mdfrc_file = file_name + ".mdfrc"
        if mden_file is None:
            mden_file = file_name + ".mden"
        if mdout_file is None:
            mdout_file = file_name + ".mdout"
        return dpdata.amber.md.read_amber_traj(
            parm7_file, nc_file, mdfrc_file, mden_file, mdout_file, use_element_symbols
        )


@Format.register("sqm/out")
class SQMOutFormat(Format):
    def from_system(self, fname, **kwargs):
        """Read from ambertools sqm.out."""
        return dpdata.amber.sqm.parse_sqm_out(fname)

    def from_labeled_system(self, fname, **kwargs):
        """Read from ambertools sqm.out."""
        data = dpdata.amber.sqm.parse_sqm_out(fname)
        assert "forces" in list(data.keys()), f"No forces in {fname}"
        return data


@Format.register("sqm/in")
class SQMINFormat(Format):
    def to_system(self, data, fname=None, frame_idx=0, **kwargs):
        """Generate input files for semi-emperical calculation in sqm software.

        Parameters
        ----------
        data : dict
            system data
        fname : str
            output file name
        frame_idx : int, default=0
            index of frame to write
        **kwargs : dict
            other parameters

        Other Parameters
        ----------------
        **kwargs : dict
            valid parameters are:
                qm_theory : str, default=dftb3
                    level of theory. Options includes AM1, RM1, MNDO, PM3-PDDG, MNDO-PDDG,
                    PM3-CARB1, MNDO/d, AM1/d, PM6, DFTB2, DFTB3
                charge : int, default=0
                    total charge in electron units
                maxcyc : int, default=0
                    maximum number of minimization cycles to allow. 0 represents a
                    single-point calculation
                mult : int, default=1
                    multiplicity. Only 1 is allowed.
        """
        return dpdata.amber.sqm.make_sqm_in(data, fname, frame_idx, **kwargs)


@Driver.register("sqm")
class SQMDriver(Driver):
    """AMBER sqm program driver.

    Parameters
    ----------
    sqm_exec : str, default=sqm
        path to sqm program
    **kwargs : dict
        other arguments to make input files. See :class:`SQMINFormat`

    Examples
    --------
    Use DFTB3 method to calculate potential energy:

    >>> labeled_system = system.predict(theory="DFTB3", driver="sqm")
    >>> labeled_system['energies'][0]
    -15.41111246
    """

    def __init__(self, sqm_exec: str = "sqm", **kwargs: dict) -> None:
        self.sqm_exec = sqm_exec
        self.kwargs = kwargs

    def label(self, data: dict) -> dict:
        ori_system = dpdata.System(data=data)
        labeled_system = dpdata.LabeledSystem()
        with tempfile.TemporaryDirectory() as d:
            for ii, ss in enumerate(ori_system):
                inp_fn = os.path.join(d, "%d.in" % ii)
                out_fn = os.path.join(d, "%d.out" % ii)
                ss.to("sqm/in", inp_fn, **self.kwargs)
                try:
                    sp.check_output(
                        [*self.sqm_exec.split(), "-O", "-i", inp_fn, "-o", out_fn]
                    )
                except sp.CalledProcessError as e:
                    with open(out_fn) as f:
                        raise RuntimeError(
                            "Run sqm failed! Output:\n" + f.read()
                        ) from e
                labeled_system.append(dpdata.LabeledSystem(out_fn, fmt="sqm/out"))
        return labeled_system.data


@Minimizer.register("sqm")
class SQMMinimizer(Minimizer):
    """SQM minimizer.

    Parameters
    ----------
    maxcyc : int, default=1000
        maximun cycle to minimize
    """

    def __init__(self, maxcyc=1000, *args, **kwargs) -> None:
        assert maxcyc > 0, "maxcyc should be more than 0 to minimize"
        self.driver = SQMDriver(maxcyc=maxcyc, **kwargs)

    def minimize(self, data: dict) -> dict:
        # sqm has minimize feature
        return self.driver.label(data)
