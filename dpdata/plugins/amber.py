from __future__ import annotations

import os
import subprocess as sp
import tempfile

import dpdata.formats.amber.md
import dpdata.formats.amber.sqm
from dpdata.driver import Driver, Minimizer
from dpdata.format import Format
from dpdata.utils import open_file


@Format.register("amber/md")
class AmberMDFormat(Format):
    """AMBER molecular-dynamics trajectory and label files.

    `AMBER <https://ambermd.org/>`_ is a suite of biomolecular simulation
    programs for molecular dynamics simulations and analysis.

    Coordinates and topology are read from ``.nc`` and ``.parm7`` files.
    Labeled loading additionally requires the ``.mdfrc`` force trajectory and
    takes energies from either the ``.mden`` or the ``.mdout`` file. The
    ``parmed`` optional dependency is required.
    """

    def from_system(
        self,
        file_name=None,
        parm7_file=None,
        nc_file=None,
        use_element_symbols=None,
        **kwargs,
    ):
        """Load an unlabeled AMBER trajectory.

        Parameters
        ----------
        file_name : str, optional
            Common prefix used to infer ``<prefix>.parm7`` and ``<prefix>.nc``.
        parm7_file : str, optional
            Explicit AMBER topology file. Overrides the inferred path.
        nc_file : str, optional
            Explicit NetCDF trajectory file. Overrides the inferred path.
        use_element_symbols : list[int] or str, optional
            Atoms whose element symbols are used instead of AMBER atom types,
            given either as a list of atom indexes or as an AMBER mask string
            selecting them.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Unlabeled trajectory data.
        """
        # assume the prefix is the same if the spefic name is not given
        if parm7_file is None:
            parm7_file = file_name + ".parm7"
        if nc_file is None:
            nc_file = file_name + ".nc"
        return dpdata.formats.amber.md.read_amber_traj(
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
        """Load a labeled AMBER trajectory.

        Parameters
        ----------
        file_name : str, optional
            Common prefix used to infer the AMBER file names.
        parm7_file : str, optional
            Explicit AMBER topology file.
        nc_file : str, optional
            Explicit NetCDF coordinate trajectory.
        mdfrc_file : str, optional
            Explicit force trajectory. Required for labeled loading; inferred
            from ``file_name`` when not given.
        mden_file : str, optional
            Explicit energy file. Used when present, otherwise ``mdout_file``
            supplies the energies.
        mdout_file : str, optional
            Explicit AMBER text output. Fallback energy source when
            ``mden_file`` is absent.
        use_element_symbols : list[int] or str, optional
            Atoms whose element symbols are used instead of AMBER atom types,
            given either as a list of atom indexes or as an AMBER mask string
            selecting them.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        dict
            Labeled trajectory data assembled from the supplied files.
        """
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
        return dpdata.formats.amber.md.read_amber_traj(
            parm7_file, nc_file, mdfrc_file, mden_file, mdout_file, use_element_symbols
        )


@Format.register("sqm/out")
class SQMOutFormat(Format):
    """AmberTools SQM output from a semiempirical calculation.

    `AmberTools <https://ambermd.org/AmberTools.php>`_ is a collection of
    complementary tools for AMBER simulations. SQM implements semiempirical
    quantum-mechanical methods.

    The same file can be loaded as an unlabeled system, or as a labeled system
    when the output contains gradients that can be converted to forces.
    """

    def from_system(self, fname, **kwargs):
        """Read coordinates from an AmberTools ``sqm.out`` file.

        Parameters
        ----------
        fname : str or os.PathLike
            SQM output file.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.
        """
        return dpdata.formats.amber.sqm.parse_sqm_out(fname)

    def from_labeled_system(self, fname, **kwargs):
        """Read coordinates, energy, and forces from ``sqm.out``.

        Parameters
        ----------
        fname : str or os.PathLike
            SQM output file containing gradients.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.
        """
        data = dpdata.formats.amber.sqm.parse_sqm_out(fname)
        assert "forces" in list(data.keys()), f"No forces in {fname}"
        return data


@Format.register("sqm/in")
class SQMINFormat(Format):
    """AmberTools SQM input for semiempirical calculations.

    `AmberTools <https://ambermd.org/AmberTools.php>`_ provides the SQM
    module for semiempirical QM calculations.

    This write-only format serializes one nonperiodic System frame with its
    charge, multiplicity, semiempirical method, and minimization-cycle limit.
    Setting ``maxcyc=0`` requests a single-point calculation; positive values
    request geometry minimization.
    """

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
        return dpdata.formats.amber.sqm.make_sqm_in(data, fname, frame_idx, **kwargs)


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

    def __init__(self, sqm_exec: str = "sqm", **kwargs) -> None:
        self.sqm_exec = sqm_exec
        self.kwargs = kwargs

    def label(self, data: dict) -> dict:
        ori_system = dpdata.System(data=data)
        labeled_system = dpdata.LabeledSystem()
        with tempfile.TemporaryDirectory() as d:
            for ii, ss in enumerate(ori_system):
                inp_fn = os.path.join(d, "%d.in" % ii)  # noqa: UP031
                out_fn = os.path.join(d, "%d.out" % ii)  # noqa: UP031
                ss.to("sqm/in", inp_fn, **self.kwargs)
                try:
                    sp.check_output(
                        [*self.sqm_exec.split(), "-O", "-i", inp_fn, "-o", out_fn]
                    )
                except sp.CalledProcessError as e:
                    with open_file(out_fn) as f:
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
