"""Command line interface for dpdata."""
import argparse
from typing import Optional

from . import __version__
from .system import LabeledSystem, MultiSystems, System


def dpdata_parser() -> argparse.ArgumentParser:
    """Returns dpdata cli parser.

    Returns
    -------
    argparse.ArgumentParser
        dpdata cli parser
    """
    parser = argparse.ArgumentParser(
        description="dpdata: Manipulating multiple atomic simulation data formats",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("from_file", type=str, help="read data from a file")
    parser.add_argument("--to_file", "-O", type=str, help="dump data to a file")
    parser.add_argument(
        "--from_format", "-i", type=str, default="auto", help="the format of from_file"
    )
    parser.add_argument("--to_format", "-o", type=str, help="the format of to_file")
    parser.add_argument(
        "--no-labeled", "-n", action="store_true", help="labels aren't provided"
    )
    parser.add_argument(
        "--multi",
        "-m",
        action="store_true",
        help="the system contains multiple directories",
    )
    parser.add_argument("--type-map", "-t", type=str, nargs="+", help="type map")

    parser.add_argument(
        "--version", action="version", version="dpdata v%s" % __version__
    )
    return parser


def dpdata_cli():
    """Dpdata cli.

    Examples
    --------
    .. code-block:: bash

        $ dpdata -iposcar POSCAR -odeepmd/npy -O data -n
    """
    parser = dpdata_parser()
    parsed_args = parser.parse_args()
    convert(**vars(parsed_args))


def convert(
    *,
    from_file: str,
    from_format: str = "auto",
    to_file: Optional[str] = None,
    to_format: Optional[str] = None,
    no_labeled: bool = False,
    multi: bool = False,
    type_map: Optional[list] = None,
    **kwargs,
):
    """Convert files from one format to another one.

    Parameters
    ----------
    from_file : str
        read data from a file
    from_format : str
        the format of from_file
    to_file : str
        dump data to a file
    to_format : str
        the format of to_file
    no_labeled : bool
        labels aren't provided
    multi : bool
        the system contains multiple directories
    type_map : list
        type map
    **kwargs : dict
        Additional arguments for the format.
    """
    if multi:
        s = MultiSystems.from_file(
            from_file, fmt=from_format, type_map=type_map, labeled=not no_labeled
        )
    elif not no_labeled:
        s = LabeledSystem(from_file, fmt=from_format, type_map=type_map)
    else:
        s = System(from_file, fmt=from_format, type_map=type_map)
    if to_format is not None:
        out = s.to(to_format, to_file)
        if isinstance(out, str):
            print(out)
    else:
        print(s)
