import csv
from typing import Any
from collections import defaultdict

# ensure all plugins are loaded!
import dpdata.plugins
from dpdata.format import Format
from dpdata.driver import Driver
from dpdata.driver import Minimizer
from dpdata.system import get_cls_name


def get_formats() -> dict:
    formats = defaultdict(list)
    for kk, ff in Format.get_formats().items():
        formats[ff].append(kk)
    return formats


def get_driver() -> dict:
    drivers = defaultdict(list)
    for kk, ff in Driver.get_drivers().items():
        drivers[ff].append(kk)
    return drivers


def get_minimizer() -> dict:
    minimizers = defaultdict(list)
    for kk, ff in Minimizer.get_minimizers().items():
        minimizers[ff].append(kk)
    return minimizers


def detect_overridden(cls: Format, method: str) -> bool:
    """Check whether a method is override

    Parameters
    ----------
    cls : Format
        a format
    method : str
        method name

    Returns
    -------
    bool
        whether a method is overridden
    """
    return method in cls.__dict__


def get_cls_link(cls: object) -> str:
    """Returns class link.

    Parameters
    ----------
    cls : object
        the class

    Returns
    -------
    str
        the link of a class
    """
    return ":class:`%s <%s>`" % (cls.__name__, ".".join([cls.__module__, cls.__name__]))


def check_supported(fmt: Format):
    methods = set()
    for mtd in [
        "from_system",
        "to_system",
        "from_labeled_system",
        "to_labeled_system",
        "from_bond_order_system",
        "to_bond_order_system",
        "from_multi_systems",
        "to_multi_systems",
    ]:
        if detect_overridden(fmt, mtd):
            methods.add(mtd)
            if mtd == "to_system":
                methods.add("to_labeled_system")
    if fmt.MultiMode != fmt.MultiModes.NotImplemented:
        methods.add("from_multi_systems")
        methods.add("to_multi_systems")
    return methods


method_links = {
    "from_system": ":func:`System() <dpdata.system.System>`",
    "to_system": ":func:`System.to() <dpdata.system.System.to>`",
    "from_labeled_system": ":func:`LabeledSystem() <dpdata.system.LabeledSystem>`",
    "to_labeled_system": ":func:`LabeledSystem.to() <dpdata.system.System.to>`",
    "from_bond_order_system": ":func:`BondOrderSystem() <dpdata.bond_order_system.BondOrderSystem>`",
    "to_bond_order_system": ":func:`BondOrderSystem.to() <dpdata.system.System.to>`",
    "from_multi_systems": ":func:`MultiSystems.load_systems_from_file() <dpdata.system.MultiSystems.load_systems_from_file>`",
    "to_multi_systems": ":func:`MultiSystems.to() <dpdata.system.MultiSystems.to>`",
}

if __name__ == "__main__":
    formats = get_formats()
    with open("formats.csv", "w", newline="") as csvfile:
        fieldnames = [
            "Class",
            "Alias",
            "Supported Functions",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for kk, vv in formats.items():
            writer.writerow(
                {
                    "Class": get_cls_link(kk),
                    "Alias": "\n".join(("``%s``" % vvv for vvv in vv)),
                    "Supported Functions": "\n".join(
                        method_links[mtd] for mtd in check_supported(kk)
                    ),
                }
            )

    drivers = get_driver()
    with open("drivers.csv", "w", newline="") as csvfile:
        fieldnames = [
            "Class",
            "Alias",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for kk, vv in drivers.items():
            writer.writerow(
                {
                    "Class": get_cls_link(kk),
                    "Alias": "\n".join(("``%s``" % vvv for vvv in vv)),
                }
            )

    minimizers = get_minimizer()
    with open("minimizers.csv", "w", newline="") as csvfile:
        fieldnames = [
            "Class",
            "Alias",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for kk, vv in minimizers.items():
            writer.writerow(
                {
                    "Class": get_cls_link(kk),
                    "Alias": "\n".join(("``%s``" % vvv for vvv in vv)),
                }
            )
