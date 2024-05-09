import csv
import os
from collections import defaultdict
from inspect import Parameter, Signature, cleandoc, signature
from typing import Literal

from numpydoc.docscrape import Parameter as numpydoc_Parameter
from numpydoc.docscrape_sphinx import SphinxDocString

from dpdata.bond_order_system import BondOrderSystem

# ensure all plugins are loaded!
from dpdata.driver import Driver, Minimizer
from dpdata.format import Format
from dpdata.system import LabeledSystem, MultiSystems, System


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
    """Check whether a method is override.

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
    return ":class:`{} <{}>`".format(
        cls.__name__, ".".join([cls.__module__, cls.__name__])
    )


def check_supported(fmt: Format):
    methods = {}
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
            methods[mtd] = None
            if mtd == "to_system":
                methods["to_labeled_system"] = None
    if fmt.MultiMode != fmt.MultiModes.NotImplemented:
        methods["from_multi_systems"] = None
        methods["to_multi_systems"] = None
    return list(methods.keys())


method_links = {
    "from_system": ":ref:`System() <{}_{}>`",
    "to_system": ":ref:`System.to() <{}_{}>`",
    "from_labeled_system": ":ref:`LabeledSystem() <{}_{}>`",
    "to_labeled_system": ":ref:`LabeledSystem.to() <{}_{}>`",
    "from_bond_order_system": ":ref:`BondOrderSystem() <{}_{}>`",
    "to_bond_order_system": ":ref:`BondOrderSystem.to() <{}_{}>`",
    "from_multi_systems": ":ref:`MultiSystems.load_systems_from_file() <{}_{}>`",
    "to_multi_systems": ":ref:`MultiSystems.to() <{}_{}>`",
}

method_classes = {
    "from_system": "System",
    "to_system": "System",
    "from_labeled_system": "LabeledSystem",
    "to_labeled_system": "LabeledSystem",
    "from_bond_order_system": "BondOrderSystem",
    "to_bond_order_system": "BondOrderSystem",
    "from_multi_systems": "MultiSystems",
    "to_multi_systems": "MultiSystems",
}

method_cls_obj = {
    "from_system": System,
    "to_system": System,
    "from_labeled_system": LabeledSystem,
    "to_labeled_system": LabeledSystem,
    "from_bond_order_system": BondOrderSystem,
    "to_bond_order_system": BondOrderSystem,
    "from_multi_systems": MultiSystems,
    "to_multi_systems": MultiSystems,
}


def generate_sub_format_pages(formats: dict):
    """Generate sub format pages."""
    os.makedirs("formats", exist_ok=True)
    for format, alias in formats.items():
        # format: Format, alias: list[str]
        buff = []
        buff.append(f".. _{format.__name__}:")
        buff.append("")
        for aa in alias:
            buff.append(f"{aa} format")
            buff.append("=" * len(buff[-1]))
            buff.append("")
        buff.append(f"Class: {get_cls_link(format)}")
        buff.append("")

        docstring = format.__doc__
        if docstring is not None:
            docstring = cleandoc(docstring)
            rst = str(SphinxDocString(docstring))
            buff.append(rst)
            buff.append("")

        buff.append("Conversions")
        buff.append("-----------")
        methods = check_supported(format)
        for method in methods:
            buff.append("")
            buff.append(f".. _{format.__name__}_{method}:")
            buff.append("")
            if method.startswith("from_"):
                buff.append(f"Convert from this format to {method_classes[method]}")
                buff.append("`" * len(buff[-1]))
            elif method.startswith("to_"):
                buff.append(f"Convert from {method_classes[method]} to this format")
                buff.append("`" * len(buff[-1]))
            buff.append("")
            method_obj = getattr(format, method)
            if (
                method == "to_labeled_system"
                and method not in format.__dict__
                and "to_system" in format.__dict__
            ):
                method_obj = getattr(format, "to_system")
            docstring = method_obj.__doc__
            if docstring is not None:
                docstring = cleandoc(docstring)
            sig = signature(method_obj)
            parameters = dict(sig.parameters)
            return_annotation = sig.return_annotation
            # del self
            del parameters[list(parameters)[0]]
            # del data
            if method.startswith("to_"):
                del parameters[list(parameters)[0]]
            if "args" in parameters:
                del parameters["args"]
            if "kwargs" in parameters:
                del parameters["kwargs"]
            if method == "to_multi_systems" or method.startswith("from_"):
                sig = Signature(
                    list(parameters.values()), return_annotation=method_cls_obj[method]
                )
            else:
                sig = Signature(
                    list(parameters.values()), return_annotation=return_annotation
                )
            sig = str(sig)
            if method.startswith("from_"):
                if method != "from_multi_systems":
                    for aa in alias:
                        parameters["fmt"] = Parameter(
                            "fmt",
                            Parameter.POSITIONAL_OR_KEYWORD,
                            default=None,
                            annotation=Literal[aa],
                        )
                        sig_fmt = Signature(
                            list(parameters.values()),
                            return_annotation=method_cls_obj[method],
                        )
                        sig_fmt = str(sig_fmt)
                        buff.append(
                            f""".. py:function:: dpdata.{method_classes[method]}{sig_fmt}"""
                        )
                        buff.append("""   :noindex:""")
                for aa in alias:
                    buff.append(
                        """.. py:function:: dpdata.{}.from_{}{}""".format(
                            method_classes[method],
                            aa.replace("/", "_").replace(".", ""),
                            sig,
                        )
                    )
                    buff.append("""   :noindex:""")
                buff.append("")
                if docstring is None or method not in format.__dict__:
                    docstring = f"""   Convert this format to :class:`{method_classes[method]}`."""
                doc_obj = SphinxDocString(docstring)
                if len(doc_obj["Parameters"]) > 0:
                    doc_obj["Parameters"] = [
                        xx
                        for xx in doc_obj["Parameters"]
                        if xx.name not in ("*args", "**kwargs")
                    ]
                else:
                    if method == "from_multi_systems":
                        doc_obj["Parameters"] = [
                            numpydoc_Parameter(
                                "directory",
                                "str",
                                ["directory of systems"],
                            )
                        ]
                doc_obj["Yields"] = []
                doc_obj["Returns"] = [
                    numpydoc_Parameter("", method_classes[method], ["converted system"])
                ]
                rst = "   " + str(doc_obj)
                buff.append(rst)
                buff.append("")
            elif method.startswith("to_"):
                for aa in alias:
                    parameters = {
                        "fmt": Parameter(
                            "fmt",
                            Parameter.POSITIONAL_OR_KEYWORD,
                            annotation=Literal[aa],
                        ),
                        **parameters,
                    }
                    if method == "to_multi_systems":
                        sig_fmt = Signature(
                            list(parameters.values()),
                            return_annotation=method_cls_obj[method],
                        )
                    else:
                        sig_fmt = Signature(
                            list(parameters.values()),
                            return_annotation=return_annotation,
                        )
                    sig_fmt = str(sig_fmt)
                    buff.append(
                        f""".. py:function:: dpdata.{method_classes[method]}.to{sig_fmt}"""
                    )
                    buff.append("""   :noindex:""")
                for aa in alias:
                    buff.append(
                        """.. py:function:: dpdata.{}.to_{}{}""".format(
                            method_classes[method],
                            aa.replace("/", "_").replace(".", ""),
                            sig,
                        )
                    )
                    buff.append("""   :noindex:""")
                buff.append("")
                if docstring is None or (
                    method not in format.__dict__
                    and not (
                        method == "to_labeled_system"
                        and method not in format.__dict__
                        and "to_system" in format.__dict__
                    )
                ):
                    docstring = (
                        f"Convert :class:`{method_classes[method]}` to this format."
                    )
                doc_obj = SphinxDocString(docstring)
                if len(doc_obj["Parameters"]) > 0:
                    doc_obj["Parameters"] = [
                        xx
                        for xx in doc_obj["Parameters"][1:]
                        if xx.name not in ("*args", "**kwargs")
                    ]
                else:
                    if method == "to_multi_systems":
                        doc_obj["Parameters"] = [
                            numpydoc_Parameter(
                                "directory",
                                "str",
                                ["directory to save systems"],
                            )
                        ]
                if method == "to_multi_systems":
                    doc_obj["Yields"] = []
                    doc_obj["Returns"] = [
                        numpydoc_Parameter("", method_classes[method], ["this system"])
                    ]
                rst = "   " + str(doc_obj)
                buff.append(rst)
                buff.append("")
            buff.append("")
            buff.append("")

        with open(f"formats/{format.__name__}.rst", "w") as rstfile:
            rstfile.write("\n".join(buff))


if __name__ == "__main__":
    formats = get_formats()
    with open("formats.csv", "w", newline="") as csvfile:
        fieldnames = [
            "Format",
            "Alias",
            "Supported Conversions",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for kk, vv in formats.items():
            writer.writerow(
                {
                    "Format": f":ref:`{kk.__name__}`",
                    "Alias": "\n".join(f"``{vvv}``" for vvv in vv),
                    "Supported Conversions": "\n".join(
                        method_links[mtd].format(kk.__name__, mtd)
                        for mtd in check_supported(kk)
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
                    "Alias": "\n".join(f"``{vvv}``" for vvv in vv),
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
                    "Alias": "\n".join(f"``{vvv}``" for vvv in vv),
                }
            )
    generate_sub_format_pages(formats)
