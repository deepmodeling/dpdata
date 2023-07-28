import csv
import os
from collections import defaultdict
from inspect import Parameter, Signature, signature
from typing import Literal

from numpydoc.docscrape_sphinx import SphinxDocString

# ensure all plugins are loaded!
from dpdata.driver import Driver, Minimizer
from dpdata.format import Format


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


def generate_sub_format_pages(formats: dict):
    """Generate sub format pages."""
    os.makedirs("formats", exist_ok=True)
    for format, alias in formats.items():
        # format: Format, alias: list[str]
        buff = []
        for aa in alias:
            buff.append("%s format" % aa)
            buff.append("=" * len(buff[-1]))
            buff.append("")
        buff.append("Class: %s" % get_cls_link(format))
        buff.append("")

        docstring = format.__doc__
        if docstring is not None:
            rst = str(SphinxDocString(docstring))
            buff.append(rst)
            buff.append("")

        buff.append("Conversation")
        buff.append("------------")
        methods = check_supported(format)
        for method in methods:
            if method.startswith("from_"):
                buff.append("Convert from this format to %s" % method_classes[method])
                buff.append("`" * len(buff[-1]))
            elif method.startswith("to_"):
                buff.append("Convert from %s to this format" % method_classes[method])
                buff.append("`" * len(buff[-1]))
            buff.append("")
            method_obj = getattr(format, method)
            docstring = method_obj.__doc__
            sig = signature(method_obj)
            parameters = dict(sig.parameters)
            # del self
            del parameters[list(parameters)[0]]
            # del data
            if method.startswith("to_"):
                del parameters[list(parameters)[0]]
            if "args" in parameters:
                del parameters["args"]
            if "kwargs" in parameters:
                del parameters["kwargs"]
            sig = Signature(list(parameters.values()))
            sig = str(sig)
            if method.startswith("from_") and method != "from_multi_systems":
                for aa in alias:
                    parameters["fmt"] = Parameter(
                        "fmt",
                        Parameter.POSITIONAL_OR_KEYWORD,
                        default=None,
                        annotation=Literal[aa],
                    )
                    sig_fmt = Signature(list(parameters.values()))
                    sig_fmt = str(sig_fmt)
                    buff.append(
                        f""".. py:function:: {method_classes[method]}.from{sig_fmt}"""
                    )
                    buff.append("""   :noindex:""")
                for aa in alias:
                    buff.append(
                        """.. py:function:: {}.from_{}{}""".format(
                            method_classes[method],
                            aa.replace("/", "_").replace(".", ""),
                            sig,
                        )
                    )
                    buff.append("""   :noindex:""")
                buff.append("")
                if docstring is not None:
                    doc_obj = SphinxDocString(docstring)
                    doc_obj["Parameters"] = [
                        xx
                        for xx in doc_obj["Parameters"]
                        if xx.name not in ("args", "kwargs")
                    ]
                    rst = str(doc_obj)
                    buff.append(rst)
                    buff.append("")
                else:
                    buff.append(
                        """   Convert this format to :func:`%s`."""
                        % (method_classes[method])
                    )
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
                    sig_fmt = Signature(list(parameters.values()))
                    sig_fmt = str(sig_fmt)
                    buff.append(
                        f""".. py:function:: {method_classes[method]}.to{sig_fmt}"""
                    )
                    buff.append("""   :noindex:""")
                for aa in alias:
                    buff.append(
                        """.. py:function:: {}.to_{}{}""".format(
                            method_classes[method],
                            aa.replace("/", "_").replace(".", ""),
                            sig,
                        )
                    )
                    buff.append("""   :noindex:""")
                buff.append("")
                if docstring is not None:
                    doc_obj = SphinxDocString(docstring)
                    doc_obj["Parameters"] = [
                        xx
                        for xx in doc_obj["Parameters"]
                        if xx.name not in ("data", "args", "kwargs")
                    ]
                    rst = str(doc_obj)
                    buff.append(rst)
                    buff.append("")
                else:
                    buff.append(
                        """   Convert :func:`%s` to this format."""
                        % (method_classes[method])
                    )
            buff.append("")
            buff.append("")

        with open("formats/%s.rst" % format.__name__, "w") as rstfile:
            rstfile.write("\n".join(buff))


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
                    "Alias": "\n".join("``%s``" % vvv for vvv in vv),
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
                    "Alias": "\n".join("``%s``" % vvv for vvv in vv),
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
                    "Alias": "\n".join("``%s``" % vvv for vvv in vv),
                }
            )
    generate_sub_format_pages(formats)
