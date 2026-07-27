from __future__ import annotations

import csv
from collections import defaultdict
from inspect import Parameter, Signature, cleandoc, signature
from pathlib import Path
from typing import Literal

from numpydoc.docscrape import Parameter as numpydoc_Parameter
from numpydoc.docscrape_sphinx import SphinxDocString

from dpdata.bond_order_system import BondOrderSystem

# ensure all plugins are loaded!
from dpdata.driver import Driver, Minimizer
from dpdata.format import Format
from dpdata.system import LabeledSystem, MultiSystems, System

DOCS_DIR = Path(__file__).resolve().parent
FORMAT_DIR = DOCS_DIR / "formats"


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


preferred_aliases = {
    "PwmatOutputFormat": "pwmat/output",
    "QuipGapXYZFormat": "extxyz",
    "VASPPoscarFormat": "vasp/poscar",
}

read_source_overrides = {
    "AmberMDFormat": '"trajectory_prefix"',
    "ASEStructureFormat": "atoms",
    "CP2KAIMDOutputFormat": '"calculation_directory"',
    "DFTBplusFormat": '("geometry.gen", "detailed.out")',
    "OPENMXFormat": '"system_name_prefix"',
    "PyMatgenMoleculeFormat": "molecule",
    "PyMatgenStructureFormat": "structure",
    "QECPTrajFormat": '"trajectory_prefix"',
    "QuipGapXYZFormat": '"data.xyz"',
}

internal_parameters = {
    "to_system": {"data"},
    "to_labeled_system": {"data"},
    "to_bond_order_system": {"data", "mol", "rdkit_mol"},
    "to_multi_systems": {"formulas"},
}

# A single ``to_system`` implementation can back several writer pages, so the
# generated docs filter the union of every dispatch-supplied name instead of
# only the current method's.
hidden_doc_parameters = {"*args"}.union(*internal_parameters.values())

location_parameters = {"file_name", "fname", "directory"}


def get_primary_alias(format_cls: type[Format], aliases: list[str]) -> str:
    """Return the alias used in headings and generated examples.

    Namespaced aliases are generally clearer than legacy short aliases. A few
    widely recognized names need an explicit preference because several
    equally valid ecosystem aliases share one implementation.
    """
    if format_cls.__name__ in preferred_aliases:
        return preferred_aliases[format_cls.__name__]
    return min(
        aliases,
        key=lambda alias: (alias.count("/") != 1, len(alias), alias),
    )


def sorted_format_items(formats: dict):
    """Return formats in deterministic, user-facing alias order."""
    return sorted(
        formats.items(),
        key=lambda item: (get_primary_alias(item[0], item[1]), item[0].__name__),
    )


def get_summary(format_cls: type[Format]) -> str:
    """Return the first docstring paragraph as a compact table summary."""
    docstring = cleandoc(format_cls.__doc__ or "")
    return " ".join(docstring.split("\n\n", maxsplit=1)[0].splitlines())


def get_method_obj(format_cls: type[Format], method: str):
    """Return the implementation that supplies a conversion's signature."""
    if (
        method == "to_labeled_system"
        and method not in format_cls.__dict__
        and "to_system" in format_cls.__dict__
    ):
        return getattr(format_cls, "to_system")
    return getattr(format_cls, method)


def get_user_parameters(
    method_obj, method: str, *, include_variadic: bool = False
) -> list[Parameter]:
    """Remove parameters supplied internally by dpdata's format dispatch."""
    hidden = internal_parameters.get(method, set())
    return [
        parameter
        for parameter in signature(method_obj).parameters.values()
        if parameter.name != "self"
        and parameter.name not in hidden
        and (
            include_variadic
            or parameter.kind not in (Parameter.VAR_POSITIONAL, Parameter.VAR_KEYWORD)
        )
    ]


def get_read_source(format_cls: type[Format], method_obj, method: str) -> str:
    """Return a meaningful placeholder for a generated loading example."""
    if method == "from_multi_systems" and format_cls.__name__ == "ASEStructureFormat":
        return '"trajectory.xyz"'
    if format_cls.__name__ in read_source_overrides:
        return read_source_overrides[format_cls.__name__]
    parameters = get_user_parameters(method_obj, method)
    if not parameters:
        return '"input_path"'
    name = parameters[0].name
    if name == "directory":
        return '"input_directory"'
    if name == "file_paths":
        return '("input_geometry", "output_file")'
    if name in {"file_name", "fname", "data"}:
        return '"input_file"'
    return name


def get_parameter_example(parameter: Parameter) -> str:
    """Return a readable placeholder for one required output parameter."""
    examples = {
        "basis": '"sto-3g"',
        "method": '"hf"',
    }
    return examples.get(parameter.name, f"{parameter.name}_value")


def get_write_arguments(method_obj, method: str) -> tuple[list[str], bool]:
    """Build user-supplied arguments for a generated ``to`` example."""
    arguments = []
    has_location = False
    for parameter in get_user_parameters(method_obj, method):
        if parameter.name in location_parameters:
            arguments.append('"output_path"')
            has_location = True
        elif parameter.default is Parameter.empty:
            arguments.append(f"{parameter.name}={get_parameter_example(parameter)}")
    return arguments, has_location


def append_quick_examples(
    buff: list[str], format_cls: type[Format], aliases: list[str], methods: list[str]
):
    """Append copyable examples for every supported dpdata object type."""
    alias = get_primary_alias(format_cls, aliases)
    lines = ["import dpdata", ""]

    read_examples = [
        ("from_system", "Geometry-only data", "system", "dpdata.System"),
        (
            "from_labeled_system",
            "Data with energies and forces",
            "labeled_system",
            "dpdata.LabeledSystem",
        ),
        (
            "from_bond_order_system",
            "Molecular graph and conformers",
            "bond_order_system",
            "dpdata.BondOrderSystem",
        ),
    ]
    for method, comment, variable, constructor in read_examples:
        if method not in methods:
            continue
        method_obj = get_method_obj(format_cls, method)
        source = get_read_source(format_cls, method_obj, method)
        lines.extend(
            [
                f"# {comment}",
                f'{variable} = {constructor}({source}, fmt="{alias}")',
                "",
            ]
        )

    if "from_multi_systems" in methods:
        method_obj = get_method_obj(format_cls, "from_multi_systems")
        source = get_read_source(format_cls, method_obj, "from_multi_systems")
        lines.extend(
            [
                "# Multiple compositions or calculation directories",
                f'systems = dpdata.MultiSystems.from_file({source}, fmt="{alias}")',
                "",
            ]
        )

    write_examples = [
        ("to_system", "system", "Write geometry-only data"),
        ("to_labeled_system", "labeled_system", "Write labeled data"),
        (
            "to_bond_order_system",
            "bond_order_system",
            "Write molecular graph data",
        ),
        ("to_multi_systems", "systems", "Write multiple systems"),
    ]
    for method, variable, comment in write_examples:
        if method not in methods:
            continue
        # ``to_labeled_system`` is automatically advertised when a System
        # writer is reused, so avoid showing the same call twice.
        if method == "to_labeled_system" and method not in format_cls.__dict__:
            continue
        method_obj = get_method_obj(format_cls, method)
        arguments, has_location = get_write_arguments(method_obj, method)
        call_args = ", ".join([f'"{alias}"', *arguments])
        call = f"{variable}.to({call_args})"
        lines.append(f"# {comment}")
        if has_location:
            lines.append(call)
        else:
            lines.append(f"converted = {call}")
        lines.append("")

    if lines[-1] == "":
        lines.pop()
    buff.extend(
        [
            "Quick examples",
            "--------------",
            "",
            f"The examples use the preferred alias ``{alias}``; any alias listed above is equivalent.",
            "",
            ".. code-block:: python",
            "",
            *(f"   {line}" for line in lines),
            "",
        ]
    )


def generate_sub_format_pages(formats: dict):
    """Generate one complete reference page per registered format class."""
    FORMAT_DIR.mkdir(exist_ok=True)
    for format_cls, aliases in sorted_format_items(formats):
        aliases = sorted(aliases)
        primary_alias = get_primary_alias(format_cls, aliases)
        buff = []
        buff.append(f".. _{format_cls.__name__}:")
        buff.append("")
        title = f"{primary_alias} format"
        buff.append(title)
        buff.append("=" * len(title))
        buff.append("")
        buff.append("Aliases")
        buff.append("-------")
        buff.append("")
        buff.append(", ".join(f"``{alias}``" for alias in aliases))
        buff.append("")
        buff.append(f"Implementation: {get_cls_link(format_cls)}")
        buff.append("")
        buff.append("Overview")
        buff.append("--------")
        buff.append("")

        docstring = format_cls.__doc__
        if docstring is not None:
            docstring = cleandoc(docstring)
            rst = str(SphinxDocString(docstring))
            buff.append(rst)
            buff.append("")

        methods = check_supported(format_cls)
        append_quick_examples(buff, format_cls, aliases, methods)

        buff.append("Conversions")
        buff.append("-----------")
        for method in methods:
            buff.append("")
            buff.append(f".. _{format_cls.__name__}_{method}:")
            buff.append("")
            if method.startswith("from_"):
                buff.append(f"Convert from this format to {method_classes[method]}")
                buff.append("`" * len(buff[-1]))
            elif method.startswith("to_"):
                buff.append(f"Convert from {method_classes[method]} to this format")
                buff.append("`" * len(buff[-1]))
            buff.append("")
            method_obj = get_method_obj(format_cls, method)
            docstring = method_obj.__doc__
            if docstring is not None:
                docstring = cleandoc(docstring)
            method_signature = signature(method_obj)
            parameters = get_user_parameters(method_obj, method, include_variadic=True)
            return_annotation = method_signature.return_annotation
            if method == "to_multi_systems" or method.startswith("from_"):
                sig = Signature(parameters, return_annotation=method_cls_obj[method])
            else:
                sig = Signature(parameters, return_annotation=return_annotation)
            sig = str(sig)
            if method.startswith("from_"):
                if method != "from_multi_systems":
                    for alias in aliases:
                        positional_parameters = [
                            parameter
                            for parameter in parameters
                            if parameter.kind
                            not in (
                                Parameter.VAR_POSITIONAL,
                                Parameter.KEYWORD_ONLY,
                                Parameter.VAR_KEYWORD,
                            )
                        ]
                        trailing_parameters = [
                            parameter
                            for parameter in parameters
                            if parameter not in positional_parameters
                        ]
                        constructor_parameters = [
                            *positional_parameters,
                            Parameter(
                                "fmt",
                                Parameter.POSITIONAL_OR_KEYWORD,
                                default=None,
                                annotation=Literal[alias],
                            ),
                            *trailing_parameters,
                        ]
                        sig_fmt = Signature(
                            constructor_parameters,
                            return_annotation=method_cls_obj[method],
                        )
                        sig_fmt = str(sig_fmt)
                        buff.append(
                            f""".. py:function:: dpdata.{method_classes[method]}{sig_fmt}"""
                        )
                        buff.append("""   :noindex:""")
                for alias in aliases:
                    buff.append(
                        """.. py:function:: dpdata.{}.from_{}{}""".format(
                            method_classes[method],
                            alias.replace("/", "_").replace(".", ""),
                            sig,
                        )
                    )
                    buff.append("""   :noindex:""")
                buff.append("")
                if docstring is None:
                    docstring = f"""   Convert this format to :class:`{method_classes[method]}`."""
                doc_obj = SphinxDocString(docstring)
                if len(doc_obj["Parameters"]) > 0:
                    doc_obj["Parameters"] = [
                        xx for xx in doc_obj["Parameters"] if xx.name != "*args"
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
                for alias in aliases:
                    to_parameters = [
                        Parameter(
                            "fmt",
                            Parameter.POSITIONAL_OR_KEYWORD,
                            annotation=Literal[alias],
                        ),
                        *parameters,
                    ]
                    if method == "to_multi_systems":
                        sig_fmt = Signature(
                            to_parameters,
                            return_annotation=method_cls_obj[method],
                        )
                    else:
                        sig_fmt = Signature(
                            to_parameters,
                            return_annotation=return_annotation,
                        )
                    sig_fmt = str(sig_fmt)
                    buff.append(
                        f""".. py:function:: dpdata.{method_classes[method]}.to{sig_fmt}"""
                    )
                    buff.append("""   :noindex:""")
                for alias in aliases:
                    buff.append(
                        """.. py:function:: dpdata.{}.to_{}{}""".format(
                            method_classes[method],
                            alias.replace("/", "_").replace(".", ""),
                            sig,
                        )
                    )
                    buff.append("""   :noindex:""")
                buff.append("")
                if docstring is None:
                    docstring = (
                        f"Convert :class:`{method_classes[method]}` to this format."
                    )
                doc_obj = SphinxDocString(docstring)
                if len(doc_obj["Parameters"]) > 0:
                    doc_obj["Parameters"] = [
                        xx
                        for xx in doc_obj["Parameters"]
                        if xx.name not in hidden_doc_parameters
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

        with (FORMAT_DIR / f"{format_cls.__name__}.rst").open(
            "w", encoding="utf-8"
        ) as rstfile:
            rstfile.write("\n".join(buff))


if __name__ == "__main__":
    formats = get_formats()
    with (DOCS_DIR / "formats.csv").open("w", encoding="utf-8", newline="") as csvfile:
        fieldnames = [
            "Format",
            "Description",
            "Alias",
            "Supported Conversions",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for kk, vv in sorted_format_items(formats):
            writer.writerow(
                {
                    "Format": f":ref:`{kk.__name__}`",
                    "Description": get_summary(kk),
                    "Alias": "\n".join(f"``{vvv}``" for vvv in sorted(vv)),
                    "Supported Conversions": "\n".join(
                        method_links[mtd].format(kk.__name__, mtd)
                        for mtd in check_supported(kk)
                    ),
                }
            )

    drivers = get_driver()
    with (DOCS_DIR / "drivers.csv").open("w", encoding="utf-8", newline="") as csvfile:
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
    with (DOCS_DIR / "minimizers.csv").open(
        "w", encoding="utf-8", newline=""
    ) as csvfile:
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
