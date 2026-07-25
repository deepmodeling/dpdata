from __future__ import annotations

import re
import unittest
from inspect import cleandoc, signature

from dpdata.format import Format

FORMAT_METHODS = (
    "from_system",
    "to_system",
    "from_labeled_system",
    "to_labeled_system",
    "from_bond_order_system",
    "to_bond_order_system",
    "from_multi_systems",
    "to_multi_systems",
)


class TestFormatDocumentation(unittest.TestCase):
    """Keep generated format reference pages useful as plugins evolve."""

    @staticmethod
    def format_classes():
        """Return each registered implementation once in deterministic order."""
        return sorted(set(Format.get_formats().values()), key=lambda cls: cls.__name__)

    def test_registered_formats_have_overviews(self):
        """Require every format page to have an implementation-owned overview."""
        for format_cls in self.format_classes():
            with self.subTest(format=format_cls.__name__):
                self.assertIsNotNone(
                    format_cls.__doc__,
                    f"{format_cls.__name__} needs a class docstring that introduces the format",
                )
                self.assertTrue(cleandoc(format_cls.__doc__))

    def test_conversion_parameters_are_documented(self):
        """Require every explicit conversion parameter in its NumPy docstring."""
        for format_cls in self.format_classes():
            for method_name in FORMAT_METHODS:
                if method_name not in format_cls.__dict__:
                    continue
                method = getattr(format_cls, method_name)
                docstring = cleandoc(method.__doc__ or "")
                with self.subTest(format=format_cls.__name__, method=method_name):
                    self.assertTrue(
                        docstring,
                        f"{format_cls.__name__}.{method_name} needs a docstring",
                    )
                    for parameter in signature(method).parameters.values():
                        if parameter.name == "self":
                            continue
                        prefix = ""
                        if parameter.kind == parameter.VAR_POSITIONAL:
                            prefix = r"\*"
                        elif parameter.kind == parameter.VAR_KEYWORD:
                            prefix = r"\*\*"
                        self.assertRegex(
                            docstring,
                            rf"(?m)^\s*{prefix}{re.escape(parameter.name)}\s*:",
                            f"{format_cls.__name__}.{method_name} must document "
                            f"the {parameter.name!r} parameter",
                        )


if __name__ == "__main__":
    unittest.main()
