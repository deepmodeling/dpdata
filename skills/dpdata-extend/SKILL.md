---
name: dpdata-extend
description: Extend dpdata for developers by adding format plugins, drivers, or minimizers. Use when implementing or packaging dpdata extensions rather than converting user data.
compatibility: Requires dpdata development knowledge; optional dependencies depend on the extension.
metadata:
  repository: https://github.com/deepmodeling/dpdata
---

# Extend dpdata

Use this skill for implementation-facing work. Choose one focused reference:

- [`format-plugin.md`](references/format-plugin.md) for a new format or an
  external `dpdata.plugins` package;
- [`driver.md`](references/driver.md) for `System.predict()` and label
  generation; or
- [`minimizer.md`](references/minimizer.md) for geometry optimization through
  `System.minimize()`.

An extension must state the data contract it produces. For a format reader,
that includes atom/type ordering, frame count, cells/PBC, labels, units, and
provenance expectations. Add a focused regression test for parser behavior and
keep producer-specific assumptions in the format reference or parser project.

This skill is deliberately separate from the user workflows. Route requests
to convert or inspect data through `dpdata-convert` or `dpdata-inspect` instead
of loading developer API references.
