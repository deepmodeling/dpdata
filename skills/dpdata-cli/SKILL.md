---
name: dpdata-cli
description: Compatibility entry point for the dpdata conversion workflow. Use dpdata-convert for user-facing conversion and data-contract validation.
compatibility: Requires uvx (uv) for the command examples.
metadata:
  repository: https://github.com/deepmodeling/dpdata
---

# dpdata CLI compatibility entry point

Conversion requests now use [`dpdata-convert`](../dpdata-convert/SKILL.md),
which combines the command reference with semantic validation and producer
contracts. This path remains so existing agents that request `dpdata-cli`
continue to resolve.

For the full command and format table, see
[`cli-reference.md`](../dpdata-convert/references/cli-reference.md).
