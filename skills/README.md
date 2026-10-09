# dpdata agent skills

Start with [`dpdata-workflows`](dpdata-workflows/SKILL.md) to route a request by
user intent.

- [`dpdata-convert`](dpdata-convert/SKILL.md) converts producer outputs and
  checks the labeled-data contract.
- [`dpdata-inspect`](dpdata-inspect/SKILL.md) inspects and validates systems,
  including the deterministic validator in `scripts/`.
- [`dpdata-extend`](dpdata-extend/SKILL.md) groups developer references for
  format plugins, drivers, and minimizers.

The historical `dpdata-cli`, `dpdata-driver`, `dpdata-minimizer`, and
`dpdata-plugin` directories remain compatibility entry points.
