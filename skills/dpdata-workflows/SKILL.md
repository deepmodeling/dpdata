# dpdata workflow routing

Use this page to choose the smallest dpdata skill that matches the user's
intent. The workflow skills own user-facing procedures; API and extension
material lives under `dpdata-extend`.

| User intent                                                           | Skill                                          | Load next                                                                                                    |
| --------------------------------------------------------------------- | ---------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| Convert a producer output, export a dataset, or prepare training data | [`dpdata-convert`](../dpdata-convert/SKILL.md) | [`semantic-validation.md`](../dpdata-convert/references/semantic-validation.md), then the producer reference |
| Inspect, validate, merge, or filter a dataset before it is consumed   | [`dpdata-inspect`](../dpdata-inspect/SKILL.md) | `validate_labeled_system.py` for deterministic checks                                                        |
| Add a format, driver, or minimizer                                    | [`dpdata-extend`](../dpdata-extend/SKILL.md)   | The matching developer reference                                                                             |

A conversion request is incomplete until the output contract is checked. Ask
whether the source is labeled, whether it contains one system or many, and
which producer/version generated it. Preserve the raw input and record the
producer contract, units, periodic-cell assumptions, atom ordering, and any
missing-label policy. Do not infer scientific conventions from a file name.

Load producer references only after the producer is known. They document the
contract to check and link to the parser or producer documentation; they do not
replace parser implementation or project-specific scientific choices.

The old `dpdata-cli`, `dpdata-driver`, `dpdata-minimizer`, and `dpdata-plugin`
paths remain compatibility entry points. New requests should use the workflow
skills above.
