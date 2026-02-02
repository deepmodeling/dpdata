# dpdata - Atomistic Data Format Manipulation

dpdata is a Python package for manipulating atomistic data from computational science software. It supports format conversion between various atomistic simulation packages including VASP, DeePMD-kit, LAMMPS, GROMACS, Gaussian, ABACUS, and many others.

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Working Effectively

- **Bootstrap and install the repository:**
  - `cd /home/runner/work/dpdata/dpdata` (or wherever the repo is cloned)
  - `uv pip install -e .` -- installs dpdata in development mode with core dependencies (numpy, scipy, h5py, monty, wcmatch)
  - Test installation: `dpdata --version` -- should show version like "dpdata v0.1.dev2+..."

- **Run tests:**
  - `cd tests && python -m unittest discover` -- runs all 1826 tests in ~10 seconds. NEVER CANCEL.
  - `cd tests && python -m unittest test_<module>.py` -- run specific test modules (individual modules take ~0.5 seconds)
  - `cd tests && coverage run --source=../dpdata -m unittest discover && coverage report` -- run tests with coverage

- **Linting and formatting:**
  - Install ruff: `uv pip install ruff`
  - `ruff check dpdata/` -- lint the main package (takes ~1 second)
  - `ruff format dpdata/` -- format code according to project style
  - `ruff check --fix dpdata/` -- auto-fix linting issues where possible

- **Pre-commit hooks:**
  - Install: `uv pip install pre-commit`
  - `pre-commit run --all-files` -- run all hooks on all files
  - Hooks include: ruff linting/formatting, trailing whitespace, end-of-file-fixer, yaml/json/toml checks

## Validation

- **Always test CLI functionality after making changes:**
  - `dpdata --help` -- ensure CLI still works
  - `dpdata --version` -- verify version is correct
  - Test a basic conversion if sample data is available

- **Always run linting before committing:**
  - `ruff check dpdata/` -- ensure no new linting errors
  - `ruff format dpdata/` -- ensure code is properly formatted

- **Run relevant tests for your changes:**
  - For format-specific changes: `cd tests && python -m unittest test_<format>*.py`
  - For core system changes: `cd tests && python -m unittest test_system*.py test_multisystems.py`
  - For CLI changes: `cd tests && python -m unittest test_cli.py` (if exists)

## Build and Documentation

- **Documentation:**
  - `cd docs && make help` -- see all available build targets
  - `cd docs && make html` -- build HTML documentation (requires additional dependencies)
  - Documentation source is in `docs/` directory using Sphinx
  - **NOTE:** Full docs build requires additional dependencies like `deepmodeling-sphinx` that may not be readily available

- **Package building:**
  - Uses setuptools with pyproject.toml configuration
  - `uv pip install build && python -m build` -- create source and wheel distributions
  - Version is managed by setuptools_scm from git tags

## Common Tasks

The following are outputs from frequently run commands. Reference them instead of re-running to save time.

### Repository structure
```
/home/runner/work/dpdata/dpdata/
├── dpdata/           # Main package code
│   ├── __init__.py
│   ├── cli.py        # Command-line interface
│   ├── system.py     # Core System classes
│   ├── format.py     # Format registry
│   ├── abacus/       # ABACUS format support
│   ├── amber/        # AMBER format support
│   ├── deepmd/       # DeePMD format support
│   ├── vasp/         # VASP format support
│   ├── xyz/          # XYZ format support
│   └── ...          # Other format modules
├── tests/            # Test suite (91 test files)
├── docs/             # Sphinx documentation
├── plugin_example/   # Example plugin
├── pyproject.toml    # Project configuration
└── README.md
```

### Key dependencies
- Core: numpy>=1.14.3, scipy, h5py, monty, wcmatch
- Optional: ase (ASE integration), parmed (AMBER), pymatgen (Materials Project), rdkit (molecular analysis)
- Testing: unittest (built-in), coverage
- Linting: ruff
- Docs: sphinx with various extensions

### Test timing expectations
- Full test suite: ~10 seconds (1826 tests). NEVER CANCEL.
- Individual test modules: ~0.5 seconds
- Linting with ruff: ~1 second
- Documentation build: ~30 seconds

### Common workflows
1. **Adding a new format:**
   - Create module in `dpdata/<format>/`
   - Implement format classes inheriting from appropriate base classes
   - Add tests in `tests/test_<format>*.py`
   - Register format in the plugin system

2. **Fixing bugs:**
   - Write test that reproduces the bug first
   - Make minimal fix to pass the test
   - Run full test suite to ensure no regressions
   - Run linting to ensure code style compliance

3. **CLI changes:**
   - Modify `dpdata/cli.py`
   - Test with `dpdata --help` and specific commands
   - Add/update tests if needed

## Troubleshooting

- **Installation timeouts:** Network timeouts during `uv pip install` are common. If this occurs, try:
  - Individual package installation: `uv pip install numpy scipy h5py monty wcmatch`
  - Use `--timeout` option: `uv pip install --timeout 300 -e .`
  - Verify existing installation works: `dpdata --version` should work even if reinstall fails

- **Optional dependency errors:** Many tests will skip or fail if optional dependencies (ase, parmed, pymatgen, rdkit) are not installed. This is expected. Core functionality will work with just the basic dependencies.

- **Documentation build failures:** The docs build requires specific dependencies like `deepmodeling-sphinx` that may not be readily available. Use `make help` to see available targets, but expect build failures without full doc dependencies.

- **Test artifacts:** The test suite generates temporary files (`tests/data_*`, `tests/tmp.*`, `tests/.coverage`). These are excluded by `.gitignore` and should not be committed.

- **Import errors:** If you see import errors for specific modules, check if the corresponding optional dependency is installed. For example, ASE functionality requires `uv pip install ase`.

## Critical Notes

- **NEVER CANCEL** test runs or builds - they complete quickly (10 seconds for tests, 30 seconds for docs)
- Always run `ruff check` and `ruff format` before committing
- Test artifacts in `tests/` directory are excluded by `.gitignore` - don't commit them
- Optional dependencies are required for some formats but core functionality works without them
- The CLI tool `dpdata` is the main user interface for format conversion

## Commit and PR Guidelines

- **Use semantic commit messages** for all commits and PR titles following the format: `type(scope): description`
  - **Types:** `feat` (new feature), `fix` (bug fix), `docs` (documentation), `style` (formatting), `refactor` (code restructuring), `test` (testing), `chore` (maintenance)
  - **Examples:**
    - `feat(vasp): add support for POSCAR format`
    - `fix(cli): resolve parsing error for multi-frame files`
    - `docs: update installation instructions`
    - `test(amber): add tests for trajectory parsing`
- **PR titles** must follow semantic commit format
- **Commit messages** should be concise but descriptive of the actual changes made
