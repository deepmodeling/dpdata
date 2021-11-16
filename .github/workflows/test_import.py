name: test Python import

on:
  - push
  - pull_request

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v2
    - uses: actions/setup-python@v2
        with:
            python-version: '3.x'
            architecture: 'x64'
    - run: python -m pip install .
    - run: python -c 'import dpdata'

