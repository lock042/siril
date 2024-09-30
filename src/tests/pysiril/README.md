# Siril Command Tests

This directory contains tests for Siril commands written in Python to simplify the testing process, using [pysiril](https://siril.org/tutorials/pysiril/) to connect to `siril-cli`.

The **Criterion C** tests can be used for unit testing. These tests focus on command testing and regression checks for important parts of large pieces of code.

Some lightweight test images can also be placed in this directory.

### Naming Conventions

- Files and test functions must start with `test_`.

### Installation

1. Install `pysiril` by following the instructions [here](https://siril.org/tutorials/pysiril/).
2. Install `pytest`.

### Running Tests

To run the tests, execute the following command in this directory:

```bash
pytest
```

### TODO List

- [ ] Make the `siril-cli` path configurable (currently works only with autotools).
- [ ] Ensure only one `siril-cli` instance is used, without affecting the main config file, and use default settings on every start.
- [ ] Write more tests.
- [ ] Add tests to the CI pipeline.
