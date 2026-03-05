# mzIdentML-polars Tests

This directory contains integration and unit tests for the Python bindings.

## Requirements
To run these tests, you'll need `pytest` and `lxml`:
```bash
pip install pytest lxml
```

## Running Tests
Run all tests using `pytest` from the project root:
```bash
pytest tests/
```

## Structure
- `conftest.py`: Common fixtures for peptide and protein DataFrames, and default metadata.
- `test_writer.py`: Integration tests for the `write_mzidentml` function, including crosslinks, looplinks, and ambiguous protein mappings.
- `data/`: Sample input/output files for testing (if any).

## Continuous Integration
The tests generate temporary `.mzid` files and validate them against the `mzIdentML1.3.0.xsd` schema found in the `context/` directory.
