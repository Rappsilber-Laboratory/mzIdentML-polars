# mzIdentML-polars

A fast Rust-based writer for mzIdentML 1.3 files using Polars DataFrames as input. This project simplifies the generation of standard-compliant proteomics identification files, with built-in support for:

- **Polars Integration**: Directly write mzIdentML from high-performance DataFrames.
- **ProForma v2**: Support for standard peptide sequence notation (e.g., `PEPT[Unimod:35]IDEK`).
- **Crosslinking**: Native encoding for crosslinked peptide matches (CSMs).
- **mzIdentML 1.3.0 Compliance**: Generates valid XML according to the latest PSI-PI standards.

## Installation

You can install the Python bindings directly from the source using `maturin`:

```bash
# Clone the repository
git clone https://github.com/Rappsilber-Laboratory/mzIdentML-polars.git
cd mzIdentML-polars

# Install via pipenv (requires a Rust toolchain and maturin)
pipenv install
pipenv run maturin develop
```

## Usage

The primary functions are `write_mzidentml` (for file output) and `serialize_mzidentml` (for string output). Both take Polars DataFrames and a dictionary for metadata. Note that `write_mzidentml` takes the output path as its first argument.

### Writing to a File (Recommended)

This method is memory-efficient as it streams the XML directly to the disk. It also supports **automatic Gzip compression** if the filename ends in `.gz`.

```python
import polars as pl
import mzidentml_polars

# ... define DataFrames ...

# Generate mzIdentML directly to a file
mzidentml_polars.write_mzidentml("output.mzid", csms, prot_seqs, spectra, metadata)

# Automatic Gzip compression
mzidentml_polars.write_mzidentml("output.mzid.gz", csms, prot_seqs, spectra, metadata)
```

### Serializing to a String
```python
# Generate mzIdentML as a string (if needed for further processing)
xml_string = mzidentml_polars.serialize_mzidentml(csms, prot_seqs, spectra, metadata)
```

## Testing

The project includes a comprehensive test suite using `pytest` that validates output against official mzIdentML XML schemas.

### Prerequisites
```bash
pip install pytest lxml
```

### Running Tests
```bash
pipenv run pytest tests/
```

## Troubleshooting

### `TypeError: ... compat_level has invalid type: 'int'`
If you see this error, it indicates a version mismatch between your Python `polars` and the `polars` Rust crate used during compilation. 

The build process now **automatically synchronizes** these versions by updating `pyproject.toml` based on `Cargo.toml`. If you encounter this after manual dependency changes, simply rebuild the project:
```bash
pipenv run maturin develop
```
This will ensure your Python environment matches the compiled extension's expected ABI.

### `No module named 'pyarrow'`
`pyo3-polars` may require `pyarrow` for internal data conversions:
```bash
pip install pyarrow
```

## Input Schemas

### `prot_seqs` (DataFrame)
| Column | Type | Description |
| :--- | :--- | :--- |
| `protein_id` | String | Unique internal ID for the protein |
| `accession` | String | Public accession (e.g., UniProt) |
| `protein_name` | String | **Optional**. Descriptive name for the protein (e.g., `BIPA_BACSU`) |
| `sequence` | String | Full amino acid sequence |
| `is_decoy` | Boolean | Whether the protein is a decoy (default: `false`) |

### `csms` (DataFrame)
| Column | Type | Description |
| :--- | :--- | :--- |
| `spectrum_id` | String | ID of the spectrum (e.g., `index=1` or `scan=123`) |
| `file_path` | String | Path to the source file to resolve duplicate IDs across files. |
| `peptide1_seq` | String | ProForma v2 sequence of the first peptide |
| `protein1_id` | String / List[Str] | ID matching `prot_seqs` |
| `peptide1_start`| UInt32 / List[U32] | Start position in protein (1-based) |
| `peptide1_end` | UInt32 / List[U32] | End position in protein (1-based) |
| `charge` | Int32 | Precursor charge state |
| `rank` | UInt32 | Identification rank (1 = top match) |
| `is_crosslink` | Boolean | Whether this is a crosslink match |
| `is_looplink`  | Boolean | Whether this is a looplink match |
| `experimental_mz`| Float64| **Recommended**. Observed precursor m/z |
| `calculated_mz`| Float64| **Recommended**. Theoretical precursor m/z |
| `score` | Float64| **Recommended**. Primary search engine score |
| `peptide1_link_pos` | Int32 | 1-based link position on peptide 1 |
| `peptide2_link_pos` | Int32 | 1-based link position on peptide 2 (or site 2 for looplink) |
| `peptide2_seq` | String | (Crosslink only) Second peptide sequence |
| `protein2_id` | String / List[Str] | (Crosslink only) Second protein ID |
| `peptide2_start`| UInt32 / List[U32] | (Crosslink only) Start position |
| `peptide2_end` | UInt32 / List[U32] | (Crosslink only) End position |
| `crosslinker_name`| String | **Recommended**. Name of the crosslinker (e.g., `DSSO`) |
| `crosslinker_accession`| String | **Recommended**. CV accession of the crosslinker (e.g., `MS:1003124`) |
| `crosslinker_mass`| Float64| **Recommended**. Mass of the crosslinker |

### `metadata` (Dictionary)
| Key | Type | Description |
| :--- | :--- | :--- |
| `software_name`| String | Name of the analysis software (default: `mzidentml-polars`) |
| `software_version`| String | Version of the software |
| `author` | String | Name of the primary researcher/author |
| `is_ppm` | Boolean | Whether tolerances are in PPM (default: `true`) |
| `parent_plus` | Float | Precursor tolerance upper bound |
| `parent_minus` | Float | Precursor tolerance lower bound |
| `frag_plus` | Float | Fragment tolerance upper bound |
| `frag_minus` | Float | Fragment tolerance lower bound |
| `enzymes` | List[Dict] | Enzymes used: `[{"name": "Trypsin", "accession": "MS:1001251"}]` |
| `modifications`| List[Dict] | Search mods: `[{"fixed": true, "mass": 57.02, "residues": "C", "name": "Carbamidomethyl", "accession": "UNIMOD:4"}]` |
| `search_params`| List[Dict] | Additional parameters: `[{"name": "xi:score", "accession": "MS:1002545", "value": "0.5"}]` |

## Protein Ambiguity

If a peptide sequence maps to multiple proteins, you can encode this using Polars **List** columns in the `csms` DataFrame. For each mapped protein, provide the corresponding ID, start, and end positions in the lists. The library will generate multiple `<PeptideEvidence>` entries for that match.

```python
csms = pl.DataFrame({
    "protein1_id": [["PROT_A", "PROT_B"], ["PROT_C"]],
    "peptide1_start": [[1, 50], [10]],
    "peptide1_end": [[10, 60], [20]],
    # ... other columns
})
```

## Development & Releases

### Version Management
This project uses **Git tags** as the single source of truth for versioning.

- **Python**: Managed by `setuptools_scm`. The version is automatically derived from the latest Git tag (e.g., `v0.1.0`). If no tag is present, it defaults to a `.dev` version.
- **Rust**: The version is hardcoded in `Cargo.toml`. To ensure consistency, always use `cargo-release` to bump versions.

### Bumping the Version
To create a new release (e.g., moving from `0.1.0` to `0.2.0`):

1. **Install cargo-release**:
   ```bash
   cargo install cargo-release
   ```

2. **Run the release command**:
   ```bash
   # Dry run to verify changes
   cargo release minor --execute --no-publish
   ```
   This will:
   - Update the version in `Cargo.toml` and `Cargo.lock`.
   - Create a Git commit and a tag (e.g., `v0.2.0`).
   - Push the commit and the tag to the remote repository.

3. **CI/CD**:
   The GitHub Action (`.github/workflows/pypi.yml`) will automatically trigger on the new tag and publish the updated wheels to PyPI.

### Syncing Polars
If you change the `polars` version in `Cargo.toml`, the build script (`build.rs`) will automatically run `sync_polars.py` to update the constraints in `pyproject.toml`.

## License

This project is licensed under the Apache-2.0 License.

## TODO
- Implementing basic Protein Grouping (ProteinDetectionList) support, even as a simple 1-to-1 mapping if full inference isn't required.