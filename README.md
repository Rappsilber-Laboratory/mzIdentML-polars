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

# Install via pip (requires a Rust toolchain and maturin)
pip install maturin
pip install .
```

## Usage

The primary function is `write_mzidentml`, which takes three Polars DataFrames and a dictionary for metadata.

```python
import polars as pl
import mzidentml_polars

# 1. Define Protein Sequences
prot_seqs = pl.DataFrame({
    "protein_id": ["PROT1", "PROT2"],
    "accession": ["P12345", "Q67890"],
    "sequence": ["MAGA...END", "MSRV...STOP"]
})

# 2. Define Identifications (CSMs)
# Supports both Linear and Crosslinked peptides
# NOTE: Columns must be cast to expected types (UInt32, Int32, Boolean)
csms = pl.DataFrame({
    "spectrum_id": ["scan=123", "scan=456"],
    "peptide1_seq": ["PEPTIDEK", "PEPT[Unimod:35]IDEK"],
    "protein1_id": ["PROT1", "PROT2"],
    "peptide1_start": [1, 10],
    "peptide1_end": [8, 18],
    "charge": [2, 3],
    "rank": [1, 1],
    "is_crosslink": [False, True],
    
    # Required for crosslinks (is_crosslink = True)
    "peptide2_seq": [None, "KLS"],
    "protein2_id": [None, "PROT1"],
    "peptide2_start": [None, 5],
    "peptide2_end": [None, 12]
}).with_columns([
    pl.col("peptide1_start").cast(pl.UInt32),
    pl.col("peptide1_end").cast(pl.UInt32),
    pl.col("charge").cast(pl.Int32),
    pl.col("rank").cast(pl.UInt32),
    pl.col("is_crosslink").cast(pl.Boolean),
    pl.col("peptide2_start").cast(pl.UInt32),
    pl.col("peptide2_end").cast(pl.UInt32),
])

# 3. Define Spectra (placeholder for now)
spectra = pl.DataFrame({
    "spectrum_id": ["scan=123", "scan=456"],
    "file_path": ["data.mzML", "data.mzML"]
})

# 4. Generate mzIdentML XML
xml_content = mzidentml_polars.write_mzidentml(csms, prot_seqs, spectra, {})

with open("output.mzid", "w") as f:
    f.write(xml_content)
```

## Troubleshooting

### `TypeError: ... compat_level has invalid type: 'int'`
If you see this error, it indicates a version mismatch between your Python `polars` and the `pyo3-polars` used during compilation. As of now, ensure you are using a compatible version of Polars:
```bash
pip install polars==1.31.0
```

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
| `sequence` | String | Full amino acid sequence |

### `csms` (DataFrame)
| Column | Type | Description |
| :--- | :--- | :--- |
| `spectrum_id` | String | ID of the spectrum in the source file |
| `peptide1_seq` | String | ProForma v2 sequence of the first peptide |
| `protein1_id` | String | ID matching `prot_seqs` |
| `peptide1_start`| UInt32 | Start position in protein (1-based) |
| `peptide1_end` | UInt32 | End position in protein (1-based) |
| `charge` | Int32 | Precursor charge state |
| `rank` | UInt32 | Identification rank (1 = top match) |
| `is_crosslink` | Boolean | Whether this is a crosslink match |
| `peptide2_seq` | String | (Crosslink only) Second peptide sequence |
| `protein2_id` | String | (Crosslink only) Second protein ID |
| `peptide2_start`| UInt32 | (Crosslink only) Start position |
| `peptide2_end` | UInt32 | (Crosslink only) End position |

## License

This project is licensed under the AGPL-3.0 License.
