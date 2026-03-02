import polars as pl
import mzidentml_polars
import sys

def test():
    print(f"Python Polars version: {pl.__version__}")
    
    # 1. Define Protein Sequences
    prot_seqs = pl.DataFrame({
        "protein_id": ["PROT1", "PROT2"],
        "accession": ["P12345", "Q67890"],
        "sequence": ["MAGA", "MSRV"]
    })

    # 2. Define Identifications (CSMs)
    csms = pl.DataFrame({
        "spectrum_id": ["scan=123", "scan=456"],
        "peptide1_seq": ["PEPTIDEK", "PEPTIDEK"],
        "protein1_id": ["PROT1", "PROT2"],
        "peptide1_start": [1, 10],
        "peptide1_end": [8, 18],
        "charge": [2, 3],
        "rank": [1, 1],
        "is_crosslink": [False, True],
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

    # 3. Define Spectra
    spectra = pl.DataFrame({
        "spectrum_id": ["scan=123", "scan=456"],
        "file_path": ["data.mzML", "data.mzML"]
    })

    print("Attempting to call write_mzidentml...")
    try:
        xml_content = mzidentml_polars.write_mzidentml(csms, prot_seqs, spectra, {})
        print("Success! XML generated.")
        with open("test_output.mzid", "w") as f:
            f.write(xml_content)
        print("XML written to test_output.mzid")
    except Exception as e:
        print(f"FAILED with error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    test()
