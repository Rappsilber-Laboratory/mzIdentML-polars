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
    # Supports Linear, Crosslinked, and Looplinked peptides
    # Standards mandate 2 SpectrumIdentificationItems for crosslinks
    csms = pl.DataFrame({
        "spectrum_id": ["scan=123", "scan=456", "scan=789"],
        "peptide1_seq": ["PEPTIDEK", "PEPT[Unimod:35]IDEK", "PEPTIDEK"],
        "protein1_id": ["PROT1", "PROT2", "PROT1"],
        "peptide1_start": [1, 10, 1],
        "peptide1_end": [8, 18, 8],
        "charge": [2, 3, 2],
        "rank": [1, 1, 1],
        "is_crosslink": [False, True, False],
        "is_looplink": [False, False, True],
        "peptide1_link_pos": [None, 8, 2],
        "peptide2_link_pos": [None, 1, 8],
        
        # Required for crosslinks (is_crosslink = True)
        "peptide2_seq": [None, "KLS", None],
        "protein2_id": [None, "PROT1", None],
        "peptide2_start": [None, 5, None],
        "peptide2_end": [None, 12, None]
    }).with_columns([
        pl.col("peptide1_start").cast(pl.UInt32),
        pl.col("peptide1_end").cast(pl.UInt32),
        pl.col("charge").cast(pl.Int32),
        pl.col("rank").cast(pl.UInt32),
        pl.col("is_crosslink").cast(pl.Boolean),
        pl.col("is_looplink").cast(pl.Boolean),
        pl.col("peptide1_link_pos").cast(pl.Int32),
        pl.col("peptide2_link_pos").cast(pl.Int32),
        pl.col("peptide2_start").cast(pl.UInt32),
        pl.col("peptide2_end").cast(pl.UInt32),
    ])

    # 3. Define Spectra
    spectra = pl.DataFrame({
        "spectrum_id": ["scan=123", "scan=456", "scan=789"],
        "file_path": ["data.mzML", "data.mzML", "data.mzML"]
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
