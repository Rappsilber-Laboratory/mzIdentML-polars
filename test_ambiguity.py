import polars as pl
import mzidentml_polars

def test_ambiguity():
    print(f"Testing ambiguity with Polars List columns...")
    
    # 1. Define Protein Sequences
    prot_seqs = pl.DataFrame({
        "protein_id": ["PROT1", "PROT2", "DECOY_PROT1"],
        "accession": ["P12345", "Q67890", "D12345"],
        "protein_name": ["BIPA_BACSU", "GCST_BACSU", None],
        "sequence": ["MAGA", "MSRV", "AGAM"],
        "is_decoy": [False, False, True]
    })

    # 2. Define Identifications (CSMs) with List columns to test ambiguity
    csms = pl.DataFrame({
        "spectrum_id": ["index=1", "index=2"],
        "peptide1_seq": ["PEPTIDEK", "PEPTIDEK"],
        # Row 1 is ambiguous for Peptide 1
        "protein1_id": [["PROT1", "PROT2"], ["PROT1"]],
        "peptide1_start": [[1, 5], [1]],
        "peptide1_end": [[8, 12], [8]],
        "charge": [2, 3],
        "rank": [1, 1],
        "is_crosslink": [True, False],
        "is_looplink": [False, False],
        "peptide1_link_pos": [8, None],
        "peptide2_link_pos": [1, None],
        "file_path": ["data1.mzML", "data1.mzML"],

        # Peptide 2 for crosslink
        "peptide2_seq": ["KLS", None],
        # Row 1 is ambiguous for Peptide 2 too
        "protein2_id": [["PROT2", "DECOY_PROT1"], [None]],
        "peptide2_start": [[10, 5], [None]],
        "peptide2_end": [[12, 7], [None]]
    }).with_columns([
        pl.col("peptide1_start").cast(pl.List(pl.UInt32)),
        pl.col("peptide1_end").cast(pl.List(pl.UInt32)),
        pl.col("peptide2_start").cast(pl.List(pl.UInt32)),
        pl.col("peptide2_end").cast(pl.List(pl.UInt32)),
        pl.col("charge").cast(pl.Int32),
        pl.col("rank").cast(pl.UInt32),
        pl.col("is_crosslink").cast(pl.Boolean),
        pl.col("is_looplink").cast(pl.Boolean),
        pl.col("peptide1_link_pos").cast(pl.Int32),
        pl.col("peptide2_link_pos").cast(pl.Int32),
    ])

    spectra = pl.DataFrame({
        "spectrum_id": ["index=1", "index=2"],
        "file_path": ["data1.mzML", "data1.mzML"]
    })

    metadata = {
        "software_name": "xi",
        "software_version": "1.7.6",
        "author": "Test Author",
        "parent_plus": 10.0,
        "parent_minus": 10.0,
        "frag_plus": 0.5,
        "frag_minus": 0.5,
    }

    try:
        xml_string = mzidentml_polars.write_mzidentml(csms, prot_seqs, spectra, metadata)
        with open("output_ambiguity_test.mzid", "w") as f:
            f.write(xml_string)
        print("Success! XML generated for ambiguity test.")
    except Exception as e:
        print(f"Error: {e}")
        exit(1)

if __name__ == "__main__":
    test_ambiguity()
