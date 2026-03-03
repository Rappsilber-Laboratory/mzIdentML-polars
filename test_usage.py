import polars as pl
import mzidentml_polars
import sys

def test():
    print(f"Python Polars version: {pl.__version__}")
    
    # 1. Define Protein Sequences
    prot_seqs = pl.DataFrame({
        "protein_id": ["PROT1", "PROT2", "DECOY_PROT1"],
        "accession": ["P12345", "Q67890", "D12345"],
        "sequence": ["MAGA", "MSRV", "AGAM"],
        "is_decoy": [False, False, True]
    })

    # 2. Define Identifications (CSMs)
    # Supports Linear, Crosslinked, and Looplinked peptides
    # Standards mandate 2 SpectrumIdentificationItems for crosslinks
    csms = pl.DataFrame({
        "spectrum_id": ["index=1", "index=2", "index=1"],
        "peptide1_seq": ["PEPTIDEK", "PEPT[Unimod:35]IDEK", "PEPTIDEK"],
        "protein1_id": ["DECOY_PROT1", "PROT2", "PROT1"],
        "peptide1_start": [1, 10, 1],
        "peptide1_end": [8, 18, 8],
        "charge": [2, 3, 2],
        "rank": [1, 1, 1],
        "is_crosslink": [False, True, False],
        "is_looplink": [False, False, True],
        "peptide1_link_pos": [None, 8, 2],
        "peptide2_link_pos": [None, 1, 8],
        
        # Explicitly link CSM to file (required for multi-file datasets)
        "file_path": ["data1.mzML", "data1.mzML", "data2.mzML"],

        # Recommended metadata (improves xiView/downstream compatibility)
        "experimental_mz": [1234.5, 678.9, 1234.5],
        "calculated_mz": [1234.4, 678.8, 1234.3],
        "score": [10.5, 20.1, 15.0],
        "crosslinker_name": ["DSSO", "DSSO", "DSSO"],
        "crosslinker_accession": ["MS:1003124", "MS:1003124", "MS:1003124"],
        "crosslinker_mass": [158.0038, 158.0038, 158.0038],

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

    # 3. Define Spectra (Linking files to IDs)
    spectra = pl.DataFrame({
        "spectrum_id": ["index=1", "index=2", "index=1"],
        "file_path": ["data1.mzML", "data1.mzML", "data2.mzML"]
    })

    # 4. Generate mzIdentML XML
    metadata = {
        "software_name": "xi",
        "parent_plus": 10.0,
        "parent_minus": 10.0,
        "frag_plus": 0.5,
        "frag_minus": 0.5,
        "is_ppm": True
    }
    
    print("Attempting to call write_mzidentml...")
    try:
        xml_content = mzidentml_polars.write_mzidentml(csms, prot_seqs, spectra, metadata)
        print("Success! XML generated.")
        output_file = "test_output.mzid"
        with open(output_file, "w") as f:
            f.write(xml_content)
        print(f"XML written to {output_file}")
    except Exception as e:
        print(f"FAILED with error: {e}")
        sys.exit(1)

    # 4. Verify with mzidentml-reader
    print("\nVerifying with mzidentml-reader...")
    from parser.process_dataset import sequences_and_residue_pairs
    import tempfile
    import os

    try:
        tmpdir = tempfile.gettempdir()
        data = sequences_and_residue_pairs(output_file, tmpdir)
        print("mzidentml-reader successfully parsed the file!")
        # Basic validation of parsed data
        if "residue_pairs" in data and len(data["residue_pairs"]) > 0:
            print(f"Found {len(data['residue_pairs'])} residue pairs (crosslinks).")
        else:
            print("Warning: No residue pairs found by mzidentml-reader.")
    except Exception as e:
        print(f"mzidentml-reader FAILED to parse the file: {e}")
        # Not exiting here so we can still see the output

    # 5. Schema Validation
    print("\nValidating against XSD...")
    try:
        from lxml import etree
        schema_file = "mzIdentML1.3.0.xsd"
        if os.path.exists(schema_file):
            with open(schema_file, 'rb') as f:
                schema_root = etree.XML(f.read())
                schema = etree.XMLSchema(schema_root)
                
            with open(output_file, 'rb') as f:
                doc = etree.XML(f.read())
                schema.assertValid(doc)
                print("Schema validation SUCCESSFUL!")
        else:
            print(f"Warning: {schema_file} not found, skipping schema validation.")
    except Exception as e:
        print(f"Schema validation FAILED: {e}")

if __name__ == "__main__":
    test()
