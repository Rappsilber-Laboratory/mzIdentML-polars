import pytest
import mzidentml_polars
import polars as pl
import os
import tempfile
from polars.testing import assert_frame_equal

def test_reader_roundtrip(default_metadata, base_protein_seqs, base_spectra):
    """Test generating a basic crosslinked identification and reading it back."""
    
    # Input DataFrame replicating test_writer
    csms_input = pl.DataFrame({
        "spectrum_id": ["index=1", "index=2", "index=3"],
        "peptide1_seq": ["PEPTIDEK", "PEPT[UNIMOD:35]IDEK", "PEPTIDEK"],
        "protein1_id": ["DECOY_PROT1", "PROT2", "PROT1"],
        "peptide1_start": [1, 10, 1],
        "peptide1_end": [8, 18, 8],
        "charge": [2, 3, 2],
        "rank": [1, 1, 1],
        "is_crosslink": [False, True, False],
        "is_looplink": [False, False, True],
        "peptide1_link_pos": [None, 8, 2],
        "peptide2_link_pos": [None, 1, 8],
        "file_path": ["data1.mzML", "data1.mzML", "data2.mzML"],
        "experimental_mz": [1234.5, 678.9, 1234.5],
        "calculated_mz": [1234.4, 678.8, 1234.3],
        "score": [10.5, 20.1, 15.0],
        "crosslinker_name": ["DSSO", "DSSO", "DSSO"],
        "crosslinker_accession": ["MS:1003124", "MS:1003124", "MS:1003124"],
        "crosslinker_mass": [158.0038, 158.0038, 158.0038],
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

    with tempfile.NamedTemporaryFile(suffix=".mzid", delete=False) as tmp:
        tmp_path = tmp.name
        
    try:
        mzidentml_polars.write_mzidentml(tmp_path, csms_input, base_protein_seqs, base_spectra, default_metadata)
        
        # Read the generated file back
        csms_read, prot_read, spectra_read = mzidentml_polars.read_mzidentml(tmp_path)
        
        # Validate output shape and data types
        assert isinstance(csms_read, pl.DataFrame)
        assert len(csms_read) == 3
        
        # Validations
        assert "peptide1_link_pos" in csms_read.columns
        assert "peptide2_link_pos" in csms_read.columns
        
        # Verify looplink extraction mapping properties directly
        looplink_row = csms_read.filter(pl.col("is_looplink")).row(0, named=True)
        assert looplink_row["peptide1_link_pos"] == 2
        assert looplink_row["peptide2_link_pos"] == 8
        
        # Verify crosslink components
        crosslink_row = csms_read.filter(pl.col("is_crosslink")).row(0, named=True)
        assert crosslink_row["peptide2_seq"] == "KLS[MS:1003124]" or crosslink_row["peptide2_seq"] == "K[MS:1003124]LS"
        assert crosslink_row["peptide1_link_pos"] == 8
        assert crosslink_row["peptide2_link_pos"] == 1
        
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
