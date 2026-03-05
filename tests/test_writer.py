import polars as pl
import pytest
import mzidentml_polars
import os
import tempfile
from lxml import etree

def validate_mzid(mzid_path, xsd_path):
    """Utility function to validate mzid against XSD."""
    if not os.path.exists(xsd_path):
        pytest.skip(f"XSD file not found at {xsd_path}")
    
    with open(xsd_path, 'rb') as f:
        schema_root = etree.XML(f.read())
        schema = etree.XMLSchema(schema_root)
        
    with open(mzid_path, 'rb') as f:
        doc = etree.XML(f.read())
        schema.assertValid(doc)

def test_basic_crosslinking(default_metadata, base_protein_seqs, base_spectra, xsd_path):
    """Test generating a basic crosslinked identification."""
    csms = pl.DataFrame({
        "spectrum_id": ["index=1", "index=2", "index=3"],
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

    xml_content = mzidentml_polars.serialize_mzidentml(csms, base_protein_seqs, base_spectra, default_metadata)
    assert xml_content is not None
    assert "mzIdentML" in xml_content
    
    with tempfile.NamedTemporaryFile(suffix=".mzid", delete=False) as tmp:
        tmp.write(xml_content.encode('utf-8'))
        tmp_path = tmp.name
        
    try:
        validate_mzid(tmp_path, xsd_path)
    finally:
        os.remove(tmp_path)

def test_ambiguity(default_metadata, base_protein_seqs, base_spectra, xsd_path):
    """Test generating mzIdentML with ambiguous protein matches."""
    csms = pl.DataFrame({
        "spectrum_id": ["index=1", "index=2"],
        "peptide1_seq": ["PEPTIDEK", "PEPTIDEK"],
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
        "peptide2_seq": ["KLS", None],
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

    xml_content = mzidentml_polars.serialize_mzidentml(csms, base_protein_seqs, base_spectra, default_metadata)
    assert xml_content is not None
    assert "PeptideEvidenceRef" in xml_content
    
    with tempfile.NamedTemporaryFile(suffix=".mzid", delete=False) as tmp:
        tmp.write(xml_content.encode('utf-8'))
        tmp_path = tmp.name
        
    try:
        validate_mzid(tmp_path, xsd_path)
    finally:
        os.remove(tmp_path)

def test_write_to_file(default_metadata, base_protein_seqs, base_spectra, xsd_path):
    """Test generating mzIdentML directly to a file."""
    csms = pl.DataFrame({
        "spectrum_id": ["index=1"],
        "peptide1_seq": ["PEPTIDEK"],
        "protein1_id": ["PROT1"],
        "peptide1_start": [1],
        "peptide1_end": [8],
        "charge": [2],
        "rank": [1],
        "is_crosslink": [False],
        "is_looplink": [False],
        "file_path": ["data1.mzML"],
        "peptide2_seq": [None],
        "protein2_id": [None],
        "peptide2_start": [None],
        "peptide2_end": [None],
        "peptide1_link_pos": [None],
        "peptide2_link_pos": [None],
    }).with_columns([
        pl.col("peptide1_start").cast(pl.UInt32),
        pl.col("peptide1_end").cast(pl.UInt32),
        pl.col("charge").cast(pl.Int32),
        pl.col("rank").cast(pl.UInt32),
        pl.col("is_crosslink").cast(pl.Boolean),
        pl.col("is_looplink").cast(pl.Boolean),
        pl.col("peptide2_start").cast(pl.UInt32),
        pl.col("peptide2_end").cast(pl.UInt32),
        pl.col("peptide1_link_pos").cast(pl.Int32),
        pl.col("peptide2_link_pos").cast(pl.Int32),
        pl.col("peptide2_seq").cast(pl.String),
        pl.col("protein2_id").cast(pl.String),
    ])

    with tempfile.NamedTemporaryFile(suffix=".mzid", delete=False) as tmp:
        tmp_path = tmp.name
        
    try:
        mzidentml_polars.write_mzidentml(csms, base_protein_seqs, base_spectra, default_metadata, tmp_path)
        assert os.path.exists(tmp_path)
        assert os.path.getsize(tmp_path) > 0
        validate_mzid(tmp_path, xsd_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
def test_filetype_derivation(default_metadata, base_protein_seqs, base_spectra, xsd_path):
    """Test that filetype is correctly derived from extension."""
    # Use .mgf for one of the spectra
    mgf_spectra = pl.DataFrame({
        "spectrum_id": ["index=1"],
        "file_path": ["test_data.mgf"]
    })
    
    csms = pl.DataFrame({
        "spectrum_id": ["index=1"],
        "peptide1_seq": ["PEPTIDEK"],
        "protein1_id": ["PROT1"],
        "peptide1_start": [1],
        "peptide1_end": [8],
        "charge": [2],
        "rank": [1],
        "is_crosslink": [False],
        "is_looplink": [False],
        "file_path": ["test_data.mgf"],
        "peptide2_seq": [None],
        "protein2_id": [None],
        "peptide2_start": [None],
        "peptide2_end": [None],
        "peptide1_link_pos": [None],
        "peptide2_link_pos": [None],
    }).with_columns([
        pl.col("peptide1_start").cast(pl.UInt32),
        pl.col("peptide1_end").cast(pl.UInt32),
        pl.col("charge").cast(pl.Int32),
        pl.col("rank").cast(pl.UInt32),
        pl.col("is_crosslink").cast(pl.Boolean),
        pl.col("is_looplink").cast(pl.Boolean),
        pl.col("peptide2_start").cast(pl.UInt32),
        pl.col("peptide2_end").cast(pl.UInt32),
        pl.col("peptide1_link_pos").cast(pl.Int32),
        pl.col("peptide2_link_pos").cast(pl.Int32),
        pl.col("peptide2_seq").cast(pl.String),
        pl.col("protein2_id").cast(pl.String),
    ])

    xml = mzidentml_polars.serialize_mzidentml(csms, base_protein_seqs, mgf_spectra, default_metadata)
    
    # Check for MGF format accession
    assert 'accession="MS:1001062"' in xml
    # Check for MGF nativeID format accession
    assert 'accession="MS:1000775"' in xml

def test_write_gzip(default_metadata, base_protein_seqs, base_spectra, xsd_path):
    """Test generating mzIdentML with Gzip compression."""
    csms = pl.DataFrame({
        "spectrum_id": ["index=1"],
        "peptide1_seq": ["PEPTIDEK"],
        "protein1_id": ["PROT1"],
        "peptide1_start": [1],
        "peptide1_end": [8],
        "charge": [2],
        "rank": [1],
        "is_crosslink": [False],
        "is_looplink": [False],
        "file_path": ["data1.mzML"],
        "peptide2_seq": [None],
        "protein2_id": [None],
        "peptide2_start": [None],
        "peptide2_end": [None],
        "peptide1_link_pos": [None],
        "peptide2_link_pos": [None],
    }).with_columns([
        pl.col("peptide1_start").cast(pl.UInt32),
        pl.col("peptide1_end").cast(pl.UInt32),
        pl.col("charge").cast(pl.Int32),
        pl.col("rank").cast(pl.UInt32),
        pl.col("is_crosslink").cast(pl.Boolean),
        pl.col("is_looplink").cast(pl.Boolean),
        pl.col("peptide2_start").cast(pl.UInt32),
        pl.col("peptide2_end").cast(pl.UInt32),
        pl.col("peptide1_link_pos").cast(pl.Int32),
        pl.col("peptide2_link_pos").cast(pl.Int32),
        pl.col("peptide2_seq").cast(pl.String),
        pl.col("protein2_id").cast(pl.String),
    ])

    with tempfile.NamedTemporaryFile(suffix=".mzid.gz", delete=False) as tmp:
        tmp_path = tmp.name
        
    try:
        mzidentml_polars.write_mzidentml(csms, base_protein_seqs, base_spectra, default_metadata, tmp_path)
        assert os.path.exists(tmp_path)
        
        # Verify it's actually a gzip file
        import gzip
        with gzip.open(tmp_path, 'rt') as f:
            content = f.read()
            assert "mzIdentML" in content
            
        # Optional: Full validation of de-compressed content
        with tempfile.NamedTemporaryFile(suffix=".mzid", delete=False) as decomp_tmp:
            decomp_path = decomp_tmp.name
            with gzip.open(tmp_path, 'rb') as f_in:
                decomp_tmp.write(f_in.read())
            
        try:
            validate_mzid(decomp_path, xsd_path)
        finally:
            os.remove(decomp_path)
            
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
