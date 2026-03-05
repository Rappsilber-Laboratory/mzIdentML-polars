import polars as pl
import mzidentml_polars
import pytest

def test_cv_lookup(default_metadata, base_protein_seqs, base_spectra):
    """Test that mods in ProForma are correctly looked up in the packaged CV."""
    csms = pl.DataFrame({
        "spectrum_id": ["index=1"],
        "peptide1_seq": ["PEPT[Unimod:35]IDEK"],
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

    xml = mzidentml_polars.serialize_mzidentml(csms, base_protein_seqs, base_spectra, default_metadata)
    
    # Check if 'Oxidation' name was retrieved from CV
    assert 'name="Oxidation"' in xml
    # Check if monoisotopic mass delta was retrieved
    assert 'monoisotopicMassDelta="15.994915"' in xml
    # Check if CV ref is correct
    assert 'cvRef="UNIMOD"' in xml

def test_cv_lookup_by_name(default_metadata, base_protein_seqs, base_spectra):
    """Test that mods in ProForma (by name) are correctly looked up in the packaged CV."""
    csms = pl.DataFrame({
        "spectrum_id": ["index=1"],
        "peptide1_seq": ["PEPT[Oxidation]IDEK"],
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

    xml = mzidentml_polars.serialize_mzidentml(csms, base_protein_seqs, base_spectra, default_metadata)
    
    # Check if 'Oxidation' name was retrieved
    assert 'name="Oxidation"' in xml
    # Check if monoisotopic mass delta was retrieved
    assert 'monoisotopicMassDelta="15.994915"' in xml
    # Check if accession was correctly mapped
    assert 'accession="UNIMOD:35"' in xml

def test_xlmod_lookup(default_metadata, base_protein_seqs, base_spectra):
    """Test that XLMOD terms are correctly looked up."""
    csms = pl.DataFrame({
        "spectrum_id": ["index=1"],
        "peptide1_seq": ["PEPT[XLMOD:02001]IDEK"],
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

    xml = mzidentml_polars.serialize_mzidentml(csms, base_protein_seqs, base_spectra, default_metadata)
    
    # Check if 'DSS' name was retrieved from XLMOD
    assert 'name="DSS"' in xml
    # Check if monoisotopic mass delta was retrieved
    assert 'monoisotopicMassDelta="138.06807961"' in xml
    # Check if CV ref is correct
    assert 'cvRef="XLMOD"' in xml

def test_xlmod_lookup_by_name(default_metadata, base_protein_seqs, base_spectra):
    """Test that XLMOD terms by name are correctly looked up."""
    csms = pl.DataFrame({
        "spectrum_id": ["index=1"],
        "peptide1_seq": ["PEPT[DSS]IDEK"],
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

    xml = mzidentml_polars.serialize_mzidentml(csms, base_protein_seqs, base_spectra, default_metadata)
    
    assert 'name="DSS"' in xml
    assert 'accession="XLMOD:02001"' in xml
    assert 'monoisotopicMassDelta="138.06807961"' in xml

def test_psimod_lookup(default_metadata, base_protein_seqs, base_spectra):
    """Test that PSI-MOD terms are correctly looked up."""
    csms = pl.DataFrame({
        "spectrum_id": ["index=1"],
        "peptide1_seq": ["PEPT[MOD:00696]IDEK"],
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

    xml = mzidentml_polars.serialize_mzidentml(csms, base_protein_seqs, base_spectra, default_metadata)
    
    # Check if name was retrieved
    assert 'name="phosphorylated residue"' in xml
    # Check if monoisotopic mass delta was retrieved
    assert 'monoisotopicMassDelta="79.966331"' in xml
    # Check if CV ref is correct
    assert 'cvRef="PSI-MOD"' in xml
