import polars as pl
import pytest
import mzidentml_polars
import os
import tempfile
import subprocess
from lxml import etree

def run_validator(mzid_path, tmpdir):
    """Run process_dataset -v -n -t on the given file."""
    env = os.environ.copy()
    try:
        import certifi
        env["SSL_CERT_FILE"] = certifi.where()
    except ImportError:
        pass
    
    result = subprocess.run(
        ["pipenv", "run", "process_dataset", "-v", mzid_path, "-n", "-t", tmpdir],
        capture_output=True,
        text=True,
        env=env
    )
    return result

def test_process_dataset_validation_mzml(default_metadata, base_protein_seqs):
    """Test that the generated mzML mzid is valid according to mzidentml-reader's process_dataset -v -n."""
    
    # Define spectra in an mzML file with various ID formats
    spectra = pl.DataFrame({
        "spectrum_id": ["controllerType=0 controllerNumber=1 scan=2231", "index=89209", "1234"],
        "file_path": ["test_data.mzML", "test_data.mzML", "test_data.mzML"]
    })
    
    csms = pl.DataFrame({
        "spectrum_id": ["controllerType=0 controllerNumber=1 scan=2231", "index=89209", "1234"],
        "peptide1_seq": ["PEPTIDEK", "PEPTIDEK", "KLS"],
        "protein1_id": ["PROT1", "PROT1", "PROT2"],
        "peptide1_start": [1, 1, 5],
        "peptide1_end": [8, 8, 7],
        "charge": [3, 2, 2],
        "rank": [1, 1, 1],
        "is_crosslink": [True, False, False],
        "is_looplink": [False, False, False],
        "file_path": ["test_data.mzML", "test_data.mzML", "test_data.mzML"],
        "peptide2_seq": ["KLS", None, None],
        "protein2_id": ["PROT2", None, None],
        "peptide2_start": [10, None, None],
        "peptide2_end": [12, None, None],
        "peptide1_link_pos": [8, None, None],
        "peptide2_link_pos": [1, None, None],
        "crosslinker_name": ["DSSO", None, None],
        "crosslinker_accession": ["MS:1003124", None, None],
        "crosslinker_mass": [158.0038, None, None],
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

    with tempfile.TemporaryDirectory() as tmpdir:
        mzid_path = os.path.join(tmpdir, "test_mzml.mzid")
        mzidentml_polars.write_mzidentml(mzid_path, csms, base_protein_seqs, spectra, default_metadata)
        
        result = run_validator(mzid_path, tmpdir)
        print("mzML STDOUT:", result.stdout)
        print("mzML STDERR:", result.stderr)
        
        assert "is schema valid" in result.stdout or "is schema valid" in result.stderr
        if result.returncode != 0 and "SSL: CERTIFICATE_VERIFY_FAILED" not in result.stderr:
            assert result.returncode == 0

def test_process_dataset_validation_mgf(default_metadata, base_protein_seqs):
    """Test that the generated MGF mzid is valid according to mzidentml-reader's process_dataset -v -n."""
    
    spectra = pl.DataFrame({
        "spectrum_id": ["TITLE=Spectrum1", "Scan 500"],
        "file_path": ["test_data.mgf", "test_data.mgf"]
    })
    
    csms = pl.DataFrame({
        "spectrum_id": ["TITLE=Spectrum1", "Scan 500"],
        "peptide1_seq": ["PEPTIDEK", "KLS"],
        "protein1_id": ["PROT1", "PROT2"],
        "peptide1_start": [1, 5],
        "peptide1_end": [8, 7],
        "charge": [2, 2],
        "rank": [1, 1],
        "is_crosslink": [False, False],
        "is_looplink": [False, False],
        "file_path": ["test_data.mgf", "test_data.mgf"],
        "peptide2_seq": [None, None],
        "protein2_id": [None, None],
        "peptide2_start": [None, None],
        "peptide2_end": [None, None],
        "peptide1_link_pos": [None, None],
        "peptide2_link_pos": [None, None],
    }).with_columns([
        pl.col("peptide1_start").cast(pl.UInt32),
        pl.col("peptide1_end").cast(pl.UInt32),
        pl.col("charge").cast(pl.Int32),
        pl.col("rank").cast(pl.UInt32),
        pl.col("is_crosslink").cast(pl.Boolean),
        pl.col("is_looplink").cast(pl.Boolean),
        pl.col("peptide1_link_pos").cast(pl.Int32),
        pl.col("peptide2_link_pos").cast(pl.Int32),
        pl.col("peptide2_seq").cast(pl.String),
        pl.col("protein2_id").cast(pl.String),
    ])

    with tempfile.TemporaryDirectory() as tmpdir:
        mzid_path = os.path.join(tmpdir, "test_mgf.mzid")
        mzidentml_polars.write_mzidentml(mzid_path, csms, base_protein_seqs, spectra, default_metadata)
        
        result = run_validator(mzid_path, tmpdir)
        print("MGF STDOUT:", result.stdout)
        print("MGF STDERR:", result.stderr)
        
        assert "is schema valid" in result.stdout or "is schema valid" in result.stderr
        if result.returncode != 0 and "SSL: CERTIFICATE_VERIFY_FAILED" not in result.stderr:
            assert result.returncode == 0
