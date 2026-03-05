import pytest
import polars as pl
import os

@pytest.fixture
def default_metadata():
    return {
        "software_name": "xi",
        "software_version": "2.0.beta",
        "author": "Max Mustermann",
        "parent_plus": 10.0,
        "parent_minus": 10.0,
        "frag_plus": 0.5,
        "frag_minus": 0.5,
        "is_ppm": True
    }

@pytest.fixture
def base_protein_seqs():
    return pl.DataFrame({
        "protein_id": ["PROT1", "PROT2", "DECOY_PROT1"],
        "accession": ["P12345", "Q67890", "D12345"],
        "protein_name": ["BIPA_BACSU", "GCST_BACSU", None],
        "sequence": ["MAGA", "MSRV", "AGAM"],
        "is_decoy": [False, False, True]
    })

@pytest.fixture
def base_spectra():
    return pl.DataFrame({
        "spectrum_id": ["index=1", "index=2", "index=3"],
        "file_path": ["data1.mzML", "data1.mzML", "data2.mzML"]
    })

@pytest.fixture
def xsd_path():
    # Return absolute path to the XSD file in context/
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(root_dir, "context", "mzIdentML1.3.0.xsd")
