import pytest
import mzidentml_polars
import polars as pl
import os

def test_reader_roundtrip(tmp_path, default_metadata, base_protein_seqs, base_spectra):
    """Test generating a basic crosslinked identification and reading it back."""
    tmp_file = str(tmp_path / "test.mzid")
    
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

    mzidentml_polars.write_mzidentml(tmp_file, csms_input, base_protein_seqs, base_spectra, default_metadata)
    csms_read, prot_read, spectra_read = mzidentml_polars.read_mzidentml(tmp_file)
    
    assert len(csms_read) == 3
    assert csms_read.filter(pl.col("is_looplink"))["peptide1_link_pos"][0] == 2
    assert csms_read.filter(pl.col("is_crosslink"))["peptide2_seq"][0] == "KLS"

def test_interaction_scores_parsing(tmp_path):
    """Test extracting interaction scores from ProteinDetectionHypothesis."""
    tmp_file = tmp_path / "interaction_test.mzid"
    xml_content = """<?xml version="1.0" encoding="UTF-8"?>
<mzIdentML id="test" xmlns="http://psidev.info/psi/pi/mzIdentML/1.3">
    <SequenceCollection>
        <DBSequence id="dbseq1" accession="P12345" />
        <Peptide id="pep1">
            <PeptideSequence>PEPTIDE</PeptideSequence>
        </Peptide>
        <PeptideEvidence id="pe1" dBSequence_ref="dbseq1" isDecoy="false" />
    </SequenceCollection>
    <AnalysisCollection>
        <SpectrumIdentification id="SI1" spectrumIdentificationProtocol_ref="SIP1" spectrumIdentificationList_ref="SIL1" />
    </AnalysisCollection>
    <DataCollection>
        <Inputs>
            <SpectraData id="SD1" location="data.mzML" />
        </Inputs>
        <AnalysisData>
            <SpectrumIdentificationList id="SIL1">
                <SpectrumIdentificationResult id="SIR1" spectrumID="index=1" spectraData_ref="SD1">
                    <SpectrumIdentificationItem id="SII1" rank="1" chargeState="2" peptide_ref="pep1" passThreshold="true">
                        <PeptideEvidenceRef peptideEvidence_ref="pe1" />
                        <cvParam cvRef="PSI-MS" accession="MS:1003344" name="residue pair ref" value="100.a" />
                    </SpectrumIdentificationItem>
                </SpectrumIdentificationResult>
            </SpectrumIdentificationList>
            <ProteinAmbiguityGroup id="PAG1">
                <ProteinDetectionHypothesis id="PDH1" dBSequence_ref="dbseq1" passThreshold="true">
                    <cvParam cvRef="PSI-MS" accession="MS:1002677" name="residue-pair-level global FDR" value="100.a:146:0.0294:true" />
                    <cvParam cvRef="PSI-MS" accession="MS:1002676" name="protein-pair-level global FDR" value="100.a:null:0.001:true" />
                </ProteinDetectionHypothesis>
            </ProteinAmbiguityGroup>
        </AnalysisData>
    </DataCollection>
</mzIdentML>
"""
    tmp_file.write_text(xml_content)
    csms_read, prot_read, spectra_read = mzidentml_polars.read_mzidentml(str(tmp_file))
    
    assert csms_read["residue_pair_fdr"][0] == 0.0294
    assert csms_read["protein_pair_fdr"][0] == 0.001

def test_gzip_roundtrip(tmp_path, default_metadata, base_protein_seqs, base_spectra):
    """Test compressed file IO."""
    tmp_file = str(tmp_path / "test.mzid.gz")
    csms_input = pl.DataFrame({
        "spectrum_id": ["index=1"], "peptide1_seq": ["PEPTIDEK"],
        "protein1_id": ["PROT1"], "peptide1_start": [1], "peptide1_end": [8],
        "charge": [2], "rank": [1], "is_crosslink": [False], "is_looplink": [False],
        "peptide1_link_pos": [None], "peptide2_link_pos": [None],
        "file_path": ["data1.mzML"], "experimental_mz": [1234.5], "calculated_mz": [1234.4],
        "score": [10.5], "crosslinker_name": ["DSSO"], "crosslinker_accession": ["MS:1003124"],
        "crosslinker_mass": [158.0038], "peptide2_seq": [None], "protein2_id": [None],
        "peptide2_start": [None], "peptide2_end": [None]
    }).with_columns([
        pl.col("peptide1_start").cast(pl.UInt32), pl.col("peptide1_end").cast(pl.UInt32),
        pl.col("charge").cast(pl.Int32), pl.col("rank").cast(pl.UInt32),
        pl.col("peptide2_seq").cast(pl.String), pl.col("protein2_id").cast(pl.String),
        pl.col("peptide1_link_pos").cast(pl.Int32), pl.col("peptide2_link_pos").cast(pl.Int32),
    ])

    mzidentml_polars.write_mzidentml(tmp_file, csms_input, base_protein_seqs, base_spectra, default_metadata)
    csms_read, prot_read, spectra_read = mzidentml_polars.read_mzidentml(tmp_file)
    
    # Detailed assertions for compressed read
    assert len(csms_read) == 1
    row = csms_read.row(0, named=True)
    assert row["peptide1_seq"] == "PEPTIDEK"
    assert row["charge"] == 2
    assert row["score"] == 10.5
    assert row["spectrum_id"] == "index=1"
    
    # Check protein and spectra DFs from compressed file
    assert len(prot_read) == 3
    assert "dbseq_PROT1" in prot_read["protein_id"].to_list()
    assert len(spectra_read) == 2
    assert "data1.mzML" in spectra_read["file_path"].to_list()



