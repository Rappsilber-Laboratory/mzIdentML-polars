use pyo3::prelude::*;
use pyo3::types::{PyDict};
use pyo3_polars::PyDataFrame;
use std::collections::HashMap;
use crate::mzidentml::psi_pi::*;
use xsd_parser::quick_xml::WithSerializer;
use xsd_parser::quick_xml::Writer;

/// Factory to manage the construction of the MzIdentML XML tree.
pub struct MzIdentMLFactory {
    pub doc: MzIdentMlType,
    peptide_map: HashMap<String, String>, // ProForma -> Peptide ID
    db_seq_map: HashMap<String, String>,  // Protein ID -> DBSequence ID
    pep_evidence_map: HashMap<(String, String), String>, // (Peptide ID, Protein ID) -> Evidence ID
    cv_map: HashMap<String, String>,      // CV ID -> URI
}

impl MzIdentMLFactory {
    pub fn new(id: String) -> Self {
        let doc = MzIdentMlType {
            id,
            name: None,
            creation_date: None, // TODO: Add current date
            version: "1.3.0".to_string(),
            cv_list: CvListType { cv: Vec::new() },
            cv_param: Vec::new(),
            analysis_software_list: Some(AnalysisSoftwareListType {
                analysis_software: Vec::new(),
            }),
            provider: None,
            audit_collection: None,
            analysis_sample_collection: None,
            sequence_collection: Some(SequenceCollectionType {
                db_sequence: Vec::new(),
                peptide: Vec::new(),
                peptide_evidence: Vec::new(),
            }),
            analysis_collection: AnalysisCollectionType {
                spectrum_identification: Vec::new(),
                protein_detection: None,
            },
            analysis_protocol_collection: AnalysisProtocolCollectionType {
                spectrum_identification_protocol: Vec::new(),
                protein_detection_protocol: None,
            },
            data_collection: DataCollectionType {
                inputs: InputsType {
                    source_file: Vec::new(),
                    search_database: Vec::new(),
                    spectra_data: Vec::new(),
                },
                analysis_data: AnalysisDataType {
                    spectrum_identification_list: Vec::new(),
                    protein_detection_list: None,
                },
            },
            bibliographic_reference: Vec::new(),
        };

        // Add standard CVs
        let mut factory = Self {
            doc,
            peptide_map: HashMap::new(),
            db_seq_map: HashMap::new(),
            pep_evidence_map: HashMap::new(),
            cv_map: HashMap::new(),
        };

        factory.add_cv("PSI-MS", "PSI-MS", "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo");
        factory.add_cv("UNIMOD", "UNIMOD", "http://www.unimod.org/obo/unimod.obo");
        factory.add_cv("UO", "Unit Ontology", "http://purl.obolibrary.org/obo/uo.obo");
        factory.add_cv("XLMOD", "PSI-XLMOD", "https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/cv/XLMOD.obo");

        factory
    }

    pub fn add_cv(&mut self, id: &str, name: &str, uri: &str) {
        if self.cv_map.contains_key(id) {
            return;
        }
        self.doc.cv_list.cv.push(CvType {
            id: id.to_string(),
            full_name: name.to_string(),
            uri: uri.to_string(),
            version: None,
        });
        self.cv_map.insert(id.to_string(), uri.to_string());
    }

    pub fn add_peptide(&mut self, proforma: &str) -> String {
        if let Some(id) = self.peptide_map.get(proforma) {
            return id.clone();
        }

        let id = format!("pep_{}", self.peptide_map.len());
        
        // For now, just store the sequence. 
        // Parsing modifications would require a more detailed rustyms integration.
        let peptide_struct = PeptideType {
            id: id.clone(),
            name: None,
            content: vec![PeptideTypeContent::PeptideSequence(proforma.to_string())],
        };

        if let Some(sc) = &mut self.doc.sequence_collection {
            sc.peptide.push(peptide_struct);
        }
        self.peptide_map.insert(proforma.to_string(), id.clone());
        id
    }

    pub fn add_db_sequence(&mut self, protein_id: &str, accession: &str, sequence: &str, db_ref: &str) -> String {
        if let Some(id) = self.db_seq_map.get(protein_id) {
            return id.clone();
        }

        let id = format!("dbseq_{}", protein_id);
        let db_seq = DbSequenceType {
            id: id.clone(),
            name: None,
            length: Some(sequence.len() as i32),
            search_database_ref: db_ref.to_string(),
            accession: accession.to_string(),
            content: vec![DbSequenceTypeContent::Seq(sequence.to_string())],
        };

        if let Some(sc) = &mut self.doc.sequence_collection {
            sc.db_sequence.push(db_seq);
        }
        self.db_seq_map.insert(protein_id.to_string(), id.clone());
        id
    }

    pub fn add_peptide_evidence(&mut self, pep_ref: &str, db_ref: &str, start: Option<u32>, end: Option<u32>, is_decoy: bool) -> String {
        let key = (pep_ref.to_string(), db_ref.to_string());
        if let Some(id) = self.pep_evidence_map.get(&key) {
            return id.clone();
        }

        let id = format!("ev_{}_{}", pep_ref, db_ref);
        let evidence = PeptideEvidenceType {
            id: id.clone(),
            name: None,
            db_sequence_ref: db_ref.to_string(),
            peptide_ref: pep_ref.to_string(),
            start: start.map(|s| s as i32),
            end: end.map(|e| e as i32),
            pre: None,
            post: None,
            translation_table_ref: None,
            frame: None,
            is_decoy,
            content: Vec::new(),
        };

        if let Some(sc) = &mut self.doc.sequence_collection {
            sc.peptide_evidence.push(evidence);
        }
        self.pep_evidence_map.insert(key, id.clone());
        id
    }

    pub fn add_sii(&mut self, result_id: &str, sii: SpectrumIdentificationItemType, spectra_data_ref: &str) {
        // Ensure the SpectrumIdentificationList exists
        let sil_id = "SIL_1";
        let sil = if let Some(sil) = self.doc.data_collection.analysis_data.spectrum_identification_list.iter_mut().find(|l| l.id == sil_id) {
            sil
        } else {
            self.doc.data_collection.analysis_data.spectrum_identification_list.push(SpectrumIdentificationListType {
                id: sil_id.to_string(),
                name: None,
                num_sequences_searched: None,
                content: Vec::new(),
            });
            self.doc.data_collection.analysis_data.spectrum_identification_list.last_mut().unwrap()
        };

        // Find or create the SpectrumIdentificationResult for result_id
        let sir = if let Some(sir) = sil.content.iter_mut().find(|c| match c { SpectrumIdentificationListTypeContent::SpectrumIdentificationResult(r) => r.id == result_id, _ => false }) {
            match sir { SpectrumIdentificationListTypeContent::SpectrumIdentificationResult(r) => r, _ => unreachable!() }
        } else {
            let new_sir = SpectrumIdentificationResultType {
                id: result_id.to_string(),
                name: None,
                spectrum_id: result_id.to_string(),
                spectra_data_ref: spectra_data_ref.to_string(),
                content: Vec::new(),
            };
            sil.content.push(SpectrumIdentificationListTypeContent::SpectrumIdentificationResult(new_sir));
            match sil.content.last_mut().unwrap() { SpectrumIdentificationListTypeContent::SpectrumIdentificationResult(r) => r, _ => unreachable!() }
        };

        sir.content.push(SpectrumIdentificationResultTypeContent::SpectrumIdentificationItem(sii));
    }

    pub fn add_spectra_data(&mut self, id: &str, location: &str) {
        self.doc.data_collection.inputs.spectra_data.push(SpectraDataType {
            id: id.to_string(),
            name: None,
            location: location.to_string(),
            external_format_documentation: None,
            file_format: FileFormatType {
                cv_param: CvParamType {
                    name: "mzML format".to_string(),
                    accession: "MS:1000584".to_string(),
                    cv_ref: "PSI-MS".to_string(),
                    value: None,
                    unit_accession: None,
                    unit_name: None,
                    unit_cv_ref: None,
                }
            },
            spectrum_id_format: SpectrumIdFormatType {
                 cv_param: CvParamType {
                    name: "mzML unique identifier".to_string(),
                    accession: "MS:1001530".to_string(),
                    cv_ref: "PSI-MS".to_string(),
                    value: None,
                    unit_accession: None,
                    unit_name: None,
                    unit_cv_ref: None,
                }
            },
        });
    }

    pub fn add_search_database(&mut self, id: &str, name: &str) {
        self.doc.data_collection.inputs.search_database.push(SearchDatabaseType {
            id: id.to_string(),
            name: Some(name.to_string()),
            location: String::new(),
            version: None,
            release_date: None,
            num_database_sequences: None,
            num_residues: None,
            external_format_documentation: None,
            file_format: FileFormatType {
                cv_param: CvParamType {
                    name: "FASTA format".to_string(),
                    accession: "MS:1001348".to_string(),
                    cv_ref: "PSI-MS".to_string(),
                    value: None,
                    unit_accession: None,
                    unit_name: None,
                    unit_cv_ref: None,
                }
            },
            database_name: ParamType::CvParam(CvParamType {
                name: name.to_string(),
                accession: "MS:1001349".to_string(),
                cv_ref: "PSI-MS".to_string(),
                value: None,
                unit_accession: None,
                unit_name: None,
                unit_cv_ref: None,
            }),
            cv_param: Vec::new(),
        });
    }

    pub fn add_software(&mut self, id: &str, name: &str, version: &str) {
        if let Some(list) = &mut self.doc.analysis_software_list {
            list.analysis_software.push(AnalysisSoftwareType {
                id: id.to_string(),
                name: Some(name.to_string()),
                version: Some(version.to_string()),
                uri: None,
                contact_role: None,
                software_name: ParamType::CvParam(CvParamType {
                    name: name.to_string(),
                    accession: "MS:1002511".to_string(), // Placeholder for now
                    cv_ref: "PSI-MS".to_string(),
                    value: None,
                    unit_accession: None,
                    unit_name: None,
                    unit_cv_ref: None,
                }),
                customizations: None,
            });
        }
    }

    pub fn add_protocol(&mut self, id: &str, software_ref: &str) {
        self.doc.analysis_protocol_collection.spectrum_identification_protocol.push(SpectrumIdentificationProtocolType {
            id: id.to_string(),
            name: None,
            analysis_software_ref: software_ref.to_string(),
            search_type: ParamType::CvParam(CvParamType {
                name: "ms-ms search".to_string(),
                accession: "MS:1001083".to_string(),
                cv_ref: "PSI-MS".to_string(),
                value: None,
                unit_accession: None,
                unit_name: None,
                unit_cv_ref: None,
            }),
            additional_search_params: None,
            modification_params: None,
            enzymes: None,
            mass_table: Vec::new(),
            fragment_tolerance: None,
            parent_tolerance: None,
            threshold: ParamListType {
                content: vec![ParamListTypeContent::CvParam(CvParamType {
                    name: "no threshold".to_string(),
                    accession: "MS:1001494".to_string(),
                    cv_ref: "PSI-MS".to_string(),
                    value: None,
                    unit_accession: None,
                    unit_name: None,
                    unit_cv_ref: None,
                })],
            },
            database_filters: None,
            database_translation: None,
        });
    }

    pub fn add_analysis(&mut self, id: &str, protocol_ref: &str, list_ref: &str, spectra_refs: Vec<String>, db_refs: Vec<String>) {
        self.doc.analysis_collection.spectrum_identification.push(SpectrumIdentificationType {
            id: id.to_string(),
            name: None,
            spectrum_identification_protocol_ref: protocol_ref.to_string(),
            spectrum_identification_list_ref: list_ref.to_string(),
            activity_date: None,
            input_spectra: spectra_refs.into_iter().map(|r| InputSpectraType { spectra_data_ref: r }).collect(),
            search_database_ref: db_refs.into_iter().map(|r| SearchDatabaseRefType { search_database_ref: r }).collect(),
        });
    }

    pub fn serialize(self) -> Result<String, String> {
        let mut writer = Writer::new_with_indent(std::io::Cursor::new(Vec::new()), b' ', 2);
        let mut serializer = self.doc.serializer(Some("MzIdentML"), true)
            .map_err(|e| format!("Serialization error: {:?}", e))?;
        
        for event in &mut serializer {
            let event = event.map_err(|e| format!("XML event error: {:?}", e))?;
            writer.write_event(event).map_err(|e| format!("XML write error: {:?}", e))?;
        }
        
        String::from_utf8(writer.into_inner().into_inner())
            .map_err(|e| format!("UTF8 error: {:?}", e))
    }
}

#[pyfunction]
pub fn write_mzidentml(
    csms: PyDataFrame,
    prot_seqs: PyDataFrame,
    spectra: PyDataFrame,
    _cvs: Bound<'_, PyDict>,
) -> PyResult<String> {
    let mut factory = MzIdentMLFactory::new("mzidentml_export".to_string());
    
    let csms_df = csms.as_ref();
    let prot_df = prot_seqs.as_ref();

    // 1. Setup metadata
    factory.add_software("AS_1", "mzidentml-polars", "0.1.0");
    factory.add_search_database("SearchDB_1", "Target Database");
    factory.add_protocol("SIP_1", "AS_1");
    
    // Process SpectraData
    let spec_df = spectra.as_ref();
    let spec_ids_col = spec_df.column("spectrum_id").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let spec_paths_col = spec_df.column("file_path").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    
    let mut spectra_data_ids = Vec::new();
    let mut path_to_sd_id = HashMap::new();

    for i in 0..spec_df.height() {
        if let Some(path) = spec_paths_col.get(i) {
            if !path_to_sd_id.contains_key(path) {
                let sd_id = format!("SD_{}", path_to_sd_id.len() + 1);
                factory.add_spectra_data(&sd_id, path);
                path_to_sd_id.insert(path.to_string(), sd_id.clone());
                spectra_data_ids.push(sd_id);
            }
        }
    }
    
    if spectra_data_ids.is_empty() {
        // Fallback SD if none provided
        factory.add_spectra_data("SD_1", "unknown_data");
        spectra_data_ids.push("SD_1".to_string());
    }

    factory.add_analysis("SI_1", "SIP_1", "SIL_1", spectra_data_ids, vec!["SearchDB_1".to_string()]);

    // 2. Process Protein Sequences
    let prot_ids = prot_df.column("protein_id").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let prot_accs = prot_df.column("accession").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let prot_seqs_col = prot_df.column("sequence").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    
    for i in 0..prot_df.height() {
        if let (Some(id), Some(acc), Some(seq)) = (prot_ids.get(i), prot_accs.get(i), prot_seqs_col.get(i)) {
            factory.add_db_sequence(id, acc, seq, "SearchDB_1");
        }
    }

    // 2. Process CSMs
    let c_spec_id = csms_df.column("spectrum_id").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_pep1 = csms_df.column("peptide1_seq").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_prot1 = csms_df.column("protein1_id").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_start1 = csms_df.column("peptide1_start").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.u32().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_end1 = csms_df.column("peptide1_end").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.u32().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_charge = csms_df.column("charge").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.i32().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_rank = csms_df.column("rank").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.u32().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    
    let is_xl = csms_df.column("is_crosslink").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.bool().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_pep2 = csms_df.column("peptide2_seq").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_prot2 = csms_df.column("protein2_id").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_start2 = csms_df.column("peptide2_start").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.u32().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_end2 = csms_df.column("peptide2_end").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.u32().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;

    // Map spectrum IDs to their SpectraData reference
    let mut spec_id_to_sd_id = HashMap::new();
    for i in 0..spec_df.height() {
        if let (Some(sid), Some(path)) = (spec_ids_col.get(i), spec_paths_col.get(i)) {
            if let Some(sd_id) = path_to_sd_id.get(path) {
                spec_id_to_sd_id.insert(sid, sd_id);
            }
        }
    }

    for i in 0..csms_df.height() {
        let spec_id = c_spec_id.get(i).unwrap();
        let pep1_ref = factory.add_peptide(c_pep1.get(i).unwrap());
        let prot1_id = c_prot1.get(i).unwrap();
        let dbseq1_ref = format!("dbseq_{}", prot1_id);
        
        let ev1_id = factory.add_peptide_evidence(
            &pep1_ref, 
            &dbseq1_ref, 
            c_start1.get(i), 
            c_end1.get(i), 
            false
        );

        let mut sii = SpectrumIdentificationItemType {
            id: format!("SII_{}_{}", spec_id, i),
            name: None,
            charge_state: c_charge.get(i).unwrap_or(2),
            experimental_mass_to_charge: 0.0, 
            calculated_mass_to_charge: None,
            calculated_pi: None,
            peptide_ref: pep1_ref,
            rank: c_rank.get(i).unwrap_or(1) as i32,
            pass_threshold: true,
            mass_table_ref: None,
            sample_ref: None,
            content: vec![SpectrumIdentificationItemTypeContent::PeptideEvidenceRef(PeptideEvidenceRefType {
                peptide_evidence_ref: ev1_id,
            })],
        };

        if is_xl.get(i).unwrap_or(false) {
            let pep2_ref = factory.add_peptide(c_pep2.get(i).unwrap());
            let prot2_id = c_prot2.get(i).unwrap();
            let dbseq2_ref = format!("dbseq_{}", prot2_id);
            
            let ev2_id = factory.add_peptide_evidence(
                &pep2_ref, 
                &dbseq2_ref, 
                c_start2.get(i), 
                c_end2.get(i), 
                false
            );
            
            sii.content.push(SpectrumIdentificationItemTypeContent::PeptideEvidenceRef(PeptideEvidenceRefType {
                peptide_evidence_ref: ev2_id,
            }));

            // Add crosslink CV param
            sii.content.push(SpectrumIdentificationItemTypeContent::CvParam(CvParamType {
                name: "cross-link spectrum match".to_string(),
                accession: "MS:1002511".to_string(),
                cv_ref: "PSI-MS".to_string(),
                value: None,
                unit_accession: None,
                unit_name: None,
                unit_cv_ref: None,
            }));
        }

        let sd_ref = spec_id_to_sd_id.get(spec_id).map(|s| s.as_str()).unwrap_or("SD_1");
        factory.add_sii(spec_id, sii, sd_ref);
    }

    factory.serialize().map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factory_creation() {
        let factory = MzIdentMLFactory::new("test_doc".to_string());
        assert_eq!(factory.doc.id, "test_doc");
        assert_eq!(factory.doc.version, "1.3.0");
        assert!(factory.doc.cv_list.cv.len() >= 4);
    }

    #[test]
    fn test_add_peptide_deduplication() {
        let mut factory = MzIdentMLFactory::new("test_doc".to_string());
        let id1 = factory.add_peptide("PEPTIDE");
        let id2 = factory.add_peptide("PEPTIDE");
        assert_eq!(id1, id2);
        assert_eq!(id1, "pep_0");
        
        let sc = factory.doc.sequence_collection.as_ref().unwrap();
        assert_eq!(sc.peptide.len(), 1);
    }

    #[test]
    fn test_add_db_sequence() {
        let mut factory = MzIdentMLFactory::new("test_doc".to_string());
        let id = factory.add_db_sequence("P12345", "ACC123", "MAGA", "DB1");
        assert_eq!(id, "dbseq_P12345");
        
        let sc = factory.doc.sequence_collection.as_ref().unwrap();
        assert_eq!(sc.db_sequence.len(), 1);
        assert_eq!(sc.db_sequence[0].accession, "ACC123");
    }

    #[test]
    fn test_serialization_basic() {
        let mut factory = MzIdentMLFactory::new("test_doc".to_string());
        factory.add_peptide("PEPTIDE");
        factory.add_db_sequence("P12", "ACC", "M", "DB");
        
        let xml = factory.serialize().unwrap();
        // println!("XML OUTPUT:\n{}", xml);
        assert!(xml.contains("test_doc"));
        assert!(xml.contains("psi-pi:Peptide"));
        assert!(xml.contains("psi-pi:id=\"pep_0\""));
        assert!(xml.contains("psi-pi:DBSequence"));
        assert!(xml.contains("psi-pi:id=\"dbseq_P12\""));
    }
}
