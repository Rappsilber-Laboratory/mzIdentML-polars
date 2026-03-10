use pyo3::prelude::*;
use pyo3::types::{PyDict};
use pyo3_polars::PyDataFrame;
use std::collections::HashMap;
use crate::mzidentml::psi_pi::*;
use polars::prelude::*;
use xsd_parser::quick_xml::WithSerializer;
use xsd_parser::quick_xml::Writer;
use flate2::write::GzEncoder;
use flate2::Compression;


mod cv_data {
    include!(concat!(env!("OUT_DIR"), "/cv_data.rs"));
}

fn get_string_list(val: AnyValue) -> Vec<String> {
    match val {
        AnyValue::List(s) => s.str().unwrap().into_iter().flatten().map(|s| s.to_string()).collect(),
        AnyValue::String(s) => vec![s.to_string()],
        AnyValue::Null => Vec::new(),
        _ => Vec::new(),
    }
}

fn get_u32_list(val: AnyValue) -> Vec<u32> {
    match val {
        AnyValue::List(s) => s.u32().unwrap().into_iter().flatten().collect(),
        AnyValue::UInt32(u) => vec![u],
        AnyValue::Null => Vec::new(),
        _ => Vec::new(),
    }
}

fn derive_spectra_data_format(location: &str) -> (CvParamType, CvParamType) {
    let mut lower = location.to_lowercase();
    if lower.ends_with(".gz") {
        lower = lower[..lower.len() - 3].to_string();
    } else if lower.ends_with(".gzip") {
        lower = lower[..lower.len() - 5].to_string();
    }
    
    if lower.ends_with(".mzml") {
        (
            CvParamType {
                name: "mzML format".to_string(),
                accession: "MS:1000584".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
            CvParamType {
                name: "mzML unique identifier".to_string(),
                accession: "MS:1001530".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
        )
    } else if lower.ends_with(".mgf") {
        (
            CvParamType {
                name: "MGF format".to_string(),
                accession: "MS:1001062".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
            CvParamType {
                name: "spectrum title".to_string(),
                accession: "MS:1000796".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
        )
    } else if lower.ends_with(".raw") {
        (
            CvParamType {
                name: "Thermo RAW format".to_string(),
                accession: "MS:1000563".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
            CvParamType {
                name: "Thermo nativeID format".to_string(),
                accession: "MS:1000768".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
        )
    } else if lower.ends_with(".mzxml") {
        (
            CvParamType {
                name: "ISB mzXML format".to_string(),
                accession: "MS:1000566".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
            CvParamType {
                name: "scan number only nativeID format".to_string(),
                accession: "MS:1000776".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
        )
    } else if lower.ends_with(".d") {
        (
            CvParamType {
                name: "Bruker BAF format".to_string(),
                accession: "MS:1000815".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
            CvParamType {
                name: "Bruker BAF nativeID format".to_string(),
                accession: "MS:1000772".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
        )
    } else if lower.ends_with(".wiff") {
        (
            CvParamType {
                name: "ABI WIFF format".to_string(),
                accession: "MS:1000562".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
            CvParamType {
                name: "WIFF nativeID format".to_string(),
                accession: "MS:1000770".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
        )
    } else if lower.ends_with(".pkl") {
        (
            CvParamType {
                name: "Micromass PKL format".to_string(),
                accession: "MS:1000565".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
            CvParamType {
                name: "Waters nativeID format".to_string(),
                accession: "MS:1000769".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
        )
    } else {
        // Default to mzML if unknown
        (
            CvParamType {
                name: "mzML format".to_string(),
                accession: "MS:1000584".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
            CvParamType {
                name: "mzML unique identifier".to_string(),
                accession: "MS:1001530".to_string(),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            },
        )
    }
}

/// Factory to manage the construction of the MzIdentML XML tree.
pub struct MzIdentMLFactory {
    pub doc: MzIdentMlType,
    peptide_map: HashMap<String, String>, // ProForma -> Peptide ID
    db_seq_map: HashMap<String, String>,  // Protein ID -> DBSequence ID
    pep_evidence_map: HashMap<(String, String), String>, // (Peptide ID, Protein ID) -> Evidence ID
    cv_map: HashMap<String, String>,      // CV ID -> URI
    decoy_map: HashMap<String, bool>,     // Protein ID -> is_decoy
    pub mod_name_overrides: HashMap<String, String>, // Accession -> Name
}

impl MzIdentMLFactory {
    pub fn new(id: String) -> Self {
        let now = chrono::Utc::now().to_rfc3339_opts(chrono::SecondsFormat::Secs, true);
        let doc = MzIdentMlType {
            id,
            name: None,
            creation_date: Some(now),
            version: "1.2.0".to_string(),
            cv_list: CvListType { cv: Vec::new() },
            cv_param: vec![CvParamType {
                name: "mzIdentML crosslinking extension document version".to_string(),
                accession: "MS:1003385".to_string(),
                cv_ref: "PSI-MS".to_string(),
                value: Some("1.0.0".to_string()),
                unit_accession: None,
                unit_name: None,
                unit_cv_ref: None,
            }],
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
            decoy_map: HashMap::new(),
            mod_name_overrides: HashMap::new(),
        };

        factory.add_cv("PSI-MS", "PSI-MS", "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo");
        factory.add_cv("UNIMOD", "UNIMOD", "http://www.unimod.org/obo/unimod.obo");
        factory.add_cv("UO", "Unit Ontology", "http://purl.obolibrary.org/obo/uo.obo");
        factory.add_cv("XLMOD", "PSI-XLMOD", "https://raw.githubusercontent.com/HUPO-PSI/xlmod-CV/refs/heads/main/XLMOD.obo");
        factory.add_cv("PSI-MOD", "PSI-MOD", "https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/master/PSI-MOD.obo");

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

    pub fn add_peptide(&mut self, proforma: &str, linkage_mods: Vec<ModificationType>) -> String {
        let (clean_seq, mut mods) = self.parse_proforma(proforma);
        
        // Add linkage modifications
        mods.extend(linkage_mods);
        // Sort mods by location for consistent keying
        mods.sort_by_key(|m| m.location.unwrap_or(0));

        // Create a unique key for this peptide variant
        let mut key = clean_seq.clone();
        for m in &mods {
            key.push('|');
            key.push_str(&m.location.unwrap_or(0).to_string());
            for cv in &m.cv_param {
                key.push(':');
                key.push_str(&cv.accession);
                if let Some(v) = &cv.value {
                    key.push('=');
                    key.push_str(v);
                }
            }
        }

        if let Some(id) = self.peptide_map.get(&key) {
            return id.clone();
        }

        let id = format!("pep_{}", self.peptide_map.len());
        let mut peptide_struct = PeptideType {
            id: id.clone(),
            name: None,
            content: vec![PeptideTypeContent::PeptideSequence(clean_seq)],
        };

        for m in mods {
            peptide_struct.content.push(PeptideTypeContent::Modification(m));
        }

        if let Some(sc) = &mut self.doc.sequence_collection {
            sc.peptide.push(peptide_struct);
        }
        self.peptide_map.insert(key, id.clone());
        id
    }

    pub fn add_db_sequence(&mut self, protein_id: &str, accession: &str, sequence: &str, db_ref: &str, is_decoy: bool, protein_name: Option<String>) -> String {
        if let Some(id) = self.db_seq_map.get(protein_id) {
            return id.clone();
        }

        let id = format!("dbseq_{}", protein_id);
        let mut content = vec![DbSequenceTypeContent::Seq(sequence.to_string())];
        let mut final_name = protein_name;
        if is_decoy {
            if final_name.is_none() {
                final_name = Some("decoy".to_string());
            }
            content.push(DbSequenceTypeContent::CvParam(CvParamType {
                name: "protein description".to_string(),
                accession: "MS:1001088".to_string(),
                cv_ref: "PSI-MS".to_string(),
                value: Some("decoy".to_string()),
                ..Default::default()
            }));
        }

        let db_seq = DbSequenceType {
            id: id.clone(),
            name: final_name,
            length: Some(sequence.len() as i32),
            search_database_ref: db_ref.to_string(),
            accession: accession.to_string(),
            content,
        };

        if let Some(sc) = &mut self.doc.sequence_collection {
            sc.db_sequence.push(db_seq);
        }
        self.db_seq_map.insert(protein_id.to_string(), id.clone());
        self.decoy_map.insert(protein_id.to_string(), is_decoy);
        id
    }

    pub fn add_peptide_evidence(&mut self, pep_ref: &str, db_ref: &str, start: Option<u32>, end: Option<u32>, is_decoy: bool) -> String {
        let key = (pep_ref.to_string(), db_ref.to_string());
        if let Some(id) = self.pep_evidence_map.get(&key) {
            return id.clone();
        }

        let id = format!("pepevid_{}_{}", pep_ref, db_ref);
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

    pub fn add_spectrum_identification_result(&mut self, sd_ref: &str, spec_id: &str, items: Vec<SpectrumIdentificationItemType>, sir_params: Vec<CvParamType>) {
        let sir_id = format!("SIR_{}_{}", sd_ref, spec_id).replace("=", "_").replace(":", "_").replace(" ", "_");
        self.add_spectrum_identification_result_with_id(&sir_id, sd_ref, spec_id, items, sir_params);
    }

    pub fn add_spectrum_identification_result_with_id(&mut self, sir_id: &str, sd_ref: &str, spec_id: &str, items: Vec<SpectrumIdentificationItemType>, sir_params: Vec<CvParamType>) {
        let sil = match self.doc.data_collection.analysis_data.spectrum_identification_list.first_mut() {
            Some(s) => s,
            None => return,
        };
        let mut content = Vec::new();
        for item in items {
            content.push(SpectrumIdentificationResultTypeContent::SpectrumIdentificationItem(item));
        }
        for p in sir_params {
            content.push(SpectrumIdentificationResultTypeContent::CvParam(p));
        }
        
        let sir = SpectrumIdentificationResultType {
            id: sir_id.to_string(),
            spectrum_id: spec_id.to_string(),
            spectra_data_ref: sd_ref.to_string(),
            content,
            ..Default::default()
        };
        sil.content.push(SpectrumIdentificationListTypeContent::SpectrumIdentificationResult(sir));
    }

    pub fn add_spectra_data(&mut self, id: &str, location: &str) {
        let (file_format, id_format) = derive_spectra_data_format(location);
        self.doc.data_collection.inputs.spectra_data.push(SpectraDataType {
            id: id.to_string(),
            name: Some(id.to_string()),
            location: location.to_string(),
            external_format_documentation: None,
            file_format: FileFormatType { cv_param: file_format },
            spectrum_id_format: SpectrumIdFormatType {
                 cv_param: id_format
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
                name: "database name".to_string(),
                accession: "MS:1001013".to_string(),
                cv_ref: "PSI-MS".to_string(),
                value: Some(name.to_string()),
                ..Default::default()
            }),
            cv_param: vec![CvParamType {
                name: "decoy DB type reverse".to_string(),
                accession: "MS:1001195".to_string(),
                cv_ref: "PSI-MS".to_string(),
                value: Some("decoy database".to_string()),
                ..Default::default()
            }],
        });
    }

    pub fn set_author(&mut self, author_name: &str) {
        let person_id = "PERSON_AUTHOR".to_string();
        
        let person = PersonType {
            id: person_id.clone(),
            name: Some(author_name.to_string()),
            last_name: None,
            first_name: None,
            mid_initials: None,
            content: Vec::new(),
        };

        if self.doc.audit_collection.is_none() {
            self.doc.audit_collection = Some(AuditCollectionType { content: Vec::new() });
        }
        
        if let Some(ac) = &mut self.doc.audit_collection {
            ac.content.push(AuditCollectionTypeContent::Person(person));
        }

        self.doc.provider = Some(ProviderType {
            id: "PROVIDER".to_string(),
            name: None,
            analysis_software_ref: None,
            contact_role: Some(ContactRoleType {
                contact_ref: person_id,
                role: RoleType {
                    cv_param: CvParamType {
                        name: "researcher".to_string(),
                        accession: "MS:1001271".to_string(),
                        cv_ref: "PSI-MS".to_string(),
                        ..Default::default()
                    }
                }
            })
        });
    }

    pub fn add_software(&mut self, id: &str, name: &str, version: &str) {
        if let Some(list) = &mut self.doc.analysis_software_list {
            let software_name = match name.to_lowercase().as_str() {
                "xi" | "xisearch" | "mzidentml-polars" => ParamType::CvParam(CvParamType {
                    name: "XI Search".to_string(),
                    accession: "MS:1001405".to_string(),
                    cv_ref: "PSI-MS".to_string(),
                    ..Default::default()
                }),
                "xifdr" => ParamType::CvParam(CvParamType {
                    name: "xiFDR".to_string(),
                    accession: "MS:1002543".to_string(),
                    cv_ref: "PSI-MS".to_string(),
                    ..Default::default()
                }),
                _ => ParamType::UserParam(UserParamType {
                    name: name.to_string(),
                    type_: None,
                    unit_accession: None,
                    unit_name: None,
                    unit_cv_ref: None,
                    value: None,
                }),
            };

            list.analysis_software.push(AnalysisSoftwareType {
                id: id.to_string(),
                name: Some(name.to_string()),
                version: Some(version.to_string()),
                uri: None,
                contact_role: None,
                software_name,
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
            additional_search_params: Some(ParamListType {
                content: vec![
                    ParamListTypeContent::CvParam(CvParamType {
                        name: "parent mass type mono".to_string(),
                        accession: "MS:1001211".to_string(),
                        cv_ref: "PSI-MS".to_string(),
                        value: None,
                        unit_accession: None,
                        unit_name: None,
                        unit_cv_ref: None,
                    }),
                    ParamListTypeContent::CvParam(CvParamType {
                        name: "fragment mass type mono".to_string(),
                        accession: "MS:1001256".to_string(),
                        cv_ref: "PSI-MS".to_string(),
                        value: None,
                        unit_accession: None,
                        unit_name: None,
                        unit_cv_ref: None,
                    }),
                    // Mandatory for crosslinking extension
                    ParamListTypeContent::CvParam(CvParamType {
                        name: "cross-linking search".to_string(),
                        accession: "MS:1002494".to_string(),
                        cv_ref: "PSI-MS".to_string(),
                        ..Default::default()
                    }),
                ],
            }),
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

        // Initialize the SIL in AnalysisData
        self.doc.data_collection.analysis_data.spectrum_identification_list.push(SpectrumIdentificationListType {
            id: list_ref.to_string(),
            ..Default::default()
        });
    }

    pub fn set_tolerances(&mut self, protocol_index: usize, parent_plus: f64, parent_minus: f64, frag_plus: f64, frag_minus: f64, is_ppm: bool) {
        if let Some(protocol) = self.doc.analysis_protocol_collection.spectrum_identification_protocol.get_mut(protocol_index) {
            let unit = if is_ppm { ("UO:0000169", "parts per million") } else { ("UO:0000221", "dalton") };
            
            protocol.parent_tolerance = Some(ToleranceType {
                cv_param: vec![
                    CvParamType {
                        name: "search tolerance plus value".to_string(),
                        accession: "MS:1001412".to_string(),
                        cv_ref: "PSI-MS".to_string(),
                        value: Some(parent_plus.to_string()),
                        unit_accession: Some(unit.0.to_string()),
                        unit_name: Some(unit.1.to_string()),
                        unit_cv_ref: Some("UO".to_string()),
                    },
                    CvParamType {
                        name: "search tolerance minus value".to_string(),
                        accession: "MS:1001413".to_string(),
                        cv_ref: "PSI-MS".to_string(),
                        value: Some(parent_minus.to_string()),
                        unit_accession: Some(unit.0.to_string()),
                        unit_name: Some(unit.1.to_string()),
                        unit_cv_ref: Some("UO".to_string()),
                    }
                ],
                ..Default::default()
            });

            protocol.fragment_tolerance = Some(ToleranceType {
                cv_param: vec![
                    CvParamType {
                        name: "search tolerance plus value".to_string(),
                        accession: "MS:1001412".to_string(),
                        cv_ref: "PSI-MS".to_string(),
                        value: Some(frag_plus.to_string()),
                        unit_accession: Some(unit.0.to_string()),
                        unit_name: Some(unit.1.to_string()),
                        unit_cv_ref: Some("UO".to_string()),
                    },
                    CvParamType {
                        name: "search tolerance minus value".to_string(),
                        accession: "MS:1001413".to_string(),
                        cv_ref: "PSI-MS".to_string(),
                        value: Some(frag_minus.to_string()),
                        unit_accession: Some(unit.0.to_string()),
                        unit_name: Some(unit.1.to_string()),
                        unit_cv_ref: Some("UO".to_string()),
                    }
                ],
                ..Default::default()
            });
        }
    }

    pub fn add_search_modification(&mut self, protocol_index: usize, fixed: bool, mass_delta: f32, residues: &str, name: &str, accession: &str) {
        if let Some(protocol) = self.doc.analysis_protocol_collection.spectrum_identification_protocol.get_mut(protocol_index) {
            if protocol.modification_params.is_none() {
                protocol.modification_params = Some(ModificationParamsType::default());
            }
            if let Some(mp) = &mut protocol.modification_params {
                mp.search_modification.push(SearchModificationType {
                    fixed_mod: fixed,
                    mass_delta,
                    residues: ListOfCharsOrAnyType::EntitiesType(crate::mzidentml::xs::EntitiesType(residues.split_whitespace().map(|s| s.to_string()).collect())),
                    cv_param: vec![CvParamType {
                        name: name.to_string(),
                        accession: accession.to_string(),
                        cv_ref: get_cv_ref(accession),
                        ..Default::default()
                    }],
                    ..Default::default()
                });
            }
        }
    }

    pub fn add_search_param(&mut self, protocol_index: usize, name: &str, accession: &str, value: Option<&str>) {
        if let Some(protocol) = self.doc.analysis_protocol_collection.spectrum_identification_protocol.get_mut(protocol_index) {
            if let Some(asp) = &mut protocol.additional_search_params {
                asp.content.push(ParamListTypeContent::CvParam(CvParamType {
                    name: name.to_string(),
                    accession: accession.to_string(),
                    cv_ref: get_cv_ref(accession),
                    value: value.map(|v| v.to_string()),
                    ..Default::default()
                }));
            }
        }
    }

    pub fn add_enzyme(&mut self, protocol_index: usize, id: &str, name: &str, accession: &str) {
        if let Some(protocol) = self.doc.analysis_protocol_collection.spectrum_identification_protocol.get_mut(protocol_index) {
            if protocol.enzymes.is_none() {
                protocol.enzymes = Some(EnzymesType::default());
            }
            if let Some(enz) = &mut protocol.enzymes {
                enz.enzyme.push(EnzymeType {
                    id: id.to_string(),
                    enzyme_name: Some(ParamListType {
                        content: vec![ParamListTypeContent::CvParam(CvParamType {
                            name: name.to_string(),
                        accession: accession.to_string(),
                        cv_ref: get_cv_ref(accession),
                            ..Default::default()
                        })]
                    }),
                    ..Default::default()
                });
            }
        }
    }

    fn serialize_internal<W: std::io::Write>(&self, writer: &mut Writer<W>) -> Result<(), String> {
        // Add XML declaration
        writer.write_event(xsd_parser::quick_xml::Event::Decl(xsd_parser::quick_xml::BytesDecl::new("1.0", Some("UTF-8"), None)))
            .map_err(|e| format!("XML write error: {:?}", e))?;

        let mut serializer = self.doc.serializer(Some("MzIdentML"), true)
            .map_err(|e| format!("Serialization error: {:?}", e))?;
        
        for event in &mut serializer {
            let event = event.map_err(|e| format!("XML event error: {:?}", e))?;
            writer.write_event(event).map_err(|e| format!("XML write error: {:?}", e))?;
        }
        Ok(())
    }

    pub fn serialize(self) -> Result<String, String> {
        let mut writer = Writer::new_with_indent(std::io::Cursor::new(Vec::new()), b' ', 2);
        self.serialize_internal(&mut writer)?;
        
        String::from_utf8(writer.into_inner().into_inner())
            .map_err(|e| format!("UTF8 error: {:?}", e))
    }

    pub fn serialize_to_file(self, path: &str) -> Result<(), String> {
        let file = std::fs::File::create(path).map_err(|e| format!("File creation error for '{}': {:?}", path, e))?;
        
        if path.ends_with(".gz") {
            let encoder = GzEncoder::new(file, Compression::default());
            let mut writer = Writer::new_with_indent(encoder, b' ', 2);
            self.serialize_internal(&mut writer)
        } else {
            let mut writer = Writer::new_with_indent(file, b' ', 2);
            self.serialize_internal(&mut writer)
        }
    }
}

pub fn prepare_factory(
    csms: PyDataFrame,
    prot_seqs: PyDataFrame,
    spectra: PyDataFrame,
    metadata: Bound<'_, PyDict>,
) -> PyResult<MzIdentMLFactory> {
    let mut factory = MzIdentMLFactory::new("mzidentml_export".to_string());
    
    let csms_df = csms.as_ref();
    let prot_df = prot_seqs.as_ref();

    // 1. Setup metadata
    let sw_name = metadata.get_item("software_name").ok().flatten().and_then(|v| v.extract::<String>().ok()).unwrap_or_else(|| "mzidentml-polars".to_string());
    let sw_version = metadata.get_item("software_version").ok().flatten().and_then(|v| v.extract::<String>().ok()).unwrap_or_else(|| "0.1.0".to_string());
    factory.add_software("AS_1", &sw_name, &sw_version);
    
    if let Some(author) = metadata.get_item("author").ok().flatten().and_then(|v| v.extract::<String>().ok()) {
        factory.set_author(&author);
    }
    factory.add_search_database("SearchDB_1", "Target Database");
    factory.add_protocol("SIP_1", "AS_1");

    // Add tolerances from metadata
    let p_plus = metadata.get_item("parent_plus").ok().flatten().and_then(|v| v.extract::<f64>().ok()).unwrap_or(5.0);
    let p_minus = metadata.get_item("parent_minus").ok().flatten().and_then(|v| v.extract::<f64>().ok()).unwrap_or(5.0);
    let f_plus = metadata.get_item("frag_plus").ok().flatten().and_then(|v| v.extract::<f64>().ok()).unwrap_or(10.0);
    let f_minus = metadata.get_item("frag_minus").ok().flatten().and_then(|v| v.extract::<f64>().ok()).unwrap_or(10.0);
    let is_ppm = metadata.get_item("is_ppm").ok().flatten().and_then(|v| v.extract::<bool>().ok()).unwrap_or(true);
    factory.set_tolerances(0, p_plus, p_minus, f_plus, f_minus, is_ppm);
    
    // Add enzymes from metadata
    if let Some(enzymes) = metadata.get_item("enzymes").ok().flatten() {
        if let Ok(enz_list) = enzymes.extract::<Vec<Bound<'_, PyDict>>>() {
            for (idx, enz_dict) in enz_list.iter().enumerate() {
                let name = enz_dict.get_item("name").ok().flatten().and_then(|v| v.extract::<String>().ok()).unwrap_or_else(|| "unknown enzyme".to_string());
                let acc = enz_dict.get_item("accession").ok().flatten().and_then(|v| v.extract::<String>().ok()).unwrap_or_else(|| "MS:1001045".to_string());
                let enz_id = format!("ENZ_{}", idx + 1);
                factory.add_enzyme(0, &enz_id, &name, &acc);
            }
        }
    }

    // Add modifications from metadata
    if let Some(mods) = metadata.get_item("modifications").ok().flatten() {
        if let Ok(mod_list) = mods.extract::<Vec<Bound<'_, PyDict>>>() {
            for mod_dict in mod_list {
                let fixed = mod_dict.get_item("fixed").ok().flatten().and_then(|v| v.extract::<bool>().ok()).unwrap_or(false);
                let mass = mod_dict.get_item("mass").ok().flatten().and_then(|v| v.extract::<f32>().ok()).unwrap_or(0.0);
                let residues = mod_dict.get_item("residues").ok().flatten().and_then(|v| v.extract::<String>().ok()).unwrap_or_else(|| "".to_string());
                let name = mod_dict.get_item("name").ok().flatten().and_then(|v| v.extract::<String>().ok());
                let acc = mod_dict.get_item("accession").ok().flatten().and_then(|v| v.extract::<String>().ok()).unwrap_or_else(|| "MS:1001460".to_string());
                
                let final_name = if let Some(n) = name {
                    factory.mod_name_overrides.insert(acc.clone(), n.clone());
                    n
                } else {
                    cv_data::lookup_mod(&acc).map(|d| d.name.to_string()).unwrap_or_else(|| "unknown modification".to_string())
                };

                factory.add_search_modification(0, fixed, mass, &residues, &final_name, &acc);
            }
        }
    }

    // Add search parameters from metadata
    if let Some(params) = metadata.get_item("search_params").ok().flatten() {
        if let Ok(param_list) = params.extract::<Vec<Bound<'_, PyDict>>>() {
            for p_dict in param_list {
                let name = p_dict.get_item("name").ok().flatten().and_then(|v| v.extract::<String>().ok()).unwrap_or_else(|| "unknown parameter".to_string());
                let acc = p_dict.get_item("accession").ok().flatten().and_then(|v| v.extract::<String>().ok()).unwrap_or_else(|| "MS:1001302".to_string());
                let value = p_dict.get_item("value").ok().flatten().and_then(|v| v.extract::<String>().ok());
                factory.add_search_param(0, &name, &acc, value.as_deref());
            }
        }
    }
    
    // Process SpectraData
    let spec_df = spectra.as_ref();
    let _spec_ids_col = spec_df.column("spectrum_id").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let spec_paths_col = spec_df.column("file_path").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    
    let mut spectra_data_ids = Vec::new();
    let mut path_to_sd_id = HashMap::new();

    for i in 0..spec_df.height() {
        if let Some(path) = spec_paths_col.get(i) {
            let path_key = path.to_lowercase();
            if !path_to_sd_id.contains_key(&path_key) {
                let filename = std::path::Path::new(path)
                    .file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or(path);
                let sd_id = format!("SD_{}", filename);
                factory.add_spectra_data(&sd_id, filename);
                path_to_sd_id.insert(path_key, sd_id.clone());
                spectra_data_ids.push(sd_id);
            }
        }
    }
    
    if spectra_data_ids.is_empty() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("No spectra files provided in the 'spectra' DataFrame."));
    }

    factory.add_analysis("SI_1", "SIP_1", "SIL_1", spectra_data_ids, vec!["SearchDB_1".to_string()]);

    // 2. Process Protein Sequences
    let prot_ids = prot_df.column("protein_id").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let prot_accs = prot_df.column("accession").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let prot_seqs_col = prot_df.column("sequence").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let prot_names_col = prot_df.column("protein_name").ok().and_then(|c| c.str().ok());
    let prot_is_decoy = prot_df.column("is_decoy").ok().and_then(|c| c.bool().ok());
    
    for i in 0..prot_df.height() {
        if let (Some(id), Some(acc), Some(seq)) = (prot_ids.get(i), prot_accs.get(i), prot_seqs_col.get(i)) {
            let is_decoy = prot_is_decoy.as_ref().and_then(|c| c.get(i)).unwrap_or(false);
            let protein_name = prot_names_col.as_ref().and_then(|c| c.get(i)).map(|s| s.to_string());
            factory.add_db_sequence(id, acc, seq, "SearchDB_1", is_decoy, protein_name);
        }
    }

    // 2. Process CSMs
    let c_spec_id = csms_df.column("spectrum_id").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_pep1 = csms_df.column("peptide1_seq").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_prot1 = csms_df.column("protein1_id").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_start1 = csms_df.column("peptide1_start").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_end1 = csms_df.column("peptide1_end").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_charge = csms_df.column("charge").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.i32().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_rank = csms_df.column("rank").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.u32().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    
    let is_xl = csms_df.column("is_crosslink").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.bool().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_pep2 = csms_df.column("peptide2_seq").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_prot2 = csms_df.column("protein2_id").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_start2 = csms_df.column("peptide2_start").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_end2 = csms_df.column("peptide2_end").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;

    let c_link_pos1 = csms_df.column("peptide1_link_pos").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.i32().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let c_link_pos2 = csms_df.column("peptide2_link_pos").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.i32().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    let is_loop_link = csms_df.column("is_looplink").map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?.bool().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;
    
    // file_path in csms_df is REQUIRED for unique linking
    let c_csm_file_path = csms_df.column("file_path")
        .map_err(|_| PyErr::new::<pyo3::exceptions::PyValueError, _>("The 'csms' DataFrame must contain a 'file_path' column to link matches to spectra files."))?
        .str().map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{}", e)))?;

    // Optional but recommended columns
    let c_exp_mz = csms_df.column("experimental_mz").ok().and_then(|c| c.f64().ok());
    let c_score = csms_df.column("score").ok().and_then(|c| c.f64().ok());
    let c_xl_name = csms_df.column("crosslinker_name").ok().and_then(|c| c.str().ok());
    let c_xl_acc = csms_df.column("crosslinker_accession").ok().and_then(|c| c.str().ok());
    let c_xl_mass = csms_df.column("crosslinker_mass").ok().and_then(|c| c.f64().ok());
    let c_calc_mz = csms_df.column("calculated_mz").ok().and_then(|c| c.f64().ok());


    // Group results by (SpectraData Ref, Spectrum ID) to avoid duplicate SIR elements
    let mut grouped_results: HashMap<(String, String), (Vec<SpectrumIdentificationItemType>, Vec<CvParamType>)> = HashMap::new();

    for i in 0..csms_df.height() {
        let spec_id = c_spec_id.get(i).unwrap();
        let path = c_csm_file_path.get(i).ok_or_else(|| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Missing 'file_path' for CSM at row {}", i))
        })?;

        let sd_ref = path_to_sd_id.get(&path.to_lowercase()).ok_or_else(|| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "The file path '{}' for spectrum '{}' was not found in the provided 'spectra' DataFrame. Please ensure all files used in identifications are registered (case-insensitive).",
                path, spec_id
            ))
        })?;

        let xl_group_id = format!("{}", i + 1);
        let calc_mz_val = c_calc_mz.as_ref().and_then(|c| c.get(i));

        let mut sir_params = Vec::new();
        let xl_name = c_xl_name.as_ref().and_then(|c| c.get(i));
        let xl_acc = c_xl_acc.as_ref().and_then(|c| c.get(i));
        
        // Add viewer-specific spectrum identifiers
        let mut path_stripped = path.to_lowercase();
        if path_stripped.ends_with(".gz") {
            path_stripped = path_stripped[..path_stripped.len() - 3].to_string();
        } else if path_stripped.ends_with(".gzip") {
            path_stripped = path_stripped[..path_stripped.len() - 5].to_string();
        }

        let is_mgf = path_stripped.ends_with(".mgf");
        
        let mut extracted_scan = None;
        if spec_id.starts_with("index=") {
            extracted_scan = Some(spec_id[6..].to_string());
        } else if spec_id.starts_with("scan=") {
            extracted_scan = Some(spec_id[5..].to_string());
        } else if let Some(idx) = spec_id.find("scan=") {
            let rest = &spec_id[idx+5..];
            let scan = rest.split_whitespace().next().unwrap_or(rest);
            extracted_scan = Some(scan.to_string());
        } else if spec_id.chars().all(|c| c.is_ascii_digit()) {
            extracted_scan = Some(spec_id.to_string());
        }

        if is_mgf {
            // For MGF, xiVIEW and others often need the title as a param even if it's the spectrumID
            let title = extracted_scan.as_deref().unwrap_or(spec_id);
            sir_params.push(CvParamType {
                name: "spectrum title".to_string(),
                accession: "MS:1000796".to_string(),
                cv_ref: "PSI-MS".to_string(),
                value: Some(title.to_string()),
                ..Default::default()
            });
        }
        
        if let Some(scan) = extracted_scan {
            // peak list scans is highly recommended for all formats to help viewers link IDs
            sir_params.push(CvParamType {
                name: "peak list scans".to_string(),
                accession: "MS:1000797".to_string(),
                cv_ref: "PSI-MS".to_string(),
                value: Some(scan),
                ..Default::default()
            });
        }


        if is_xl.get(i).unwrap_or(false) {
            // CROSS-LINK MATCH
            
            // 1. Add Peptide 1 with Linker DONOR modification
            let mut linkage1 = Vec::new();
            if let Some(pos1) = c_link_pos1.get(i) {
                let mut params = Vec::new();
                let mut xl_name_param = "unknown modification".to_string();
                let mut xl_acc_param = "MS:1001460".to_string();
                let mut xl_mass_val = c_xl_mass.as_ref().and_then(|c| c.get(i)).unwrap_or(0.0);

                if let Some(acc) = xl_acc {
                    xl_acc_param = acc.to_string();
                    if let Some(override_name) = factory.mod_name_overrides.get(acc) {
                        xl_name_param = override_name.clone();
                        if let Some(data) = cv_data::lookup_mod(acc) {
                            if xl_mass_val == 0.0 {
                                xl_mass_val = data.mono_mass;
                            }
                        }
                    } else if let Some(data) = cv_data::lookup_mod(acc) {
                        xl_name_param = data.name.to_string();
                        if xl_mass_val == 0.0 {
                            xl_mass_val = data.mono_mass;
                        }
                    } else if let Some(name) = xl_name {
                        xl_name_param = name.to_string();
                    }
                } else if let Some(name) = xl_name {
                    xl_name_param = name.to_string();
                    if let Some((acc, mass)) = cv_data::lookup_mod_by_name(name) {
                        xl_acc_param = acc.to_string();
                        // Check override for the accession found by name
                        if let Some(override_name) = factory.mod_name_overrides.get(&xl_acc_param) {
                            xl_name_param = override_name.clone();
                        }
                        if xl_mass_val == 0.0 {
                            xl_mass_val = mass;
                        }
                    }
                }

                params.push(CvParamType {
                    name: "unknown modification".to_string(),
                    accession: "MS:1001460".to_string(),
                    cv_ref: "PSI-MS".to_string(),
                    ..Default::default()
                });
                params.push(CvParamType {
                    name: "cross-link donor".to_string(),
                    accession: "MS:1002509".to_string(),
                    cv_ref: "PSI-MS".to_string(),
                    value: Some(xl_group_id.clone()),
                    ..Default::default()
                });
                linkage1.push(ModificationType {
                    location: Some(pos1),
                    cv_param: params,
                    monoisotopic_mass_delta: Some(xl_mass_val),
                    ..Default::default()
                });
            }
            let pep1_id = factory.add_peptide(c_pep1.get(i).unwrap(), linkage1);
            let _pep1_id_plain = pep1_id.replace("ident_", "").replace("pep_", "");

            // 2. Add Peptide 2 with Linker ACCEPTOR modification
            let mut linkage2 = Vec::new();
            if let Some(pos2) = c_link_pos2.get(i) {
                let mut params = Vec::new();
                // For acceptor, we use the same identified name/acc but usually 0 mass
                let mut xl_name_param = "unknown modification".to_string();
                let mut xl_acc_param = "MS:1001460".to_string();

                if let Some(acc) = xl_acc {
                    xl_acc_param = acc.to_string();
                    if let Some(override_name) = factory.mod_name_overrides.get(acc) {
                        xl_name_param = override_name.clone();
                    } else if let Some(data) = cv_data::lookup_mod(acc) {
                        xl_name_param = data.name.to_string();
                    } else if let Some(name) = xl_name {
                        xl_name_param = name.to_string();
                    }
                } else if let Some(name) = xl_name {
                    xl_name_param = name.to_string();
                    if let Some((acc, _)) = cv_data::lookup_mod_by_name(name) {
                        xl_acc_param = acc.to_string();
                        if let Some(override_name) = factory.mod_name_overrides.get(&xl_acc_param) {
                            xl_name_param = override_name.clone();
                        }
                    }
                }

                params.push(CvParamType {
                    name: "unknown modification".to_string(),
                    accession: "MS:1001460".to_string(),
                    cv_ref: "PSI-MS".to_string(),
                    ..Default::default()
                });
                params.push(CvParamType {
                    name: "cross-link acceptor".to_string(),
                    accession: "MS:1002510".to_string(),
                    cv_ref: "PSI-MS".to_string(),
                    value: Some(xl_group_id.clone()),
                    ..Default::default()
                });
                linkage2.push(ModificationType {
                    location: Some(pos2),
                    cv_param: params,
                    monoisotopic_mass_delta: Some(0.0), // Acceptor always 0.0 in most standards
                    ..Default::default()
                });
            }
            let pep2_id = factory.add_peptide(c_pep2.get(i).unwrap(), linkage2);
            let _pep2_id_plain = pep2_id.replace("ident_", "").replace("pep_", "");

            let prot1_list = get_string_list(c_prot1.get(i).unwrap());
            let start1_list = get_u32_list(c_start1.get(i).unwrap());
            let end1_list = get_u32_list(c_end1.get(i).unwrap());

            let mut ev1_ids = Vec::new();
            for j in 0..prot1_list.len() {
                let p1 = &prot1_list[j];
                let dbseq1_ref = format!("dbseq_{}", p1);
                let is_decoy1 = factory.decoy_map.get(p1).cloned().unwrap_or(false);
                let start = start1_list.get(j).cloned();
                let end = end1_list.get(j).cloned();
                ev1_ids.push(factory.add_peptide_evidence(&pep1_id, &dbseq1_ref, start, end, is_decoy1));
            }

            let prot2_list = get_string_list(c_prot2.get(i).unwrap());
            let start2_list = get_u32_list(c_start2.get(i).unwrap());
            let end2_list = get_u32_list(c_end2.get(i).unwrap());

            let mut ev2_ids = Vec::new();
            for j in 0..prot2_list.len() {
                let p2 = &prot2_list[j];
                let dbseq2_ref = format!("dbseq_{}", p2);
                let is_decoy2 = factory.decoy_map.get(p2).cloned().unwrap_or(false);
                let start = start2_list.get(j).cloned();
                let end = end2_list.get(j).cloned();
                ev2_ids.push(factory.add_peptide_evidence(&pep2_id, &dbseq2_ref, start, end, is_decoy2));
            }

            let mut content1 = Vec::new();
            for ev_id in ev1_ids {
                content1.push(SpectrumIdentificationItemTypeContent::PeptideEvidenceRef(PeptideEvidenceRefType {
                    peptide_evidence_ref: ev_id,
                }));
            }
            content1.push(SpectrumIdentificationItemTypeContent::CvParam(CvParamType {
                name: "cross-link spectrum identification item".to_string(),
                accession: "MS:1002511".to_string(),
                value: Some(xl_group_id.clone()),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            }));
            let sii1_id = format!("SII_{}_{}_p1", i + 1, 1);
            let sii1 = SpectrumIdentificationItemType {
                id: sii1_id,
                charge_state: c_charge.get(i).unwrap(),
                experimental_mass_to_charge: c_exp_mz.as_ref().and_then(|c| c.get(i)).unwrap_or(0.0),
                calculated_mass_to_charge: calc_mz_val,
                rank: c_rank.get(i).unwrap() as i32,
                pass_threshold: true,
                peptide_ref: pep1_id,
                content: content1,
                ..Default::default()
            };

            let mut content2 = Vec::new();
            for ev_id in ev2_ids {
                content2.push(SpectrumIdentificationItemTypeContent::PeptideEvidenceRef(PeptideEvidenceRefType {
                    peptide_evidence_ref: ev_id,
                }));
            }
            content2.push(SpectrumIdentificationItemTypeContent::CvParam(CvParamType {
                name: "cross-link spectrum identification item".to_string(),
                accession: "MS:1002511".to_string(),
                value: Some(xl_group_id.clone()),
                cv_ref: "PSI-MS".to_string(),
                ..Default::default()
            }));

            let sii2_id = format!("SII_{}_{}_p2", i + 1, 2);
            let sii2 = SpectrumIdentificationItemType {
                id: sii2_id,
                charge_state: c_charge.get(i).unwrap(),
                experimental_mass_to_charge: c_exp_mz.as_ref().and_then(|c| c.get(i)).unwrap_or(0.0),
                calculated_mass_to_charge: calc_mz_val,
                rank: c_rank.get(i).unwrap() as i32,
                pass_threshold: true,
                peptide_ref: pep2_id,
                content: content2,
                ..Default::default()
            };

            let mut sii_list = vec![sii1, sii2];
            if let Some(scores) = c_score {
                if let Some(s) = scores.get(i) {
                    for sii in &mut sii_list {
                        sii.content.push(SpectrumIdentificationItemTypeContent::CvParam(CvParamType {
                            name: "xi:score".to_string(),
                            accession: "MS:1002545".to_string(),
                            cv_ref: "PSI-MS".to_string(),
                            value: Some(s.to_string()),
                            ..Default::default()
                        }));
                    }
                }
            }

            let list = &mut grouped_results.entry((sd_ref.to_string(), spec_id.to_string())).or_insert((Vec::new(), sir_params)).0;
            list.extend(sii_list);
        } else {
            // LINEAR OR LOOP-LINK
            let mut linkage = Vec::new();
            let pep1_seq = c_pep1.get(i).unwrap();
            if is_loop_link.get(i).unwrap_or(false) {
                // Add Donor and Acceptor modifications to the SAME peptide
                if let Some(pos1) = c_link_pos1.get(i) {
                    let mut params = Vec::new();
                    let xl_name = c_xl_name.as_ref().and_then(|c| c.get(i));
                    let xl_acc = c_xl_acc.as_ref().and_then(|c| c.get(i));

                    let mut xl_name_param = "unknown modification".to_string();
                    let mut xl_acc_param = "MS:1001460".to_string();
                    let mut xl_mass_val = c_xl_mass.as_ref().and_then(|c| c.get(i)).unwrap_or(0.0);

                    if let Some(acc) = xl_acc {
                        xl_acc_param = acc.to_string();
                        if let Some(override_name) = factory.mod_name_overrides.get(acc) {
                            xl_name_param = override_name.clone();
                            if let Some(data) = cv_data::lookup_mod(acc) {
                                if xl_mass_val == 0.0 {
                                    xl_mass_val = data.mono_mass;
                                }
                            }
                        } else if let Some(data) = cv_data::lookup_mod(acc) {
                            xl_name_param = data.name.to_string();
                            if xl_mass_val == 0.0 {
                                xl_mass_val = data.mono_mass;
                            }
                        } else if let Some(name) = xl_name {
                            xl_name_param = name.to_string();
                        }
                    } else if let Some(name) = xl_name {
                        xl_name_param = name.to_string();
                        if let Some((acc, mass)) = cv_data::lookup_mod_by_name(name) {
                            xl_acc_param = acc.to_string();
                            if let Some(override_name) = factory.mod_name_overrides.get(&xl_acc_param) {
                                xl_name_param = override_name.clone();
                            }
                            if xl_mass_val == 0.0 {
                                xl_mass_val = mass;
                            }
                        }
                    }

                    params.push(CvParamType {
                        name: xl_name_param.clone(),
                        accession: xl_acc_param.clone(),
                        cv_ref: get_cv_ref(&xl_acc_param),
                        ..Default::default()
                    });
                    params.push(CvParamType {
                        name: "crosslink donor".to_string(),
                        accession: "MS:1002509".to_string(),
                        cv_ref: "PSI-MS".to_string(),
                        value: Some(xl_group_id.clone()),
                        ..Default::default()
                    });
                    linkage.push(ModificationType {
                        location: Some(pos1),
                        cv_param: params,
                        monoisotopic_mass_delta: Some(xl_mass_val),
                        ..Default::default()
                    });
                }
                if let Some(pos2) = c_link_pos2.get(i) {
                    let mut params = Vec::new();
                    let xl_name = c_xl_name.as_ref().and_then(|c| c.get(i));
                    let xl_acc = c_xl_acc.as_ref().and_then(|c| c.get(i));

                    let mut xl_name_param = "unknown modification".to_string();
                    let mut xl_acc_param = "MS:1001460".to_string();

                    if let Some(acc) = xl_acc {
                        xl_acc_param = acc.to_string();
                        if let Some(override_name) = factory.mod_name_overrides.get(acc) {
                            xl_name_param = override_name.clone();
                        } else if let Some(data) = cv_data::lookup_mod(acc) {
                            xl_name_param = data.name.to_string();
                        } else if let Some(name) = xl_name {
                            xl_name_param = name.to_string();
                        }
                    } else if let Some(name) = xl_name {
                        xl_name_param = name.to_string();
                        if let Some((acc, _)) = cv_data::lookup_mod_by_name(name) {
                            xl_acc_param = acc.to_string();
                            if let Some(override_name) = factory.mod_name_overrides.get(&xl_acc_param) {
                                xl_name_param = override_name.clone();
                            }
                        }
                    }

                    params.push(CvParamType {
                        name: xl_name_param.clone(),
                        accession: xl_acc_param.clone(),
                        cv_ref: get_cv_ref(&xl_acc_param),
                        ..Default::default()
                    });
                    params.push(CvParamType {
                        name: "crosslink acceptor".to_string(),
                        accession: "MS:1002510".to_string(),
                        cv_ref: "PSI-MS".to_string(),
                        value: Some(xl_group_id.clone()),
                        ..Default::default()
                    });
                    linkage.push(ModificationType {
                        location: Some(pos2),
                        cv_param: params,
                        monoisotopic_mass_delta: Some(0.0),
                        ..Default::default()
                    });
                }
            }

            let pep_id = factory.add_peptide(pep1_seq, linkage);
            let _pep_id_plain = pep_id.replace("ident_", "").replace("pep_", "");
            
            let prot1_list = get_string_list(c_prot1.get(i).unwrap());
            let start1_list = get_u32_list(c_start1.get(i).unwrap());
            let end1_list = get_u32_list(c_end1.get(i).unwrap());

            let mut ev_ids = Vec::new();
            for j in 0..prot1_list.len() {
                let p1 = &prot1_list[j];
                let dbseq1_ref = format!("dbseq_{}", p1);
                let is_decoy1 = factory.decoy_map.get(p1).cloned().unwrap_or(false);
                let start = start1_list.get(j).cloned();
                let end = end1_list.get(j).cloned();
                ev_ids.push(factory.add_peptide_evidence(&pep_id, &dbseq1_ref, start, end, is_decoy1));
            }

            let mut content = Vec::new();
            for ev_id in ev_ids {
                content.push(SpectrumIdentificationItemTypeContent::PeptideEvidenceRef(PeptideEvidenceRefType {
                    peptide_evidence_ref: ev_id,
                }));
            }

            let mut sii = SpectrumIdentificationItemType {
                id: format!("SII_{}_{}", i + 1, 1),
                charge_state: c_charge.get(i).unwrap(),
                experimental_mass_to_charge: c_exp_mz.as_ref().and_then(|c| c.get(i)).unwrap_or(0.0),
                calculated_mass_to_charge: calc_mz_val,
                rank: c_rank.get(i).unwrap() as i32,
                pass_threshold: true,
                peptide_ref: pep_id,
                content,
                ..Default::default()
            };
            
            if is_loop_link.get(i).unwrap_or(false) {
                sii.content.push(SpectrumIdentificationItemTypeContent::CvParam(CvParamType {
                    name: "looplink spectrum identification item".to_string(),
                    accession: "MS:1003329".to_string(),
                    cv_ref: "PSI-MS".to_string(),
                    ..Default::default()
                }));
            }

            if let Some(scores) = c_score {
                if let Some(s) = scores.get(i) {
                    sii.content.push(SpectrumIdentificationItemTypeContent::CvParam(CvParamType {
                        name: "xi:score".to_string(),
                        accession: "MS:1002545".to_string(),
                        cv_ref: "PSI-MS".to_string(),
                        value: Some(s.to_string()),
                        ..Default::default()
                    }));
                }
            }

            grouped_results.entry((sd_ref.to_string(), spec_id.to_string())).or_insert((Vec::new(), sir_params)).0.push(sii);
        }
    }

    // 4. Add grouped results to factory
    let mut sir_idx = 1;
    for ((sd_ref, spec_id), (items, sir_params)) in grouped_results {
        let sir_id = format!("SIR_{}", sir_idx);
        sir_idx += 1;
        factory.add_spectrum_identification_result_with_id(&sir_id, &sd_ref, &spec_id, items, sir_params);
    }

    Ok(factory)
}

#[pyfunction]
pub fn serialize_mzidentml(
    csms: PyDataFrame,
    prot_seqs: PyDataFrame,
    spectra: PyDataFrame,
    metadata: Bound<'_, PyDict>,
) -> PyResult<String> {
    let factory = prepare_factory(csms, prot_seqs, spectra, metadata)?;
    factory.serialize().map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))
}

#[pyfunction]
pub fn write_mzidentml(
    path: String,
    csms: PyDataFrame,
    prot_seqs: PyDataFrame,
    spectra: PyDataFrame,
    metadata: Bound<'_, PyDict>,
) -> PyResult<()> {
    let factory = prepare_factory(csms, prot_seqs, spectra, metadata)?;
    factory.serialize_to_file(&path).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))
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
        let id1 = factory.add_peptide("PEPTIDE", Vec::new());
        let id2 = factory.add_peptide("PEPTIDE", Vec::new());
        assert_eq!(id1, id2);
        assert_eq!(id1, "pep_0");
        
        let sc = factory.doc.sequence_collection.as_ref().unwrap();
        assert_eq!(sc.peptide.len(), 1);
    }

    #[test]
    fn test_add_db_sequence() {
        let mut factory = MzIdentMLFactory::new("test_doc".to_string());
        let id = factory.add_db_sequence("P12345", "ACC123", "MAGA", "DB1", false, None);
        assert_eq!(id, "dbseq_P12345");
        
        let sc = factory.doc.sequence_collection.as_ref().unwrap();
        assert_eq!(sc.db_sequence.len(), 1);
        assert_eq!(sc.db_sequence[0].accession, "ACC123");
    }

    #[test]
    fn test_serialization_basic() {
        let mut factory = MzIdentMLFactory::new("test_doc".to_string());
        factory.add_peptide("PEPTIDE", Vec::new());
        factory.add_db_sequence("P12", "ACC", "M", "DB", false, None);
        
        let xml = factory.serialize().unwrap();
        // println!("XML OUTPUT:\n{}", xml);
        assert!(xml.starts_with("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"));
        assert!(xml.contains("test_doc"));
        assert!(xml.contains("<Peptide"));
        assert!(xml.contains("id=\"pep_0\""));
        assert!(xml.contains("<DBSequence"));
        assert!(xml.contains("id=\"dbseq_P12\""));
    }
}
impl MzIdentMLFactory {
    pub fn parse_proforma(&self, proforma: &str) -> (String, Vec<ModificationType>) {
        let mut clean_seq = String::new();
        let mut mods = Vec::new();
        let chars: Vec<char> = proforma.chars().collect();
        let mut i = 0;
        while i < chars.len() {
            if chars[i] == '[' {
                let start = i + 1;
                let mut end = start;
                while end < chars.len() && chars[end] != ']' {
                    end += 1;
                }
                if end < chars.len() {
                    let mod_str = &proforma[start..end];
                    let mut cv_param = CvParamType {
                        cv_ref: "PSI-MS".to_string(),
                        ..Default::default()
                    };
                    let mut mass_delta = 0.0;
                    
                    if mod_str.contains(':') {
                        let parts: Vec<&str> = mod_str.split(':').collect();
                        let cv_id = parts[0].to_uppercase();
                        cv_param.cv_ref = match cv_id.as_str() {
                            "UNIMOD" => "UNIMOD",
                            "XLMOD" => "XLMOD",
                            "MOD" => "PSI-MOD",
                            _ => "PSI-MS",
                        }.to_string();
                        cv_param.accession = mod_str.to_uppercase();

                        // Priority: 1. Metadata Overrides, 2. CV lookup, 3. fall back to ID part
                        if let Some(override_name) = self.mod_name_overrides.get(&cv_param.accession) {
                            cv_param.name = override_name.clone();
                            if let Some(data) = cv_data::lookup_mod(&cv_param.accession) {
                                mass_delta = data.mono_mass;
                            }
                        } else if let Some(data) = cv_data::lookup_mod(&cv_param.accession) {
                            cv_param.name = data.name.to_string();
                            mass_delta = data.mono_mass;
                        } else {
                            cv_param.name = parts[1].to_string();
                        }
                    } else {
                        if let Some((acc, mass)) = cv_data::lookup_mod_by_name(mod_str) {
                            cv_param.name = mod_str.to_string();
                            cv_param.accession = acc.to_string();
                            // If we have an override for this accession, use it
                            if let Some(override_name) = self.mod_name_overrides.get(&cv_param.accession) {
                                cv_param.name = override_name.clone();
                            }
                            
                            cv_param.cv_ref = get_cv_ref(&acc);
                            mass_delta = mass;
                        } else {
                            cv_param.name = mod_str.to_string();
                            cv_param.accession = "MS:1001460".to_string();
                            cv_param.cv_ref = get_cv_ref(&cv_param.accession);
                            cv_param.value = Some(mod_str.to_string());
                        }
                    }
                    
                    mods.push(ModificationType {
                        location: Some(clean_seq.len() as i32),
                        cv_param: vec![cv_param],
                        monoisotopic_mass_delta: Some(mass_delta),
                        ..Default::default()
                    });
                    i = end + 1;
                } else {
                    i += 1;
                }
            } else {
                clean_seq.push(chars[i]);
                i += 1;
            }
        }
        (clean_seq, mods)
    }
}

pub fn get_cv_ref(accession: &str) -> String {
    if accession.starts_with("MS:") { 
        "PSI-MS".to_string() 
    } else if accession.starts_with("UNIMOD:") { 
        "UNIMOD".to_string() 
    } else if accession.starts_with("XLMOD:") { 
        "XLMOD".to_string() 
    } else if accession.starts_with("MOD:") { 
        "PSI-MOD".to_string() 
    } else if accession.starts_with("UO:") { 
        "UO".to_string() 
    } else { 
        "PSI-MS".to_string() 
    }
}
