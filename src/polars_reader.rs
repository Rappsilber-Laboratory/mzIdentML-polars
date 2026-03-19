use quick_xml::events::Event;
use quick_xml::Reader;
use std::collections::HashMap;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3_polars::PyDataFrame;

#[derive(Debug)]
struct DBSequence {
    accession: String,
}

#[derive(Debug)]
struct PeptideMod {
    location: i32,
    mass_delta: Option<f64>,
    accession: Option<String>,
    name: Option<String>,
}

#[derive(Debug)]
struct Peptide {
    sequence: String,
    mods: Vec<PeptideMod>,
}

impl Peptide {
    fn to_proforma(&self) -> String {
        if self.mods.is_empty() {
            return self.sequence.clone();
        }

        let mut proforma = String::new();
        let mut mods_by_loc: HashMap<i32, Vec<&PeptideMod>> = HashMap::new();
        for m in &self.mods {
            mods_by_loc.entry(m.location).or_default().push(m);
        }

        // N-term mods (location 0)
        if let Some(n_mods) = mods_by_loc.get(&0) {
            for m in n_mods {
                if let Some(acc) = &m.accession {
                    proforma.push_str(&format!("[{}]-", acc));
                } else if let Some(mass) = m.mass_delta {
                    proforma.push_str(&format!("[{:+}]-", mass));
                }
            }
        }

        for (i, c) in self.sequence.chars().enumerate() {
            proforma.push(c);
            let loc = (i + 1) as i32;
            if let Some(aa_mods) = mods_by_loc.get(&loc) {
                for m in aa_mods {
                    if let Some(acc) = &m.accession {
                        proforma.push_str(&format!("[{}]", acc));
                    } else if let Some(mass) = m.mass_delta {
                        proforma.push_str(&format!("[{:+}]", mass));
                    }
                }
            }
        }

        // C-term mods (location = seq.len() + 1)
        let c_loc = self.sequence.len() as i32 + 1;
        if let Some(c_mods) = mods_by_loc.get(&c_loc) {
            for m in c_mods {
                if let Some(acc) = &m.accession {
                    proforma.push_str(&format!("-[{}]", acc));
                } else if let Some(mass) = m.mass_delta {
                    proforma.push_str(&format!("-[{:+}]", mass));
                }
            }
        }

        proforma
    }
}

#[derive(Debug)]
struct PeptideEvidence {
    dbseq_ref: String,
    start: Option<u32>,
    end: Option<u32>,
    is_decoy: bool,
}

#[derive(Debug, Clone)]
struct SpectrumIdentificationItem {
    id: String,
    rank: u32,
    charge_state: i32,
    calc_mz: Option<f64>,
    exp_mz: Option<f64>,
    peptide_ref: String,
    peptide_evidence_refs: Vec<String>,
    score: Option<f64>,
    crosslinker_donor: bool,
    crosslinker_acceptor: bool,
    cross_link_ref: Option<String>,
}

pub fn parse_mzidentml_to_dfs(path: &str) -> std::result::Result<(DataFrame, DataFrame, DataFrame), String> {
    let xml = std::fs::read_to_string(path)
        .map_err(|e| format!("Failed to open file: {}", e))?;
    let mut reader = Reader::from_str(&xml);
    reader.config_mut().trim_text(true);

    let mut buf = Vec::new();
    
    let mut db_sequences: HashMap<String, DBSequence> = HashMap::new();
    let mut peptides: HashMap<String, Peptide> = HashMap::new();
    let mut peptide_evidences: HashMap<String, PeptideEvidence> = HashMap::new();
    let mut spectra_data: HashMap<String, String> = HashMap::new(); // id -> location
    
    struct SirData {
        spectrum_id: String,
        spectra_data_ref: String,
        items: Vec<SpectrumIdentificationItem>,
    }
    let mut sir_list: Vec<SirData> = Vec::new();

    let mut current_peptide_seq = String::new();
    let mut current_peptide_mods = Vec::new();
    let mut in_peptide_seq = false;
    
    let mut current_peptide_id = String::new();
    let mut current_mod = PeptideMod { location: 0, mass_delta: None, accession: None, name: None };
    let mut in_mod = false;

    let mut current_sir = SirData {
        spectrum_id: String::new(),
        spectra_data_ref: String::new(),
        items: Vec::new(),
    };
    let mut current_sii = SpectrumIdentificationItem {
        id: String::new(), rank: 0, charge_state: 0, calc_mz: None, exp_mz: None, peptide_ref: String::new(), peptide_evidence_refs: Vec::new(), score: None, crosslinker_donor: false, crosslinker_acceptor: false, cross_link_ref: None
    };
    
    let mut in_sir = false;
    let mut in_sii = false;

    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(ref e)) | Ok(Event::Empty(ref e)) => {
                let local_name = e.name().local_name();
                let name_str = std::str::from_utf8(local_name.as_ref()).unwrap_or("");
                
                match name_str {
                    "DBSequence" => {
                        let mut id = String::new();
                        let mut acc = String::new();
                        for attr in e.attributes().filter_map(|a| a.ok()) {
                            match attr.key.as_ref() {
                                b"id" => id = String::from_utf8_lossy(&attr.value).into_owned(),
                                b"accession" => acc = String::from_utf8_lossy(&attr.value).into_owned(),
                                _ => ()
                            }
                        }
                        db_sequences.insert(id, DBSequence { accession: acc });
                    },
                    "Peptide" => {
                        for attr in e.attributes().filter_map(|a| a.ok()) {
                            if attr.key.as_ref() == b"id" {
                                current_peptide_id = String::from_utf8_lossy(&attr.value).into_owned();
                            }
                        }
                        current_peptide_seq.clear();
                        current_peptide_mods.clear();
                    },
                    "PeptideSequence" => {
                        in_peptide_seq = true;
                    },
                    "Modification" => {
                        in_mod = true;
                        current_mod = PeptideMod { location: 0, mass_delta: None, accession: None, name: None };
                        for attr in e.attributes().filter_map(|a| a.ok()) {
                            match attr.key.as_ref() {
                                b"location" => current_mod.location = String::from_utf8_lossy(&attr.value).parse().unwrap_or(0),
                                b"monoisotopicMassDelta" => current_mod.mass_delta = String::from_utf8_lossy(&attr.value).parse().ok(),
                                _ => ()
                            }
                        }
                    },
                    "cvParam" if in_mod => {
                        for attr in e.attributes().filter_map(|a| a.ok()) {
                            match attr.key.as_ref() {
                                b"accession" => {
                                    let acc = String::from_utf8_lossy(&attr.value).into_owned();
                                    // if it's UNIMOD or XLMOD etc
                                    if acc.starts_with("UNIMOD:") || acc.starts_with("MOD:") || acc.starts_with("XLMOD:") {
                                        current_mod.accession = Some(acc);
                                    }
                                },
                                b"name" => {
                                    current_mod.name = Some(String::from_utf8_lossy(&attr.value).into_owned());
                                }
                                _ => ()
                            }
                        }
                    },
                    "PeptideEvidence" => {
                        let mut id = String::new();
                        let mut dbseq_ref = String::new();
                        let mut start = None;
                        let mut end = None;
                        let mut is_decoy = false;
                        for attr in e.attributes().filter_map(|a| a.ok()) {
                            match attr.key.as_ref() {
                                b"id" => id = String::from_utf8_lossy(&attr.value).into_owned(),
                                b"dBSequence_ref" => dbseq_ref = String::from_utf8_lossy(&attr.value).into_owned(),
                                b"start" => start = String::from_utf8_lossy(&attr.value).parse::<u32>().ok(),
                                b"end" => end = String::from_utf8_lossy(&attr.value).parse::<u32>().ok(),
                                b"isDecoy" => is_decoy = match String::from_utf8_lossy(&attr.value).as_ref() {
                                    "true" | "1" => true,
                                    _ => false,
                                },
                                _ => ()
                            }
                        }
                        peptide_evidences.insert(id, PeptideEvidence {
                            dbseq_ref, start, end, is_decoy
                        });
                    },
                    "SpectraData" => {
                        let mut id = String::new();
                        let mut location = String::new();
                        for attr in e.attributes().filter_map(|a| a.ok()) {
                            match attr.key.as_ref() {
                                b"id" => id = String::from_utf8_lossy(&attr.value).into_owned(),
                                b"location" => location = String::from_utf8_lossy(&attr.value).into_owned(),
                                _ => ()
                            }
                        }
                        spectra_data.insert(id, location);
                    },
                    "SpectrumIdentificationResult" => {
                        in_sir = true;
                        current_sir = SirData {
                            spectrum_id: String::new(),
                            spectra_data_ref: String::new(),
                            items: Vec::new(),
                        };
                        for attr in e.attributes().filter_map(|a| a.ok()) {
                            match attr.key.as_ref() {
                                b"spectrumID" => current_sir.spectrum_id = String::from_utf8_lossy(&attr.value).into_owned(),
                                b"spectraData_ref" => current_sir.spectra_data_ref = String::from_utf8_lossy(&attr.value).into_owned(),
                                _ => ()
                            }
                        }
                    },
                    "SpectrumIdentificationItem" => {
                        in_sii = true;
                        current_sii = SpectrumIdentificationItem {
                            id: String::new(), rank: 0, charge_state: 0, calc_mz: None, exp_mz: None, peptide_ref: String::new(), peptide_evidence_refs: Vec::new(), score: None, crosslinker_donor: false, crosslinker_acceptor: false, cross_link_ref: None
                        };
                        for attr in e.attributes().filter_map(|a| a.ok()) {
                            match attr.key.as_ref() {
                                b"id" => current_sii.id = String::from_utf8_lossy(&attr.value).into_owned(),
                                b"rank" => current_sii.rank = String::from_utf8_lossy(&attr.value).parse().unwrap_or(0),
                                b"chargeState" => current_sii.charge_state = String::from_utf8_lossy(&attr.value).parse().unwrap_or(0),
                                b"calculatedMassToCharge" => current_sii.calc_mz = String::from_utf8_lossy(&attr.value).parse().ok(),
                                b"experimentalMassToCharge" => current_sii.exp_mz = String::from_utf8_lossy(&attr.value).parse().ok(),
                                b"peptide_ref" => current_sii.peptide_ref = String::from_utf8_lossy(&attr.value).into_owned(),
                                _ => ()
                            }
                        }
                    },
                    "PeptideEvidenceRef" if in_sii => {
                        for attr in e.attributes().filter_map(|a| a.ok()) {
                            if attr.key.as_ref() == b"peptideEvidence_ref" {
                                current_sii.peptide_evidence_refs.push(String::from_utf8_lossy(&attr.value).into_owned());
                            }
                        }
                    },
                    "cvParam" if in_sii => {
                        let mut name = String::new();
                        let mut acc = String::new();
                        let mut val = String::new();
                        for attr in e.attributes().filter_map(|a| a.ok()) {
                            match attr.key.as_ref() {
                                b"name" => name = String::from_utf8_lossy(&attr.value).into_owned(),
                                b"accession" => acc = String::from_utf8_lossy(&attr.value).into_owned(),
                                b"value" => val = String::from_utf8_lossy(&attr.value).into_owned(),
                                _ => ()
                            }
                        }
                        
                        let lower_name = name.to_lowercase();
                        if lower_name.contains("score") {
                            if current_sii.score.is_none() {
                                current_sii.score = val.parse::<f64>().ok();
                            }
                        }
                        
                        if lower_name == "cross-link donor" { current_sii.crosslinker_donor = true; }
                        if lower_name == "cross-link acceptor" { current_sii.crosslinker_acceptor = true; }
                        if name == "cross-link spectrum identification item" || acc == "MS:1002511" { current_sii.cross_link_ref = Some(val); }
                    },
                    _ => ()
                }
            },
            Ok(Event::End(ref e)) => {
                let local_name = e.name().local_name();
                let name = std::str::from_utf8(local_name.as_ref()).unwrap_or("");
                match name {
                    "Modification" => {
                        in_mod = false;
                        let m = std::mem::replace(&mut current_mod, PeptideMod { location: 0, mass_delta: None, accession: None, name: None });
                        current_peptide_mods.push(m);
                    },
                    "Peptide" => {
                        peptides.insert(current_peptide_id.clone(), Peptide {
                            sequence: current_peptide_seq.clone(),
                            mods: std::mem::replace(&mut current_peptide_mods, Vec::new()),
                        });
                    },
                    "PeptideSequence" => in_peptide_seq = false,
                    "SpectrumIdentificationResult" => {
                        in_sir = false;
                        let s = std::mem::replace(&mut current_sir, SirData { spectrum_id: String::new(), spectra_data_ref: String::new(), items: Vec::new() });
                        sir_list.push(s);
                    },
                    "SpectrumIdentificationItem" => {
                        in_sii = false;
                        let s = std::mem::replace(&mut current_sii, SpectrumIdentificationItem { id: String::new(), rank: 0, charge_state: 0, calc_mz: None, exp_mz: None, peptide_ref: String::new(), peptide_evidence_refs: Vec::new(), score: None, crosslinker_donor: false, crosslinker_acceptor: false, cross_link_ref: None });
                        current_sir.items.push(s);
                    },
                    _ => ()
                }
            },
            Ok(Event::Text(ref e)) => {
                if in_peptide_seq {
                    current_peptide_seq = String::from_utf8_lossy(e.as_ref()).into_owned();
                }
            },
            Ok(Event::Eof) => break,
            Err(_) => break,
            _ => (),
        }
        buf.clear();
    }

    // Now build DataFrames
    // prot_seqs
    let mut prot_ids = Vec::new();
    let mut prot_accs = Vec::new();
    let mut prot_seqs = Vec::new();
    let mut prot_names = Vec::new();
    let mut prot_decoys = Vec::new();
    
    for (id, db) in &db_sequences {
        prot_ids.push(id.clone());
        prot_accs.push(db.accession.clone());
        prot_seqs.push(String::new());
        prot_names.push(None::<String>);
        prot_decoys.push(false); // will be inferred
    }

    let prot_df = df!(
        "protein_id" => prot_ids,
        "accession" => prot_accs,
        "sequence" => prot_seqs,
        "protein_name" => prot_names,
        "is_decoy" => prot_decoys
    ).unwrap();

    // spectra
    let mut spec_ids = Vec::new();
    let mut spec_locs = Vec::new();
    for (id, loc) in &spectra_data {
        spec_ids.push(id.clone());
        spec_locs.push(loc.clone());
    }
    let spectra_df = df!(
        "spectrum_id" => spec_ids,
        "file_path" => spec_locs
    ).unwrap();

    // csms
    let mut csm_spectrum_id = Vec::new();
    let mut csm_file_path = Vec::new();
    let mut csm_pep1_seq = Vec::new();
    let mut csm_charge = Vec::new();
    let mut csm_rank = Vec::new();
    let mut csm_is_crosslink = Vec::new();
    let mut csm_is_looplink = Vec::new();
    let mut csm_exp_mz = Vec::new();
    let mut csm_calc_mz = Vec::new();
    let mut csm_score = Vec::new();
    
    // crosslink peptide 2 features
    let mut csm_pep2_seq = Vec::new();

    let mut mapped_prots_builder = polars::chunked_array::builder::ListStringChunkedBuilder::new("protein1_id".into(), sir_list.len(), 5);
    let mut mapped_starts_builder = polars::chunked_array::builder::ListPrimitiveChunkedBuilder::<UInt32Type>::new("peptide1_start".into(), sir_list.len(), 5, polars::datatypes::DataType::UInt32);
    let mut mapped_ends_builder = polars::chunked_array::builder::ListPrimitiveChunkedBuilder::<UInt32Type>::new("peptide1_end".into(), sir_list.len(), 5, polars::datatypes::DataType::UInt32);

    let mut mapped_prots2_builder = polars::chunked_array::builder::ListStringChunkedBuilder::new("protein2_id".into(), sir_list.len(), 5);
    let mut mapped_starts2_builder = polars::chunked_array::builder::ListPrimitiveChunkedBuilder::<UInt32Type>::new("peptide2_start".into(), sir_list.len(), 5, polars::datatypes::DataType::UInt32);
    let mut mapped_ends2_builder = polars::chunked_array::builder::ListPrimitiveChunkedBuilder::<UInt32Type>::new("peptide2_end".into(), sir_list.len(), 5, polars::datatypes::DataType::UInt32);

    // Mapping items for crosslinks
    let mut sii_map: HashMap<String, SpectrumIdentificationItem> = HashMap::new();
    for sir in &sir_list {
        for sii in &sir.items {
            sii_map.insert(sii.id.clone(), sii.clone());
        }
    }

    for sir in &sir_list {
        let file_path = spectra_data.get(&sir.spectra_data_ref).cloned().unwrap_or_default();
        for sii in &sir.items {
            let mut is_donor = false;
            let mut is_acceptor = false;
            let mut is_looplink = false;
            let mut current_proforma = String::new();

            if let Some(pep) = peptides.get(&sii.peptide_ref) {
                current_proforma = pep.to_proforma();
                
                let mut donor_count = 0;
                for m in &pep.mods {
                    if let Some(acc) = &m.accession {
                        if acc == "MS:1002509" { // cross-link donor
                            is_donor = true;
                            donor_count += 1;
                        }
                        if acc == "MS:1002510" { // cross-link acceptor
                            is_acceptor = true;
                        }
                    } else if let Some(name) = &m.name {
                        let ln = name.to_lowercase();
                        if ln == "cross-link donor" {
                            is_donor = true;
                            donor_count += 1;
                        }
                        if ln == "cross-link acceptor" {
                            is_acceptor = true;
                        }
                    }
                }
                if is_donor && is_acceptor {
                    is_looplink = true;
                }
            }

            if is_acceptor && !is_donor {
                continue;
            }

            csm_pep1_seq.push(current_proforma);
            csm_spectrum_id.push(sir.spectrum_id.clone());
            csm_file_path.push(file_path.clone());
            csm_charge.push(sii.charge_state);
            csm_rank.push(sii.rank);
            csm_exp_mz.push(sii.exp_mz);
            csm_calc_mz.push(sii.calc_mz);
            csm_score.push(sii.score);
            csm_is_crosslink.push(is_donor && !is_looplink);
            csm_is_looplink.push(is_looplink);

            let mut builder_prots = Vec::new();
            let mut builder_starts = Vec::new();
            let mut builder_ends = Vec::new();

            for ev_ref in &sii.peptide_evidence_refs {
                if let Some(ev) = peptide_evidences.get(ev_ref) {
                    builder_prots.push(ev.dbseq_ref.clone());
                    builder_starts.push(ev.start);
                    builder_ends.push(ev.end);
                }
            }
            mapped_prots_builder.append_series(&Series::new("".into(), &builder_prots)).unwrap_or_default();
            mapped_starts_builder.append_series(&Series::new("".into(), &builder_starts)).unwrap_or_default();
            mapped_ends_builder.append_series(&Series::new("".into(), &builder_ends)).unwrap_or_default();

            if is_donor && !is_looplink {
                let mut p2_seq = String::new();
                let mut builder_prots2 = Vec::new();
                let mut builder_starts2 = Vec::new();
                let mut builder_ends2 = Vec::new();

                if let Some(cross_ref) = &sii.cross_link_ref {
                    if let Some(sii_acceptor) = sii_map.get(cross_ref) {
                        if let Some(pep2) = peptides.get(&sii_acceptor.peptide_ref) {
                            p2_seq = pep2.to_proforma();
                        }
                        
                        for ev_ref in &sii_acceptor.peptide_evidence_refs {
                            if let Some(ev) = peptide_evidences.get(ev_ref) {
                                builder_prots2.push(ev.dbseq_ref.clone());
                                builder_starts2.push(ev.start);
                                builder_ends2.push(ev.end);
                            }
                        }
                    }
                }
                csm_pep2_seq.push(Some(p2_seq));
                mapped_prots2_builder.append_series(&Series::new("".into(), &builder_prots2)).unwrap_or_default();
                mapped_starts2_builder.append_series(&Series::new("".into(), &builder_starts2)).unwrap_or_default();
                mapped_ends2_builder.append_series(&Series::new("".into(), &builder_ends2)).unwrap_or_default();
            } else {
                csm_pep2_seq.push(None);
                mapped_prots2_builder.append_null();
                mapped_starts2_builder.append_null();
                mapped_ends2_builder.append_null();
            }
        }
    }

    let mut prots_series = mapped_prots_builder.finish().into_series();
    let mut starts_series = mapped_starts_builder.finish().into_series();
    let mut ends_series = mapped_ends_builder.finish().into_series();

    let mut prots2_series = mapped_prots2_builder.finish().into_series();
    let mut starts2_series = mapped_starts2_builder.finish().into_series();
    let mut ends2_series = mapped_ends2_builder.finish().into_series();

    let csms_df = df!(
        "spectrum_id" => csm_spectrum_id,
        "file_path" => csm_file_path,
        "peptide1_seq" => csm_pep1_seq,
        "protein1_id" => prots_series,
        "peptide1_start" => starts_series,
        "peptide1_end" => ends_series,
        "charge" => csm_charge,
        "rank" => csm_rank,
        "is_crosslink" => csm_is_crosslink,
        "is_looplink" => csm_is_looplink,
        "experimental_mz" => csm_exp_mz,
        "calculated_mz" => csm_calc_mz,
        "score" => csm_score,
        "peptide2_seq" => csm_pep2_seq,
        "protein2_id" => prots2_series,
        "peptide2_start" => starts2_series,
        "peptide2_end" => ends2_series
    ).map_err(|e| format!("{}", e))?;

    Ok((csms_df, prot_df, spectra_df))
}

#[pyfunction]
pub fn read_mzidentml(path: String) -> PyResult<(PyDataFrame, PyDataFrame, PyDataFrame)> {
    let (csms_df, prot_df, spectra_df) = parse_mzidentml_to_dfs(&path)
        .map_err(|e| PyErr::new::<PyValueError, _>(e))?;
    Ok((PyDataFrame(csms_df), PyDataFrame(prot_df), PyDataFrame(spectra_df)))
}


