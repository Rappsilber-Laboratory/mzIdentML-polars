use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

fn main() {
    let out_dir = env::var("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("cv_data.rs");
    let mut f = File::create(&dest_path).unwrap();

    let cv_dir = Path::new("cv");
    
    // Parse Unimod
    let unimod_path = cv_dir.join("unimod.obo");
    let mods = parse_unimod(&unimod_path);

    writeln!(f, "pub struct ModData {{").unwrap();
    writeln!(f, "    pub name: &'static str,").unwrap();
    writeln!(f, "    pub mono_mass: f64,").unwrap();
    writeln!(f, "}}").unwrap();
    writeln!(f, "").unwrap();

    writeln!(f, "pub fn lookup_mod(accession: &str) -> Option<ModData> {{").unwrap();
    writeln!(f, "    match accession {{").unwrap();
    for (acc, name, mass) in &mods {
        let mut mass_str = mass.to_string();
        if !mass_str.contains('.') && !mass_str.contains('e') {
            mass_str.push_str(".0");
        }
        writeln!(f, "        {:?} => Some(ModData {{ name: {:?}, mono_mass: {} }}),", acc, name, mass_str).unwrap();
    }
    writeln!(f, "        _ => None,").unwrap();
    writeln!(f, "    }}").unwrap();
    writeln!(f, "}}").unwrap();
    writeln!(f, "").unwrap();

    writeln!(f, "pub fn lookup_mod_by_name(name: &str) -> Option<(&'static str, f64)> {{").unwrap();
    writeln!(f, "    match name {{").unwrap();
    // Use a HashMap to avoid duplicate names in the match statement (if any)
    let mut seen_names = std::collections::HashSet::new();
    for (acc, name, mass) in &mods {
        if seen_names.contains(name) {
            continue;
        }
        seen_names.insert(name.clone());
        let mut mass_str = mass.to_string();
        if !mass_str.contains('.') && !mass_str.contains('e') {
            mass_str.push_str(".0");
        }
        writeln!(f, "        {:?} => Some(({:?}, {})),", name, acc, mass_str).unwrap();
    }
    writeln!(f, "        _ => None,").unwrap();
    writeln!(f, "    }}").unwrap();
    writeln!(f, "}}").unwrap();

    println!("cargo:rerun-if-changed=cv/unimod.obo");
    println!("cargo:rerun-if-changed=build.rs");
}

fn parse_unimod(path: &Path) -> Vec<(String, String, f64)> {
    let mut results = Vec::new();
    if !path.exists() {
        return results;
    }

    let file = File::open(path).expect("Could not open unimod.obo");
    let reader = BufReader::new(file);

    let mut current_id = None;
    let mut current_name = None;
    let mut current_mass = None;

    for line in reader.lines() {
        let line = line.unwrap();
        let line = line.trim();

        if line == "[Term]" {
            if let (Some(id), Some(name), Some(mass)) = (current_id.take(), current_name.take(), current_mass.take()) {
                results.push((id, name, mass));
            }
            continue;
        }

        if line.starts_with("id: ") {
            current_id = Some(line[4..].to_string());
        } else if line.starts_with("name: ") {
            current_name = Some(line[6..].to_string());
        } else if line.starts_with("xref: delta_mono_mass \"") {
            let start = "xref: delta_mono_mass \"".len();
            let end = line.rfind('"').unwrap();
            if let Ok(m) = line[start..end].parse::<f64>() {
                current_mass = Some(m);
            }
        }
    }

    // Push last term
    if let (Some(id), Some(name), Some(mass)) = (current_id, current_name, current_mass) {
        results.push((id, name, mass));
    }

    results
}
