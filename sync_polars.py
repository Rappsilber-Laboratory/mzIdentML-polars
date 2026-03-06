import toml
import os

def sync_versions():
    cargo_toml_path = "Cargo.toml"
    pyproject_toml_path = "pyproject.toml"
    
    if not os.path.exists(cargo_toml_path) or not os.path.exists(pyproject_toml_path):
        print("Required files not found.")
        return

    cargo_data = toml.load(cargo_toml_path)
    polars_entry = cargo_data.get("dependencies", {}).get("polars", "")
    polars_rust_ver = polars_entry.get("version", polars_entry) if isinstance(polars_entry, dict) else polars_entry

    import re
    match = re.search(r"0\.(\d+)", polars_rust_ver)
    if match:
        rust_minor = int(match.group(1))
        # Current mapping: Rust 0.51 -> Python 1.31 (offset is -20)
        py_minor = rust_minor - 20
        python_polars_constraint = f"polars>=0.20.0,<1.{py_minor + 1}.0"
    else:
        print(f"Could not parse Polars Rust version: {polars_rust_ver}. Skipping sync.")
        return

    pyproject_data = toml.load(pyproject_toml_path)
    current_deps = pyproject_data.get("project", {}).get("dependencies", [])
    
    new_deps = []
    found = False
    for dep in current_deps:
        if dep.startswith("polars"):
            new_deps.append(python_polars_constraint)
            found = True
        else:
            new_deps.append(dep)
            
    if not found:
        new_deps.append(python_polars_constraint)
        
    pyproject_data["project"]["dependencies"] = new_deps
    
    with open(pyproject_toml_path, "w") as f:
        toml.dump(pyproject_data, f)
    
    print(f"Synced pyproject.toml Polars dependency to {python_polars_constraint} based on Rust version {polars_rust_ver}")

if __name__ == "__main__":
    sync_versions()
