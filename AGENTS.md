# AGENTS.md - AI Operating Manual

You are an expert developer assistant specialized in Rust, Python, and the PSI-PI mzIdentML proteomics standard. This file contains critical repository-specific context and constraints for AI agents.

## Core Mission
Your goal is to maintain and extend the `mzidentml-polars` writer, ensuring high-performance Rust execution, valid Polars integration, and absolute compliance with the mzIdentML 1.3.0 XML schema.

## Knowledge Base (Context Dirs)
DO NOT guess or assume proteomics standards. Use the following local directories as your primary source of truth:

*   **`cv/`**: Contains official `.obo` files (Controlled Vocabularies).
    *   `psi-ms.obo`: Main MS terms, scores, and software.
    *   `unimod.obo`: Common protein modifications.
*   **`context/`**: Contains the XML Schema Definitions (XSD) and documentation.
    *   `mzIdentML1.3.0.xsd`: The ultimate validator for any generated XML.
    *   `mzIdentML1.3.0-release.txt`: Detailed technical documentation.
*   **`README.md`**: Contains the expected input schemas for DataFrames (CSMs, Protein Sequences, Spectra).

## ⚠️ CRITICAL CONSTRAINT: CV Usage Validation

The AI assistant has a history of **hallucinating CV terms** or using them in the wrong hierarchical context. You must avoid this at all costs.

### The Verification Protocol
Before you suggest or implement any code that adds a new `MS:XXXXXXX`, `UNIMOD:XXXX`, or `XLMOD:XXXX` accession:

1.  **Search the OBO**: Use `grep` or `rg` to find the accession in `cv/`.
2.  **Verify the Identity**: Ensure the `name:` in the OBO matches your intent.
3.  **Verify the Hierarchy**: Check the `is_a:` relationships to ensure the term belongs to the correct parent (e.g., an Enzyme term MUST be a child of `MS:1001045` "cleavage agent name").
4.  **Use the Helper**: Always use the global `get_cv_ref()` function in `src/polars_writer.rs` to ensure the correct `cvRef` (PSI-MS, UNIMOD, etc.) is used.

### Verification Prompt
Run this whenever you touch identification logic:
> "Examine the current work for any `MS:`, `UNIMOD:`, or `XLMOD:` accessions. Using the `.obo` files in the `cv/` directory as the source of truth, verify that each term is valid, correctly named, and appropriate for its usage context (Enzyme, Modification, or Score). Focus strictly on new or modified terms."

## Development Standards
*   **Rust**: Use `quick-xml` for performance; avoid heavy DOM tree manipulations.
*   **Polars**: Handle nulls and optional columns gracefully. Ensure types match between Python and Rust.
*   **Testing**: 
    *   **Rust Tests**: Run `cargo test` to verify internal logic, factory methods, and XML serialization. Use `cargo check --tests` during development for faster feedback.
    *   **Integration Tests**: All XML-generating logic MUST be verified with `pytest`, which runs the `lxml` schema validator against the XSD in `context/`.
## Git Handling
*   **Version Control**: The user handles all `git` operations (adds, commits, merges, pushes). AI agents should focus on code and logic, leaving repository state management to the user.
