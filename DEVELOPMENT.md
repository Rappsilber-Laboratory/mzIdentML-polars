# Development Guide

This document provides information for developers working on the `mzIdentML-polars` project.

## AI-Assisted Development

This project has been developed with heavy usage of **AI assistance (Antigravity)**. While this has significantly accelerated development, it requires specific attention to detail during code reviews and future contributions.

### LLM Context Directories

If you are continuing development using an LLM, the following directories are crucial for providing necessary context:

*   **`cv/`**: Contains the Controlled Vocabulary (CV) files (`.obo`) used for proteomics standards. These are essential for looking up correct accessions and names.
*   **`context/`**: Contains documentation, schemas (XSD), and example files that define the expected structure of mzIdentML files.

## Critical: CV Usage Validation

**IMPORTANT WARNING:** During previous development iterations, the AI has multiple times "hallucinated" CV terms or accessions.

### Recommended Verification Prompt

When adding new features or modifying identification logic, use the following prompt to ensure CV integrity without scope creep:

> "Examine the `git diff` for any new `MS:`, `UNIMOD:`, or `XLMOD:` accessions. Using the `.obo` files in the `cv/` directory as the source of truth, verify that each new term is valid, correctly named, and appropriate for its usage context (e.g., as an Enzyme, Modification, or Score). 
> 
> **Important**: Focus only on the new CV usages. Do not optimize code, refactor existing logic, or re-verify unchanged terms unless explicitly requested."

## Build and Test

Refer to the [README.md](README.md) for basic installation and testing instructions. Always ensure that `pytest` passes before submitting a pull request, as it performs schema validation against the official mzIdentML XSD.
