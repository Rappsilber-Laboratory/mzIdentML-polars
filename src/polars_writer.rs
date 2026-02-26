use polars::prelude::*;
use pyo3::prelude::*;
use pyo3::types::{PyDict};
use pyo3_polars::PyDataFrame;

#[pyfunction]
pub fn write_mzidentml(
    csms: &PyDataFrame,
    prot_seqs: &PyDataFrame,
    spectra: &PyDataFrame,
    cvs: &PyDict,
) -> Result<String, PolarsError> {
    Ok("TODO".parse()?)
}
