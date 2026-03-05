use pyo3::prelude::*;

pub mod mzidentml;
pub mod polars_writer;

#[pymodule]
fn mzidentml_polars(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(polars_writer::serialize_mzidentml, m)?)?;
    m.add_function(wrap_pyfunction!(polars_writer::write_mzidentml, m)?)?;
    Ok(())
}