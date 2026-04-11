use pyo3::prelude::*;

pub mod mzidentml;
pub mod polars_writer;

#[pymodule(gil_used = false)]
fn _mzidentml_polars(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(polars_writer::serialize_mzidentml, m)?)?;
    m.add_function(wrap_pyfunction!(polars_writer::write_mzidentml, m)?)?;
    Ok(())
}