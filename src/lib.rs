//! # SCARAB
//!
//! **S**atellite **C**onstellation **A**utonomous **R**elative-motion **A**nalysis & **B**udgeting
//!
//! A toolkit for constellation station-keeping and maintenance.
//! Provides orbital element handling, TLE parsing, numerical propagation,
//! relative orbital element computation, slot-keeping box logic,
//! and maintenance maneuver planning.

pub mod constants;
pub mod elements;
pub mod tle;
pub mod propagator;
pub mod slotbox;
pub mod maneuver;

#[cfg(feature = "python")]
mod pybridge;

#[cfg(feature = "python")]
use pyo3::prelude::*;

#[cfg(feature = "python")]
#[pymodule]
fn scarab(m: &Bound<'_, PyModule>) -> PyResult<()> {
    pybridge::register(m)?;
    Ok(())
}
