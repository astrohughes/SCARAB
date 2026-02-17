//! Python bindings via PyO3 for SCARAB.
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::constants;
use crate::elements::{KeplerianElements, MeanElements, RelativeOrbitalElements};
use crate::tle::Tle;
use crate::propagator::{StateVector, ForceConfig, NumericalPropagator, IntegratorConfig};
use crate::slotbox::{SlotBox, SlotStatus};
use crate::maneuver;


// MeanElements
#[pyclass(name = "MeanElements")]
#[derive(Clone)]
pub struct PyMeanElements {
    pub(crate) inner: MeanElements,
}

#[pymethods]
impl PyMeanElements {
    #[new]
    fn new(a: f64, e: f64, i_deg: f64, raan_deg: f64, aop_deg: f64, ma_deg: f64, epoch: f64) -> Self {
        let kep = KeplerianElements::from_degrees(a, e, i_deg, raan_deg, aop_deg, ma_deg);
        PyMeanElements { inner: kep.to_mean(epoch) }
    }

    fn propagate(&self, dt: f64) -> PyMeanElements {
        PyMeanElements { inner: self.inner.propagate(dt) }
    }

    fn raan_rate_deg_day(&self) -> f64 {
        self.inner.raan_rate() * constants::RAD2DEG * constants::SOLAR_DAY
    }

    fn aop_rate_deg_day(&self) -> f64 {
        self.inner.aop_rate() * constants::RAD2DEG * constants::SOLAR_DAY
    }

    fn period(&self) -> f64 {
        constants::TAU / (constants::MU_EARTH / self.inner.a.powi(3)).sqrt()
    }

    /// Convert to Cartesian state vector [x,y,z,vx,vy,vz] in km and km/s.
    fn to_state_vector(&self) -> Vec<f64> {
        let sv = StateVector::from_mean_elements(&self.inner);
        vec![sv.r[0], sv.r[1], sv.r[2], sv.v[0], sv.v[1], sv.v[2]]
    }

    #[getter] fn a(&self) -> f64 { self.inner.a }
    #[getter] fn e(&self) -> f64 { self.inner.e }
    #[getter] fn i_deg(&self) -> f64 { self.inner.i * constants::RAD2DEG }
    #[getter] fn raan_deg(&self) -> f64 { self.inner.raan * constants::RAD2DEG }
    #[getter] fn aop_deg(&self) -> f64 { self.inner.aop * constants::RAD2DEG }
    #[getter] fn ma_deg(&self) -> f64 { self.inner.ma * constants::RAD2DEG }
    #[getter] fn epoch(&self) -> f64 { self.inner.epoch }

    fn __repr__(&self) -> String {
        format!(
            "MeanElements(a={:.3} km, e={:.6}, i={:.4}°, RAAN={:.4}°, AoP={:.4}°, MA={:.4}°)",
            self.inner.a, self.inner.e,
            self.inner.i * constants::RAD2DEG,
            self.inner.raan * constants::RAD2DEG,
            self.inner.aop * constants::RAD2DEG,
            self.inner.ma * constants::RAD2DEG,
        )
    }
}

// ROE
#[pyclass(name = "ROE")]
#[derive(Clone)]
pub struct PyROE {
    inner: RelativeOrbitalElements,
}

#[pymethods]
impl PyROE {
    #[staticmethod]
    fn from_elements(chief: &PyMeanElements, deputy: &PyMeanElements) -> Self {
        PyROE { inner: RelativeOrbitalElements::from_absolute(&chief.inner, &deputy.inner) }
    }

    fn da_meters(&self, a_chief: f64) -> f64 { self.inner.da_meters(a_chief) }
    fn in_track_km(&self, a_chief: f64) -> f64 { self.inner.in_track_km(a_chief) }
    fn cross_track_amplitude_km(&self, a_chief: f64) -> f64 { self.inner.cross_track_amplitude_km(a_chief) }
    fn de_magnitude(&self) -> f64 { self.inner.de_magnitude() }
    fn di_arcsec(&self) -> f64 { self.inner.di_magnitude() * 180.0 * 3600.0 / std::f64::consts::PI }

    #[getter] fn da(&self) -> f64 { self.inner.da }
    #[getter] fn dlambda(&self) -> f64 { self.inner.dlambda }
    #[getter] fn dex(&self) -> f64 { self.inner.dex }
    #[getter] fn dey(&self) -> f64 { self.inner.dey }
    #[getter] fn dix(&self) -> f64 { self.inner.dix }
    #[getter] fn diy(&self) -> f64 { self.inner.diy }

    fn __repr__(&self) -> String {
        format!(
            "ROE(δa={:.2e}, δλ={:.2e}, δex={:.2e}, δey={:.2e}, δix={:.2e}, δiy={:.2e})",
            self.inner.da, self.inner.dlambda,
            self.inner.dex, self.inner.dey,
            self.inner.dix, self.inner.diy,
        )
    }
}

// SlotBox
#[pyclass(name = "SlotBox")]
#[derive(Clone)]
pub struct PySlotBox {
    inner: SlotBox,
}

#[pymethods]
impl PySlotBox {
    #[new]
    fn new(da_meters: f64, dlambda_km: f64, de: f64, di_arcsec: f64, a_chief: f64) -> Self {
        PySlotBox {
            inner: SlotBox::from_physical_units(da_meters, dlambda_km, de, di_arcsec, a_chief),
        }
    }

    fn check(&self, chief: &PyMeanElements, roe: &PyROE, py: Python<'_>) -> PyResult<Py<PyDict>> {
        let status = self.inner.check(&roe.inner, &chief.inner);
        let dict = PyDict::new(py);
        dict.set_item("violated", status.violation.any())?;
        dict.set_item("da_margin", status.da_margin)?;
        dict.set_item("dlambda_margin", status.dlambda_margin)?;
        dict.set_item("de_margin", status.de_margin)?;
        dict.set_item("di_margin", status.di_margin)?;
        dict.set_item("time_to_violation_s", status.time_to_violation)?;
        dict.set_item("da_violated", status.violation.da_violated)?;
        dict.set_item("dlambda_violated", status.violation.dlambda_violated)?;
        dict.set_item("de_violated", status.violation.de_violated)?;
        dict.set_item("di_violated", status.violation.di_violated)?;
        Ok(dict.into())
    }
}

// TLE
#[pyclass(name = "TLE")]
#[derive(Clone)]
pub struct PyTle {
    inner: Tle,
}

#[pymethods]
impl PyTle {
    /// Parse a TLE from two lines.
    #[staticmethod]
    fn parse(line1: &str, line2: &str) -> PyResult<Self> {
        Tle::parse(line1, line2)
            .map(|t| PyTle { inner: t })
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Parse a TLE from three lines (name + line1 + line2).
    #[staticmethod]
    fn parse_3line(name: &str, line1: &str, line2: &str) -> PyResult<Self> {
        Tle::parse_3line(name, line1, line2)
            .map(|t| PyTle { inner: t })
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Parse a batch of TLEs from a multi-line string.
    #[staticmethod]
    fn parse_batch(text: &str) -> PyResult<Vec<PyTle>> {
        Tle::parse_batch(text)
            .map(|tles| tles.into_iter().map(|t| PyTle { inner: t }).collect())
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Convert to MeanElements.
    fn to_mean_elements(&self) -> PyMeanElements {
        PyMeanElements { inner: self.inner.to_mean_elements() }
    }

    /// Semi-major axis (km).
    fn semi_major_axis(&self) -> f64 { self.inner.semi_major_axis() }

    /// Altitude (km).
    fn altitude(&self) -> f64 { self.inner.altitude() }

    /// Orbital period (seconds).
    fn period(&self) -> f64 { self.inner.period() }

    /// Whether this is a sun-synchronous orbit.
    fn is_sun_synchronous(&self) -> bool { self.inner.is_sun_synchronous() }

    /// Epoch as seconds since J2000.
    fn epoch_j2000(&self) -> f64 { self.inner.epoch_j2000_seconds() }

    #[getter] fn name(&self) -> Option<String> { self.inner.name.clone() }
    #[getter] fn norad_id(&self) -> u32 { self.inner.norad_id }
    #[getter] fn inclination_deg(&self) -> f64 { self.inner.inclination_deg }
    #[getter] fn raan_deg(&self) -> f64 { self.inner.raan_deg }
    #[getter] fn eccentricity(&self) -> f64 { self.inner.eccentricity }
    #[getter] fn arg_perigee_deg(&self) -> f64 { self.inner.arg_perigee_deg }
    #[getter] fn mean_anomaly_deg(&self) -> f64 { self.inner.mean_anomaly_deg }
    #[getter] fn mean_motion(&self) -> f64 { self.inner.mean_motion_rev_day }
    #[getter] fn bstar(&self) -> f64 { self.inner.bstar }
    #[getter] fn epoch_year(&self) -> u16 { self.inner.epoch_year }
    #[getter] fn epoch_day(&self) -> f64 { self.inner.epoch_day }
    #[getter] fn classification(&self) -> char { self.inner.classification }
    #[getter] fn intl_designator(&self) -> String { self.inner.intl_designator.clone() }

    fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }
}

// Numerical Propagator
#[pyclass(name = "Propagator")]
pub struct PyPropagator {
    inner: NumericalPropagator,
}

#[pymethods]
impl PyPropagator {
    /// Create a numerical propagator.
    ///
    /// Args:
    ///     j2: Include J2 perturbation (default: True)
    ///     drag: Include atmospheric drag (default: False)
    ///     ballistic_coefficient: BC in kg/m² for drag (default: 50.0)
    #[new]
    #[pyo3(signature = (j2=true, drag=false, ballistic_coefficient=50.0))]
    fn new(j2: bool, drag: bool, ballistic_coefficient: f64) -> Self {
        let config = ForceConfig {
            j2,
            drag,
            ballistic_coefficient,
        };
        PyPropagator {
            inner: NumericalPropagator::new(config),
        }
    }

    /// Propagate from MeanElements for duration_s seconds.
    ///
    /// Returns list of [epoch, x, y, z, vx, vy, vz] arrays.
    fn propagate_elements(&self, elements: &PyMeanElements, duration_s: f64) -> Vec<Vec<f64>> {
        let sv0 = StateVector::from_mean_elements(&elements.inner);
        let states = self.inner.propagate(&sv0, duration_s);
        states
            .iter()
            .map(|s| vec![s.epoch, s.r[0], s.r[1], s.r[2], s.v[0], s.v[1], s.v[2]])
            .collect()
    }

    /// Propagate from a state vector [x,y,z,vx,vy,vz] at given epoch.
    ///
    /// Returns list of [epoch, x, y, z, vx, vy, vz] arrays.
    fn propagate_state(
        &self,
        state: Vec<f64>,
        epoch: f64,
        duration_s: f64,
    ) -> PyResult<Vec<Vec<f64>>> {
        if state.len() != 6 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "State vector must have 6 elements [x,y,z,vx,vy,vz]",
            ));
        }
        let sv0 = StateVector {
            r: [state[0], state[1], state[2]],
            v: [state[3], state[4], state[5]],
            epoch,
        };
        let states = self.inner.propagate(&sv0, duration_s);
        Ok(states
            .iter()
            .map(|s| vec![s.epoch, s.r[0], s.r[1], s.r[2], s.v[0], s.v[1], s.v[2]])
            .collect())
    }
}

// Free functions
#[pyfunction]
fn annual_drag_dv(a: f64, bc: f64, rho: f64) -> f64 {
    maneuver::annual_drag_dv(a, bc, rho)
}

#[pyfunction]
fn correct_sma(chief: &PyMeanElements, roe: &PyROE) -> (f64, f64, f64) {
    let m = maneuver::correct_sma(&chief.inner, &roe.inner, chief.inner.epoch);
    (m.dv_r * 1000.0, m.dv_t * 1000.0, m.dv_n * 1000.0)
}

#[pyfunction]
fn correct_inclination(chief: &PyMeanElements, roe: &PyROE) -> (f64, f64, f64) {
    let m = maneuver::correct_inclination(&chief.inner, &roe.inner, chief.inner.epoch);
    (m.dv_r * 1000.0, m.dv_t * 1000.0, m.dv_n * 1000.0)
}

// Module registration
pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMeanElements>()?;
    m.add_class::<PyROE>()?;
    m.add_class::<PySlotBox>()?;
    m.add_class::<PyTle>()?;
    m.add_class::<PyPropagator>()?;
    m.add_function(wrap_pyfunction!(annual_drag_dv, m)?)?;
    m.add_function(wrap_pyfunction!(correct_sma, m)?)?;
    m.add_function(wrap_pyfunction!(correct_inclination, m)?)?;
    Ok(())
}
