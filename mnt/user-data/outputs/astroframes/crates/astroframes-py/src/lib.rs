use pyo3::prelude::*;

use astroframes_core::frames::FrameId as RustFrameId;
use astroframes_core::time::TimeScale as RustTimeScale;
use astroframes_core::time::Epoch as RustEpoch;

// ── TimeScale ──────────────────────────────────────────────────────

#[pyclass(name = "TimeScale", eq)]
#[derive(Clone, PartialEq)]
pub struct PyTimeScale(RustTimeScale);

#[pymethods]
impl PyTimeScale {
    #[classattr] fn TAI() -> Self { PyTimeScale(RustTimeScale::TAI) }
    #[classattr] fn UTC() -> Self { PyTimeScale(RustTimeScale::UTC) }
    #[classattr] fn TT()  -> Self { PyTimeScale(RustTimeScale::TT) }
    #[classattr] fn TDB() -> Self { PyTimeScale(RustTimeScale::TDB) }
    #[classattr] fn GPS() -> Self { PyTimeScale(RustTimeScale::GPS) }
    #[classattr] fn UT1() -> Self { PyTimeScale(RustTimeScale::UT1) }

    fn __repr__(&self) -> String {
        format!("TimeScale.{}", self.0)
    }
}

// ── FrameId ────────────────────────────────────────────────────────

#[pyclass(name = "Frame", eq)]
#[derive(Clone, PartialEq)]
pub struct PyFrame(RustFrameId);

#[pymethods]
impl PyFrame {
    #[classattr] fn GCRF() -> Self { PyFrame(RustFrameId::GCRF) }
    #[classattr] fn ITRF() -> Self { PyFrame(RustFrameId::ITRF) }
    #[classattr] fn TEME() -> Self { PyFrame(RustFrameId::TEME) }
    #[classattr] fn LVLH() -> Self { PyFrame(RustFrameId::LVLH) }
    #[classattr] fn RSW()  -> Self { PyFrame(RustFrameId::RSW) }
    #[classattr] fn NED()  -> Self { PyFrame(RustFrameId::NED) }
    #[classattr] fn ENU()  -> Self { PyFrame(RustFrameId::ENU) }
    #[classattr] fn Body() -> Self { PyFrame(RustFrameId::Body) }

    fn __repr__(&self) -> String {
        format!("Frame.{}", self.0)
    }
}

// ── Epoch ──────────────────────────────────────────────────────────

#[pyclass(name = "Epoch")]
#[derive(Clone)]
pub struct PyEpoch(RustEpoch);

#[pymethods]
impl PyEpoch {
    /// Create from a UTC calendar date.
    ///
    /// >>> t = Epoch.from_utc(2024, 6, 15, 12, 0, 0.0)
    #[staticmethod]
    fn from_utc(year: i32, month: u8, day: u8, hour: u8, min: u8, sec: f64) -> Self {
        PyEpoch(RustEpoch::from_utc(year, month, day, hour, min, sec))
    }

    /// Create from TAI seconds since J2000.
    #[staticmethod]
    fn from_tai_seconds(secs: f64) -> Self {
        PyEpoch(RustEpoch::from_tai_seconds(secs))
    }

    /// Create from GPS seconds since GPS epoch.
    #[staticmethod]
    fn from_gps_seconds(secs: f64) -> Self {
        PyEpoch(RustEpoch::from_gps_seconds(secs))
    }

    /// TAI seconds since J2000.
    fn tai_seconds(&self) -> f64 {
        self.0.to_tai_seconds()
    }

    /// Seconds since J2000 in the given timescale.
    fn to_timescale(&self, ts: &PyTimeScale) -> f64 {
        self.0.to_timescale(ts.0)
    }

    /// Julian Date in the given timescale.
    fn jd(&self, ts: &PyTimeScale) -> f64 {
        self.0.to_jd(ts.0)
    }

    /// Modified Julian Date in the given timescale.
    fn mjd(&self, ts: &PyTimeScale) -> f64 {
        self.0.to_mjd(ts.0)
    }

    /// Add seconds, returning a new Epoch.
    fn add_seconds(&self, secs: f64) -> Self {
        PyEpoch(self.0.add_seconds(secs))
    }

    /// Duration in seconds between two epochs.
    fn __sub__(&self, other: &PyEpoch) -> f64 {
        self.0.duration_since(&other.0)
    }

    fn __repr__(&self) -> String {
        format!("{}", self.0)
    }
}

// ── CartesianState ─────────────────────────────────────────────────

#[pyclass(name = "CartesianState")]
#[derive(Clone)]
pub struct PyCartesianState {
    #[pyo3(get)]
    position: [f64; 3],
    #[pyo3(get)]
    velocity: [f64; 3],
    #[pyo3(get)]
    epoch: PyEpoch,
    #[pyo3(get)]
    frame: PyFrame,
}

#[pymethods]
impl PyCartesianState {
    #[new]
    fn new(position: [f64; 3], velocity: [f64; 3], epoch: PyEpoch, frame: PyFrame) -> Self {
        Self { position, velocity, epoch, frame }
    }

    /// Transform into a different reference frame.
    ///
    /// >>> state_ecef = state_eci.in_frame(Frame.ITRF)
    fn in_frame(&self, _target: &PyFrame) -> PyResult<PyCartesianState> {
        // TODO: wire up to FrameGraph
        Err(pyo3::exceptions::PyNotImplementedError::new_err(
            "Frame conversion not yet wired up — coming in v0.2",
        ))
    }

    /// Position magnitude [km].
    fn radius(&self) -> f64 {
        let [x, y, z] = self.position;
        (x * x + y * y + z * z).sqrt()
    }

    /// Velocity magnitude [km/s].
    fn speed(&self) -> f64 {
        let [vx, vy, vz] = self.velocity;
        (vx * vx + vy * vy + vz * vz).sqrt()
    }

    fn __repr__(&self) -> String {
        format!(
            "CartesianState(r=[{:.3}, {:.3}, {:.3}] km, v=[{:.6}, {:.6}, {:.6}] km/s, {}, {})",
            self.position[0], self.position[1], self.position[2],
            self.velocity[0], self.velocity[1], self.velocity[2],
            self.frame.__repr__(), self.epoch.__repr__(),
        )
    }
}

// ── Module ─────────────────────────────────────────────────────────

/// astroframes: fast, correct aerospace reference frames and time systems.
#[pymodule]
fn _astroframes(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTimeScale>()?;
    m.add_class::<PyFrame>()?;
    m.add_class::<PyEpoch>()?;
    m.add_class::<PyCartesianState>()?;
    Ok(())
}
