//! Orbital element representations and conversions.
//!
//! Implements classical (Keplerian) orbital elements, mean elements
//! with J2 secular rates, and quasi-nonsingular relative orbital elements
//! (ROEs) following D'Amico's formulation.

use serde::{Deserialize, Serialize};
use crate::constants::*;

/// Classical (osculating) Keplerian orbital elements.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct KeplerianElements {
    /// Semi-major axis (km)
    pub a: f64,
    /// Eccentricity (dimensionless)
    pub e: f64,
    /// Inclination (rad)
    pub i: f64,
    /// Right ascension of ascending node (rad)
    pub raan: f64,
    /// Argument of perigee (rad)
    pub aop: f64,
    /// Mean anomaly (rad)
    pub ma: f64,
}

/// Mean orbital elements with J2 secular drift rates.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct MeanElements {
    /// Mean semi-major axis (km)
    pub a: f64,
    /// Mean eccentricity
    pub e: f64,
    /// Mean inclination (rad)
    pub i: f64,
    /// Mean RAAN (rad)
    pub raan: f64,
    /// Mean argument of perigee (rad)
    pub aop: f64,
    /// Mean anomaly (rad)
    pub ma: f64,
    /// Epoch (seconds since some reference — typically J2000 TDB or mission epoch)
    pub epoch: f64,
}

/// Quasi-nonsingular relative orbital elements (dimensionless, normalized by a_chief).
///
/// Following D'Amico (2010) formulation:
///   δa  = (a_dep - a_chief) / a_chief
///   δλ  = (u_dep - u_chief) + (Ω_dep - Ω_chief) cos(i_chief)
///   δex = e_dep cos(ω_dep) - e_chief cos(ω_chief)
///   δey = e_dep sin(ω_dep) - e_chief sin(ω_chief)
///   δix = i_dep - i_chief
///   δiy = (Ω_dep - Ω_chief) sin(i_chief)
///
/// where u = ω + M is the mean argument of latitude.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct RelativeOrbitalElements {
    /// Relative semi-major axis (δa/a)
    pub da: f64,
    /// Relative mean longitude (rad)
    pub dlambda: f64,
    /// Relative eccentricity vector x-component
    pub dex: f64,
    /// Relative eccentricity vector y-component
    pub dey: f64,
    /// Relative inclination vector x-component (rad)
    pub dix: f64,
    /// Relative inclination vector y-component (rad)
    pub diy: f64,
}

impl KeplerianElements {
    /// Create new Keplerian elements from values in degrees (convenience constructor).
    pub fn from_degrees(a: f64, e: f64, i_deg: f64, raan_deg: f64, aop_deg: f64, ma_deg: f64) -> Self {
        Self {
            a,
            e,
            i: i_deg * DEG2RAD,
            raan: raan_deg * DEG2RAD,
            aop: aop_deg * DEG2RAD,
            ma: ma_deg * DEG2RAD,
        }
    }

    /// Mean motion (rad/s).
    pub fn mean_motion(&self) -> f64 {
        (MU_EARTH / self.a.powi(3)).sqrt()
    }

    /// Orbital period (seconds).
    pub fn period(&self) -> f64 {
        TAU / self.mean_motion()
    }

    /// Mean argument of latitude u = ω + M (rad).
    pub fn mean_arg_latitude(&self) -> f64 {
        normalize_angle(self.aop + self.ma)
    }

    /// Convert to MeanElements at a given epoch.
    pub fn to_mean(self, epoch: f64) -> MeanElements {
        MeanElements {
            a: self.a,
            e: self.e,
            i: self.i,
            raan: self.raan,
            aop: self.aop,
            ma: self.ma,
            epoch,
        }
    }
}

impl MeanElements {
    /// J2 secular drift rate of RAAN (rad/s).
    ///
    /// dΩ/dt = -3/2 * n * J2 * (R_E/a)² * cos(i) / (1-e²)²
    pub fn raan_rate(&self) -> f64 {
        let n = (MU_EARTH / self.a.powi(3)).sqrt();
        let p = self.a * (1.0 - self.e.powi(2));
        -1.5 * n * J2 * (R_EARTH / p).powi(2) * self.i.cos()
    }

    /// J2 secular drift rate of argument of perigee (rad/s).
    ///
    /// dω/dt = 3/2 * n * J2 * (R_E/a)² * (2 - 5/2 sin²i) / (1-e²)²
    pub fn aop_rate(&self) -> f64 {
        let n = (MU_EARTH / self.a.powi(3)).sqrt();
        let p = self.a * (1.0 - self.e.powi(2));
        1.5 * n * J2 * (R_EARTH / p).powi(2) * (2.0 - 2.5 * self.i.sin().powi(2))
    }

    /// J2-perturbed mean motion (rad/s).
    ///
    /// n_J2 = n * [1 + 3/2 * J2 * (R_E/a)² * (1 - 3/2 sin²i) / (1-e²)^(3/2)]
    pub fn mean_motion_j2(&self) -> f64 {
        let n = (MU_EARTH / self.a.powi(3)).sqrt();
        let eta = (1.0 - self.e.powi(2)).sqrt();
        let ratio = R_EARTH / self.a;
        n * (1.0 + 1.5 * J2 * ratio.powi(2) * (1.0 - 1.5 * self.i.sin().powi(2)) / eta.powi(3))
    }

    /// Propagate mean elements forward by dt seconds (J2 secular only).
    pub fn propagate(&self, dt: f64) -> MeanElements {
        let n_j2 = self.mean_motion_j2();
        let raan_dot = self.raan_rate();
        let aop_dot = self.aop_rate();

        MeanElements {
            a: self.a,          // No secular change in a from J2
            e: self.e,          // No secular change in e from J2
            i: self.i,          // No secular change in i from J2
            raan: normalize_angle(self.raan + raan_dot * dt),
            aop: normalize_angle(self.aop + aop_dot * dt),
            ma: normalize_angle(self.ma + n_j2 * dt),
            epoch: self.epoch + dt,
        }
    }
}

impl RelativeOrbitalElements {
    /// Compute quasi-nonsingular ROEs of `deputy` relative to `chief`.
    pub fn from_absolute(chief: &MeanElements, deputy: &MeanElements) -> Self {
        let da = (deputy.a - chief.a) / chief.a;

        let u_chief = chief.aop + chief.ma;
        let u_deputy = deputy.aop + deputy.ma;
        let d_raan = deputy.raan - chief.raan;

        let dlambda = normalize_angle(u_deputy - u_chief + d_raan * chief.i.cos());

        let dex = deputy.e * deputy.aop.cos() - chief.e * chief.aop.cos();
        let dey = deputy.e * deputy.aop.sin() - chief.e * chief.aop.sin();

        let dix = deputy.i - chief.i;
        let diy = d_raan * chief.i.sin();

        RelativeOrbitalElements { da, dlambda, dex, dey, dix, diy }
    }

    /// Relative semi-major axis difference in meters.
    pub fn da_meters(&self, a_chief: f64) -> f64 {
        self.da * a_chief * 1000.0
    }

    /// Magnitude of the relative eccentricity vector.
    pub fn de_magnitude(&self) -> f64 {
        (self.dex.powi(2) + self.dey.powi(2)).sqrt()
    }

    /// Magnitude of the relative inclination vector (rad).
    pub fn di_magnitude(&self) -> f64 {
        (self.dix.powi(2) + self.diy.powi(2)).sqrt()
    }

    /// In-track separation at mean argument of latitude u=0 (km).
    /// Approximate: δλ * a_chief
    pub fn in_track_km(&self, a_chief: f64) -> f64 {
        self.dlambda * a_chief
    }

    /// Cross-track amplitude (km).
    /// Approximate: di_magnitude * a_chief
    pub fn cross_track_amplitude_km(&self, a_chief: f64) -> f64 {
        self.di_magnitude() * a_chief
    }

    /// Along-track drift rate (rad/s) due to relative semi-major axis.
    ///
    /// dδλ/dt ≈ -3/2 * n * δa/a  (from the HCW equations / secular J2)
    pub fn along_track_drift_rate(&self, n: f64) -> f64 {
        -1.5 * n * self.da
    }

    /// Time until in-track drift exceeds a given limit (seconds).
    /// Returns None if da ≈ 0 (no drift).
    pub fn time_to_drift_limit(&self, limit_km: f64, a_chief: f64, n: f64) -> Option<f64> {
        let drift_rate_km_s = self.along_track_drift_rate(n).abs() * a_chief;
        if drift_rate_km_s < 1e-15 {
            return None;
        }
        let remaining = limit_km - self.in_track_km(a_chief).abs();
        if remaining <= 0.0 {
            return Some(0.0);
        }
        Some(remaining / drift_rate_km_s)
    }
}

/// Normalize angle to [0, 2π).
pub fn normalize_angle(angle: f64) -> f64 {
    let a = angle % TAU;
    if a < 0.0 { a + TAU } else { a }
}

/// Normalize angle to [-π, π).
pub fn normalize_angle_pm(angle: f64) -> f64 {
    let a = normalize_angle(angle);
    if a >= std::f64::consts::PI { a - TAU } else { a }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_iss_like_orbit_raan_rate() {
        // ISS-like: 420 km altitude, 51.6° inc, circular
        let elems = MeanElements {
            a: R_EARTH + 420.0,
            e: 0.0001,
            i: 51.6 * DEG2RAD,
            raan: 0.0,
            aop: 0.0,
            ma: 0.0,
            epoch: 0.0,
        };

        let raan_rate_deg_day = elems.raan_rate() * RAD2DEG * SOLAR_DAY;
        // Expected: approximately -5°/day for ISS
        assert_relative_eq!(raan_rate_deg_day, -5.0, epsilon = 0.5);
    }

    #[test]
    fn test_sun_sync_raan_rate() {
        // Sun-synchronous: ~600 km, ~97.8° inc
        let elems = MeanElements {
            a: R_EARTH + 600.0,
            e: 0.001,
            i: 97.8 * DEG2RAD,
            raan: 0.0,
            aop: 0.0,
            ma: 0.0,
            epoch: 0.0,
        };

        let raan_rate_deg_day = elems.raan_rate() * RAD2DEG * SOLAR_DAY;
        // Sun-sync should be ~+0.9856°/day
        assert_relative_eq!(raan_rate_deg_day, 0.9856, epsilon = 0.05);
    }

    #[test]
    fn test_roe_zero_for_identical() {
        let chief = MeanElements {
            a: R_EARTH + 500.0,
            e: 0.001,
            i: 97.4 * DEG2RAD,
            raan: 30.0 * DEG2RAD,
            aop: 45.0 * DEG2RAD,
            ma: 120.0 * DEG2RAD,
            epoch: 0.0,
        };

        let roe = RelativeOrbitalElements::from_absolute(&chief, &chief);
        assert_relative_eq!(roe.da, 0.0, epsilon = 1e-14);
        assert_relative_eq!(roe.dlambda, 0.0, epsilon = 1e-14);
        assert_relative_eq!(roe.dex, 0.0, epsilon = 1e-14);
        assert_relative_eq!(roe.dey, 0.0, epsilon = 1e-14);
        assert_relative_eq!(roe.dix, 0.0, epsilon = 1e-14);
        assert_relative_eq!(roe.diy, 0.0, epsilon = 1e-14);
    }

    #[test]
    fn test_roe_altitude_difference() {
        let chief = MeanElements {
            a: R_EARTH + 500.0,
            e: 0.0,
            i: 97.4 * DEG2RAD,
            raan: 0.0,
            aop: 0.0,
            ma: 0.0,
            epoch: 0.0,
        };

        // Deputy 1 km higher
        let mut deputy = chief;
        deputy.a += 1.0;

        let roe = RelativeOrbitalElements::from_absolute(&chief, &deputy);
        let da_m = roe.da_meters(chief.a);
        assert_relative_eq!(da_m, 1000.0, epsilon = 1.0);
    }

    #[test]
    fn test_propagation_preserves_period() {
        let elems = MeanElements {
            a: R_EARTH + 550.0,
            e: 0.001,
            i: 53.0 * DEG2RAD,
            raan: 0.0,
            aop: 0.0,
            ma: 0.0,
            epoch: 0.0,
        };

        let one_period = TAU / (MU_EARTH / elems.a.powi(3)).sqrt();
        let prop = elems.propagate(one_period);

        // SMA, ecc, inc should be unchanged
        assert_relative_eq!(prop.a, elems.a, epsilon = 1e-10);
        assert_relative_eq!(prop.e, elems.e, epsilon = 1e-10);
        assert_relative_eq!(prop.i, elems.i, epsilon = 1e-10);
    }
}
