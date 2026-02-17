//! Maneuver planning for constellation station-keeping.
//!
//! Computes delta-v for common maintenance maneuvers:
//! - In-plane (along-track) for SMA correction / drift control
//! - Cross-track for inclination correction
//! - Combined RAAN phasing via differential drag or altitude offset

use serde::{Deserialize, Serialize};
use crate::constants::*;
use crate::elements::{MeanElements, RelativeOrbitalElements};

/// An impulsive maneuver.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Maneuver {
    /// Delta-v in the radial direction (km/s). Positive = outward.
    pub dv_r: f64,
    /// Delta-v in the along-track direction (km/s). Positive = prograde.
    pub dv_t: f64,
    /// Delta-v in the cross-track direction (km/s). Positive = toward north.
    pub dv_n: f64,
    /// Maneuver epoch (seconds since reference).
    pub epoch: f64,
    /// Description of the maneuver purpose.
    pub description: &'static str,
}

impl Maneuver {
    /// Total delta-v magnitude (km/s).
    pub fn magnitude(&self) -> f64 {
        (self.dv_r.powi(2) + self.dv_t.powi(2) + self.dv_n.powi(2)).sqrt()
    }

    /// Total delta-v in m/s.
    pub fn magnitude_ms(&self) -> f64 {
        self.magnitude() * 1000.0
    }
}

/// Compute an along-track maneuver to correct relative semi-major axis.
///
/// For a Hohmann-like correction of δa:
///   Δv_t ≈ (n * a / 2) * δa
///
/// This is a tangential burn that changes the semi-major axis.
pub fn correct_sma(chief: &MeanElements, roe: &RelativeOrbitalElements, epoch: f64) -> Maneuver {
    let n = chief.mean_motion_j2();
    let v_circ = n * chief.a; // circular velocity approximation

    // δa is normalized by a_chief, so actual Δa = δa * a
    // Δv ≈ v/2 * Δa/a = v/2 * δa
    let dv_t = v_circ / 2.0 * (-roe.da); // negative to correct the offset

    Maneuver {
        dv_r: 0.0,
        dv_t,
        dv_n: 0.0,
        epoch,
        description: "SMA correction (along-track)",
    }
}

/// Compute a cross-track maneuver to correct relative inclination.
///
/// Single-impulse inclination change at the ascending/descending node:
///   Δv_n ≈ v * Δi  (for small Δi)
pub fn correct_inclination(
    chief: &MeanElements,
    roe: &RelativeOrbitalElements,
    epoch: f64,
) -> Maneuver {
    let v_circ = chief.mean_motion_j2() * chief.a;

    // Cross-track impulse to cancel δix (applied at node)
    let dv_n = v_circ * (-roe.dix);

    Maneuver {
        dv_r: 0.0,
        dv_t: 0.0,
        dv_n,
        epoch,
        description: "Inclination correction (cross-track)",
    }
}

/// Estimate the annual delta-v budget for drag makeup maneuvers.
///
/// Very simplified: given an estimated ballistic coefficient and
/// atmospheric density at altitude, compute the drag acceleration
/// and the delta-v needed per year to maintain altitude.
///
/// # Arguments
/// * `a` - Semi-major axis (km)
/// * `bc` - Ballistic coefficient m/(C_D * A) in kg/m² (typical: 20-100)
/// * `rho` - Atmospheric density at altitude (kg/m³)
///           (e.g., ~1e-13 at 500km, ~1e-12 at 400km — varies hugely with solar activity)
pub fn annual_drag_dv(a: f64, bc: f64, rho: f64) -> f64 {
    let v = (MU_EARTH / a).sqrt() * 1000.0; // m/s
    let a_drag = rho * v.powi(2) / (2.0 * bc); // m/s²
    let dv_per_year = a_drag * 365.25 * SOLAR_DAY; // m/s/year
    dv_per_year
}

/// Plan a complete station-keeping cycle for one satellite.
///
/// Given current ROEs and a slot box, determine which maneuvers are needed.
/// Returns a prioritized list of maneuvers.
pub fn plan_sk_cycle(
    chief: &MeanElements,
    roe: &RelativeOrbitalElements,
    epoch: f64,
) -> Vec<Maneuver> {
    let mut maneuvers = Vec::new();

    // Threshold: only maneuver if offset is significant
    let da_threshold = 1e-6;  // ~7m for LEO
    let di_threshold = 1e-6;  // ~0.2 arcsec

    if roe.da.abs() > da_threshold {
        maneuvers.push(correct_sma(chief, roe, epoch));
    }

    if roe.dix.abs() > di_threshold {
        maneuvers.push(correct_inclination(chief, roe, epoch));
    }

    maneuvers
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_sma_correction_direction() {
        let chief = MeanElements {
            a: R_EARTH + 550.0,
            e: 0.0,
            i: 53.0 * DEG2RAD,
            raan: 0.0, aop: 0.0, ma: 0.0, epoch: 0.0,
        };

        // Deputy is 100m too high — need a retrograde burn to lower
        let roe = RelativeOrbitalElements {
            da: 0.1 / chief.a, // 100m
            dlambda: 0.0, dex: 0.0, dey: 0.0, dix: 0.0, diy: 0.0,
        };

        let mnvr = correct_sma(&chief, &roe, 0.0);
        assert!(mnvr.dv_t < 0.0, "Should be retrograde to lower orbit");
        // For 100m at ~550km, expect ~few mm/s
        assert!(mnvr.magnitude_ms() < 100.0);
        assert!(mnvr.magnitude_ms() > 0.001);
    }

    #[test]
    fn test_annual_drag_dv_reasonable() {
        // 500km altitude, moderate solar activity
        let dv = annual_drag_dv(R_EARTH + 500.0, 50.0, 3e-13);
        // Should be on the order of a few m/s/year for this altitude
        assert!(dv > 0.1 && dv < 50.0, "Annual drag dv={dv} m/s seems wrong");
    }
}
