//! Slot-keeping control box logic.
//!
//! Defines keep-in boxes in ROE space and detects violations.
//! When a satellite drifts outside its box, a maneuver is triggered.
use serde::{Deserialize, Serialize};
use crate::elements::{MeanElements, RelativeOrbitalElements};

/// Control box definition in ROE space.
///
/// Each limit is a half-width: the satellite must stay within ±limit
/// of the nominal (zero) ROE value. Units match the ROE convention
/// (da is dimensionless, dlambda/dix/diy in rad, dex/dey dimensionless).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SlotBox {
    /// Max |δa/a| (dimensionless). Typical: ~1e-5 to 1e-4.
    pub da_limit: f64,
    /// Max |δλ| (rad). Controls in-track deadband.
    pub dlambda_limit: f64,
    /// Max |δe| magnitude (dimensionless). Relative eccentricity limit.
    pub de_limit: f64,
    /// Max |δi| magnitude (rad). Relative inclination limit.
    pub di_limit: f64,
}

/// Which ROE components are currently violating the control box.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct BoxViolation {
    pub da_violated: bool,
    pub dlambda_violated: bool,
    pub de_violated: bool,
    pub di_violated: bool,
}

/// Status of a satellite relative to its slot box, with margins.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SlotStatus {
    /// Current ROE values.
    pub roe: RelativeOrbitalElements,
    /// Which limits are violated.
    pub violation: BoxViolation,
    /// Fraction of da limit used (0.0 = centered, 1.0 = at boundary).
    pub da_margin: f64,
    /// Fraction of dlambda limit used.
    pub dlambda_margin: f64,
    /// Fraction of de limit used.
    pub de_margin: f64,
    /// Fraction of di limit used.
    pub di_margin: f64,
    /// Estimated seconds until next violation (most restrictive axis).
    /// None if no drift detected or already violated.
    pub time_to_violation: Option<f64>,
}

impl SlotBox {
    /// Create a slot box from limits specified in more intuitive units.
    ///
    /// - `da_meters`: half-width in meters (converted using chief SMA)
    /// - `dlambda_km`: in-track deadband half-width in km
    /// - `de`: relative eccentricity magnitude limit
    /// - `di_arcsec`: relative inclination limit in arcseconds
    /// - `a_chief`: chief semi-major axis in km (for unit conversion)
    pub fn from_physical_units(
        da_meters: f64,
        dlambda_km: f64,
        de: f64,
        di_arcsec: f64,
        a_chief: f64,
    ) -> Self {
        SlotBox {
            da_limit: da_meters / (a_chief * 1000.0),
            dlambda_limit: dlambda_km / a_chief,
            de_limit: de,
            di_limit: di_arcsec * std::f64::consts::PI / (180.0 * 3600.0),
        }
    }

    /// Check a satellite's ROE state against this box.
    pub fn check(&self, roe: &RelativeOrbitalElements, chief: &MeanElements) -> SlotStatus {
        let de_mag = roe.de_magnitude();
        let di_mag = roe.di_magnitude();

        let violation = BoxViolation {
            da_violated: roe.da.abs() > self.da_limit,
            dlambda_violated: roe.dlambda.abs() > self.dlambda_limit,
            de_violated: de_mag > self.de_limit,
            di_violated: di_mag > self.di_limit,
        };

        let da_margin = if self.da_limit > 0.0 { roe.da.abs() / self.da_limit } else { 0.0 };
        let dlambda_margin = if self.dlambda_limit > 0.0 { roe.dlambda.abs() / self.dlambda_limit } else { 0.0 };
        let de_margin = if self.de_limit > 0.0 { de_mag / self.de_limit } else { 0.0 };
        let di_margin = if self.di_limit > 0.0 { di_mag / self.di_limit } else { 0.0 };

        // Estimate time to violation from in-track drift (dominant effect)
        let n = chief.mean_motion_j2();
        let time_to_violation = if !violation.dlambda_violated && roe.da.abs() > 1e-15 {
            roe.time_to_drift_limit(self.dlambda_limit * chief.a, chief.a, n)
        } else {
            None
        };

        SlotStatus {
            roe: *roe,
            violation,
            da_margin,
            dlambda_margin,
            de_margin,
            di_margin,
            time_to_violation,
        }
    }

    /// Check whether any limit is violated.
    pub fn is_violated(&self, roe: &RelativeOrbitalElements) -> bool {
        roe.da.abs() > self.da_limit
            || roe.dlambda.abs() > self.dlambda_limit
            || roe.de_magnitude() > self.de_limit
            || roe.di_magnitude() > self.di_limit
    }
}

impl BoxViolation {
    pub fn any(&self) -> bool {
        self.da_violated || self.dlambda_violated || self.de_violated || self.di_violated
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::*;

    fn make_chief() -> MeanElements {
        MeanElements {
            a: R_EARTH + 550.0,
            e: 0.0001,
            i: 53.0 * DEG2RAD,
            raan: 0.0,
            aop: 0.0,
            ma: 0.0,
            epoch: 0.0,
        }
    }

    #[test]
    fn test_nominal_satellite_passes() {
        let chief = make_chief();
        let sbox = SlotBox::from_physical_units(200.0, 20.0, 1e-4, 50.0, chief.a);

        // Satellite at nominal slot — zero ROE
        let roe = RelativeOrbitalElements {
            da: 0.0, dlambda: 0.0, dex: 0.0, dey: 0.0, dix: 0.0, diy: 0.0,
        };

        let status = sbox.check(&roe, &chief);
        assert!(!status.violation.any());
    }

    #[test]
    fn test_drifted_satellite_fails() {
        let chief = make_chief();
        let sbox = SlotBox::from_physical_units(200.0, 20.0, 1e-4, 50.0, chief.a);

        // Satellite 500m off in SMA — should violate da
        let roe = RelativeOrbitalElements {
            da: 500.0 / (chief.a * 1000.0),
            dlambda: 0.0,
            dex: 0.0, dey: 0.0, dix: 0.0, diy: 0.0,
        };

        let status = sbox.check(&roe, &chief);
        assert!(status.violation.da_violated);
        assert!(!status.violation.dlambda_violated);
    }
}
