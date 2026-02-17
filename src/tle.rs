//! Two-Line Element (TLE) set parser.
//!
//! Parses standard NORAD/Space-Track TLE format (2-line and 3-line with name).
//! Supports batch parsing of multi-TLE files and conversion to mean Keplerian elements.
//!
//! # TLE Format Reference
//! ```text
//! Line 0 (optional): Satellite Name (up to 24 chars)
//! Line 1: 1 NNNNNC NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN
//! Line 2: 2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN
//! ```
//!
//! # Example
//! ```
//! use scarab::tle::Tle;
//!
//! let line1 = "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9003";
//! let line2 = "2 25544  51.6400 208.5000 0007417  68.0000 292.1000 15.49560000400000";
//!
//! let tle = Tle::parse(line1, line2).unwrap();
//! assert_eq!(tle.norad_id, 25544);
//! ```

use serde::{Deserialize, Serialize};
use thiserror::Error;
use crate::constants::*;
use crate::elements::{KeplerianElements, MeanElements};

/// TLE parsing errors.
#[derive(Error, Debug)]
pub enum TleError {
    #[error("Line 1 must start with '1', got '{0}'")]
    InvalidLine1Start(char),

    #[error("Line 2 must start with '2', got '{0}'")]
    InvalidLine2Start(char),

    #[error("Line 1 length must be 69 characters, got {0}")]
    InvalidLine1Length(usize),

    #[error("Line 2 length must be 69 characters, got {0}")]
    InvalidLine2Length(usize),

    #[error("NORAD IDs don't match between lines: {0} vs {1}")]
    NoradIdMismatch(u32, u32),

    #[error("Checksum failed on line {line}: expected {expected}, computed {computed}")]
    ChecksumFailed {
        line: u8,
        expected: u8,
        computed: u8,
    },

    #[error("Failed to parse field '{field}': {source}")]
    ParseField {
        field: &'static str,
        source: std::num::ParseFloatError,
    },

    #[error("Failed to parse integer field '{field}': {source}")]
    ParseIntField {
        field: &'static str,
        source: std::num::ParseIntError,
    },

    #[error("Failed to parse implied-decimal field '{0}'")]
    ImpliedDecimal(String),

    #[error("No TLEs found in input")]
    Empty,
}

/// A parsed Two-Line Element set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tle {
    /// Satellite name (from line 0, if present).
    pub name: Option<String>,
    /// NORAD catalog number.
    pub norad_id: u32,
    /// International designator (launch year, launch number, piece).
    pub intl_designator: String,
    /// Classification (U=unclassified, C=classified, S=secret).
    pub classification: char,
    /// Epoch year (full 4-digit year).
    pub epoch_year: u16,
    /// Epoch day of year (fractional).
    pub epoch_day: f64,
    /// First derivative of mean motion (rev/day²) / 2.
    pub mean_motion_dot: f64,
    /// Second derivative of mean motion (rev/day³) / 6.
    pub mean_motion_ddot: f64,
    /// B* drag term (1/Earth radii).
    pub bstar: f64,
    /// Ephemeris type (usually 0).
    pub ephemeris_type: u8,
    /// Element set number.
    pub element_set: u16,
    /// Inclination (degrees).
    pub inclination_deg: f64,
    /// Right ascension of ascending node (degrees).
    pub raan_deg: f64,
    /// Eccentricity (dimensionless).
    pub eccentricity: f64,
    /// Argument of perigee (degrees).
    pub arg_perigee_deg: f64,
    /// Mean anomaly (degrees).
    pub mean_anomaly_deg: f64,
    /// Mean motion (revolutions per day).
    pub mean_motion_rev_day: f64,
    /// Revolution number at epoch.
    pub rev_number: u32,
}

impl Tle {
    /// Parse a TLE from two lines (without satellite name).
    pub fn parse(line1: &str, line2: &str) -> Result<Self, TleError> {
        Self::parse_with_name(None, line1, line2)
    }

    /// Parse a TLE from three lines (with satellite name on line 0).
    pub fn parse_3line(line0: &str, line1: &str, line2: &str) -> Result<Self, TleError> {
        let name = line0.trim().to_string();
        Self::parse_with_name(Some(name), line1, line2)
    }

    /// Parse with optional name.
    fn parse_with_name(name: Option<String>, line1: &str, line2: &str) -> Result<Self, TleError> {
        let line1 = line1.trim_end();
        let line2 = line2.trim_end();

        // Validate lengths (allow trailing spaces, pad if short)
        let l1: String = format!("{:<69}", line1);
        let l2: String = format!("{:<69}", line2);

        if l1.len() < 69 {
            return Err(TleError::InvalidLine1Length(line1.len()));
        }
        if l2.len() < 69 {
            return Err(TleError::InvalidLine2Length(line2.len()));
        }

        // Validate line numbers
        let c1 = l1.chars().next().unwrap();
        let c2 = l2.chars().next().unwrap();
        if c1 != '1' {
            return Err(TleError::InvalidLine1Start(c1));
        }
        if c2 != '2' {
            return Err(TleError::InvalidLine2Start(c2));
        }

        // Verify checksums
        let cs1_expected = parse_digit(l1.as_bytes()[68])?;
        let cs1_computed = compute_checksum(&l1[..68]);
        if cs1_expected != cs1_computed {
            return Err(TleError::ChecksumFailed {
                line: 1,
                expected: cs1_expected,
                computed: cs1_computed,
            });
        }

        let cs2_expected = parse_digit(l2.as_bytes()[68])?;
        let cs2_computed = compute_checksum(&l2[..68]);
        if cs2_expected != cs2_computed {
            return Err(TleError::ChecksumFailed {
                line: 2,
                expected: cs2_expected,
                computed: cs2_computed,
            });
        }

        // ── Parse Line 1 ──
        let norad_id_1 = l1[2..7].trim().parse::<u32>().map_err(|e| TleError::ParseIntField {
            field: "norad_id (line 1)",
            source: e,
        })?;

        let classification = l1.as_bytes()[7] as char;
        let intl_designator = l1[9..17].trim().to_string();

        let epoch_year_2d = l1[18..20].trim().parse::<u16>().map_err(|e| TleError::ParseIntField {
            field: "epoch_year",
            source: e,
        })?;
        let epoch_year = if epoch_year_2d >= 57 {
            1900 + epoch_year_2d
        } else {
            2000 + epoch_year_2d
        };

        let epoch_day = l1[20..32].trim().parse::<f64>().map_err(|e| TleError::ParseField {
            field: "epoch_day",
            source: e,
        })?;

        let mean_motion_dot = l1[33..43].trim().parse::<f64>().map_err(|e| TleError::ParseField {
            field: "mean_motion_dot",
            source: e,
        })?;

        let mean_motion_ddot = parse_implied_decimal(&l1[44..52])?;
        let bstar = parse_implied_decimal(&l1[53..61])?;

        let ephemeris_type = l1[62..63].trim().parse::<u8>().unwrap_or(0);
        let element_set = l1[64..68].trim().parse::<u16>().unwrap_or(0);

        // ── Parse Line 2 ──
        let norad_id_2 = l2[2..7].trim().parse::<u32>().map_err(|e| TleError::ParseIntField {
            field: "norad_id (line 2)",
            source: e,
        })?;

        if norad_id_1 != norad_id_2 {
            return Err(TleError::NoradIdMismatch(norad_id_1, norad_id_2));
        }

        let inclination_deg = l2[8..16].trim().parse::<f64>().map_err(|e| TleError::ParseField {
            field: "inclination",
            source: e,
        })?;

        let raan_deg = l2[17..25].trim().parse::<f64>().map_err(|e| TleError::ParseField {
            field: "raan",
            source: e,
        })?;

        // Eccentricity has implied leading decimal point
        let ecc_str = format!("0.{}", l2[26..33].trim());
        let eccentricity = ecc_str.parse::<f64>().map_err(|e| TleError::ParseField {
            field: "eccentricity",
            source: e,
        })?;

        let arg_perigee_deg = l2[34..42].trim().parse::<f64>().map_err(|e| TleError::ParseField {
            field: "arg_perigee",
            source: e,
        })?;

        let mean_anomaly_deg = l2[42..51].trim().parse::<f64>().map_err(|e| TleError::ParseField {
            field: "mean_anomaly",
            source: e,
        })?;

        let mean_motion_rev_day = l2[52..63].trim().parse::<f64>().map_err(|e| TleError::ParseField {
            field: "mean_motion",
            source: e,
        })?;

        let rev_number = l2[63..68].trim().parse::<u32>().unwrap_or(0);

        Ok(Tle {
            name,
            norad_id: norad_id_1,
            intl_designator,
            classification,
            epoch_year,
            epoch_day,
            mean_motion_dot,
            mean_motion_ddot,
            bstar,
            ephemeris_type,
            element_set,
            inclination_deg,
            raan_deg,
            eccentricity,
            arg_perigee_deg,
            mean_anomaly_deg,
            mean_motion_rev_day,
            rev_number,
        })
    }

    /// Parse a string containing multiple TLEs (2-line or 3-line format).
    ///
    /// Handles mixed formats: lines starting with '1' begin a 2-line TLE,
    /// other non-empty lines are treated as satellite names (line 0).
    pub fn parse_batch(input: &str) -> Result<Vec<Self>, TleError> {
        let lines: Vec<&str> = input
            .lines()
            .map(|l| l.trim_end())
            .filter(|l| !l.is_empty())
            .collect();

        if lines.is_empty() {
            return Err(TleError::Empty);
        }

        let mut tles = Vec::new();
        let mut i = 0;

        while i < lines.len() {
            if lines[i].starts_with('1') && i + 1 < lines.len() && lines[i + 1].starts_with('2') {
                // 2-line TLE
                tles.push(Tle::parse(lines[i], lines[i + 1])?);
                i += 2;
            } else if i + 2 < lines.len()
                && lines[i + 1].starts_with('1')
                && lines[i + 2].starts_with('2')
            {
                // 3-line TLE (line 0 is name)
                tles.push(Tle::parse_3line(lines[i], lines[i + 1], lines[i + 2])?);
                i += 3;
            } else {
                // Skip unrecognized lines
                i += 1;
            }
        }

        if tles.is_empty() {
            return Err(TleError::Empty);
        }

        Ok(tles)
    }

    /// Semi-major axis derived from mean motion (km).
    ///
    /// Uses Kepler's third law: a = (μ / n²)^(1/3)
    /// where n is in rad/s.
    pub fn semi_major_axis(&self) -> f64 {
        let n_rad_s = self.mean_motion_rev_day * TAU / SOLAR_DAY;
        (MU_EARTH / n_rad_s.powi(2)).powf(1.0 / 3.0)
    }

    /// Altitude above Earth's surface (km), assuming circular orbit.
    pub fn altitude(&self) -> f64 {
        self.semi_major_axis() - R_EARTH
    }

    /// Orbital period (seconds).
    pub fn period(&self) -> f64 {
        SOLAR_DAY / self.mean_motion_rev_day
    }

    /// Convert to Keplerian elements (osculating, from TLE mean motion).
    pub fn to_keplerian(&self) -> KeplerianElements {
        KeplerianElements::from_degrees(
            self.semi_major_axis(),
            self.eccentricity,
            self.inclination_deg,
            self.raan_deg,
            self.arg_perigee_deg,
            self.mean_anomaly_deg,
        )
    }

    /// Convert to mean elements at TLE epoch.
    ///
    /// Epoch is expressed as seconds since J2000 (2000-01-01 12:00:00 TT).
    pub fn to_mean_elements(&self) -> MeanElements {
        let epoch_j2000 = self.epoch_j2000_seconds();
        self.to_keplerian().to_mean(epoch_j2000)
    }

    /// TLE epoch as seconds since J2000.
    pub fn epoch_j2000_seconds(&self) -> f64 {
        // Days from J2000 (2000-01-01 12:00 TT) to start of epoch year
        let year = self.epoch_year as i32;
        let days_from_j2000 = days_since_j2000_jan1(year) + self.epoch_day - 1.0;
        // -0.5 because J2000 is at noon, Jan 1 day-of-year starts at midnight
        (days_from_j2000 - 0.5) * SOLAR_DAY
    }

    /// Check if this is a sun-synchronous orbit (RAAN drift ≈ +0.9856°/day).
    pub fn is_sun_synchronous(&self) -> bool {
        let mean_elems = self.to_mean_elements();
        let raan_rate = mean_elems.raan_rate() * RAD2DEG * SOLAR_DAY;
        (raan_rate - 0.9856).abs() < 0.05
    }

    /// Ballistic coefficient derived from B* (kg/m²).
    ///
    /// B* = (C_D * A) / (2 * m) * ρ₀
    /// where ρ₀ = 0.15696615 kg/m²/R_E (reference density).
    /// BC = m / (C_D * A) = ρ₀ / (2 * B*)
    pub fn ballistic_coefficient(&self) -> Option<f64> {
        if self.bstar.abs() < 1e-20 {
            return None;
        }
        let rho0 = 2.461e-5; // kg/m²/km — reference atmospheric density * R_E
        Some(rho0 / (2.0 * self.bstar.abs()))
    }
}

impl std::fmt::Display for Tle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} (NORAD {}) — {:.1} km, {:.1}° inc, {:.4} ecc, {:.2} rev/day",
            self.name.as_deref().unwrap_or("UNKNOWN"),
            self.norad_id,
            self.altitude(),
            self.inclination_deg,
            self.eccentricity,
            self.mean_motion_rev_day,
        )
    }
}

/// Parse the TLE "implied decimal" format: " NNNNN-N" → float.
///
/// Examples: " 16538-4" → 0.16538e-4, "-11606-4" → -0.11606e-4
fn parse_implied_decimal(s: &str) -> Result<f64, TleError> {
    let s = s.trim();
    if s.is_empty() || s == "00000-0" || s == " 00000-0" || s == "00000+0" {
        return Ok(0.0);
    }

    // Find the exponent sign (last + or - that isn't the leading sign)
    let bytes = s.as_bytes();
    let mut exp_pos = None;

    for i in (1..bytes.len()).rev() {
        if bytes[i] == b'+' || bytes[i] == b'-' {
            exp_pos = Some(i);
            break;
        }
    }

    match exp_pos {
        Some(pos) => {
            let mantissa_str = &s[..pos];
            let exp_str = &s[pos..];

            // Add implied leading "0."
            let sign = if mantissa_str.starts_with('-') { "-" } else { "" };
            let digits = mantissa_str.trim_start_matches(['+', '-', ' ']);

            let full = format!("{}0.{}e{}", sign, digits, exp_str);
            full.parse::<f64>()
                .map_err(|_| TleError::ImpliedDecimal(s.to_string()))
        }
        None => {
            // No exponent — try parsing with implied leading "0."
            let sign = if s.starts_with('-') { "-" } else { "" };
            let digits = s.trim_start_matches(['+', '-', ' ']);
            let full = format!("{}0.{}", sign, digits);
            full.parse::<f64>()
                .map_err(|_| TleError::ImpliedDecimal(s.to_string()))
        }
    }
}

/// Parse a single ASCII digit, or return 0 for space.
fn parse_digit(b: u8) -> Result<u8, TleError> {
    match b {
        b'0'..=b'9' => Ok(b - b'0'),
        b' ' => Ok(0),
        _ => Ok(0), // Be lenient
    }
}

/// Compute TLE checksum (mod-10 of sum of digits, '-' counts as 1).
fn compute_checksum(line: &str) -> u8 {
    let sum: u32 = line
        .bytes()
        .map(|b| match b {
            b'0'..=b'9' => (b - b'0') as u32,
            b'-' => 1,
            _ => 0,
        })
        .sum();
    (sum % 10) as u8
}

/// Days from J2000 epoch (2000-01-01) to January 1 of the given year.
fn days_since_j2000_jan1(year: i32) -> f64 {
    // Simple calculation — good enough for TLE epoch conversion
    let y = year - 2000;
    let leap_days = if y > 0 {
        (y - 1) / 4 - (y - 1) / 100 + (y - 1) / 400 + 1
    } else if y < 0 {
        y / 4 - y / 100 + y / 400
    } else {
        0
    };
    (y * 365 + leap_days) as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const ISS_LINE1: &str = "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9003";
    const ISS_LINE2: &str = "2 25544  51.6400 208.5000 0007417  68.0000 292.1000 15.49560000400000";

    #[test]
    fn test_parse_iss() {
        let tle = Tle::parse(ISS_LINE1, ISS_LINE2).unwrap();
        assert_eq!(tle.norad_id, 25544);
        assert_eq!(tle.epoch_year, 2024);
        assert_relative_eq!(tle.epoch_day, 1.5, epsilon = 1e-8);
        assert_relative_eq!(tle.inclination_deg, 51.64, epsilon = 1e-4);
        assert_relative_eq!(tle.raan_deg, 208.5, epsilon = 1e-4);
        assert_relative_eq!(tle.eccentricity, 0.0007417, epsilon = 1e-8);
        assert_relative_eq!(tle.mean_motion_rev_day, 15.4956, epsilon = 1e-4);
    }

    #[test]
    fn test_iss_altitude() {
        let tle = Tle::parse(ISS_LINE1, ISS_LINE2).unwrap();
        let alt = tle.altitude();
        // ISS is roughly 410-420 km
        assert!(alt > 400.0 && alt < 430.0, "ISS altitude={alt} km");
    }

    #[test]
    fn test_parse_3line() {
        let name = "ISS (ZARYA)";
        let tle = Tle::parse_3line(name, ISS_LINE1, ISS_LINE2).unwrap();
        assert_eq!(tle.name.as_deref(), Some("ISS (ZARYA)"));
        assert_eq!(tle.norad_id, 25544);
    }

    #[test]
    fn test_parse_batch() {
        let input = format!(
            "ISS (ZARYA)\n{}\n{}\nHUBBLE\n1 20580U 90037B   24001.50000000  .00000764  00000-0  34340-4 0  9998\n2 20580  28.4700 100.2000 0002500 300.0000  60.0000 15.09000000400000\n",
            ISS_LINE1, ISS_LINE2
        );
        let tles = Tle::parse_batch(&input).unwrap();
        assert_eq!(tles.len(), 2);
        assert_eq!(tles[0].name.as_deref(), Some("ISS (ZARYA)"));
        assert_eq!(tles[1].name.as_deref(), Some("HUBBLE"));
    }

    #[test]
    fn test_implied_decimal() {
        assert_relative_eq!(parse_implied_decimal("10270-3").unwrap(), 0.10270e-3, epsilon = 1e-12);
        assert_relative_eq!(parse_implied_decimal("00000-0").unwrap(), 0.0, epsilon = 1e-15);
        assert_relative_eq!(parse_implied_decimal("-11606-4").unwrap(), -0.11606e-4, epsilon = 1e-12);
        assert_relative_eq!(parse_implied_decimal("16538-4").unwrap(), 0.16538e-4, epsilon = 1e-12);
    }

    #[test]
    fn test_checksum() {
        // ISS line 1 checksum should be 3
        let cs = compute_checksum(&ISS_LINE1[..68]);
        assert_eq!(cs, 3);
    }

    #[test]
    fn test_to_mean_elements() {
        let tle = Tle::parse(ISS_LINE1, ISS_LINE2).unwrap();
        let mean = tle.to_mean_elements();
        assert_relative_eq!(mean.i * RAD2DEG, 51.64, epsilon = 0.01);
        assert!(mean.a > R_EARTH + 400.0);
        assert!(mean.a < R_EARTH + 430.0);
    }

    #[test]
    fn test_sun_sync_detection() {
        // ~600km sun-sync orbit
        let line1 = "1 99999U 24001A   24001.50000000  .00000000  00000-0  00000-0 0  9991";
        let line2 = "2 99999  97.8000  45.0000 0012000  90.0000 270.0000 14.95000000 00014";
        let tle = Tle::parse(line1, line2).unwrap();
        assert!(tle.is_sun_synchronous());
    }

    #[test]
    fn test_non_sun_sync() {
        let tle = Tle::parse(ISS_LINE1, ISS_LINE2).unwrap();
        assert!(!tle.is_sun_synchronous());
    }
}
