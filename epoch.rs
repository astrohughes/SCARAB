use serde::{Deserialize, Serialize};

use super::leapseconds::tai_utc_offset;
use super::timescale::TimeScale;
use super::{GPS_TAI_OFFSET_NS, TT_TAI_OFFSET_NS};

const NANOS_PER_SEC: i128 = 1_000_000_000;
const SECS_PER_DAY: i128 = 86_400;
const NANOS_PER_DAY: i128 = SECS_PER_DAY * NANOS_PER_SEC;

/// An instant in time, internally stored as TAI nanoseconds since J2000.
///
/// J2000 epoch = 2000-01-01T12:00:00.000 TAI
///
/// # Examples
/// ```
/// use astroframes_core::time::{Epoch, TimeScale};
///
/// let t = Epoch::from_utc(2024, 6, 15, 12, 0, 0.0);
/// let tai = t.to_tai_seconds();
/// let gps = t.to_timescale(TimeScale::GPS);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Epoch {
    /// TAI nanoseconds since J2000 (2000-01-01T12:00:00 TAI)
    tai_ns: i128,
}

impl Epoch {
    // ── Constructors ──────────────────────────────────────────────

    /// Create an Epoch from raw TAI nanoseconds since J2000.
    pub const fn from_tai_ns(tai_ns: i128) -> Self {
        Self { tai_ns }
    }

    /// Create an Epoch from TAI seconds since J2000.
    pub fn from_tai_seconds(secs: f64) -> Self {
        Self {
            tai_ns: (secs * NANOS_PER_SEC as f64) as i128,
        }
    }

    /// Create an Epoch from a UTC calendar date and time.
    ///
    /// Handles leap second offset lookup automatically.
    pub fn from_utc(year: i32, month: u8, day: u8, hour: u8, min: u8, sec: f64) -> Self {
        let days_since_j2000 = utc_calendar_to_j2000_days(year, month, day);
        let time_of_day_ns = (hour as i128) * 3_600 * NANOS_PER_SEC
            + (min as i128) * 60 * NANOS_PER_SEC
            + (sec * NANOS_PER_SEC as f64) as i128;

        // UTC instant in "UTC nanoseconds since J2000 noon"
        let utc_ns = days_since_j2000 * NANOS_PER_DAY - 12 * 3_600 * NANOS_PER_SEC + time_of_day_ns;

        // TAI = UTC + leap_seconds
        let leap_s = tai_utc_offset(year, month, day) as i128;
        let tai_ns = utc_ns + leap_s * NANOS_PER_SEC;

        Self { tai_ns }
    }

    /// Create an Epoch from GPS seconds since the GPS epoch (1980-01-06T00:00:00 UTC).
    pub fn from_gps_seconds(gps_secs: f64) -> Self {
        // GPS epoch in TAI ns since J2000:
        // 1980-01-06T00:00:00 UTC = 1980-01-06T00:00:19 TAI
        // Days from J2000 (2000-01-01) to 1980-01-06 = -7300 days (approx)
        let gps_epoch_tai_ns: i128 = {
            let days = utc_calendar_to_j2000_days(1980, 1, 6);
            days * NANOS_PER_DAY - 12 * 3_600 * NANOS_PER_SEC + GPS_TAI_OFFSET_NS
        };
        let tai_ns = gps_epoch_tai_ns + (gps_secs * NANOS_PER_SEC as f64) as i128;
        Self { tai_ns }
    }

    // ── Accessors ─────────────────────────────────────────────────

    /// Raw TAI nanoseconds since J2000.
    pub const fn as_tai_ns(&self) -> i128 {
        self.tai_ns
    }

    /// TAI seconds since J2000 (floating point).
    pub fn to_tai_seconds(&self) -> f64 {
        self.tai_ns as f64 / NANOS_PER_SEC as f64
    }

    /// Convert to seconds in the given timescale since J2000.
    pub fn to_timescale(&self, ts: TimeScale) -> f64 {
        match ts {
            TimeScale::TAI => self.to_tai_seconds(),
            TimeScale::TT => (self.tai_ns + TT_TAI_OFFSET_NS) as f64 / NANOS_PER_SEC as f64,
            TimeScale::GPS => (self.tai_ns - GPS_TAI_OFFSET_NS) as f64 / NANOS_PER_SEC as f64,
            TimeScale::UTC => {
                // TODO: proper inverse leap-second lookup for continuous conversion
                // For now, approximate using current leap second count
                let approx_secs = self.to_tai_seconds();
                // This needs a proper iterative approach for edge cases
                approx_secs - 37.0 // placeholder — will be replaced with proper lookup
            }
            TimeScale::TDB => {
                // TDB ≈ TT + periodic relativistic terms (Fairhead & Bretagnon)
                // For sub-millisecond work, TDB ≈ TT is sufficient
                // TODO: implement full Fairhead-Bretagnon series
                let tt_secs = (self.tai_ns + TT_TAI_OFFSET_NS) as f64 / NANOS_PER_SEC as f64;
                tt_secs // placeholder — periodic terms < 1.7ms
            }
            TimeScale::UT1 => {
                // UT1 = UTC + (UT1-UTC) from EOP data
                // TODO: requires EOP provider
                unimplemented!("UT1 requires EOP data — use Epoch::to_ut1_with_eop()")
            }
        }
    }

    /// Julian date in the given timescale.
    pub fn to_jd(&self, ts: TimeScale) -> f64 {
        let secs = self.to_timescale(ts);
        // J2000 = JD 2451545.0
        2_451_545.0 + secs / (SECS_PER_DAY as f64)
    }

    /// Modified Julian Date in the given timescale.
    pub fn to_mjd(&self, ts: TimeScale) -> f64 {
        self.to_jd(ts) - 2_400_000.5
    }

    // ── Arithmetic ────────────────────────────────────────────────

    /// Duration between two epochs in seconds.
    pub fn duration_since(&self, other: &Epoch) -> f64 {
        (self.tai_ns - other.tai_ns) as f64 / NANOS_PER_SEC as f64
    }

    /// Add seconds to this epoch.
    pub fn add_seconds(&self, secs: f64) -> Self {
        Self {
            tai_ns: self.tai_ns + (secs * NANOS_PER_SEC as f64) as i128,
        }
    }
}

impl std::fmt::Display for Epoch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Display as TAI seconds since J2000 for now
        // TODO: calendar date formatting
        write!(f, "Epoch({:.9} TAI s past J2000)", self.to_tai_seconds())
    }
}

impl std::ops::Sub for Epoch {
    type Output = f64;
    /// Returns duration in seconds (TAI).
    fn sub(self, rhs: Self) -> f64 {
        self.duration_since(&rhs)
    }
}

// ── Helper: Gregorian calendar to J2000 days ───────────────────────

/// Convert a UTC calendar date to days since J2000 (2000-01-01).
/// Uses the standard algorithm valid for dates after 1582-10-15.
fn utc_calendar_to_j2000_days(year: i32, month: u8, day: u8) -> i128 {
    // Algorithm: Meeus, Astronomical Algorithms, Ch.7
    let y = if month <= 2 { year - 1 } else { year } as i128;
    let m = if month <= 2 {
        month as i128 + 12
    } else {
        month as i128
    };
    let d = day as i128;

    let a = y / 100;
    let b = 2 - a + a / 4;

    let jd_noon = (365.25 * (y + 4716) as f64).floor() as i128
        + (30.6001 * (m + 1) as f64).floor() as i128
        + d
        + b
        - 1524;

    // J2000 = JD 2451545 (2000-01-01 12:00:00)
    jd_noon - 2_451_545
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_j2000_epoch_is_zero() {
        // 2000-01-01 12:00:00 TAI should be 0 + 32s of leap seconds at that date
        let t = Epoch::from_utc(2000, 1, 1, 12, 0, 0.0);
        // TAI = UTC + 32s (leap seconds at 2000-01-01)
        let expected_tai_ns = 32 * NANOS_PER_SEC;
        assert_eq!(t.as_tai_ns(), expected_tai_ns);
    }

    #[test]
    fn test_tai_tt_offset() {
        let t = Epoch::from_tai_ns(0); // J2000 in TAI
        let tt_secs = t.to_timescale(TimeScale::TT);
        // TT = TAI + 32.184s
        assert!((tt_secs - 32.184).abs() < 1e-9);
    }

    #[test]
    fn test_epoch_arithmetic() {
        let t1 = Epoch::from_tai_ns(0);
        let t2 = t1.add_seconds(100.0);
        assert!((t2 - t1 - 100.0).abs() < 1e-9);
    }

    #[test]
    fn test_julian_date() {
        let t = Epoch::from_tai_ns(0);
        let jd = t.to_jd(TimeScale::TAI);
        assert!((jd - 2_451_545.0).abs() < 1e-6);
    }
}
