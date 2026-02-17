//! Physical and astrodynamic constants.
//! TODO: Confirm and add more as needed 
//  (e.g., atmospheric model parameters, solar radiation pressure constants, etc.)


/// Earth gravitational parameter (km³/s²) — WGS84
pub const MU_EARTH: f64 = 398600.4418;

/// Earth equatorial radius (km) — WGS84
pub const R_EARTH: f64 = 6378.137;

/// Earth J2 zonal harmonic — WGS84/EGM96
pub const J2: f64 = 1.08262668e-3;

/// Earth rotation rate (rad/s)
pub const OMEGA_EARTH: f64 = 7.2921159e-5;

/// Seconds per sidereal day
pub const SIDEREAL_DAY: f64 = 86164.0905;

/// Seconds per solar day
pub const SOLAR_DAY: f64 = 86400.0;

/// Two pi
pub const TAU: f64 = std::f64::consts::TAU;

/// Degrees to radians
pub const DEG2RAD: f64 = std::f64::consts::PI / 180.0;

/// Radians to degrees
pub const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;
