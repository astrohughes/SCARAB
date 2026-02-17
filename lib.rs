//! # astroframes-core
//!
//! Fast, correct aerospace reference frames and time systems.
//!
//! ## Architecture
//!
//! - **time**: Epoch representation with nanosecond precision, multi-timescale support
//!   (TAI, UTC, TT, TDB, GPS, UT1), leap second handling, and EOP-aware conversions.
//!
//! - **frames**: Reference frame definitions and a frame graph that computes rotation
//!   chains automatically. Supports inertial (GCRF/J2000), Earth-fixed (ITRF),
//!   orbital (LVLH/RSW/VNC), topocentric (NED/ENU), and body-fixed frames.
//!
//! - **states**: CartesianState (pos + vel + epoch + frame) with ergonomic `.in_frame()`
//!   conversions.
//!
//! - **eop**: Earth Orientation Parameter loading and interpolation (IERS finals2000A).
pub mod eop;
pub mod frames;
pub mod states;
pub mod time;

// Re-exports for convenience
pub use frames::FrameId;
pub use states::CartesianState;
pub use time::{Epoch, TimeScale};
