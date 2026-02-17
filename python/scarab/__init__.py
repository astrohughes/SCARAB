"""
SCARAB — Satellite Constellation Autonomous Relative-motion Analysis & Budgeting

A Rust-powered toolkit for constellation station-keeping and maintenance.

Modules:
    - TLE: Parse Two-Line Element sets from Space-Track or other sources
    - MeanElements: Mean orbital elements with J2 secular propagation
    - ROE: Quasi-nonsingular relative orbital elements (D'Amico formulation)
    - SlotBox: Control box violation detection
    - Propagator: Numerical orbit propagation (RK4 adaptive with J2 + drag)
    - Maneuver planning: SMA correction, inclination correction, drag budgets

Quick start:
    >>> from scarab import TLE, MeanElements, ROE, SlotBox, Propagator
    >>>
    >>> # Parse a TLE
    >>> tle = TLE.parse(line1, line2)
    >>> print(f"{tle.name} at {tle.altitude():.1f} km")
    >>>
    >>> # Compute relative state
    >>> chief = tle.to_mean_elements()
    >>> roe = ROE.from_elements(chief, deputy)
    >>>
    >>> # Numerical propagation with J2 + drag
    >>> prop = Propagator(j2=True, drag=True, ballistic_coefficient=40.0)
    >>> states = prop.propagate_elements(chief, duration_s=86400.0)
"""

from scarab.scarab import (
    MeanElements,
    ROE,
    SlotBox,
    TLE,
    Propagator,
    annual_drag_dv,
    correct_sma,
    correct_inclination,
)

__all__ = [
    "MeanElements",
    "ROE",
    "SlotBox",
    "TLE",
    "Propagator",
    "annual_drag_dv",
    "correct_sma",
    "correct_inclination",
]

__version__ = "0.1.0"
