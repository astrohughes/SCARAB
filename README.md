# 🪲 SCARAB

**S**atellite **C**onstellation **A**utonomous **R**elative-motion **A**nalysis & **B**udgeting

A toolkit for satellite constellation station-keeping and maintenance — Rust core with Python bindings.

**Author:** Kyle Hughes ([@](https://github.com/huqhesy)astrohughes) — kyle.evan.hughes@gmail.com

---

## What this does

If you operate a constellation of satellites, you need to keep them in their assigned orbital slots. SCARAB handles the math:

- **TLE Parsing** — Ingest Two-Line Element sets from Space-Track, batch-parse entire catalogs, and convert to mean Keplerian elements.
- **Numerical Propagation** — Adaptive RK4 integrator with J2 and atmospheric drag force models. Propagate single satellites or entire constellations in parallel.
- **Relative Orbital Elements** — D'Amico quasi-nonsingular ROE formulation. Know exactly how far each satellite has drifted from its nominal slot.
- **Slot-Keeping Boxes** — Define control deadbands in ROE space (SMA, in-track, eccentricity, inclination). Detect violations and estimate time to boundary crossing.
- **Maneuver Planning** — Impulsive Δv for SMA correction, inclination maintenance, and annual drag makeup budgets.
- **Walker Constellations** — Generate Walker Delta and Star patterns with proper phasing.

The Rust core handles the numerics (fast enough for 1000+ satellite constellations). Python bindings via PyO3 make it usable from Jupyter notebooks and ops scripts.

## Quick Start

### Prerequisites

- [Rust toolchain](https://rustup.rs/) (1.75+)
- Python 3.9+
- [maturin](https://www.maturin.rs/)

### Build & Install

```bash
git clone https://github.com/YOUR_USERNAME/scarab.git
cd scarab

# Build Rust + install Python package
pip install maturin
maturin develop --release

# Run tests
cargo test
pytest tests/ -v

# Run the example
python examples/fleet_health.py
```

### Python Usage

```python
from scarab import TLE, MeanElements, ROE, SlotBox, Propagator

# ── Parse TLEs ──
tle = TLE.parse(line1, line2)
print(f"{tle.name} — {tle.altitude():.1f} km, {tle.inclination_deg:.2f}°")

# Batch parse from Space-Track download
tles = TLE.parse_batch(open("catalog.tle").read())

# ── Convert to mean elements ──
chief = tle.to_mean_elements()
print(f"RAAN drift: {chief.raan_rate_deg_day():.4f} °/day")

# ── Compute relative state ──
deputy = tles[1].to_mean_elements()
roe = ROE.from_elements(chief, deputy)
print(f"Relative SMA: {roe.da_meters(chief.a):.1f} m")
print(f"In-track: {roe.in_track_km(chief.a):.3f} km")

# ── Check slot box ──
box = SlotBox(da_meters=100.0, dlambda_km=10.0, de=5e-5,
              di_arcsec=20.0, a_chief=chief.a)
status = box.check(chief, roe)
print(f"Violated: {status['violated']}")

# ── Numerical propagation (J2 + drag) ──
prop = Propagator(j2=True, drag=True, ballistic_coefficient=40.0)
states = prop.propagate_elements(chief, duration_s=86400.0)
# states: list of [epoch, x, y, z, vx, vy, vz]
```

### Rust Usage

```rust
use scarab::tle::Tle;
use scarab::propagator::{NumericalPropagator, ForceConfig, StateVector};

let tle = Tle::parse(line1, line2).unwrap();
let elements = tle.to_mean_elements();

let prop = NumericalPropagator::new(ForceConfig::j2_drag(50.0));
let sv = StateVector::from_mean_elements(&elements);
let states = prop.propagate(&sv, 86400.0);
```

## Architecture

```
scarab/
├── src/
│   ├── lib.rs           # Crate root + PyO3 module entry
│   ├── constants.rs     # WGS84, J2, physical constants
│   ├── elements.rs      # Keplerian, mean, and relative orbital elements
│   ├── tle.rs           # TLE parser with validation and batch support
│   ├── propagator.rs    # J2 secular (analytical) + numerical (adaptive RK4)
│   ├── slotbox.rs       # Control box definition and violation detection
│   ├── maneuver.rs      # Δv computation for station-keeping
│   └── pybridge.rs      # PyO3 bindings
├── python/scarab/
│   ├── __init__.py      # Python package re-exports
│   └── walker.py        # Walker Delta/Star generators
├── examples/
│   └── fleet_health.py  # Full demo: TLE → propagate → ROE → maneuver
├── tests/
│   └── test_scarab.py   # Python integration tests
├── .github/workflows/
│   └── ci.yml           # GitHub Actions: cargo test + pytest
├── Cargo.toml
└── pyproject.toml
```

## Roadmap

- [X] **Phase 1:** Mean elements, ROEs, slot boxes, impulsive maneuvers
- [X] **Phase 2:** TLE parser, numerical propagation with J2 + drag
- [ ] **Phase 3:** Low-thrust maneuver planning (electric propulsion)
- [ ] **Phase 4:** RAAN phasing / constellation build-up planner
- [ ] **Phase 5:** CCSDS CDM/OEM import/export
- [ ] Dormand-Prince 7(8) with dense output (upgrade from RK4 adaptive)
- [ ] Osculating ↔ mean element conversion (Brouwer)
- [ ] Benchmarks & validation against Vallado test cases
- [ ] Jupyter notebook examples with plotly visualization
- [ ] CLI tool (`scarab check`, `scarab plan`, `scarab propagate`)

## References

- D'Amico, S. (2010). *Autonomous Formation Flying in Low Earth Orbit.* PhD Dissertation, TU Delft.
- Vallado, D. (2013). *Fundamentals of Astrodynamics and Applications*, 4th ed.
- Schaub, H. & Junkins, J. (2018). *Analytical Mechanics of Space Systems*, 4th ed.
- Kelso, T.S. — [CelesTrak TLE format documentation](https://celestrak.org/columns/v04n03/)

## License

MIT OR Apache-2.0
