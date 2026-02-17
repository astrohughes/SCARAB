"""
SCARAB Example: TLE Ingestion → Numerical Propagation → Fleet Health Check

Demonstrates the full workflow:
1. Parse TLEs from a batch string (simulating Space-Track data)
2. Convert to mean elements
3. Numerically propagate with J2 + drag
4. Compute ROEs relative to constellation nominal slots
5. Check slot box violations and plan corrections
"""
from scarab import TLE, MeanElements, ROE, SlotBox, Propagator, correct_sma, annual_drag_dv
from scarab.walker import walker_delta

# ═══════════════════════════════════════════════════════════════
# 1. PARSE TLEs
# ═══════════════════════════════════════════════════════════════

# Example: a few satellites from a hypothetical constellation
tle_data = """SCARAB-01
1 55001U 23001A   24180.50000000  .00001200  00000-0  62000-4 0  9991
2 55001  53.0000   0.0000 0001500  90.0000   0.0000 15.05800000  0018
SCARAB-02
1 55002U 23001B   24180.50000000  .00001180  00000-0  61000-4 0  9996
2 55002  53.0100   0.0500 0001600  91.0000  90.0000 15.05850000  0017
SCARAB-03
1 55003U 23001C   24180.50000000  .00001250  00000-0  63000-4 0  9993
2 55003  52.9900  59.9800 0001400  89.0000 180.0000 15.05750000  0012
SCARAB-04
1 55004U 23001D   24180.50000000  .00001150  00000-0  60000-4 0  9998
2 55004  53.0200  60.0200 0001700  88.0000 270.0000 15.05900000  0015
"""

print("=" * 65)
print("  SCARAB — Satellite Constellation Autonomous")
print("  Relative-motion Analysis & Budgeting")
print("=" * 65)

tles = TLE.parse_batch(tle_data)
print(f"\nParsed {len(tles)} TLEs:")
print("-" * 65)

for tle in tles:
    sun_sync = " [SSO]" if tle.is_sun_synchronous() else ""
    print(f"  {tle.name:12s}  NORAD {tle.norad_id}  "
          f"alt={tle.altitude():6.1f} km  "
          f"inc={tle.inclination_deg:6.2f}°  "
          f"e={tle.eccentricity:.5f}{sun_sync}")

# ═══════════════════════════════════════════════════════════════
# 2. NUMERICAL PROPAGATION
# ═══════════════════════════════════════════════════════════════

print(f"\n{'=' * 65}")
print("NUMERICAL PROPAGATION (J2 + Drag)")
print(f"{'=' * 65}")

prop = Propagator(j2=True, drag=True, ballistic_coefficient=40.0)

# Propagate first satellite for 1 day
tle0 = tles[0]
elem0 = tle0.to_mean_elements()
states = prop.propagate_elements(elem0, duration_s=86400.0)

print(f"\n  Propagated {tle0.name} for 1 day ({len(states)} integration steps)")
sv0 = states[0]
svf = states[-1]

r0_mag = (sv0[1]**2 + sv0[2]**2 + sv0[3]**2)**0.5
rf_mag = (svf[1]**2 + svf[2]**2 + svf[3]**2)**0.5

print(f"  Initial altitude: {r0_mag - 6378.137:.2f} km")
print(f"  Final altitude:   {rf_mag - 6378.137:.2f} km")
print(f"  Altitude change:  {(rf_mag - r0_mag)*1000:.1f} m (drag decay)")

# ═══════════════════════════════════════════════════════════════
# 3. FLEET HEALTH — ROE Analysis
# ═══════════════════════════════════════════════════════════════

print(f"\n{'=' * 65}")
print("FLEET HEALTH — Relative Orbital Elements")
print(f"{'=' * 65}")

# Use first satellite as the chief (reference)
chief = tles[0].to_mean_elements()
a_chief = chief.a

# Define control box
slot_box = SlotBox(
    da_meters=150.0,     # ±150m SMA
    dlambda_km=15.0,     # ±15 km in-track
    de=1e-4,             # relative eccentricity
    di_arcsec=30.0,      # ±30 arcsec relative inclination
    a_chief=a_chief,
)

print(f"\n  Chief: {tles[0].name} (a = {a_chief:.3f} km)")
print(f"  Slot box: ±150m SMA, ±15km in-track, ±30 arcsec inc")
print()
print(f"  {'SAT':12s} {'δa (m)':>8} {'δλ (km)':>9} {'δe':>10} {'δi (as)':>9} {'STATUS':>8}")
print(f"  {'-'*56}")

for tle in tles[1:]:
    deputy = tle.to_mean_elements()
    roe = ROE.from_elements(chief, deputy)
    status = slot_box.check(chief, roe)

    flag = "⚠ VIOL" if status["violated"] else "✓ OK"

    print(
        f"  {tle.name:12s} "
        f"{roe.da_meters(a_chief):>+8.1f} "
        f"{roe.in_track_km(a_chief):>+9.3f} "
        f"{roe.de_magnitude():>10.2e} "
        f"{roe.di_arcsec():>9.2f} "
        f"{flag:>8}"
    )

    if status["violated"]:
        dv_r, dv_t, dv_n = correct_sma(chief, roe)
        dv_mag = (dv_r**2 + dv_t**2 + dv_n**2)**0.5
        print(f"  {'':12s} └─ Correction: Δv_t = {dv_t:+.3f} m/s (total {dv_mag:.3f} m/s)")

# ═══════════════════════════════════════════════════════════════
# 4. ANNUAL BUDGET
# ═══════════════════════════════════════════════════════════════

print(f"\n{'=' * 65}")
print("ANNUAL Δv BUDGET ESTIMATE")
print(f"{'=' * 65}")

rho = 2e-13    # kg/m³ at ~550 km, moderate solar activity
bc = 40.0      # kg/m²

dv = annual_drag_dv(a_chief, bc, rho)
print(f"\n  Drag makeup:  {dv:.2f} m/s/year")
print(f"  (BC = {bc} kg/m², ρ = {rho:.0e} kg/m³)")
print(f"\n  Note: Actual budget varies 5-50x with solar cycle (F10.7)")
print(f"  Low activity: ~{dv*0.3:.1f} m/s/yr  |  High activity: ~{dv*5:.1f} m/s/yr")
