"""
Tests for SCARAB Python bindings.
"""
import pytest
import math


def test_tle_parse():
    from scarab import TLE

    line1 = "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9003"
    line2 = "2 25544  51.6400 208.5000 0007417  68.0000 292.1000 15.49560000400000"

    tle = TLE.parse(line1, line2)
    assert tle.norad_id == 25544
    assert abs(tle.inclination_deg - 51.64) < 0.01
    assert abs(tle.eccentricity - 0.0007417) < 1e-7
    assert 400 < tle.altitude() < 430


def test_tle_parse_3line():
    from scarab import TLE

    name = "ISS (ZARYA)"
    line1 = "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9003"
    line2 = "2 25544  51.6400 208.5000 0007417  68.0000 292.1000 15.49560000400000"

    tle = TLE.parse_3line(name, line1, line2)
    assert tle.name == "ISS (ZARYA)"
    assert tle.norad_id == 25544


def test_tle_batch():
    from scarab import TLE

    batch = """ISS (ZARYA)
1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9003
2 25544  51.6400 208.5000 0007417  68.0000 292.1000 15.49560000400000
"""
    tles = TLE.parse_batch(batch)
    assert len(tles) == 1
    assert tles[0].name == "ISS (ZARYA)"


def test_tle_to_mean_elements():
    from scarab import TLE

    line1 = "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9003"
    line2 = "2 25544  51.6400 208.5000 0007417  68.0000 292.1000 15.49560000400000"

    tle = TLE.parse(line1, line2)
    elem = tle.to_mean_elements()
    assert abs(elem.i_deg - 51.64) < 0.01
    assert elem.a > 6378.137 + 400


def test_mean_elements_propagation():
    from scarab import MeanElements

    elem = MeanElements(a=6928.137, e=0.001, i_deg=97.6,
                        raan_deg=0.0, aop_deg=0.0, ma_deg=0.0, epoch=0.0)
    prop = elem.propagate(86400.0)

    # Sun-sync RAAN should drift ~0.9856 deg/day
    assert abs(prop.raan_deg - 0.9856) < 0.1


def test_roe_identical():
    from scarab import MeanElements, ROE

    chief = MeanElements(a=6928.137, e=0.001, i_deg=53.0,
                         raan_deg=0.0, aop_deg=0.0, ma_deg=0.0, epoch=0.0)
    roe = ROE.from_elements(chief, chief)

    assert abs(roe.da) < 1e-14
    assert abs(roe.dlambda) < 1e-14


def test_roe_altitude_difference():
    from scarab import MeanElements, ROE

    chief = MeanElements(a=6928.137, e=0.0, i_deg=53.0,
                         raan_deg=0.0, aop_deg=0.0, ma_deg=0.0, epoch=0.0)
    deputy = MeanElements(a=6929.137, e=0.0, i_deg=53.0,
                          raan_deg=0.0, aop_deg=0.0, ma_deg=0.0, epoch=0.0)

    roe = ROE.from_elements(chief, deputy)
    assert abs(roe.da_meters(chief.a) - 1000.0) < 1.0


def test_slot_box():
    from scarab import MeanElements, ROE, SlotBox

    chief = MeanElements(a=6928.137, e=0.0, i_deg=53.0,
                         raan_deg=0.0, aop_deg=0.0, ma_deg=0.0, epoch=0.0)

    # Deputy 200m higher — should violate ±100m box
    deputy = MeanElements(a=6928.337, e=0.0, i_deg=53.0,
                          raan_deg=0.0, aop_deg=0.0, ma_deg=0.0, epoch=0.0)

    roe = ROE.from_elements(chief, deputy)
    box = SlotBox(da_meters=100.0, dlambda_km=10.0, de=1e-4,
                  di_arcsec=20.0, a_chief=chief.a)
    status = box.check(chief, roe)
    assert status["violated"] is True
    assert status["da_violated"] is True


def test_numerical_propagator_two_body():
    from scarab import MeanElements, Propagator

    elem = MeanElements(a=6878.137, e=0.001, i_deg=51.6,
                        raan_deg=0.0, aop_deg=0.0, ma_deg=0.0, epoch=0.0)

    # Two-body propagation should conserve energy
    prop = Propagator(j2=False, drag=False)
    states = prop.propagate_elements(elem, duration_s=5400.0)  # ~1 orbit

    assert len(states) > 10

    # Check first and last altitudes are similar (circular orbit)
    r0 = (states[0][1]**2 + states[0][2]**2 + states[0][3]**2)**0.5
    rf = (states[-1][1]**2 + states[-1][2]**2 + states[-1][3]**2)**0.5
    assert abs(r0 - rf) < 1.0  # Within 1 km


def test_numerical_propagator_drag():
    from scarab import MeanElements, Propagator

    elem = MeanElements(a=6798.137, e=0.0001, i_deg=51.6,
                        raan_deg=0.0, aop_deg=0.0, ma_deg=0.0, epoch=0.0)

    prop = Propagator(j2=True, drag=True, ballistic_coefficient=40.0)
    states = prop.propagate_elements(elem, duration_s=86400.0)

    # Drag should lower the orbit
    r0 = (states[0][1]**2 + states[0][2]**2 + states[0][3]**2)**0.5
    rf = (states[-1][1]**2 + states[-1][2]**2 + states[-1][3]**2)**0.5
    # Final mean altitude should be lower (check SMA via energy)
    mu = 398600.4418
    v0 = (states[0][4]**2 + states[0][5]**2 + states[0][6]**2)**0.5
    vf = (states[-1][4]**2 + states[-1][5]**2 + states[-1][6]**2)**0.5
    e0 = v0**2 / 2 - mu / r0
    ef = vf**2 / 2 - mu / rf
    assert ef < e0, "Drag should decrease orbital energy"


def test_correction_maneuver():
    from scarab import MeanElements, ROE, correct_sma

    chief = MeanElements(a=6928.137, e=0.0, i_deg=53.0,
                         raan_deg=0.0, aop_deg=0.0, ma_deg=0.0, epoch=0.0)
    deputy = MeanElements(a=6928.237, e=0.0, i_deg=53.0,
                          raan_deg=0.0, aop_deg=0.0, ma_deg=0.0, epoch=0.0)

    roe = ROE.from_elements(chief, deputy)
    dv_r, dv_t, dv_n = correct_sma(chief, roe)

    # Should be an along-track burn
    assert abs(dv_r) < 1e-6
    assert abs(dv_n) < 1e-6
    assert abs(dv_t) > 0.001  # Some tangential dv needed
    assert dv_t < 0  # Retrograde to lower orbit
