"""Walker constellation generator for SCARAB.
"""
from scarab import MeanElements


def walker_delta(
    total_sats: int,
    num_planes: int,
    phasing: int,
    altitude_km: float,
    inclination_deg: float,
    raan_base_deg: float = 0.0,
    epoch: float = 0.0,
) -> list[tuple[int, MeanElements]]:
    """
    Generate a Walker Delta constellation (i:T/P/F notation).

    Args:
        total_sats: Total number of satellites (T)
        num_planes: Number of orbital planes (P)
        phasing: Phasing parameter (F), 0 <= F < P
        altitude_km: Orbital altitude above Earth's surface (km)
        inclination_deg: Orbital inclination (degrees)
        raan_base_deg: RAAN of the first plane (degrees)
        epoch: Epoch in seconds

    Returns:
        List of (satellite_id, MeanElements) tuples.
    """
    R_EARTH = 6378.137
    a = R_EARTH + altitude_km
    sats_per_plane = total_sats // num_planes

    satellites = []
    sat_id = 0

    for plane in range(num_planes):
        raan = raan_base_deg + plane * (360.0 / num_planes)

        for sat_in_plane in range(sats_per_plane):
            ma = (
                sat_in_plane * (360.0 / sats_per_plane)
                + plane * phasing * (360.0 / total_sats)
            ) % 360.0

            elem = MeanElements(
                a=a,
                e=0.0,
                i_deg=inclination_deg,
                raan_deg=raan % 360.0,
                aop_deg=0.0,
                ma_deg=ma,
                epoch=epoch,
            )
            satellites.append((sat_id, elem))
            sat_id += 1

    return satellites


def walker_star(
    total_sats: int,
    num_planes: int,
    altitude_km: float,
    inclination_deg: float,
    raan_base_deg: float = 0.0,
    epoch: float = 0.0,
) -> list[tuple[int, MeanElements]]:
    """
    Generate a Walker Star constellation (RAAN spacing = 180°/P).
    """
    R_EARTH = 6378.137
    a = R_EARTH + altitude_km
    sats_per_plane = total_sats // num_planes

    satellites = []
    sat_id = 0

    for plane in range(num_planes):
        raan = raan_base_deg + plane * (180.0 / num_planes)

        for sat_in_plane in range(sats_per_plane):
            ma = sat_in_plane * (360.0 / sats_per_plane)

            elem = MeanElements(
                a=a,
                e=0.0,
                i_deg=inclination_deg,
                raan_deg=raan % 360.0,
                aop_deg=0.0,
                ma_deg=ma,
                epoch=epoch,
            )
            satellites.append((sat_id, elem))
            sat_id += 1

    return satellites
