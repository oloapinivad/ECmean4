"""Tests for bounds and area functions"""

import numpy as np
import pytest
import xarray as xr

from ecmean.libs.areas import AreaCalculator, guess_bounds

# pressure levels of the ECmean climatologies (Pa), stored descending
PLEV = np.array(
    [100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000],
    dtype=float,
)


def _widths(bounds):
    return np.abs(bounds[:, 0] - bounds[:, 1])


@pytest.mark.parametrize("plev", [PLEV, PLEV[::-1]])
def test_guess_bounds_plev(plev):
    """Vertical weights must cover the column from the surface to the top, whatever the orientation"""
    bounds = guess_bounds(plev, name="plev")
    widths = _widths(bounds)
    assert np.isclose(widths.sum(), 100000)
    assert np.isclose(widths[np.argmax(plev)], 3750)
    assert np.isclose(widths[np.argmin(plev)], 1500)


def test_area_lat_orientation():
    """Cell areas do not depend on latitude orientation and sum to the Earth surface (longitude is never clamped)"""
    lon = np.arange(0.5, 360, 1.0)
    lat = np.arange(-89.5, 90, 1.0)
    calc = AreaCalculator()
    areas = []
    for lats in [lat, lat[::-1]]:
        xfield = xr.Dataset(
            {"dummy": (("lat", "lon"), np.zeros((len(lats), len(lon))))},
            coords={"lat": lats, "lon": lon},
        )
        areas.append(calc.calculate_area(xfield, gridtype="lonlat"))
    earth = 4 * np.pi * calc.earth_radius**2
    assert np.isclose(areas[0].sum(), earth, rtol=1e-6)
    np.testing.assert_allclose(areas[0].values, areas[1].values[::-1], rtol=1e-10)
