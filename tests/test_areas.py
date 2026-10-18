"""Tests for the vertical bounds/weights used by the 3D performance indices.

Regression test for issue #150: guess_bounds() assumed an ascending axis and
clamped the outer bounds on fixed ends, so a descending pressure axis (the
convention used by all ECmean climatologies: 100000 -> 1000 Pa) got its first
and last level weighted by ~96 kPa each, ~67% of the total vertical weight
sitting on two levels.
"""

import numpy as np
import xarray as xr

from ecmean.libs.areas import guess_bounds

# EC24 ta climatology levels (Pa), as in the issue report
EC24_PLEV = np.array(
    [
        100000,
        92500,
        85000,
        77500,
        70000,
        60000,
        50000,
        40000,
        30000,
        25000,
        22500,
        20000,
        17500,
        15000,
        12500,
        10000,
        1000,
    ],
    dtype=float,
)


def test_guess_bounds_descending_plev_weights():
    """Descending plev: weights must sum to 1e5 Pa with no huge edge level."""
    bounds = guess_bounds(EC24_PLEV, name="plev")
    weights = np.abs(bounds[:, 0] - bounds[:, 1])

    assert np.isclose(weights.sum(), 100000), f"weights sum {weights.sum()} != 1e5 Pa"
    # surface half-level: 3750 Pa for the EC24 axis (issue #150 expects this)
    assert np.isclose(weights[0], 3750), f"surface weight {weights[0]} != 3750 Pa"
    # top level 1000 Pa below a 10000 Pa level: bounds (0, 5500), not (5500, -3500)
    assert np.isclose(weights[-1], 5500), f"top weight {weights[-1]} != 5500 Pa"
    assert weights.max() < 20000, f"no level should carry ~96 kPa, max is {weights.max()}"


def test_guess_bounds_ascending_plev_unchanged():
    """Ascending plev keeps the previous behaviour (issue #150 fix must not regress)."""
    ascending = EC24_PLEV[::-1].copy()
    bounds = guess_bounds(ascending, name="plev")
    weights = np.abs(bounds[:, 0] - bounds[:, 1])

    assert np.isclose(weights.sum(), 100000)
    # ascending: the 1000 Pa level below a 10000 Pa neighbour spans (0, 5500)
    assert np.isclose(weights[0], 5500)
    # surface level 100000 with a 92500 neighbour spans (96250, 1e5)
    assert np.isclose(weights[-1], 3750)
    # ascending: first lower bound clamped to 0, last upper bound to 1e5
    assert bounds[0, 0] == 0
    assert bounds[-1, 1] == 100000


def test_guess_bounds_descending_plev_clamped_on_right_end():
    """Descending plev: the 1e5 Pa clamp must land on the surface end."""
    bounds = guess_bounds(EC24_PLEV, name="plev")
    # descending: bounds[:,0] holds the surface-side edge of each level
    assert bounds[0, 0] == 100000  # surface side of the first level
    assert bounds[-1, 1] == 0  # top side of the last level


def test_guess_bounds_descending_plev_surface_level_not_dominant():
    """A field that is 1 only on the 1000 hPa level must contribute its ~3.75%
    pressure thickness to the vertical mean, not the ~34% a 96 kPa clamped
    weight would give it (issue #150: ~67% of the weight on two levels)."""
    bounds = guess_bounds(EC24_PLEV, name="plev")
    weights = np.abs(bounds[:, 0] - bounds[:, 1])
    www = xr.DataArray(weights, coords=[EC24_PLEV], dims=["plev"])
    data = np.zeros((len(EC24_PLEV), 3, 4))
    data[0, :, :] = 1.0  # 1 only on the 1000 hPa (surface) level
    field = xr.DataArray(
        data,
        coords=[EC24_PLEV, np.arange(3), np.arange(4)],
        dims=["plev", "lat", "lon"],
    )
    out = field.weighted(www).mean(dim="plev")  # -> (lat, lon) field
    surface_level_mean = float(out.mean())
    assert np.isclose(surface_level_mean, weights[0] / 1e5), (
        f"surface-level fraction {surface_level_mean} != {weights[0] / 1e5}"
    )
    assert surface_level_mean < 0.1, f"surface level carries {surface_level_mean:.2f} of the mean"


def test_guess_bounds_lat_orientation():
    """Latitude bounds: clamps at +-90 land on the right ends for both orientations."""
    asc_lat = np.linspace(-90, 90, 19)
    bounds = guess_bounds(asc_lat, name="lat")
    assert bounds[0, 0] == -90
    assert bounds[-1, 1] == 90

    desc_lat = np.linspace(90, -90, 19)
    bounds = guess_bounds(desc_lat, name="lat")
    # for a descending axis min_bounds/max_bounds hold the swapped physical
    # edges, so the clamps land on bounds[:, 0] at the start and bounds[:, 1]
    # at the end (pre-fix they overwrote the interior edges with +-90)
    assert bounds[0, 0] == 90
    assert bounds[-1, 1] == -90
    # clamped interior weights must stay sane: sum of |widths| == 180
    weights = np.abs(bounds[:, 0] - bounds[:, 1])
    assert np.isclose(weights.sum(), 180)


def test_guess_bounds_lon_unaffected():
    """Longitude (default name) is untouched by the clamping."""
    lon = np.arange(0, 360, 5.0)
    bounds = guess_bounds(lon, name="lon")
    assert np.allclose(np.abs(bounds[:, 0] - bounds[:, 1]), 5.0)
