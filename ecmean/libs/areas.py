#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

import logging
import numpy as np
import xarray as xr

##################################
# AREA-WEIGHT FUNCTIONS #
##################################

# def cdo_identify_grid(file):
#     cdogrid = cdo.sinfo(input=os.path.join(dir, file))
#     for gg in ['unstructured', 'lonlat', 'curvilinear', 'gaussian_reduced', 'gaussian']:
#         grid_type = [gg for i in cdogrid if gg in i]
#         if grid_type :
#             break
#     return grid_type[0]

loggy = logging.getLogger(__name__)

def guess_bounds(axis, name='lon'):
    """Basic function that estimates the boundaries for lon and lat if they are not
    available. Works only with regular grids.
    It also avoid having values larger than 90N/90S as well as 0 and 10^5 Pa"""
    # inspired by
    # https://gist.github.com/dennissergeev/60bf7b03443f1b2c8eb96ce0b1880150

    # this define the proportion of the bounds, assumed to half of the levels
    bound_position = 0.5
    diffs = np.diff(axis)
    diffs = np.insert(diffs, 0, diffs[0])
    diffs = np.append(diffs, diffs[-1])

    # pair of dounds
    min_bounds = axis - diffs[:-1] * bound_position
    max_bounds = axis + diffs[1:] * (1 - bound_position)

    # safety check, to be generalized
    if name in 'lat':
        max_bounds[-1] = 90
        min_bounds[0] = -90
    if name in 'plev':
        min_bounds[0] = 0
        max_bounds[-1] = 100000

    # should we use a xarray object instead of a numpy?
    bounds = np.array([min_bounds, max_bounds]).transpose()
    return bounds


def _lonlat_to_sphere(lon, lat):
    """Convert from lon lat coordinates to a 3d sphere of unity radius"""

    vec = np.array([
        np.cos(np.deg2rad(lon)) * np.cos(np.deg2rad(lat)),
        np.sin(np.deg2rad(lon)) * np.cos(np.deg2rad(lat)),
        np.sin(np.deg2rad(lat))
    ])
    return vec


def _huilier(a, b, c):
    """Apply the L'Huilier theorem from the three side of the spherical triangle
    obtaining the spherical excess, i.e. the solid angle of the triangle,
    i.e. the area of the spherical surface"""
    # More info at https://mathworld.wolfram.com/LHuiliersTheorem.html

    s = (a + b + c) * 0.5
    t = np.tan(s * 0.5) * np.tan((s - a) * 0.5) * \
        np.tan((s - b) * 0.5) * np.tan((s - c) * 0.5)
    area = abs(4. * np.arctan(np.sqrt(abs(t))))
    return area


def _vector_spherical_triangle(p1, p2, p3):
    """Given the coordinates of three points on a sphere, estimate the length
    of their vectors a,b,c connecting to the centre of the spere.
    Then using L'Huilier formula derive the solid angle among the three
    vectors, which is multiplied by squared Earth Radius is
    exactly the surface of the corresponding spherical triangle. """
    # This is inspired by CDO code found at
    # https://code.mpimet.mpg.de/projects/cdo/repository/cdo/revisions/331ab3f7fd18295cf6a433fb799034c7589a4a61/entry/src/grid_area.cc

    a = np.arcsin(np.linalg.norm(np.cross(p1, p2), axis=1))
    b = np.arcsin(np.linalg.norm(np.cross(p1, p3), axis=1))
    c = np.arcsin(np.linalg.norm(np.cross(p3, p2), axis=1))

    area = _huilier(a, b, c)
    return area


def area_cell(xfield, gridtype=None, formula='triangles'):
    """
    Function which estimate the area cell from bounds. This is done assuming
    making use of spherical triangels.
    Working also on regular grids which does not have lon/lat bounds
    via the guess_bounds function. Curvilinear/unstructured grids are not supported,
    especially if with more with more than 4 vertices are not supported.

    Args:
    xfield: a generic xarray dataset
    formula: 'squares' or 'trapezoids' or 'triangles' equation for the area cell
        'triangles' is the default, uses the spherical triangles - same as
        used by CDO - and it is very accurate

    Returns:
    An xarray dataarray with the area for each grid point
    """

    # safe check to operate only on single timeframe
    if 'time' in xfield.dims:
        xfield = xfield.isel(time=0)

    if all(x in xfield.data_vars for x in ['lon_bnds', 'lat_bnds']):
        loggy.debug('cmor lon/lat_bounds found...')
        cmor_bounds = True
    else:
        cmor_bounds = False

    # this is a nightmare, so far working only for ECE4 gaussian reduced
    if gridtype in ['unstructured', 'gaussian_reduced', 'curvilinear']:

        loggy.info('Curvilinear/Unstructured grid, tryin to get grid info...')

        # use next to assign the first occurence in data vars for the boundaries
        blondim = next((t for t in xfield.data_vars if t in ['lon_bnds', 'bounds_lon']), None)
        blatdim = next((t for t in xfield.data_vars if t in ['lat_bnds', 'bounds_lat']), None)

        # checking
        if blondim is None and blatdim is None:
            raise ValueError(
                "ERROR: Can't find any lon/lat boundaries and grid is unstructured, need some help!")

        loggy.info('Unstructured grid, special ECE4 treatment...')
        # ATTENTION: this is a very specific ECE4 definition, it will not work
        # with other unstructured grids. The assumption of the vertex position
        # is absolutely random. Needs to be generalized.
        bounds_lon = np.column_stack((xfield[blondim].isel(nvertex=1),
                                      xfield[blondim].isel(nvertex=2)))
        bounds_lat = np.column_stack((xfield[blatdim].isel(nvertex=2),
                                      xfield[blatdim].isel(nvertex=3)))
        area_dims = 'cell'
        # set full lat
        full_lat = xfield['lat'].data

    # if we are dealing with a regular grid
    elif gridtype in ['gaussian', 'lonlat']:

        # if we have bounds, just check they have the right dimensions names
        if cmor_bounds:

            # if dimension is not called bnds, rename it
            if 'bnds' not in xfield.dims:
                loggy.info('bnds not found, trying to rename it...')
                bdim = [g for g in xfield.dims if g not in ['lon', 'lat', 'time']][0]
                xfield = xfield.rename_dims({bdim: "bnds"})

        # else use guess_bounds() and expand the xarray dataset including them
        if not cmor_bounds:

            loggy.debug('Bounds estimation from lon/lat...')
            # create and xarray dataset which the boundaries
            xfield['lat_bnds'] = (('lat', 'bnds'), guess_bounds(xfield['lat'], name='lat'))
            xfield['lon_bnds'] = (('lon', 'bnds'), guess_bounds(xfield['lon'], name='lon'))

        # create numpy array
        blon = np.column_stack((xfield['lon_bnds'].isel(bnds=0),
                                xfield['lon_bnds'].isel(bnds=1)))
        blat = np.column_stack((xfield['lat_bnds'].isel(bnds=0),
                                xfield['lat_bnds'].isel(bnds=1)))
        full_lat = np.repeat(xfield['lat'].data, len(xfield['lon']), axis=0)

        # 2d matrix of bounds
        # expansion and speed up using numpy instead of list comprehension
        bounds_lon = np.stack([blon] * len(blat), axis=0)
        bounds_lat = np.stack([blat] * len(blon), axis=1)
        bounds_lon = bounds_lon.reshape(-1, 2)
        bounds_lat = bounds_lat.reshape(-1, 2)

        area_dims = ('lat', 'lon')

    else:
        raise ValueError('Gridtype undefined!')

    # compute the area
    area = _area_computation(bounds_lon, bounds_lat, formula=formula, full_lat=full_lat)

    if gridtype in ['gaussian', 'lonlat']:
        area = area.reshape([len(xfield['lat']), len(xfield['lon'])])

    # since we are using numpy need to bring them back into xarray dataset
    outfield = xr.DataArray(area, dims=area_dims, coords=xfield.coords, name='area')

    # check the total area
    loggy.debug('Total Earth Surface: %s Km2',
                 str(outfield.sum().values / 10**6))

    return outfield


def _area_computation(bounds_lon, bounds_lat, formula='triangles', full_lat=None):
    """Inner function for the area computation"""

    earth_radius = 6371000.

    # cell dimension
    if formula == "triangles":

        p1 = _lonlat_to_sphere(bounds_lon[:, 0], bounds_lat[:, 0]).transpose()
        p2 = _lonlat_to_sphere(bounds_lon[:, 0], bounds_lat[:, 1]).transpose()
        p3 = _lonlat_to_sphere(bounds_lon[:, 1], bounds_lat[:, 1]).transpose()
        p4 = _lonlat_to_sphere(bounds_lon[:, 1], bounds_lat[:, 0]).transpose()
        area = _vector_spherical_triangle(
            p1, p2, p3) + _vector_spherical_triangle(p1, p4, p3)
        area = area * earth_radius**2
    else:

        # cell dimension
        dlon = abs(bounds_lon[:, 0] - bounds_lon[:, 1])
        dlat = abs(bounds_lat[:, 0] - bounds_lat[:, 1])

        # safe check on cosine of 90 included:
        # assume a trapezoid or a squared cell
        if formula == 'trapezoids':
            arclon1 = earth_radius * \
                abs(np.cos(abs(np.deg2rad(bounds_lat[:, 0])))) * np.deg2rad(dlon)
            arclon2 = earth_radius * \
                abs(np.cos(abs(np.deg2rad(bounds_lat[:, 1])))) * np.deg2rad(dlon)
        if formula == 'squares':
            # need to provide full_lat
            arclon1 = arclon2 = earth_radius * \
                abs(np.cos(abs(np.deg2rad(full_lat)))) * np.deg2rad(dlon)

        arclat = earth_radius * np.deg2rad(dlat)

        # trapezoid area
        area = (arclon1 + arclon2) * arclat / 2

    return area
