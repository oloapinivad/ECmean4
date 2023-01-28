#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

import numpy as np
import logging
import sys
import xarray as xr
from ecmean.libs.ncfixers import xr_preproc
from ecmean.libs.files import inifiles_priority

##################################
# AREA-WEIGHT FUNCTIONS #
##################################


def areas_dictionary(component, atmdict, ocedict):
    """Create a dictionary with atmospheric and oceanic area weights"""

    # get atmospheric areas
    atmareafile = inifiles_priority(atmdict)
    atm_areas = _make_atm_areas(component['atm'], atmareafile)

    # get oceanic areas, assuming AMIP if no oceanic area is found
    oceareafile = inifiles_priority(ocedict)
    if oceareafile:
        oce_areas = _make_oce_areas(component['oce'], oceareafile)
    else:
        logging.warning("Ocereafile cannot be found, assuming this is an AMIP run")
        oce_areas = None

    areas = {
        'atm_areas': atm_areas,
        'oce_areas': oce_areas
    }

    return areas


def _make_atm_areas(component, atmareafile):
    """Create atmospheric weights for area operations. Load once defined to allow
    for parallel computation."""

    logging.info('Atmareafile is ' + atmareafile)
    if not atmareafile:
        sys.exit("ERROR: Atmareafile cannot be found")

    # this might be universal, but keep this as for supported components only
    if component in ['oifs', 'cmoratm', 'globo']:
        xfield = xr.open_mfdataset(atmareafile, preprocess=xr_preproc)
        area = _area_cell(xfield)
    else:
        sys.exit("ERROR: Area for this configuration cannot be handled!")
    return area.load()


def _make_oce_areas(component, oceareafile):
    """Create atmospheric weights for area operations. Load once defined to allow
    for parallel computation
    """

    logging.info('Oceareafile is ' + oceareafile)

    xfield = xr.open_mfdataset(oceareafile, preprocess=xr_preproc)

    # this might be universal, but keep this as for supported components only
    if component in ['nemo', 'cmoroce']:

        if 'areacello' in xfield.data_vars:  # CMOR case
            area = xfield['areacello']
        elif 'cell_area' in xfield.data_vars:  # ECE4 NEMO case for nemo-initial-state.nc
            area = xfield['cell_area']
        elif 'e1t' in xfield.data_vars:  # ECE4 NEMO case for domaing_cfg.nc, deprecated
            area = xfield['e1t'] * xfield['e2t']
        else:  # automatic solution!
            area = _area_cell(xfield)

    else:
        sys.exit("ERROR: Area for this configuration cannot be handled!")

    return area.load()


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
        min_bounds[0] = (-90)
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


def _area_cell(xfield, formula='triangles'):
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

    earth_radius = 6371000.

    # some check to starts
    if all(x in xfield.dims for x in ['lon', 'lat']):
        logging.debug('Regulard grid recognized..')
        regular_grid = True
    else:
        regular_grid = False

    if all(x in xfield.data_vars for x in ['lon_bnds', 'lat_bnds']):
        logging.debug('cmor lon/lat_bounds found...')
        cmor_bounds = True
    else:
        cmor_bounds = False

    # this is a nightmare, so far working only for ECE4 gaussian reduced
    if not regular_grid:

        logging.debug('Curvilinear/Unstructured grid, tryin to get grid info...')

        blondim = None
        blatdim = None
        # trying to find bounderies
        for g in (
            t for t in list(
                xfield.data_vars) if t in [
                'lon_bnds',
                'bounds_lon']):
            blondim = g
        for g in (
            t for t in list(
                xfield.data_vars) if t in [
                'lat_bnds',
                'bounds_lat']):
            blatdim = g

        # checking
        if blondim is None and blatdim is None:
            sys.exit(
                "ERROR: Can't find any lon/lat boundaries and grid is unstructured, need some help!")

        logging.debug('Unstructured grid, special ECE4 treatment...')
        # ATTENTION: this is a very specific ECE4 definition, it will not work
        # with other unstructured grids. The assumption of the vertex position
        # is absolutely random. Needs to be generalized.
        bounds_lon = np.column_stack((xfield[blondim].isel(nvertex=1),
                                      xfield[blondim].isel(nvertex=2)))
        bounds_lat = np.column_stack((xfield[blatdim].isel(nvertex=2),
                                      xfield[blatdim].isel(nvertex=3)))
        area_dims = 'cell'
        # set full lat
        full_lat = xfield['lat'].values

    # if we are dealing with a regular grid
    if regular_grid:

        # if we have bounds, just check they have the right dimensions names
        if cmor_bounds:

            # if dimension is not called bnds, rename it
            if 'bnds' not in list(xfield.dims):
                logging.debug('bnds not found, trying to rename it...')
                for g in (
                    t for t in list(
                        xfield.dims) if t not in [
                        'lon',
                        'lat',
                        'time']):
                    bdim = g
                xfield = xfield.rename_dims({bdim: "bnds"})

        # else use guess_bounds() and expand the xarray dataset including them
        if not cmor_bounds:

            logging.debug('Bounds estimation from lon/lat...')
            # create and xarray dataset which the boundaries
            xbounds = xr.Dataset(
                data_vars=dict(
                    lat_bnds=(
                        ('lat', 'bnds'), guess_bounds(
                            xfield['lat'], name='lat')), lon_bnds=(
                        ('lon', 'bnds'), guess_bounds(
                            xfield['lon'], name='lon'))), coords=dict(
                    lat=(
                        'lat', xfield['lat'].values), lon=(
                        'lon', xfield['lon'].values)))
            xfield = xfield.merge(xbounds)

        # create numpy array
        blon = np.column_stack((xfield['lon_bnds'].isel(bnds=0),
                                xfield['lon_bnds'].isel(bnds=1)))
        blat = np.column_stack((xfield['lat_bnds'].isel(bnds=0),
                                xfield['lat_bnds'].isel(bnds=1)))
        full_lat = np.repeat(xfield['lat'].values, len(xfield['lon']), axis=0)

        # 2d matrix of bounds
        expansion = np.array([(y, x) for y in blat for x in blon])
        bounds_lon = expansion[:, 1, :]
        bounds_lat = expansion[:, 0, :]
        area_dims = ('lat', 'lon')

    # cell dimension
    if formula == "triangles":
        p1 = _lonlat_to_sphere(bounds_lon[:, 0], bounds_lat[:, 0]).transpose()
        p2 = _lonlat_to_sphere(bounds_lon[:, 0], bounds_lat[:, 1]).transpose()
        p3 = _lonlat_to_sphere(bounds_lon[:, 1], bounds_lat[:, 1]).transpose()
        p4 = _lonlat_to_sphere(bounds_lon[:, 1], bounds_lat[:, 0]).transpose()
        area_cell = _vector_spherical_triangle(
            p1, p2, p3) + _vector_spherical_triangle(p1, p4, p3)
        area_cell = area_cell * earth_radius**2
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
            full_lat = np.repeat(
                xfield['lat'].values, len(
                    xfield['lon']), axis=0)
            arclon1 = arclon2 = earth_radius * \
                abs(np.cos(abs(np.deg2rad(full_lat)))) * np.deg2rad(dlon)

        arclat = earth_radius * np.deg2rad(dlat)

        # trapezoid area
        area_cell = (arclon1 + arclon2) * arclat / 2

    if regular_grid:
        area_cell = area_cell.reshape([len(xfield['lat']), len(xfield['lon'])])

    # since we are using numpy need to bring them back into xarray dataset
    xfield['area'] = (area_dims, area_cell)

    # check the total area
    logging.debug('Total Earth Surface: %s Km2',
                  str(xfield['area'].sum().values / 10**6))

    return xfield['area']
