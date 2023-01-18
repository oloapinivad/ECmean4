#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

import logging
import xarray as xr
import xesmf as xe
import sys
from ecmean.ncfixers import xr_preproc

###########################
# INTERPOLATION FUNCTIONS #
###########################


def remap_dictionary(component, atmareafile, oceareafile, target_grid):
    """Create a dicitionary with atmospheric and oceanic weights for
    interpolation. There is an option of fix grid before the real
    interpolation: this is used for Gaussian reduced grids"""

    atmfix, atmremap = _make_atm_interp_weights(
        component['atm'], atmareafile, target_grid)
    if oceareafile:
        ocefix, oceremap = _make_oce_interp_weights(
            component['oce'], oceareafile, target_grid)
    else:
        ocefix = None
        oceremap = None

    remap = {
        'atm_fix': atmfix,
        'atm_remap': atmremap,
        'oce_fix': ocefix,
        'oce_remap': oceremap,
    }

    return remap


def _make_atm_interp_weights(component, atmareafile, target_grid):
    """Create atmospheric interpolator"""

    logging.debug('Atmareafile is ' + atmareafile)
    if not atmareafile:
        sys.exit("ERROR: Atmareafile cannot be found")

    if component == 'oifs':

        # this is to get lon and lat from the Equator
        xfield = xr.open_mfdataset(atmareafile, preprocess=xr_preproc).load()
        xname = list(xfield.data_vars)[-1]
        m = xfield[xname].isel(time=0).load()
        g = sorted(list(set(m.lat.values)))
        f = sorted(list(m.sel(cell=m.lat == g[int(len(g) / 2)]).lon.values))

        # this creates a a gaussian non reduced grid
        ds_out = xr.Dataset({"lon": (["lon"], f), "lat": (["lat"], g)})

        # use nearest neighbour to remap to gaussian regular
        fix = xe.Regridder(
            xfield[xname],
            ds_out,
            method="nearest_s2d",
            locstream_in=True,
            periodic=True)

        # create bilinear interpolator
        interp = xe.Regridder(
            fix(xfield[xname]), target_grid, periodic=True, method="bilinear")

    elif component == 'cmoratm':

        fix = None
        xfield = xr.open_mfdataset(atmareafile, preprocess=xr_preproc).load()
        interp = xe.Regridder(
            xfield,
            target_grid,
            periodic=True,
            method="bilinear")

    elif component == 'globo':

        fix = None
        xfield = xr.open_mfdataset(atmareafile, preprocess=xr_preproc).load()
        interp = xe.Regridder(
            xfield,
            target_grid,
            periodic=True,
            ignore_degenerate=True,
            method="bilinear")

    else:
        sys.exit(
            "ERROR: Atm weights not defined for this component, this cannot be handled!")

    return fix, interp


def _make_oce_interp_weights(component, oceareafile, target_grid):
    """Create oceanic interpolator weights"""

    logging.debug('Oceareafile is ' + oceareafile)
    if not oceareafile:
        sys.exit("ERROR: Oceareafile cannot be found")

    if component == 'nemo':
        fix = None
        xfield = xr.open_mfdataset(oceareafile, preprocess=xr_preproc).load()

        # set coordinates which are missing
        for cl in ['nav_lon', 'nav_lat', 'nav_lev', 'time_counter', 'x', 'y']:
            if cl in xfield.data_vars:
                xfield = xfield.set_coords([cl])

        # rename dimensions and coordinates
        xfield = xfield.rename(
            {"nav_lon": "lon", "nav_lat": "lat", "nav_lev": "deptht"})

        # use grid distance as generic variable
        interp = xe.Regridder(
            xfield['e1t'],
            target_grid,
            method="bilinear",
            periodic=True,
            ignore_degenerate=True)

    elif component == 'cmoroce':

        fix = None
        xfield = xr.open_mfdataset(oceareafile, preprocess=xr_preproc)
        xname = list(xfield.data_vars)[-1]
        # print(len(xfield.coords['lon'].shape))

        # check if oceanic grid is regular: lon/lat dims should be 1d
        # if not all(x in xfield.dims for x in ['lon', 'lat']) and (len(xfield.dims) < 3) :
        if len(xfield.coords['lon'].shape) == 1 and len(xfield.coords['lat'].shape) == 1:

            print("Detecting a unstructured grid, using nearest neighbour!")
            interp = xe.Regridder(
                xfield[xname].load(),
                target_grid,
                method="nearest_s2d",
                locstream_in=True,
                periodic=True)
        else:
            print("Detecting regular or curvilinear grid, using bilinear!")
            interp = xe.Regridder(
                xfield[xname].load(),
                target_grid,
                method="bilinear",
                ignore_degenerate=True,
                periodic=True)

    else:
        sys.exit(
            "ERROR: Oce weights not defined for this component, this cannot be handled!")

    return fix, interp