#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

import logging
import xarray as xr
import xesmf as xe
import sys
from ecmean.libs.ncfixers import xr_preproc
from ecmean.libs.files import inifiles_priority
from ecmean.libs.areas import identify_grid

###########################
# INTERPOLATION FUNCTIONS #
###########################


def remap_dictionary(component, atmdict, ocedict, target_grid):
    """Create a dicitionary with atmospheric and oceanic weights for
    interpolation. There is an option of fix grid before the real
    interpolation: this is used for Gaussian reduced grids"""

    atmareafile = inifiles_priority(atmdict)
    atmfix, atmremap = _make_atm_interp_weights(
        component['atm'], atmareafile, target_grid)

    # get oceanic areas, assuming AMIP if no oceanic area is found
    oceareafile = inifiles_priority(ocedict)
    if oceareafile:
        ocefix, oceremap = _make_oce_interp_weights(
            component['oce'], oceareafile, target_grid)
    else:
        logging.warning("Ocereafile cannot be found, assuming this is an AMIP run")
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

    xfield = xr.open_mfdataset(atmareafile, preprocess=xr_preproc).load()
    gridtype = identify_grid(xfield)
    logging.warning(f'Atmosphere grid is is a {gridtype} grid!')

    if component == 'oifs':

        # this is to get lon and lat from the Equator
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

    elif component in ['cmoratm', 'globo']:

        fix = None
        interp = xe.Regridder(
            xfield,
            target_grid,
            periodic=True,
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

    xfield = xr.open_mfdataset(oceareafile, preprocess=xr_preproc).load()
    gridtype = identify_grid(xfield)
    logging.warning(f'Ocean grid is is a {gridtype} grid!')

    if component == 'nemo':
        xname = 'cell_area'
    elif component == 'cmoroce':
        xname = list(xfield.data_vars)[-1]
    else:
        sys.exit(
            "ERROR: Oce weights not defined for this component, this cannot be handled!")

    if gridtype in ['unstructured'] : 
        #print("Detecting a unstructured grid, using nearest neighbour!")
        fix = None
        interp = xe.Regridder(
            xfield[xname],
            target_grid,
            method="nearest_s2d",
            locstream_in=True,
            periodic=True)
    
    else: 
       #print("Detecting regular or curvilinear grid, using bilinear!")
        fix = None
        interp = xe.Regridder(
            xfield[xname],
            target_grid,
            method="bilinear",
            periodic=True,
            ignore_degenerate=True)

    return fix, interp
