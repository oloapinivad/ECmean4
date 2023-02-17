#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

import logging
import sys
import xarray as xr
import numpy as np
from ecmean.libs.ncfixers import xr_preproc
from ecmean.libs.units import units_converter

##########################
# MASK-RELATED FUNCTIONS #
##########################


def masks_dictionary(component, maskatmfile, maskocefile, remap_dictionary=None):
    """
    Create a dictionary with atmospheric land-sea mask
    Oceanic mask are not used yet, since it is assumed that oceanic variables
    are computed by definition on the oceanic grid.
    """

    # default mask: mandatory
    matm = _make_atm_masks(component['atm'],
                           maskatmfile, remap_dictionary=remap_dictionary)

    # if I have the oceanic mask, use it
    if maskocefile:
        moce = _make_oce_masks(component['oce'],
                               maskocefile, remap_dictionary=remap_dictionary)
    else:
        # if it is missing, for PI I can use the atmospheric one
        if remap_dictionary:
            moce = matm
        # for global mean, I have no solution: I will not use it over the ocean
        else:
            moce = None

    mask = {
        'atm_mask': matm,
        'oce_mask': moce,
    }

    return mask


def _make_atm_masks(component, maskatmfile, remap_dictionary=None):
    """Create land-sea masks for atmosphere model"""

    # prepare ATM LSM
    logging.info('maskatmfile is' + maskatmfile)
    if not maskatmfile:
        sys.exit("ERROR: maskatmfile cannot be found")

    if component == 'oifs':
        # create mask: opening a grib and loading only lsm to avoid
        # inconsistencies # in the grib structure ->
        # see here https://github.com/ecmwf/cfgrib/issues/13
        mask = xr.open_mfdataset(
            maskatmfile,
            engine="cfgrib",
            indexpath=None,
            filter_by_keys={
                'shortName': 'lsm'},
            preprocess=xr_preproc)['lsm']

    elif component in ['cmoratm', 'globo']:
        dmask = xr.open_mfdataset(maskatmfile, preprocess=xr_preproc)
        if 'sftlf' in dmask.data_vars:
            mask = dmask['sftlf']
            mask = mask / 100  # cmor mask are %
        elif 'lsm' in dmask.data_vars:
            mask = dmask['lsm']

        # globo has a reversed mask
        if component == 'globo':
            mask = abs(1 - mask)
    else:
        sys.exit("ERROR: _make_atm_masks -> Mask undefined yet mismatch, this cannot be handled!")

    # safe check to operate only on single timeframe
    if 'time' in mask.dims:
        mask = mask.isel(time=0)

    # interp the mask if required
    if remap_dictionary is not None:
        if remap_dictionary['atm_fix']:
            mask = remap_dictionary['atm_fix'](mask, keep_attrs=True)
        mask = remap_dictionary['atm_remap'](mask, keep_attrs=True)

    return mask


def _make_oce_masks(component, maskocefile, remap_dictionary=None):
    """Create land-sea masks for oceanic model. This is used only for CMIP"""

    # prepare ocean LSM:
    logging.info('maskocefile is' + maskocefile)
    if not maskocefile:
        sys.exit("ERROR: maskocefile cannot be found")

    if component == 'cmoroce':
        dmask = xr.open_mfdataset(maskocefile, preprocess=xr_preproc)
        if 'sftof' in dmask.data_vars:
            # use fillna to have a binary max (0 land, 1 sea)
            mask = dmask['sftof'].fillna(0)

        # check if we need to convert from % to fraction
        # offset should not count!
        if mask.units:
            offset, factor = units_converter(mask.units, 'frac')
            mask = (mask * factor) + offset

        # to keep coehrence in other operations, oceanic mask is set to be
        # identical to atmospheric mask, i.e. 1 over land and 0 over ocean
        mask = abs(1 - mask)

    else:
        mask = None
        sys.exit("ERROR: _make_oce_masks -> Mask undefined yet mismatch, this cannot be handled!")

    # safe check to operate only on single timeframe
    if 'time' in mask.dims:
        mask = mask.isel(time=0)

    # interp the mask if required
    if remap_dictionary is not None:
        if remap_dictionary['oce_fix']:
            mask = remap_dictionary['oce_fix'](mask, keep_attrs=True)
        mask = remap_dictionary['oce_remap'](mask, keep_attrs=True)

    return mask


def masked_meansum(xfield, weights, mask, operation, domain, mask_type):
    """For global variables evaluate the weighted averaged
    or weighted integral when required by the variable properties"""

    # call the mask_field to mask where necessary
    # the mask field is area
    masked = mask_field(xfield, mask_type, domain, mask)

    # no time dimensions
    notimedim = [dim for dim in xfield.dims if dim != 'time']

    # global mean
    if operation in ['average', 'mean']:
        out = masked.weighted(weights.fillna(0)).mean(dim=notimedim).data
    # global integrals
    elif operation in ['integral', 'sum']:
        out = masked.weighted(weights.fillna(0)).sum(dim=notimedim).data
    else:
        sys.exit("ERROR: masked_meansum-> mask undefined, this cannot be handled!")

    return out


# def _merge_mask(var, xfield, mask):

#     # check that we are receiving a dataset and not a datarray
#     if isinstance(xfield, xr.DataArray):
#         xfield = xfield.to_dataset(name=var)

#     # convert from datarray to dataset and merge
#     mask = mask.to_dataset(name='mask')

#     # the compat='override' option forces the merging. some CMIP6 data might
#     # have different float type, this simplies the handling
#     bfield = xr.merge([xfield, mask], compat='override')

#     return bfield

def mask_field(xfield, mask_type, dom, mask):
    """Apply a land/sea mask on a xarray variable var"""

    # nothing to be done
    if mask_type == 'global':
        out = xfield
    # northern and southern hemisphere
    elif mask_type == 'north':
        out = xfield.where(xfield['lat'] > 0)
    elif mask_type == 'south':
        out = xfield.where(xfield['lat'] < 0)
    else:
        # if oceanic, apply the ocenanic mask (if it exists!)
        if dom == 'oce':
            if isinstance(mask, xr.DataArray):
                out = xfield.where(mask.data < 0.5)
            else:
                logging.warning('No oceanic mask available for oceanic vars, this might lead to inconsistent results...')
                out = xfield
        elif dom == 'atm':

            # conditions
            if mask_type == 'land':
                out = xfield.where(mask.data >= 0.5)
            elif mask_type in ['sea', 'ocean']:
                out = xfield.where(mask.data < 0.5)
            else:
                sys.exit("ERROR: mask_field -> Mask undefined, this cannot be handled!")

    return out

# def mask_field_old(xfield, var, mask_type, dom, mask):
#     """Apply a land/sea mask on a xarray variable var"""

#     # if oceanic, apply the ocenanic mask (if it exists!)
#     if dom == 'oce' and isinstance(mask, xr.DataArray):
#         bfield = _merge_mask(var, xfield, mask)
#         xfield = bfield[var].where(bfield['mask'] < 0.5)
#         # xfield.to_netcdf(var+'okmasked.nc')

#     # nothing to be done
#     if mask_type == 'global':
#         out = xfield
#     elif mask_type == 'north':
#         out = xfield.where(xfield['lat'] > 0)
#     elif mask_type == 'south':
#         out = xfield.where(xfield['lat'] < 0)
#     else:

#         bfield = _merge_mask(var, xfield, mask)
#         # conditions
#         if mask_type == 'land':
#             out = bfield[var].where(bfield['mask'] >= 0.5)
#         elif mask_type in ['sea', 'ocean']:
#             out = bfield[var].where(bfield['mask'] < 0.5)
#         else:
#             sys.exit("ERROR: mask_field -> Mask undefined, this cannot be handled!")

#     # out.to_netcdf(var+'masked.nc')

#     return out


def select_region(xfield, region):
    """Trivial function to convert region definition to xarray
    sliced array to compute the PIs or global means on selected regions"""

    # fixed for the order of latitudes
    #xfield = xfield.sortby('lat')

    if region == 'Global':
        return xfield
    else: 
        if region == 'North Midlat':
            lat_min, lat_max = 30.0, 90.0
        elif region == 'South Midlat':
            lat_min, lat_max = -90.0, -30.0 
        elif region == 'Tropical':
            lat_min, lat_max = -30.0, 30.0
        else:
            sys.exit(region + "region not supported!!!")
    
        return xfield.where((xfield.lat >= lat_min) & (xfield.lat <= lat_max))
         

def select_region_old(xfield, region):
    """Trivial function to convert region definition to xarray
    sliced array to compute the PIs or global means on selected regions"""

    # fixed for the order of latitudes
    #xfield = xfield.sortby('lat')


    if region == 'Global':
        slicearray = xfield
    elif region == 'North Midlat':
        slicearray = xfield.sel(lat=slice(30, 90))
    elif region == 'South Midlat':
        slicearray = xfield.sel(lat=slice(-90, -30))
    elif region == 'Tropical':
        slicearray = xfield.sel(lat=slice(-30, 30))
    else:
        sys.exit(region + "region not supported!!!")

    return slicearray