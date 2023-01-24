#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

import logging
import sys
import xarray as xr
from ecmean.libs.ncfixers import xr_preproc

##########################
# MASK-RELATED FUNCTIONS #
##########################


def masks_dictionary(component, maskatmfile, remap_dictionary=None):
    """Create a dictionary with atmospheric land-sea mask"""

    mask = {
        'atm_mask': _make_atm_masks(
            component['atm'],
            maskatmfile,
            remap_dictionary=remap_dictionary),
    }

    return mask


def _make_atm_masks(component, maskatmfile, remap_dictionary=None):
    """Create land-sea masks for atmosphere model"""

    # prepare ATM LSM: this needs to be improved, since it is clearly model
    # dependent
    logging.debug('maskatmfile is' + maskatmfile)
    if not maskatmfile:
        sys.exit("ERROR: maskatmfile cannot be found")

    if component == 'oifs':
        # create mask: opening a grib and loading only lsm to avoid
        # inconsistencies # in the grib structure ->
        # see here https://github.com/ecmwf/cfgrib/issues/13
        mask = xr.open_mfdataset(
            maskatmfile,
            engine="cfgrib",
            filter_by_keys={
                'shortName': 'lsm'},
            preprocess=xr_preproc)['lsm']
        
    elif component in ['cmoratm', 'globo']:
        dmask = xr.open_mfdataset(maskatmfile, preprocess=xr_preproc)
        if 'sftlf' in dmask.data_vars : 
            mask = dmask['sftlf']
        elif 'lsm' in dmask.data_vars : 
            mask = dmask['lsm']

        #globo has a reversed mask
        if component == 'globo':
            mask = abs(1-mask)
    else:
        sys.exit("ERROR: _make_atm_masks -> Mask undefined yet mismatch, this cannot be handled!")

    if remap_dictionary is not None:
        if remap_dictionary['atm_fix']:
            mask = remap_dictionary['atm_fix'](mask, keep_attrs=True)
        mask = remap_dictionary['atm_remap'](mask, keep_attrs=True)

    return mask


def masked_meansum(xfield, var, weights, mask_type, mask):
    """For global variables evaluate the weighted averaged
    or weighted integral when required by the variable properties"""

    # call the mask_field to mask where necessary
    # the mask field is area
    masked = mask_field(xfield, var, mask_type, mask)

    # global mean
    if mask_type in ['global']:
        out = masked.weighted(weights.fillna(0)).mean().values
    # global integrals
    elif mask_type in ['land', 'ocean', 'sea', 'north', 'south']:
        out = masked.weighted(weights.fillna(0)).sum().values
    else:
        sys.exit("ERROR: masked_meansum-> mask undefined, this cannot be handled!")

    return float(out)


def mask_field(xfield, var, mask_type, mask):
    """Apply a land/sea mask on a xarray variable var"""

    # nothing to be done
    if mask_type == 'global':
        out = xfield
    elif mask_type == 'north':
        out = xfield.where(xfield['lat'] > 0)
    elif mask_type == 'south':
        out = xfield.where(xfield['lat'] < 0)
    else:
        # check that we are receiving a dataset and not a datarray
        if isinstance(xfield, xr.DataArray):
            xfield = xfield.to_dataset(name=var)

        # convert from datarray to dataset and merge
        mask = mask.to_dataset(name='mask')

        # the compat='override' option forces the merging. some CMIP6 data might
        # have different float type, this simplies the handling
        bfield = xr.merge([xfield, mask], compat='override')

        # conditions
        if mask_type == 'land':
            out = bfield[var].where(bfield['mask'] >= 0.5)
        elif mask_type in ['sea', 'ocean']:
            out = bfield[var].where(bfield['mask'] < 0.5)
        else:
            sys.exit("ERROR: mask_field -> Mask undefined, this cannot be handled!")

    return out


def select_region(xfield, region):
    """Trivial function to convert region definition to xarray
    sliced array to compute the PIs or global means on selected regions"""

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


