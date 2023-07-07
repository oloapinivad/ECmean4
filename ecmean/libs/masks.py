#!/usr/bin/env python3
'''
Shared functions for XArray ECmean4
'''

import logging
import xarray as xr

##########################
# MASK-RELATED FUNCTIONS #
##########################

loggy = logging.getLogger(__name__)

def masked_meansum(xfield, weights, mask, operation, domain, mask_type):
    """For global variables evaluate the weighted averaged
    or weighted integral when required by the variable properties"""

    # call the mask_field to mask where necessary
    # the mask field is area
    masked = mask_field(xfield, mask_type, domain, mask)

    # no time dimensions
    notimedim = [dim for dim in xfield.dims if dim not in ['time', 'year', 'month']]

    # global mean
    if operation in ['average', 'mean']:
        out = masked.weighted(weights.fillna(0)).mean(dim=notimedim).data
    # global integrals
    elif operation in ['integral', 'sum']:
        out = masked.weighted(weights.fillna(0)).sum(dim=notimedim).data
    else:
        raise ValueError("ERROR: masked_meansum-> mask undefined, this cannot be handled!")

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
                out = xfield
        elif dom == 'atm':

            # conditions
            if mask_type == 'land':
                out = xfield.where(mask.data >= 0.5)
            elif mask_type == 'land-no-antarctica':
                out = xfield.where(mask.lat>(-60)).where(mask.data >= 0.5)
            elif mask_type in ['sea', 'ocean']:
                out = xfield.where(mask.data < 0.5)
            else:
                raise ValueError("ERROR: mask_field -> Mask undefined, this cannot be handled!")

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
            raise KeyError(region + "region not supported!!!")

        # new version more flexible than the slice one
        return xfield.where((xfield.lat >= lat_min) & (xfield.lat <= lat_max))


# def select_region_old(xfield, region):
#     """Trivial function to convert region definition to xarray
#     sliced array to compute the PIs or global means on selected regions"""

#     # fixed for the order of latitudes
#     #xfield = xfield.sortby('lat')


#     if region == 'Global':
#         slicearray = xfield
#     elif region == 'North Midlat':
#         slicearray = xfield.sel(lat=slice(30, 90))
#     elif region == 'South Midlat':
#         slicearray = xfield.sel(lat=slice(-90, -30))
#     elif region == 'Tropical':
#         slicearray = xfield.sel(lat=slice(-30, 30))
#     else:
#         sys.exit(region + "region not supported!!!")

#     return slicearray
