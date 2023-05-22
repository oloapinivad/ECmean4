#!/usr/bin/env python3
'''
Shared functions for Support class for ECmean4
'''

import logging
import xarray as xr
import xesmf as xe
import numpy as np
from ecmean.libs.ncfixers import xr_preproc
from ecmean.libs.files import inifiles_priority
from ecmean.libs.areas import area_cell
from ecmean.libs.units import UnitsHandler

loggy = logging.getLogger(__name__)

class Supporter():

    """
    Support class for ECmean4, including areas and masks to be used
    in global mean and performance indices
    """

    def __init__(self, component, atmdict, ocedict, areas=True, remap=False, targetgrid=None):
        """Class for masks, areas and interpolation (xESMF-based)
        for both atmospheric and oceanic component"""

        # define the basics
        self.atmareafile = inifiles_priority(atmdict)
        self.oceareafile = inifiles_priority(ocedict)
        self.atmmaskfile = atmdict['maskfile']
        self.ocemaskfile = ocedict['maskfile']
        self.atmcomponent = component['atm']
        self.ocecomponent = component['oce']
        self.targetgrid = targetgrid

        # remapping default
        self.ocefix, self.oceremap = None, None
        self.atmfix, self.atmremap = None, None

        # areas and mask for amip case
        self.ocemask = None
        self.ocearea = None

        # loading and examining atmospheric file
        self.atmfield = self.load_field(self.atmareafile, comp='atm')
        self.atmgridtype = identify_grid(self.atmfield)
        loggy.warning('Atmosphere grid is is a %s grid!', self.atmgridtype)

        # compute atmopheric area
        if areas:
            self.atmarea = self.make_areas(self.atmgridtype, self.atmfield)

        # initialize the interpolation for atmosphere
        if self.targetgrid and remap:
            self.atmfix, self.atmremap = self.make_atm_interp_weights(self.atmfield)

        # init the land-sea mask for atm (mandatory)
        self.atmmask = self.make_atm_masks()

        # do the same if oceanic file is found
        if self.oceareafile:
            self.ocefield = self.load_field(self.oceareafile, comp='oce')
            self.ocegridtype = identify_grid(self.ocefield)
            loggy.warning('Oceanic grid is is a %s grid!', self.ocegridtype)

            # compute oceanic area
            if areas:
                self.ocearea = self.make_areas(self.ocegridtype, self.ocefield)

            # init the ocean interpolation
            if self.targetgrid and remap:
                self.ocefix, self.oceremap = self.make_oce_interp_weights(self.ocefield)

            # ocean mask
            if self.ocemaskfile:
                self.ocemask = self.make_oce_masks()
            else:
                # if it is missing, when remapping I can use the atmospheric one
                if self.targetgrid and remap:
                    self.ocemask = self.atmmask
                # otherwise, no solution!
                else:
                    loggy.warning('Oceanic mask not found!')

        else:
            loggy.warning("Ocereafile cannot be found, assuming this is an AMIP run")

    def make_atm_masks(self):
        """Create land-sea masks for atmosphere model"""

        # prepare ATM LSM
        loggy.info('maskatmfile is %s', self.atmmaskfile)
        if not self.atmmaskfile:
            raise KeyError("ERROR: maskatmfile cannot be found")

        if self.atmcomponent == 'oifs':
            # create mask: opening a grib and loading only lsm to avoid
            # inconsistencies # in the grib structure ->
            # see here https://github.com/ecmwf/cfgrib/issues/13
            mask = xr.open_mfdataset(
                self.atmmaskfile,
                engine="cfgrib",
                indexpath=None,
                filter_by_keys={
                    'shortName': 'lsm'},
                preprocess=xr_preproc)['lsm']

        elif self.atmcomponent in ['cmoratm', 'globo']:
            dmask = xr.open_mfdataset(self.atmmaskfile, preprocess=xr_preproc)
            if 'sftlf' in dmask.data_vars:
                mask = dmask['sftlf']
                mask = mask / 100  # cmor mask are %
            elif 'lsm' in dmask.data_vars:
                mask = dmask['lsm']

            # globo has a reversed mask
            if self.atmcomponent == 'globo':
                mask = abs(1 - mask)
        else:
            raise KeyError("ERROR: _make_atm_masks -> Mask undefined yet mismatch, this cannot be handled!")

        # safe check to operate only on single timeframe
        if 'time' in mask.dims:
            mask = mask.isel(time=0)

        # interp the mask if required
        if self.atmremap is not None:
            if self.atmfix:
                mask = self.atmfix(mask, keep_attrs=True)
            mask = self.atmremap(mask, keep_attrs=True)

        return mask

    def make_oce_masks(self):
        """Create land-sea masks for oceanic model. This is used only for CMIP"""

        # prepare ocean LSM:
        loggy.info('maskocefile is %s', self.ocemaskfile)
        if not self.ocemaskfile:
            raise KeyError("ERROR: maskocefile cannot be found")

        if self.ocecomponent == 'cmoroce':
            dmask = xr.open_mfdataset(self.ocemaskfile, preprocess=xr_preproc)
            if 'sftof' in dmask.data_vars:
                # use fillna to have a binary max (0 land, 1 sea)
                mask = dmask['sftof'].fillna(0)

            # check if we need to convert from % to fraction
            # offset should not count!
            if mask.units:
                units_handler = UnitsHandler(org_units=mask.units, tgt_units='frac')
                offset, factor = units_handler.offset, units_handler.factor
                # offset, factor = units_converter(mask.units, 'frac')
                mask = (mask * factor) + offset

            # to keep coehrence in other operations, oceanic mask is set to be
            # identical to atmospheric mask, i.e. 1 over land and 0 over ocean
            mask = abs(1 - mask)

        else:
            raise KeyError("ERROR: _make_oce_masks -> Mask undefined yet mismatch, this cannot be handled!")

        # safe check to operate only on single timeframe
        if 'time' in mask.dims:
            mask = mask.isel(time=0)

        # interp the mask if required
        if self.oceremap is not None:
            if self.ocefix:
                mask = self.ocefix(mask, keep_attrs=True)
            mask = self.oceremap(mask, keep_attrs=True)

        return mask

    def load_field(self, areafile, comp):
        """Loading files for area and interpolation"""

        # loading and examining atmospheric file
        if areafile:
            loggy.info(f'{comp}mareafile is ' + areafile)
            if not areafile:
                raise FileExistsError(f'ERROR: {comp}reafile cannot be found')
            return xr.open_mfdataset(areafile, preprocess=xr_preproc).load()
        else:
            raise FileExistsError('ERROR: Cannot find any file to load! Does your experiment exit?')

    def make_areas(self, gridtype, xfield):
        """Create weights for area operations.
        Minimal structure."""

        # this might be universal, but keep this as for supported components only
        if 'areacello' in xfield.data_vars:  # as oceanic CMOR case
            area = xfield['areacello']
        elif 'cell_area' in xfield.data_vars:  # as ECE4 NEMO case for nemo-initial-state.nc
            area = xfield['cell_area']
        elif 'e1t' in xfield.data_vars:  # ECE4 NEMO case for domaing_cfg.nc
            area = xfield['e1t'] * xfield['e2t']
        else:  # automatic solution, wish you luck!
            area = area_cell(xfield, gridtype)

        return area

    def make_atm_interp_weights(self, xfield):
        """Create atmospheric interpolator weights"""

        if self.atmcomponent == 'oifs':

            # this is to get lon and lat from the Equator
            xname = list(xfield.data_vars)[-1]
            m = xfield[xname].isel(time=0).load()
            # use numpy since it is faster
            g = np.unique(m.lat.data)
            f = np.unique(m.sel(cell=m.lat == g[int(len(g) / 2)]).lon.data)

            # this creates a a gaussian non reduced grid
            gaussian_regular = xr.Dataset({"lon": (["lon"], f), "lat": (["lat"], g)})

            # use nearest neighbour to remap to gaussian regular
            fix = xe.Regridder(
                xfield[xname], gaussian_regular,
                method="nearest_s2d", locstream_in=True,
                periodic=True)

            # create bilinear interpolator
            remap = xe.Regridder(
                fix(xfield[xname]), self.targetgrid,
                periodic=True, method="bilinear")

        elif self.atmcomponent in ['cmoratm', 'globo']:

            fix = None
            remap = xe.Regridder(
                xfield, self.targetgrid, periodic=True,
                method="bilinear")

        else:
            raise KeyError(
                "ERROR: Atm weights not defined for this component, this cannot be handled!")

        return fix, remap

    def make_oce_interp_weights(self, xfield):
        """Create oceanic interpolator weights"""

        if self.ocecomponent in ['nemo', 'cmoroce']:
            if 'areacello' in xfield.data_vars:  # CMOR case
                xname = 'areacello'
            elif 'cell_area' in xfield.data_vars:  # ECE4 NEMO case for nemo-initial-state.nc
                xname = 'cell_area'
            else:
                # tentative extraction
                xname = list(xfield.data_vars)[-1]
        else:
            raise KeyError(
                "ERROR: Oce weights not defined for this component, this cannot be handled!")

        if self.ocegridtype in ['unstructured']:
            # print("Detecting a unstructured grid, using nearest neighbour!")
            fix = None
            remap = xe.Regridder(
                xfield[xname],
                self.targetgrid,
                method="nearest_s2d",
                locstream_in=True,
                periodic=True)

        else:
            # print("Detecting regular or curvilinear grid, using bilinear!")
            fix = None
            remap = xe.Regridder(
                xfield[xname],
                self.targetgrid,
                method="bilinear",
                periodic=True,
                ignore_degenerate=True)

        return fix, remap


def identify_grid(xfield):
    """Receiveng an xarray object (DataArray or Dataset) investigates its coordinates
    and dimensions and provide the grid type (regular, gaussian, curvilinear,
    gaussian reduced, unstructured). It assumes that data is defined by 'lon' and 'lat'
    dimensions

    Args :
        xfield:

    Returns
        string with the grid type

    """

    # if the coordinates are lon/lat proceed
    if all(x in xfield.coords for x in ['lon', 'lat']):

        # if dimensions are lon/lat as well, this is a regular grid
        if all(x in xfield.dims for x in ['lon', 'lat']):
            lat = xfield.coords['lat']

            # if lat grid spacing is equal, is regular lonlat, otherwise gaussian
            if (lat[3] - lat[2]) == (lat[1] - lat[0]):
                gridtype = 'lonlat'
            else:
                gridtype = 'gaussian'
        else:
            # if the coords are 2D, we are curvilinear
            if xfield.coords['lon'].ndim == 2 and xfield.coords['lon'].ndim == 2:
                gridtype = 'curvilinear'
            else:
                # check the first four elements of the grid (arbitrary)
                lat = xfield.coords['lat'].values[0:5]

                # if they are all the same, we have a gaussian reduced, else unstructured
                if (lat == lat[0]).all():
                    gridtype = 'gaussian_reduced'
                else:
                    gridtype = 'unstructured'
    else:
        raise ValueError("Cannot find any lon/lat dimension, aborting...")

    return gridtype
