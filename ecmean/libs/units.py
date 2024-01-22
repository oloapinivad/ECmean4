#!/usr/bin/env python3
'''
Shared class for ECmean4 units
'''

#############################
# UNIT ADJUSTMENT FUNCTIONS #
#############################

import logging
from metpy.units import units

loggy = logging.getLogger(__name__)

class UnitsHandler():
    """
    Class for dealing with units format and conversion within ECmean4
    """

    def __init__(self, var=None, org_units=None, tgt_units=None,
                 org_direction='down', tgt_direction='down',
                 operation='mean', cumulation_time=None,
                 clim=None, face=None, convert=True):
        """
        Class for handling units and their conversion

        Args:
            var (string): variable name, necessary for ECmean4 calls
            org_units (string): original units
            tgt_units (string): target units
            org_direction (string): down or up, direction for fluxes
            tgt_direction (string): down or up, direction for fluxes
            operation (string): mean or integral operation for target units
            cumulation_time (integer): expressed in seconds, cumulation time for fluxes
            clim (dict): climatology dictionary for ECmean4 calls
            face (dict): interface dictionary for Ecmean4 calss
            convert (bool): if you want to perform conversion immediately

        Returns:
            offset (float): units offset of conversion
            factor (float): units factor of conversion
        """

        # generic object initialization
        self.var = var
        self.org_units = org_units
        self.tgt_units = tgt_units
        self.org_direction = org_direction
        self.tgt_direction = tgt_direction
        self.cumulation_time = cumulation_time
        self.operation = operation

        # init in the format from ecmean
        if face and clim:
            self.init_ecmean(clim, face)

        loggy.debug(vars(self))
        # if units are defined and convert called
        if convert and self.org_units and self.tgt_units:

            self.parse_units()
            self.units_are_integrals()
            self.offset, self.factor = self.units_converter()

    def init_ecmean(self, clim, face):
        """Specific initialization used within ECmean"""

        if not self.org_units:
            loggy.warning('Source unit for %s undefined, assuming fraction', self.var)
            self.org_units = 'frac'

        if clim[self.var]['units']:
            self.tgt_units = clim[self.var]['units']
        else:
            loggy.warning('Target unit undefined, assuming fraction')
            self.tgt_units = 'frac'

        self.org_direction = face['variables'][self.var].get('direction', 'down')
        self.tgt_direction = clim[self.var].get('direction', 'down')
        self.operation = clim[self.var].get('operation', 'mean')
        self.cumulation_time = face['variables'][self.var].get('cumulation_time', None)

    def parse_units(self):

        """Parse units once defined with metpy.units (i.e. pint)"""

        self.org_units = units(self.org_units)
        self.tgt_units = units(self.tgt_units)


    def units_are_integrals(self):
        """Check functions for spatially integrated variables"""

        if self.operation in ['integral', 'sum']:
            self.org_units = self.org_units * units('m^2')

    def units_converter(self):
        """Units conversion using metpy and pint.
        From a org_units convert to tgt_units providing offset and factor.
        Some assumptions are done for water fluxes and cumulated fluxes.
        It will not work if BOTH factor and offset are required"""

        units_relation = (self.org_units / self.tgt_units).to_base_units()
        loggy.debug('Original vs Target units relation is %s', units_relation)

        if units_relation.units == units('dimensionless'):

            if units_relation.magnitude != 1:
                loggy.debug('Unit conversion required...')
                offset_standard = 0 * self.org_units
                offset = offset_standard.to(self.tgt_units).magnitude
                if offset == 0:
                    factor_standard = 1 * self.org_units
                    factor = factor_standard.to(self.tgt_units).magnitude
                else:
                    factor = 1.
            else:
                offset = 0.
                factor = 1.

        elif units_relation.units == units('kg / m^3'):
            loggy.debug("Assuming this as a water flux! Am I correct?")
            loggy.debug("Dividing by water density...")
            density_water = units('kg / m^3') * 1000
            offset = 0.
            factor_standard = 1 * self.org_units
            factor = (factor_standard / density_water).to(self.tgt_units).magnitude

        elif units_relation.units == units('s'):
            loggy.debug("Assuming this is a cumulated flux...")
            loggy.debug("Dividing by cumulation time (expressed in seconds)...")
            if self.cumulation_time:
                cumtime = units('s') * self.cumulation_time
                offset = 0.
                factor_standard = 1 * self.org_units
                factor = (factor_standard / cumtime).to(self.tgt_units).magnitude
            else:
                loggy.error('This variable seems cumulated over time but has no cumulation time defined')
                raise ValueError("ERROR: Units mismatch, this cannot be handled!")

        else:
            loggy.error(units_relation)
            raise ValueError("ERROR: Units mismatch, this cannot be handled!")

        if self.org_direction != self.tgt_direction:
            factor = -1. * factor

        loggy.debug('Offset is %s', str(offset))
        loggy.debug('Factor is %s', str(factor))
        return offset, factor



def units_extra_definition():
    """Add units to the pint registry required by ECMean4"""

    # special units definition
    # needed to work with metpy 1.4.0 see
    # https://github.com/Unidata/MetPy/issues/2884
    units._on_redefinition = 'ignore'
    units.define('fraction = [] = frac')
    units.define('Fraction = [] = frac')
    units.define('psu = 1e-3 frac')
    units.define('PSU = 1e-3 frac')
    units.define('million = 1e6 = M')
    units.define('Sv = 1e+6 m^3/s')  # Replace Sievert with Sverdrup
