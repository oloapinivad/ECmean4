#!/usr/bin/env python3
'''
Shared functions for ECmean4
'''

#############################
# UNIT ADJUSTMENT FUNCTIONS #
#############################

from metpy.units import units
import logging
import sys

def units_extra_definition():
    """Add units to the pint registry required by ECMean4"""

    # special units definition
    #needed to work with metpy 1.4.0 see
    #https://github.com/Unidata/MetPy/issues/2884
    units._on_redefinition = 'warn' 
    units.define('fraction = [] = frac')
    units.define('psu = 1e-3 frac')
    units.define('PSU = 1e-3 frac')
    units.define('Sv = 1e+6 m^3/s')  # Replace Sievert with Sverdrup


def units_converter(org_units, tgt_units):
    """Units conversion using metpy and pint.
    From a org_units convert to tgt_units providing offset and factor.
    Some assumptions are done for precipitation field: must be extended
    to other vars. It will not work if BOTH factor and offset are required"""

    units_relation = (units(org_units) / units(tgt_units)).to_base_units()
    logging.debug(units_relation)
    if units_relation.magnitude != 1:
        logging.info('Unit conversion required...')
        offset_standard = 0 * units(org_units)
        factor_standard = 1 * units(org_units)
        if units_relation.units == units('dimensionless'):
            offset = offset_standard.to(tgt_units).magnitude
            if offset == 0:
                factor = factor_standard.to(tgt_units).magnitude
            else:
                factor = 1.

        elif units_relation.units == units('kg / m^3'):
            logging.info("Assuming this as a water flux! Am I correct?")
            logging.info("Dividing by water density...")
            density_water = units('kg / m^3') * 1000
            offset = 0.
            factor = (factor_standard / density_water).to(tgt_units).magnitude

        else:
            logging.error(units_relation)
            sys.exit("ERROR: Units mismatch, this cannot be handled!")
    else:
        offset = 0.
        factor = 1.

    logging.info('Offset is ' + str(offset))
    logging.info('Factor is ' + str(factor))
    return offset, factor


def units_are_integrals(org_units, ref_var):
    """Check functions for spatially integrated variables"""

    if 'total' in ref_var.keys():
        new_units = str((units(org_units) * units('m^2')).units)
    else:
        new_units = org_units
    return new_units


def directions_match(org, dst):
    """Check function for fluxes direction: they should match. Default is down"""

    direction_org = org.get('direction', 'down')
    direction_dst = dst.get('direction', 'down')
    if direction_org != direction_dst:
        factor = -1.
    else:
        factor = 1.
    return factor