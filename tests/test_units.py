"""Tests for the units conversion helpers"""

import glob

import pytest

from ecmean.libs.files import load_yaml
from ecmean.libs.units import UnitsHandler, units_extra_definition

INTERFACES = sorted(glob.glob("ecmean/interfaces/interface_*.yml"))
REFERENCES = sorted(glob.glob("ecmean/reference/gm_reference_*.yml"))


@pytest.mark.parametrize("reffile", REFERENCES)
@pytest.mark.parametrize("interface", INTERFACES)
def test_offset_and_factor_exclusive(interface, reffile):
    """units_converter never returns a non-zero offset together with a factor != 1,
    which is what makes the order of operations in apply_units_conversion immaterial"""
    units_extra_definition()
    face = load_yaml(interface)
    ref = load_yaml(reffile)
    for var, info in face["variables"].items():
        if var not in ref:
            continue
        try:
            offset, factor = UnitsHandler(var, org_units=info.get("units"), clim=ref, face=face).units_converter()
        except ValueError:
            continue  # unsupported pair, never converted by ecmean
        assert offset == 0 or abs(factor) == 1, f"{var}: offset {offset} and factor {factor}"
