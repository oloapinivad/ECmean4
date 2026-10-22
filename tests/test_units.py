"""Tests for the units conversion helpers"""

import glob
import logging

import numpy as np
import pytest

from ecmean.libs.files import load_yaml
from ecmean.libs.general import Timer
from ecmean.libs.units import UnitsHandler, apply_units_conversion, units_extra_definition

INTERFACES = sorted(glob.glob("ecmean/interfaces/interface_*.yml"))
REFERENCES = sorted(glob.glob("ecmean/reference/gm_reference_*.yml"))


def test_apply_units_conversion():
    """Order is data * factor + offset"""
    np.testing.assert_allclose(apply_units_conversion(np.array([1.0, 2.0]), 273.15, 1.0), [274.15, 275.15])
    np.testing.assert_allclose(apply_units_conversion(np.array([1.0, 2.0]), 0.0, -2.0), [-2.0, -4.0])


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


def test_timer(caplog):
    """Timer logs elapsed and total time"""
    timer = Timer(logging.getLogger("timer-test"))
    with caplog.at_level(logging.INFO, logger="timer-test"):
        timer.toc("Step")
        timer.final_toc()
    assert "Step time:" in caplog.text
    assert "Total execution time:" in caplog.text
