"""Tests for performance indices functions"""

# test for PIs: run on ECE4 test data, both amip and coupled, and on CMIP6 EC-Earth3 data.
# all run for both EC23 and RK08 climatologies
import os
import subprocess
from filecmp import cmp
import pytest
import xarray as xr
from ecmean.performance_indices import performance_indices
from ecmean.libs.ncfixers import xr_preproc

# set up coverage env var
env = {**os.environ, "COVERAGE_PROCESS_START": ".coveragerc"}

# test on coupled
@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_performance_indices_cpld(clim):

    performance_indices('cpld', 1990, 1990, numproc=4, climatology=clim, config='tests/config.yml')
    assert cmp('tests/table/PI4_' + clim + '_cpld_EC-Earth4_r1i1p1f1_1990_1990.yml',
               'tests/table/PI4_' + clim + '_cpld_1990_1990.ref')

# test on amip
@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_performance_indices_amip(clim):
    performance_indices('amip', 1990, 1990, numproc=1, climatology=clim, config='tests/config.yml', outputdir='pluto')

    assert cmp('pluto/YAML/PI4_' + clim + '_amip_EC-Earth4_r1i1p1f1_1990_1990.yml',
               'tests/table/PI4_' + clim + '_amip_1990_1990.ref')

# test performance_indices from commnand line + debug
@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_cmd_performance_indices_CMIP6(clim):
    subprocess.run(['performance_indices', 'historical', '1990', '1990', '-j', '2', '-c',
                   'tests/config_CMIP6.yml', '-k', clim, '-m', 'EC-Earth3', '-v', 'debug'],
                    env=env, check=True)
    assert cmp('tests/table/PI4_' + clim + '_historical_EC-Earth3_r1i1p1f1_1990_1990.yml',
               'tests/table/PI4_' + clim + '_CMIP6_1990_1990.ref')

# test on amip but with access from xarray dataset


@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_performance_indices_amip_xdataset(clim):
    thefile = 'tests/table/PI4_' + clim + '_amip_EC-Earth4_r1i1p1f1_1990_1990.yml'
    if os.path.isfile(thefile):
        os.remove(thefile)
    xfield = xr.open_mfdataset('tests/data/amip/output/oifs/*.nc', preprocess=xr_preproc)
    performance_indices('amip', 1990, 1990, numproc=2, climatology=clim,
                        config='tests/config.yml', xdataset=xfield)
    assert cmp(thefile, 'tests/table/PI4_' + clim + '_amip_1990_1990.ref')

# @pytest.mark.parametrize("clim", ['RK08', 'EC23'])
# def test_performance_indices_CMIP6(clim):
#    performance_indices('historical', 1990, 1990, numproc = 2, climatology = clim, config = 'tests/config_CMIP6.yml')
#    assert cmp('tests/table/PI4_' + clim + '_historical_EC-Earth3_r1i1p1f1_1990_1990.yml',
#               'tests/table/PI4_' + clim + '_CMIP6_1990_1990.ref')
