"""Tests for performance indices functions"""

# test for PIs: run on ECE4 test data, both amip and coupled, and on CMIP6 EC-Earth3 data.
# all run for both EC23 and RK08 climatologies
import os
import subprocess
import pytest
import xarray as xr
import yaml
from ecmean.performance_indices import performance_indices
from ecmean.libs.ncfixers import xr_preproc
from ecmean.libs.general import are_dicts_equal

# set TOLERANCE
TOLERANCE = 1e-3

# set up coverage env var
env = {**os.environ, "COVERAGE_PROCESS_START": ".coveragerc"}

# test on coupled
@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_performance_indices_cpld(clim):
    performance_indices('cpld', 1990, 1990, numproc=4, climatology=clim, config='tests/config.yml')
    file1 = 'tests/table/PI4_' + clim + '_cpld_EC-Earth4_r1i1p1f1_1990_1990.yml'
    file2 = 'tests/table/PI4_' + clim + '_cpld_1990_1990.ref'
    with open(file1, 'r', encoding='utf8') as f1, open(file2, 'r', encoding='utf8') as f2:
        data1 = yaml.safe_load(f1)
        data2 = yaml.safe_load(f2)

    assert are_dicts_equal(data1, data2, TOLERANCE), "YAML files are not identical."


# test on amip
@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_performance_indices_amip(clim):
    performance_indices('amip', 1990, 1990, numproc=1, climatology=clim, config='tests/config.yml', outputdir='tests/pluto')
    file1 = 'tests/pluto/YAML/PI4_' + clim + '_amip_EC-Earth4_r1i1p1f1_1990_1990.yml'
    file2 = 'tests/table/PI4_' + clim + '_amip_1990_1990.ref'
    with open(file1, 'r', encoding='utf8') as f1, open(file2, 'r', encoding='utf8') as f2:
        data1 = yaml.safe_load(f1)
        data2 = yaml.safe_load(f2)

    assert are_dicts_equal(data1, data2, TOLERANCE), "YAML files are not identical."

# test performance_indices from commnand line + debug
@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_cmd_performance_indices_CMIP6(clim):
    subprocess.run(['performance_indices', 'historical', '1990', '1990', '-j', '2', '-c',
                   'tests/config_CMIP6.yml', '-k', clim, '-m', 'EC-Earth3', '-v', 'debug'],
                    env=env, check=True)
    file1 = 'tests/table/PI4_' + clim + '_historical_EC-Earth3_r1i1p1f1_1990_1990.yml'
    file2 = 'tests/table/PI4_' + clim + '_CMIP6_1990_1990.ref'
    with open(file1, 'r', encoding='utf8') as f1, open(file2, 'r', encoding='utf8') as f2:
        data1 = yaml.safe_load(f1)
        data2 = yaml.safe_load(f2)

    assert are_dicts_equal(data1, data2, TOLERANCE), "YAML files are not identical."

# test on amip but with access from xarray dataset
@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_performance_indices_amip_xdataset(clim):
    file1 = 'tests/table/PI4_' + clim + '_amip_EC-Earth4_r1i1p1f1_1990_1990.yml'
    file2 = 'tests/table/PI4_' + clim + '_amip_1990_1990.ref'
    if os.path.isfile(file1):
        os.remove(file1)
    xfield = xr.open_mfdataset('tests/data/amip/output/oifs/*.nc', preprocess=xr_preproc)
    performance_indices('amip', 1990, 1990, numproc=2, climatology=clim,
                        config='tests/config.yml', xdataset=xfield)
    with open(file1, 'r', encoding='utf8') as f1, open(file2, 'r', encoding='utf8') as f2:
        data1 = yaml.safe_load(f1)
        data2 = yaml.safe_load(f2)

    assert are_dicts_equal(data1, data2, TOLERANCE), "YAML files are not identical."



