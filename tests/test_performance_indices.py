# test for PIs: run on ECE4 test data, both amip and coupled, and on CMIP6 EC-Earth3 data.
# all run for both EC23 and RK08 climatologies

import pytest
from filecmp import cmp
import subprocess
from ecmean.performance_indices import performance_indices


@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_performance_indices_cpld(clim):

    performance_indices('cpld', 1990, 1990, numproc=2, climatology=clim, config='tests/config.yml')
    assert cmp('tests/table/PI4_' + clim + '_cpld_EC-Earth4_r1i1p1f1_1990_1990.yml',
               'tests/table/PI4_' + clim + '_cpld_1990_1990.ref')


@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_performance_indices_amip(clim):
    performance_indices('amip', 1990, 1990, numproc=2, climatology=clim, config='tests/config.yml')

    assert cmp('tests/table/PI4_' + clim + '_amip_EC-Earth4_r1i1p1f1_1990_1990.yml',
               'tests/table/PI4_' + clim + '_amip_1990_1990.ref')

# test performance_indices from commnand line + debug


@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_cmd_performance_indices_CMIP6(clim):
    subprocess.run(['performance_indices', 'historical', '1990', '1990', '-j', '2', '-c', 'tests/config_CMIP6.yml',
                    '-k', clim, '-m', 'EC-Earth4', '-v', 'debug'])
    assert cmp('tests/table/PI4_' + clim + '_historical_EC-Earth3_r1i1p1f1_1990_1990.yml',
               'tests/table/PI4_' + clim + '_CMIP6_1990_1990.ref')


# @pytest.mark.parametrize("clim", ['RK08', 'EC23'])
# def test_performance_indices_CMIP6(clim):
#    performance_indices('historical', 1990, 1990, numproc = 2, climatology = clim, config = 'tests/config_CMIP6.yml')
#    assert cmp('tests/table/PI4_' + clim + '_historical_EC-Earth3_r1i1p1f1_1990_1990.yml',
#               'tests/table/PI4_' + clim + '_CMIP6_1990_1990.ref')
