# test for PIs: run on ECE4 test data, both amip and coupled, and on CMIP6 EC-Earth3 data.
# all run for both EC23 and RK08 climatologies

import pytest
from filecmp import cmp
from performance_indices import pi_main


@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_performance_indices_cpld(clim):
    pi_main(['cpld', '1990', '1990', '-j', '2', '-k', clim, '-c', 'test/config.yml'])
    assert cmp('test/table/PI4_' + clim + '_cpld_EC-Earth4_r1i1p1f1_1990_1990.yml',
               'test/table/PI4_' + clim + '_cpld_1990_1990.ref')


@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_performance_indices_amip(clim):
    pi_main(['amip', '1990', '1990', '-j', '2', '-k', clim, '-c', 'test/config.yml'])
    assert cmp('test/table/PI4_' + clim + '_amip_EC-Earth4_r1i1p1f1_1990_1990.yml',
               'test/table/PI4_' + clim + '_amip_1990_1990.ref')


@pytest.mark.parametrize("clim", ['RK08', 'EC23'])
def test_performance_indices_CMIP6(clim):
    pi_main(['historical', '1990', '1990', '-j', '2', '-k', clim, '-c', 'test/config_CMIP6.yml'])
    assert cmp('test/table/PI4_' + clim + '_historical_EC-Earth3_r1i1p1f1_1990_1990.yml',
               'test/table/PI4_' + clim + '_CMIP6_1990_1990.ref')
