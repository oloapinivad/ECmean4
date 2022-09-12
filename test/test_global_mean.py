import pytest
from global_mean import main
from filecmp import cmp


#def test_global_mean_coupled():
#    main(['cpld', '1990', '1990', '-j', '2', '-c', 'test/config.yml'])
#    assert cmp('test/table/global_mean_cpld_EC-Earth4_r1i1p1f1_1990_1990.txt', 'test/table/global_mean_cpld_1990_1990.ref')


#def test_global_mean_amip():
#    main(['amip', '1990', '1990', '-j', '2', '-c', 'test/config.yml'])
#    assert cmp('test/table/global_mean_amip_EC-Earth4_r1i1p1f1_1990_1990.txt', 'test/table/global_mean_amip_1990_1990.ref')


def test_global_mean_CMIP6():
    main(['historical', '1990', '1990', '-j', '2', '-c', 'test/config_CMIP6.yml'])
    assert cmp('test/table/global_mean_historical_EC-Earth3_r1i1p1f1_1990_1990.txt', 'test/table/global_mean_CMIP6_1990_1990.ref')
