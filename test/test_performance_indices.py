from filecmp import cmp
from performance_indices import pi_main


def test_performance_indices_cpld():
    pi_main(['cpld', '1990', '1990', '-j', '2', '-c', 'test/config.yml'])
    assert cmp('test/table/PI4_RK08_cpld_EC-Earth4_r1i1p1f1_1990_1990.txt', 'test/table/PI4_RK08_cpld_1990_1990.ref')


def test_performance_indices_amip():
    pi_main(['amip', '1990', '1990', '-j', '2', '-c', 'test/config.yml'])
    assert cmp('test/table/PI4_RK08_amip_EC-Earth4_r1i1p1f1_1990_1990.txt', 'test/table/PI4_RK08_amip_1990_1990.ref')


def test_performance_indices_CMIP6():
    pi_main(['historical', '1990', '1990', '-j', '2', '-c', 'test/config_CMIP6.yml'])
    assert cmp('test/table/PI4_RK08_historical_EC-Earth3_r1i1p1f1_1990_1990.txt', 'test/table/PI4_RK08_CMIP6_1990_1990.ref')
