from filecmp import cmp
from global_mean import global_mean


def test_global_mean_coupled():
    #global_mean(['cpld', '1990', '1990', '-j', '2', '-c', 'test/config.yml'])
    global_mean(exp = 'cpld', year1 = 1990, year2 = 1990, numproc = 2, config = 'test/config.yml')
    assert cmp('test/table/global_mean_cpld_EC-Earth4_r1i1p1f1_1990_1990.txt', 'test/table/global_mean_cpld_1990_1990.ref')


def test_global_mean_amip():
    #global_mean(['amip', '1990', '1990', '-j', '2', '-c', 'test/config.yml'])
    global_mean(exp = 'amip', year1 = 1990, year2 = 1990, numproc = 2, config = 'test/config.yml')
    assert cmp('test/table/global_mean_amip_EC-Earth4_r1i1p1f1_1990_1990.txt', 'test/table/global_mean_amip_1990_1990.ref')


def test_global_mean_CMIP6():
    #global_mean(['historical', '1990', '1990', '-j', '2', '-c', 'test/config_CMIP6.yml'])
    global_mean(exp = 'historical', year1 = 1990, year2 = 1990, numproc = 2, config = 'test/config_CMIP6.yml')
    assert cmp('test/table/global_mean_historical_EC-Earth3_r1i1p1f1_1990_1990.txt',
               'test/table/global_mean_CMIP6_1990_1990.ref')
