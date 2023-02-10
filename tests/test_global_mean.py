from filecmp import cmp
from ecmean.global_mean import global_mean
import subprocess
import os

# set up coverage env var
env = {**os.environ, "COVERAGE_PROCESS_START": ".coveragerc"}


# call on coupled ECE using parser and debug mode
def test_cmd_global_mean_coupled():
    thefile = 'tests/table/global_mean_cpld_EC-Earth4_r1i1p1f1_1990_1990.txt'
    if os.path.isfile(thefile):
        os.remove(thefile)
    subprocess.run(['global_mean', 'cpld', '1990', '1990', '-j', '2', '-c', 'tests/config.yml', '-t', '-v', 'debug'], env=env)
    assert cmp(thefile, 'tests/table/global_mean_cpld_1990_1990.ref')

# call on amip ECE
def test_global_mean_amip():
    thefile = 'tests/table/global_mean_amip_EC-Earth4_r1i1p1f1_1990_1990.txt'
    if os.path.isfile(thefile):
        os.remove(thefile)
    global_mean(exp='amip', year1=1990, year2=1990, numproc=2, config='tests/config.yml')
    assert cmp(thefile, 'tests/table/global_mean_amip_1990_1990.ref')

# call on historical CMIP6
def test_global_mean_CMIP6():
    thefile = 'tests/table/global_mean_historical_EC-Earth3_r1i1p1f1_1990_1991.txt'
    if os.path.isfile(thefile):
        os.remove(thefile)
    global_mean(exp='historical', year1=1990, year2=1991, numproc=2, config='tests/config_CMIP6.yml', trend=True)
    assert cmp(thefile, 'tests/table/global_mean_CMIP6_1990_1991.ref')
