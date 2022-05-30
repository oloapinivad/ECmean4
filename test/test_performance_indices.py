from performance_indices import main
from filecmp import cmp


def test_performance_indices_cpld():
    main(['cpld', '1990', '1990', '-j', '2', '-c', 'test/config.yml'])
    assert cmp('test/table/PI4_RK08_cpld_EC-Earth4_r1i1p1f1_1990_1990.txt', 'test/table/PI4_RK08_cpld_1990_1990.ref')


def test_performance_indices_amip():
    main(['amip', '1990', '1990', '-j', '2', '-c', 'test/config.yml'])
    assert cmp('test/table/PI4_RK08_amip_EC-Earth4_r1i1p1f1_1990_1990.txt', 'test/table/PI4_RK08_amip_1990_1990.ref')
