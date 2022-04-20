from performance_indices import main
from filecmp import cmp

def test_performance_indices():
    main(['mock', '1990', '1990', '-j', '2', '-c', 'test/config.yml'])
    assert cmp('test/table/global_mean_mock_1990_1990.txt', 'test/table/global_mean_mock_1990_1990.ref')
