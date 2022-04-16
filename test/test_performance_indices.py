from performance_indices import main

def test_performance_indices():
    main(['BETA', '1990', '1990', '-j', '2'])
    assert True
