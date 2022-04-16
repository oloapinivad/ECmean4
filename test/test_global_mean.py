from global_mean import main

def test_global_mean():
    main(['BETA', '1990', '1990', '-j', '2'])
    assert True
