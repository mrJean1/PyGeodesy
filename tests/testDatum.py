
# -*- coding: utf-8 -*-

# Test datums, ellipsoids and transforms.

__version__ = '16.10.13'

if __name__ == '__main__':

    from tests import Tests

    import geodesy
    import geodesy.datum as D

    t = Tests(__file__, __version__, D)
    t.testDatum(geodesy)
    t.results()

    # Typical test results (on MacOS X):

    # testing geodesy.datum version 16.10.10
    # test 1 ellipsoid: True
    # test 2 transform: True
    # test 3 datum: True
    # test 4 meanR: 6371008.77142
    # all geodesy.datum tests passed (Python 2.7.10)

    # testing geodesy.datum version 16.10.10
    # test 1 ellipsoid: True
    # test 2 transform: True
    # test 3 datum: True
    # test 4 meanR: 6371008.77142
    # all geodesy.datum tests passed (Python 3.5.1)
