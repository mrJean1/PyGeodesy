
# -*- coding: utf-8 -*-

# Test datums, ellipsoids and transforms.

__version__ = '16.09.14'

if __name__ == '__main__':

    from tests import Tests

    import geodesy
    import geodesy.datum as D

    t = Tests(__file__, __version__, D)
    t.testDatum(geodesy)
    t.results()

    # Typical test results (on MacOS X):

    # testing datum.py version 16.09.14
    # test 1 ellipsoid: True
    # test 2 transform: True
    # test 3 datum: True
    # test 4 meanR: 6371008.77142
    # all datum.py tests passed (Python 2.7.10)

    # testing datum.py version 16.09.14
    # test 1 ellipsoid: True
    # test 2 transform: True
    # test 3 datum: True
    # test 4 meanR: 6371008.77142
    # all datum.py tests passed (Python 3.5.1)
