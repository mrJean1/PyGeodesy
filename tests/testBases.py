
# -*- coding: utf-8 -*-

# Test base classes.

__version__ = '16.12.07'

if __name__ == '__main__':

    from tests import Tests

    from geodesy import bases as B
    t = Tests(__file__, __version__, B)
    t.testBases(B._LatLonHeightBase)
    t.results()
    t.exit()

    # Typical test results (on MacOS X):

    # testing geodesy.bases version 16.11.11
    # test 1 lat, lon: 50.06632°N, 005.71475°W
    # test 2 lat, lon: 50.066389°N, 005.714722°W
    # test 3 equals: True
    # test 4 toStr: 51°28′40″N, 000°00′06″W
    # test 5 toStr: 51.4778°N, 000.0016°W
    # test 6 toStr: 51°28′40″N, 000°00′06″W, +42.00m
    # all geodesy.bases tests passed (Python 2.7.10)

    # testing geodesy.bases version 16.11.11
    # test 1 lat, lon: 50.06632°N, 005.71475°W
    # test 2 lat, lon: 50.066389°N, 005.714722°W
    # test 3 equals: True
    # test 4 toStr: 51°28′40″N, 000°00′06″W
    # test 5 toStr: 51.4778°N, 000.0016°W
    # test 6 toStr: 51°28′40″N, 000°00′06″W, +42.00m
    # all geodesy.bases tests passed (Python 3.5.2)
