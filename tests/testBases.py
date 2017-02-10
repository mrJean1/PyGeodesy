
# -*- coding: utf-8 -*-

# Test base classes.

__version__ = '17.02.09'

if __name__ == '__main__':

    from tests import Tests

    from geodesy import bases  # private

    t = Tests(__file__, __version__, bases)
    t.testBases(bases.LatLonHeightBase)
    t.results()
    t.exit()

    # Typical test results (on MacOS 10.12.2):

    # testing geodesy.bases version 17.02.09
    # test 1 lat, lon: 50.06632°N, 005.71475°W
    # test 2 lat, lon: 50.066389°N, 005.714722°W
    # test 3 equals: True
    # test 4 toStr: 51°28′40″N, 000°00′06″W
    # test 5 toStr: 51.4778°N, 000.0016°W
    # test 6 precision: 0
    # test 7 toStr: 51°28′40″N, 000°00′06″W, +42.00m
    # all geodesy.bases tests passed (Python 2.7.13 64bit)

    # testing geodesy.bases version 17.02.09
    # test 1 lat, lon: 50.06632°N, 005.71475°W
    # test 2 lat, lon: 50.066389°N, 005.714722°W
    # test 3 equals: True
    # test 4 toStr: 51°28′40″N, 000°00′06″W
    # test 5 toStr: 51.4778°N, 000.0016°W
    # test 6 precision: 0
    # test 7 toStr: 51°28′40″N, 000°00′06″W, +42.00m
    # all geodesy.bases tests passed (Python 3.6.0 64bit)
