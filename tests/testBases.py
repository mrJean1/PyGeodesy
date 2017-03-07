
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '17.03.07'

from tests import Tests as _Tests

from geodesy import F_D, F_DMS, precision


class Tests(_Tests):

    def testBases(self, LatLon):

        p = LatLon(50.06632, -5.71475)
        self.test('lat, lon', p, '50.06632°N, 005.71475°W')
        q = LatLon('50°03′59″N', """005°42'53"W""")
        self.test('lat, lon', q, '50.066389°N, 005.714722°W')

        p = LatLon(52.205, 0.119)
        q = LatLon(52.205, 0.119)
        self.test('equals', p.equals(q), 'True')

        p = LatLon(51.4778, -0.0016)
        precision(F_DMS, 0)
        self.test('toStr', p.toStr(), '''51°28'40"N, 000°00'06"W''')
        self.test('toStr', p.toStr(F_D), '51.4778°N, 000.0016°W')
        p = LatLon(51.4778, -0.0016, 42)
        self.test('precision', precision(F_DMS), '0')
        self.test('toStr', p.toStr(), '''51°28'40"N, 000°00'06"W, +42.00m''')


if __name__ == '__main__':

    from geodesy import bases  # private

    t = Tests(__file__, __version__, bases)
    t.testBases(bases.LatLonHeightBase)
    t.results()
    t.exit()

    # Typical test results (on MacOS 10.12.3):

    # testing geodesy.bases version 17.03.07
    # test 1 lat, lon: 50.06632°N, 005.71475°W
    # test 2 lat, lon: 50.066389°N, 005.714722°W
    # test 3 equals: True
    # test 4 toStr: 51°28′40″N, 000°00′06″W
    # test 5 toStr: 51.4778°N, 000.0016°W
    # test 6 precision: 0
    # test 7 toStr: 51°28′40″N, 000°00′06″W, +42.00m
    # all geodesy.bases tests passed (Python 2.7.13 64bit)

    # testing geodesy.bases version 17.03.07
    # test 1 lat, lon: 50.06632°N, 005.71475°W
    # test 2 lat, lon: 50.066389°N, 005.714722°W
    # test 3 equals: True
    # test 4 toStr: 51°28′40″N, 000°00′06″W
    # test 5 toStr: 51.4778°N, 000.0016°W
    # test 6 precision: 0
    # test 7 toStr: 51°28′40″N, 000°00′06″W, +42.00m
    # all geodesy.bases tests passed (Python 3.6.0 64bit)
