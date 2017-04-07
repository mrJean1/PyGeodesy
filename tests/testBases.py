
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '17.04.07'

from tests import Tests as _Tests

from pygeodesy import F_D, F_DMS, precision


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

    from pygeodesy import bases  # private

    t = Tests(__file__, __version__, bases)
    t.testBases(bases.LatLonHeightBase)
    t.results()
    t.exit()
