
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '18.10.20'

from base import TestsBase

from pygeodesy import F_D, F_DMS, precision


class Tests(TestsBase):

    def testBases(self, LatLon):

        p = LatLon(50.06632, -5.71475)
        self.test('lat, lon', p, '50.06632°N, 005.71475°W')
        q = LatLon('50°03′59″N', """005°42'53"W""")
        self.test('lat, lon', q, '50.066389°N, 005.714722°W')

        p = LatLon(52.205, 0.119)
        q = LatLon(52.205, 0.119)
        self.test('isequalTo', p.isequalTo(q), True)

        t = precision(F_DMS, 0)
        p = LatLon(51.4778, -0.0016)
        self.test('toStr', p.toStr(), '''51°28'40"N, 000°00'06"W''')
        self.test('toStr', p.toStr(F_D), '51.4778°N, 000.0016°W')

        self.test('precision', precision(F_DMS), 0)
        p = LatLon(51.4778, -0.0016, 42)
        self.test('toStr', p.toStr(), '''51°28'40"N, 000°00'06"W, +42.00m''')
        precision(F_DMS, t)  # restore

        q.latlon = p.latlon
        self.test('isequalTo', p.isequalTo(q), True)
        self.test('isequalTo3', p.isequalTo3(q), False)
        q.latlon = p.latlon + (p.height,)
        self.test('isequalTo3', p.isequalTo3(q), True)


if __name__ == '__main__':

    from pygeodesy import ellipsoidalNvector, ellipsoidalVincenty, \
                          sphericalNvector, sphericalTrigonometry

    t = Tests(__file__, __version__)
    t.testBases(ellipsoidalNvector.LatLon)
    t.testBases(ellipsoidalVincenty.LatLon)
    t.testBases(sphericalNvector.LatLon)
    t.testBases(sphericalTrigonometry.LatLon)
    t.results()
    t.exit()
