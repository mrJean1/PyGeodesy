
# -*- coding: utf-8 -*-

u'''Test projection L{ExactTransverseMercator}.
'''

__all__ = ('Tests',)
__version__ = '20.04.04'

from base import isiOS, isNix, isWindows, TestsBase

from pygeodesy import etm, ExactTransverseMercator


class Tests(TestsBase):

    def testEtm(self, LatLon):
        self.subtitle(etm, LatLon.__name__)

        e = etm.toEtm8(ellipsoidalVincenty.LatLon(-2, 88))  # coverage
        t = '45 S -20297797 5336899'
        self.test('toEtm8', e, t)
        u = e.toUtm()
        self.test('toUtm', u, t)
        self.test('toETM5', etm.parseETM5(t), e)

        self.testCopy(e)

    def testExactTM(self, extendp):
        self.subtitle(etm, ExactTransverseMercator.__name__)

        xtm = ExactTransverseMercator(extendp=extendp)
        t = xtm.toStr()
        self.test('toStr', t, t)  # for coverage

        e, n, g, k = xtm.forward(40.4, -3.7, lon0=-3)  # Madrid , UTM zone 30?
        self.test('easting',  e,  '-59401.921148', fmt='%.6f')
        self.test('northing', n, '4472390.031129', fmt='%.6f')
        self.test('gamma',    g,      '-0.453697', fmt='%.6f')
        self.test('scale',    k,       '0.999643', fmt='%.6f')

        lat, lon, g, k = xtm.reverse(e, n, lon0=-3)
        self.test('lat', lat, '40.400000', fmt='%.6f')
        self.test('lon', lon, '-3.700000', fmt='%.6f')
        self.test('gamma', g, '-0.453697', fmt='%.6f')
        self.test('scale', k,  '0.999643', fmt='%.6f')

        e, n, g, k = xtm.forward(lat, lon, lon0=-3)
        self.test('easting',  e,  '-59401.921148', fmt='%.6f')
        self.test('northing', n, '4472390.031129', fmt='%.6f')
        self.test('gamma',    g,      '-0.453697', fmt='%.6f')
        self.test('scale',    k,       '0.999643', fmt='%.6f')

        e, n, g, k = xtm.forward(40.3, -74.7, lon0=-75)
        self.test('easting',  e,   '25495.511523', fmt='%.6f')
        self.test('northing', n, '4461098.320889', fmt='%.6f')
        self.test('gamma',    g,       '0.194038', fmt='%.6f')
        self.test('scale',    k,       '0.999608', fmt='%.6f')

        lat, lon, g, k = xtm.reverse(e, n, lon0=-75)
        self.test('lat', lat,  '40.300000', fmt='%.6f')
        self.test('lon', lon, '-74.700000', fmt='%.6f')
        self.test('gamma', g,   '0.194038', fmt='%.6f')
        self.test('scale', k,   '0.999608', fmt='%.6f')

        e, n, g, k = xtm.forward(lat, lon, lon0=-75)
        self.test('easting',  e,   '25495.511523', fmt='%.6f')
        self.test('northing', n, '4461098.320889', fmt='%.6f')
        self.test('gamma',    g,       '0.194038', fmt='%.6f')
        self.test('scale',    k,       '0.999608', fmt='%.6f')

        # <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1TransverseMercatorExact.html>
        # lat, lon, _, _ = xtm.reverse(22051449.037349, -7131237.022729)
        lat, lon, g, k = xtm.reverse(29735142.378357, 4235043.607933)
        self.test('lat', lat, '-2.00000000', fmt='%.8f')
        self.test('lon', lon, '88.00000000', fmt='%.8f')
        self.test('gamma', g, '67.63332900', fmt='%.8f')
        self.test('scale', k, '26.33699547', fmt='%.8f')

        if extendp:
            e, n, g, k = xtm.forward(lat, lon)
            self.test('easting',  e, '29735142.37835701', fmt='%.8f', known=isNix or isWindows)
            self.test('northing', n,  '4235043.60793304', fmt='%.8f', known=isNix or isWindows or isiOS)
            self.test('gamma',    g,       '67.63332900', fmt='%.8f')
            self.test('scale',    k,       '26.33699547', fmt='%.8f')

        else:
            e, n, g, k = xtm.forward(-90, -120)  # coverage
            self.test('easting',  e,       '-0.000', fmt='%.3f', known=True)
            self.test('northing', n, '-9997964.943', fmt='%.3f', known=True)
            self.test('gamma',    g,      '120.000', fmt='%.3f')
            self.test('scale',    k,        '1.000', fmt='%.3f')

            lat, lon, g, k = xtm.reverse(e, n)
            self.test('lat', lat,  '-90.000', fmt='%.3f')
            self.test('lon', lon,    '0.000', fmt='%.3f', known=True)
            self.test('gamma', g,    '0.000', fmt='%.3f')
            self.test('scale', k,    '1.000', fmt='%.3f')

        self.testCopy(xtm)


if __name__ == '__main__':

    from pygeodesy import ellipsoidalNvector, ellipsoidalVincenty

    t = Tests(__file__, __version__, etm)
    t.testExactTM(True)
    t.testExactTM(False)
    t.testEtm(ellipsoidalNvector.LatLon)
    t.testEtm(ellipsoidalVincenty.LatLon)
    t.results()
    t.exit()
