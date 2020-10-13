
# -*- coding: utf-8 -*-

u'''Test projection L{ExactTransverseMercator}.
'''

__all__ = ('Tests',)
__version__ = '20.10.13'

from base import isiOS, isNix, isWindows, TestsBase

from pygeodesy import etm, ExactTransverseMercator


class Tests(TestsBase):

    def testEtm(self, LatLon):
        self.subtitle(etm, LatLon.__name__)

        e = etm.toEtm8(ellipsoidalVincenty.LatLon(-2, 88), name='test')  # coverage
        t = '45 S -20297797 5336899'
        self.test('toEtm8', e, t)
        self.test('name', e.name, 'test')
        u = e.toUtm()
        self.test('toUtm', u, t)
        self.test('name', u.name, 'test')
        self.test('toETM5', etm.parseETM5(t), e)

        self.testCopy(e, 'name')

        e = e.parse('31 N 448251, 5411932', name='parse')  # coverage
        t = '31 N 448251 5411932'
        self.test('parse', e, t)
        self.test('name', e.name, 'parse')
        e = e.parseETM(t.replace(' ', ','))
        self.test('parseETM', e, t)
        self.test('name', e.name, 'parse')

    def testExactTM(self, extendp):
        self.subtitle(etm, ExactTransverseMercator.__name__)

        xtm = ExactTransverseMercator(extendp=extendp, name='test')
        self.test('name', xtm.name, 'test')  # coverage
        t = xtm.toStr()
        self.test('toStr', t, t)  # coverage

        e, n, g, k = xtm.forward(40.4, -3.7, lon0=-3)  # Madrid , UTM zone 30?
        self.test('easting',  e,  '-59401.921148', prec=6)
        self.test('northing', n, '4472390.031129', prec=6)
        self.test('gamma',    g,      '-0.453697', prec=6)
        self.test('scale',    k,       '0.999643', prec=6)

        lat, lon, g, k = xtm.reverse(e, n, lon0=-3)
        self.test('lat', lat, '40.400000', prec=6)
        self.test('lon', lon, '-3.700000', prec=6)
        self.test('gamma', g, '-0.453697', prec=6)
        self.test('scale', k,  '0.999643', prec=6)

        e, n, g, k = xtm.forward(lat, lon, lon0=-3)
        self.test('easting',  e,  '-59401.921148', prec=6)
        self.test('northing', n, '4472390.031129', prec=6)
        self.test('gamma',    g,      '-0.453697', prec=6)
        self.test('scale',    k,       '0.999643', prec=6)

        e, n, g, k = xtm.forward(40.3, -74.7, lon0=-75)
        self.test('easting',  e,   '25495.511523', prec=6)
        self.test('northing', n, '4461098.320889', prec=6)
        self.test('gamma',    g,       '0.194038', prec=6)
        self.test('scale',    k,       '0.999608', prec=6)

        lat, lon, g, k = xtm.reverse(e, n, lon0=-75)
        self.test('lat', lat,  '40.300000', prec=6)
        self.test('lon', lon, '-74.700000', prec=6)
        self.test('gamma', g,   '0.194038', prec=6)
        self.test('scale', k,   '0.999608', prec=6)

        e, n, g, k = xtm.forward(lat, lon, lon0=-75)
        self.test('easting',  e,   '25495.511523', prec=6)
        self.test('northing', n, '4461098.320889', prec=6)
        self.test('gamma',    g,       '0.194038', prec=6)
        self.test('scale',    k,       '0.999608', prec=6)

        # <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1TransverseMercatorExact.html>
        # lat, lon, _, _ = xtm.reverse(22051449.037349, -7131237.022729)
        lat, lon, g, k = xtm.reverse(29735142.378357, 4235043.607933)
        self.test('lat', lat, '-2.00000000', prec=8)
        self.test('lon', lon, '88.00000000', prec=8)
        self.test('gamma', g, '67.63332900', prec=8)
        self.test('scale', k, '26.33699547', prec=8)

        if extendp:
            e, n, g, k = xtm.forward(lat, lon)
            self.test('easting',  e, '29735142.37835703', prec=8, known=isNix or isWindows or isiOS)
            self.test('northing', n,  '4235043.60793304', prec=8, known=isNix or isWindows or isiOS)
            self.test('gamma',    g,       '67.63332900', prec=8)
            self.test('scale',    k,       '26.33699547', prec=8)

        else:
            e, n, g, k = xtm.forward(-90, -120)  # coverage
            self.test('easting',  e,       '-0.000', prec=3, known=True)
            self.test('northing', n, '-9997964.943', prec=3, known=True)
            self.test('gamma',    g,      '120.000', prec=3)
            self.test('scale',    k,        '1.000', prec=3)

            lat, lon, g, k = xtm.reverse(e, n)
            self.test('lat', lat,  '-90.000', prec=3)
            self.test('lon', lon,    '0.000', prec=3, known=True)
            self.test('gamma', g,    '0.000', prec=3, known=abs(g) == 0)
            self.test('scale', k,    '1.000', prec=3)

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
