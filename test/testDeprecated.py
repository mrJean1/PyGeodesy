
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '21.02.11'

from base import TestsBase

from pygeodesy import R_MA, map2, \
                      HeightIDW, HeightIDW2, HeightIDW3, \
                      HeightIDWequirectangular, HeightIDWeuclidean, \
                      HeightIDWhaversine, \
                      anStr, areaof, bounds, clipStr, decodeEPSG2, encodeEPSG, \
                      equirectangular3, fStr, hypot3, isenclosedby, \
                      nearestOn3, nearestOn4, \
                      parseUTM, perimeterof, polygon,\
                      simplify2, toUtm, utmZoneBand2
from testRoutes import RdpFFI

from math import sqrt


class Tests(TestsBase):

    def testDeprecated(self, LatLon):

        c = HeightIDW  # == HeightIDWeuclidean in Python 3.7+
        self.test(c.__name__, issubclass(c, HeightIDWeuclidean), True)
        c = HeightIDW2  # == HeightIDWequirectangular in Python 3.7+
        self.test(c.__name__, issubclass(c, HeightIDWequirectangular), True)
        c = HeightIDW3  # == HeightIDWhaversine in Python 3.7+
        self.test(c.__name__, issubclass(c, HeightIDWhaversine), True)

        p = LatLon(0, 0), LatLon(1, 0), LatLon(0, 1)
        self.test('areaof', areaof(p, radius=R_MA), '7.086883e+09', fmt='%.6e')

        p = LatLon(85, 90), LatLon(-85, 0), LatLon(85, -90), LatLon(85, -180)
        b = map2(float, bounds(p))
        self.test('bounds', b, '(-85.0, -180.0, 85.0, 90.0)')

        self.test('anStr', anStr('a-b?_'), 'a-b__')

        self.test('clipStr', clipStr('test/testBasics.py', limit=12), 'test/t....ics.py')

        self.test('decodeEPSG2', decodeEPSG2(32712), "(12, 'S')")
        self.test('encodeEPSG', encodeEPSG(12, hemipole='S'), '32712')

        t = equirectangular3(0, 2, 3, 4)
        self.test('equirectangular3', len(t), 3)
        self.test('equirectangular3', t[0], 12.997, fmt='%.3f')

        self.test('fStr', fStr(0.123, prec=-6), '0.123000')
        self.test('fStr', fStr(0.123, prec=+6), '0.123')
        self.test('fStr', fStr((0.123, 456.789), prec=+6), '0.123, 456.789')
        self.test('fStr', fStr(0.123, prec=-5, fmt='%.*e'), '1.23000e-01')
        self.test('fStr', fStr(0.123, prec=+5, fmt='%.*e'), '1.23e-01')
        self.test('fStr', fStr(0.123, prec=+6, fmt='%.*f'), '0.123')

        h = hypot3(3000, 200, 10)
        s = sqrt(3000**2 + 200**2 + 10**2)
        self.test('hypot3', h, s, fmt='%.6f')

        b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
        self.test('isenclosedby', isenclosedby(LatLon(45.5, 1.5), b), True)

        p = LatLon(45, 2)
        b = LatLon(45, 1), LatLon(47, 3)
        t = nearestOn3(p, b, adjust=False)
        self.test('nearestOn3', len(t), 3)
        self.test('nearestOn3', t[:2], (45.5, 1.5))
        t = nearestOn4(p, b, adjust=False)
        self.test('nearestOn4', len(t), 4)
        self.test('nearestOn4', t[:2], (45.5, 1.5))

        t = parseUTM('18 N 516620 4574500', Utm=None)  # Milford, PA
        self.test('parseUTM', t, "(18, 'N', 516620.0, 4574500.0)")

        p = LatLon(0, 0), LatLon(1, 0), LatLon(0, 1)
        self.test('perimeterof', perimeterof(p, radius=R_MA), '2.687460e+05', fmt='%.6e')

        p = LatLon(0, 0), LatLon(1, 0), LatLon(0, 1)
        self.test('polygon', polygon(p)[0], 3)

        t = simplify2(RdpFFI, 16, adjust=True, shortest=False)
        self.test('simplify2', len(t), 4)

        t = toUtm('50°52′10″N', '115°39′03″W', Utm=None, name='Mt Assiniboine')
        self.test('toUtm', len(t), 6)

        t = utmZoneBand2('50°52′10″N', '115°39′03″W')
        self.test('utmZoneBand2', t, "(11, 'U')")


if __name__ == '__main__':

    from pygeodesy import deprecated, LatLon_, \
                          ellipsoidalVincenty, sphericalTrigonometry

    t = Tests(__file__, __version__, deprecated)
    try:
        t.testDeprecated(LatLon_)
        t.testDeprecated(ellipsoidalVincenty.LatLon)
        t.testDeprecated(sphericalTrigonometry.LatLon)
    except DeprecationWarning:
        t.skip(DeprecationWarning.__name__, n=87 - t.total)
    t.results()
    t.exit()
