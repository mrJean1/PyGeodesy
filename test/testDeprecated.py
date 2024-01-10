
# -*- coding: utf-8 -*-

# Test L{deprecated} classes, functions and methods.

__all__ = ('Tests',)
__version__ = '24.01.02'

from bases import TestsBase

from pygeodesy import R_MA, deprecated, isDEPRECATED, map2, unstr, \
                      HeightIDWequirectangular, HeightIDWeuclidean, HeightIDWhaversine
try:
    from pygeodesy.deprecated import __star__
    __star__(__name__)
except ImportError:
    from pygeodesy.deprecated import *  # PYCHOK expected
from testRoutes import RdpFFI

from math import sqrt


class Tests(TestsBase):

    def testDeprecated(self, LatLon):

        ks = tuple(LatLon(k, k, height=k) for k in range(17))
        c = HeightIDW(ks)  # PYCHOK == HeightIDWeuclidean in Python 3.7+
        self.test(c.__class__.__name__, isinstance(c, HeightIDWeuclidean), True, nl=1)
        c = HeightIDW2(ks)  # PYCHOK == HeightIDWequirectangular in Python 3.7+
        self.test(c.__class__.__name__, isinstance(c, HeightIDWequirectangular), True)
        c = HeightIDW3(ks)  # PYCHOK == HeightIDWhaversine in Python 3.7+
        self.test(c.__class__.__name__, isinstance(c, HeightIDWhaversine), True, nt=1)
        del ks

        p = LatLon(0, 0), LatLon(1, 0), LatLon(0, 1)
        self.test('areaof', areaof(p, radius=R_MA), '7.086883e+09', fmt='%.6e')  # PYCHOK DEPRECATED

        p = LatLon(85, 90), LatLon(-85, 0), LatLon(85, -90), LatLon(85, -180)
        b = map2(float, bounds(p))  # PYCHOK DEPRECATED
        self.test('bounds', b, '(-85.0, -180.0, 85.0, 90.0)')

        self.test('anStr', anStr('a-b?_'), 'a-b__')  # PYCHOK DEPRECATED

        self.test('clipStr', clipStr('test/testBasics.py', limit=12), 'test/t....ics.py')  # PYCHOK DEPRECATED

        self.test('decodeEPSG2', decodeEPSG2(32712), "(12, 'S')")  # PYCHOK DEPRECATED
        self.test('encodeEPSG', encodeEPSG(12, hemipole='S'), '32712')  # PYCHOK DEPRECATED

        t = equirectangular3(0, 2, 3, 4)  # PYCHOK DEPRECATED
        self.test('equirectangular3', len(t), 3)
        self.test('equirectangular3', t[0], 12.997, fmt='%.3f')

        self.test('fStr', fStr(0.123, prec=-6), '0.123000')  # PYCHOK DEPRECATED
        self.test('fStr', fStr(0.123, prec=+6), '0.123')  # PYCHOK DEPRECATED
        self.test('fStr', fStr((0.123, 456.789), prec=+6), '0.123, 456.789')  # PYCHOK DEPRECATED
        self.test('fStr', fStr(0.123, prec=-5, fmt='%.*e'), '1.23000e-01')  # PYCHOK DEPRECATED
        self.test('fStr', fStr(0.123, prec=+5, fmt='%.*e'), '1.23e-01')  # PYCHOK DEPRECATED
        self.test('fStr', fStr(0.123, prec=+6, fmt='%.*f'), '0.123')  # PYCHOK DEPRECATED

        h = hypot3(3000, 200, 10)  # PYCHOK DEPRECATED
        s = sqrt(3000**2 + 200**2 + 10**2)
        self.test('hypot3', h, s, fmt='%.6f')

        b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
        self.test('isenclosedby', isenclosedby(LatLon(45.5, 1.5), b), True)  # PYCHOK DEPRECATED

        p = LatLon(45, 2)
        b = LatLon(45, 1), LatLon(47, 3)
        t = nearestOn3(p, b, adjust=False)  # PYCHOK DEPRECATED
        self.test('nearestOn3', len(t), 3)
        self.test('nearestOn3', t[:2], (45.5, 1.5))
        t = nearestOn4(p, b, adjust=False)  # PYCHOK DEPRECATED
        self.test('nearestOn4', len(t), 4)
        self.test('nearestOn4', t[:2], (45.5, 1.5))

        t = parseUTM('18 N 516620 4574500', Utm=None)  # PYCHOK Milford, PA
        self.test('parseUTM', t, "(18, 'N', 516620.0, 4574500.0)")

        p = LatLon(0, 0), LatLon(1, 0), LatLon(0, 1)
        self.test('perimeterof', perimeterof(p, radius=R_MA), '2.687460e+05', fmt='%.6e')  # PYCHOK DEPRECATED

        p = LatLon(0, 0), LatLon(1, 0), LatLon(0, 1)
        self.test('polygon', polygon(p)[0], 3)  # PYCHOK DEPRECATED

        t = simplify2(RdpFFI, 16, adjust=True, shortest=False)  # PYCHOK DEPRECATED
        self.test('simplify2', len(t), 4)

        t = toUtm('50°52′10″N', '115°39′03″W', Utm=None, name='Mt Assiniboine')  # PYCHOK DEPRECATED
        self.test('toUtm', len(t), 6)

        t = utmZoneBand2('50°52′10″N', '115°39′03″W')  # PYCHOK DEPRECATED
        self.test('utmZoneBand2', t, "(11, 'U')")

    def testDEPRECATED(self, *known):
        # from pygeodesy.interns import _DOT_, _UNDER_  # from .lazily
        from pygeodesy.lazily import _ALL_MODS as _MODS, _ALL_DEPRECATED,  \
                                     _attrof, _DOT_, _headof, _UNDER_

        for m, t in _ALL_DEPRECATED.items():
            if _headof(m) == 'deprecated':
                self.test(m, len(t), len(t), nl=1)
                m = _MODS.getmodule(m.replace(_UNDER_, _DOT_))
                for a in t:
                    a = _attrof(a)  # a_ or a
                    n =  unstr(isDEPRECATED, a)
                    try:
                        v = getattr(m, a)
                        self.test(n, isDEPRECATED(v), True, known=v in known)
                    except Exception as x:
                        self.test(n, str(x), True)


if __name__ == '__main__':

    from pygeodesy import LatLon_, ellipsoidalVincenty, sphericalTrigonometry

    t = Tests(__file__, __version__, deprecated)
    try:
        t.testDeprecated(LatLon_)
        t.testDeprecated(ellipsoidalVincenty.LatLon)
        t.testDeprecated(sphericalTrigonometry.LatLon)
        t.testDEPRECATED(isDEPRECATED)
    except DeprecationWarning:
        t.skip(DeprecationWarning.__name__, n=87 - t.total)
    t.results()
    t.exit()
