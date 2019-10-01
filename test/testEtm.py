
# -*- coding: utf-8 -*-

u'''Test projection L{ExactTransverseMercator}.
'''

__all__ = ('Tests',)
__version__ = '19.10.01'

from base import isiOS, isNix, isWindows, TestsBase

from pygeodesy import etm, ExactTransverseMercator, wrap180

_ETM = ExactTransverseMercator()


class Tests(TestsBase):

    def testExactTM(self):

        e, n, g, k = _ETM.forward(40.4, -3.7, lon0=-3)  # Madrid , UTM zone 30?
        self.test('easting',  e,  '-59401.921148', fmt='%.6f')
        self.test('northing', n, '4472390.031129', fmt='%.6f')
        self.test('gamma',    g,      '-0.453697', fmt='%.6f')
        self.test('scale',    k,       '0.999643', fmt='%.6f')

        lat, lon, g, k = _ETM.reverse(e, n, lon0=-3)
        self.test('lat', lat, '40.400000', fmt='%.6f')
        self.test('lon', lon, '-3.700000', fmt='%.6f')
        self.test('gamma', g, '-0.453697', fmt='%.6f')
        self.test('scale', k,  '0.999643', fmt='%.6f')

        e, n, g, k = _ETM.forward(lat, lon, lon0=-3)
        self.test('easting',  e,  '-59401.921148', fmt='%.6f')
        self.test('northing', n, '4472390.031129', fmt='%.6f')
        self.test('gamma',    g,      '-0.453697', fmt='%.6f')
        self.test('scale',    k,       '0.999643', fmt='%.6f')

        e, n, g, k = _ETM.forward(40.3, -74.7, lon0=-75)
        self.test('easting',  e,   '25495.511523', fmt='%.6f')
        self.test('northing', n, '4461098.320889', fmt='%.6f')
        self.test('gamma',    g,       '0.194038', fmt='%.6f')
        self.test('scale',    k,       '0.999608', fmt='%.6f')

        lat, lon, g, k = _ETM.reverse(e, n, lon0=-75)
        self.test('lat', lat,  '40.300000', fmt='%.6f')
        self.test('lon', lon, '-74.700000', fmt='%.6f')
        self.test('gamma', g,   '0.194038', fmt='%.6f')
        self.test('scale', k,   '0.999608', fmt='%.6f')

        e, n, g, k = _ETM.forward(lat, lon, lon0=-75)
        self.test('easting',  e,   '25495.511523', fmt='%.6f')
        self.test('northing', n, '4461098.320889', fmt='%.6f')
        self.test('gamma',    g,       '0.194038', fmt='%.6f')
        self.test('scale',    k,       '0.999608', fmt='%.6f')

        # <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1TransverseMercatorExact.html>
        # lat, lon, _, _ = _ETM.reverse(22051449.037349, -7131237.022729)
        lat, lon, g, k = _ETM.reverse(29735142.378357, 4235043.607933)
        self.test('lat', lat, '-2.00000000', fmt='%.8f')
        self.test('lon', lon, '88.00000000', fmt='%.8f')
        self.test('gamma', g, '67.63332900', fmt='%.8f')
        self.test('scale', k, '26.33699547', fmt='%.8f')

        e, n, g, k = _ETM.forward(lat, lon)
        self.test('easting',  e, '29735142.37835701', fmt='%.8f', known=isNix or isWindows)
        self.test('northing', n,  '4235043.60793304', fmt='%.8f', known=isNix or isWindows or isiOS)
        self.test('gamma',    g,       '67.63332900', fmt='%.8f')
        self.test('scale',    k,       '26.33699547', fmt='%.8f')

        # compare geomath.Math.AngDiff with mimicked etm._diff182
        _diff = etm._diff182
        n = '%s.%s(%%d, %%d)' % (_diff.__module__, _diff.__name__)
        for a in range(-180, 181, 90):
            for b in range(-180, 181, 90):
                d, _ = _diff(a, b)
                x, _ = _diff(b, a)
                self.test(n % (a, b), d, -x, known=d in (0, 180))

        # compare geomath.Math.AngNormalize with mimicked etm._wrap180
        _wrap = etm._wrap180
        n = '%s.%s(%%s)' % (_wrap.__module__, _wrap.__name__)
        for a, x in zip((-361, -360, -180,   -90,   -0,   0,   90,   180,   360, 361),
                        (-1.0, -0.0,  180.0, -90.0,  0.0, 0.0, 90.0, 180.0, 0.0, 1.0)):
            w = _wrap(a)
            self.test(n % (a), float(w), float(x), known=w in (0, -180))
            w = wrap180(a)
            self.test('pygeodesy.wrap180(%s)' % (a), float(w), float(x), known=w in (0, -180))

        # compare geomath.Math.sum with mimicked etm._sum2
        _sum2 = etm._sum2
        n = '%s.%s' % (_sum2.__module__, _sum2.__name__)
        s = t = 0  # see test.Fmath.py
        for x in (7, 1e100, -7, -1e100, 9e-20, -8e-20):
            s, x = _sum2(s, x)
            t, _ = _sum2(t, x)
        self.test(n, s, '1.0e-20', fmt='%.1e')
        self.test(n, t, '0.0e+00', fmt='%.1e')


if __name__ == '__main__':

    t = Tests(__file__, __version__, etm)
    t.testExactTM()
    t.results()
    t.exit()
