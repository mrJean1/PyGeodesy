
# -*- coding: utf-8 -*-

u'''Test Karney wrappers.
'''

__all__ = ('Tests',)
__version__ = '20.06.05'

from base import geographiclib, TestsBase

from pygeodesy import karney, NEG0, wrap180

# some tests from <https://PyPI.org/project/geographiclib>
_testCases = ((35.60777, -139.44815, 111.098748429560326,
              -11.17491,  -69.95921, 129.289270889708762,
               8935244.5604818305, 80.50729714281974,
               6273170.2055303837,  0.16606318447386067,
               0.16479116945612937, 12841384694976.432),
              (55.52454, 106.05087,  22.020059880982801,
               77.03196, 197.18234, 109.112041110671519,
               4105086.1713924406, 36.892740690445894,
               3828869.3344387607,  0.80076349608092607,
               0.80101006984201008, 61674961290615.615),
             (-21.97856, 142.59065, -32.44456876433189,
               41.84138,  98.56635, -41.84359951440466,
               8394328.894657671, 75.62930491011522,
               6161154.5773110616, 0.24816339233950381,
               0.24930251203627892, -6637997720646.717),
             (-66.99028, 112.2363, 173.73491240878403,
              -12.70631, 285.90344,  2.512956620913668,
               11150344.2312080241, 100.278634181155759,
               6289939.5670446687, -0.17199490274700385,
              -0.17722569526345708, -121287239862139.744),
             (-17.42761, 173.34268, -159.033557661192928,
              -15.84784,   5.93557,  -20.787484651536988,
               16076603.1631180673, 144.640108810286253,
               3732902.1583877189, -0.81273638700070476,
              -0.81299800519154474, 97825992354058.708))


class Tests(TestsBase):

    def testDelta(self, n, x, v, e):
        self.test('dict.' + n, v, x, known=abs(v - x) <= e)

    def testDirect(self, g):
        self.subtitle(karney, 'Direct')

        # from <https://PyPI.org/project/geographiclib>
        m = g.ALL | g.LONG_UNROLL
        for (lat1, lon1, azi1, lat2, lon2, azi2,
             s12, a12, m12, M12, M21, S12) in _testCases:
            d = g.Direct(lat1, lon1, azi1, s12, m)
            self.testDelta('lat2', lat2, d.lat2, 1e-13)
            self.testDelta('lon2', lon2, d.lon2, 1e-13)
            self.testDelta('azi2', azi2, d.azi2, 1e-13)
            self.testDelta('a12',  a12,  d.a12,  1e-13)
            self.testDelta('m12',  m12,  d.m12,  1e-8)
            self.testDelta('M12',  M12,  d.M12,  1e-15)
            self.testDelta('M21',  M21,  d.M21,  1e-15)
            self.testDelta('S12',  S12,  d.S12,  0.1)

    def testInverse(self, g):
        self.subtitle(karney, 'Inverse')

        # from <https://PyPI.org/project/geographiclib>
        m = g.ALL | g.LONG_UNROLL
        for (lat1, lon1, azi1, lat2, lon2, azi2,
             s12, a12, m12, M12, M21, S12) in _testCases:
            d = g.Inverse(lat1, lon1, lat2, lon2, m)
            self.testDelta('lat2', lat2, d.lat2, 1e-13)
            self.testDelta('lon2', lon2, d.lon2, 1e-13)
            self.testDelta('azi1', azi1, d.azi1, 1e-13)
            self.testDelta('azi2', azi2, d.azi2, 1e-13)
            self.testDelta('s12',  s12,  d.s12,  1e-8)
            self.testDelta('a12',  a12,  d.a12,  1e-13)
            self.testDelta('m12',  m12,  d.m12,  1e-8)
            self.testDelta('M12',  M12,  d.M12,  1e-15)
            self.testDelta('M21',  M21,  d.M21,  1e-15)
            self.testDelta('S12',  S12,  d.S12,  0.1)

    def testMask(self):
        self.subtitle(karney, 'Mask')

        g = karney._wrapped.Geodesic(1, 0)  # .WGS84
        for n in ('EMPTY', 'LATITUDE', 'LONGITUDE', 'AZIMUTH',
                  'DISTANCE', 'STANDARD', 'DISTANCE_IN',
                  'REDUCEDLENGTH', 'GEODESICSCALE', 'AREA',
                  'ALL', 'LONG_UNROLL'):
            m = getattr(g, n)
            self.test('Geodesic.' + n, m, m)

    def testMath(self):
        self.subtitle(karney, 'Math')

        # compare geomath.Math.AngDiff with mimicked _diff182
        _diff = karney._diff182
        n = '%s(%%d, %%d)' % (_diff.__name__,)
        for a in range(-180, 181, 90):
            for b in range(-180, 181, 90):
                d, _ = _diff(a, b)
                x, _ = _diff(b, a)
                self.test(n % (a, b), d, -x, known=d in (0, 180))

        # compare geomath.Math.AngNormalize with mimicked _norm180
        _norm180 = karney._norm180
        n = '%s(%%s)' % (_norm180.__name__,)
        for a, x in zip((-361, -360, -180,   -90,   -0,   0,   90,   180,   360, 361),
                        (-1.0, -0.0,  180.0, -90.0,  0.0, 0.0, 90.0, 180.0, 0.0, 1.0)):
            w = _norm180(a)
            self.test(n % (a), float(w), float(x), known=w in (0, -180))
            w = wrap180(a)
            self.test(' wrap180(%s)' % (a), float(w), float(x), known=w in (0, -180))

        # compare geomath.Math.sum with mimicked _sum2
        _fsum2_ = karney._fsum2_
        n = _fsum2_.__name__
        s, t = _fsum2_(7, 1e100, -7, -1e100, 9e-20, -8e-20)
        self.test(n, s, '1.0e-20', fmt='%.1e')
        self.test(n, t, 0.0, known=t in (NEG0, 0.0))
        # however, summation fails after some shuffling ...
        s, _ = _fsum2_(7, 1e100, 9e-20, -7, -1e100, -8e-20)
        self.test(n, s, '1.0e-20', fmt='%.1e', known=True)
        # ... but works after some other shuffling
        s, _ = _fsum2_(7, 1e100, -7, 9e-20, -1e100, -8e-20)
        self.test(n, s, '1.0e-20', fmt='%.1e', known=True)

        # Knuth/Kulisch, TAOCP, vol 2, p 245, sec 4.2.2, item 31, see also .testFmath.py
        # <https://SeriousComputerist.Atariverse.com/media/pdf/book/
        #          Art%20of%20Computer%20Programming%20-%20Volume%202%20(Seminumerical%20Algorithms).pdf>
        x = 408855776.0
        y = 708158977.0
        s, t = _fsum2_(2*y**2, 9*x**4, -(y**4))
        self.test(n, s, 1.0, known=True)  # -3.589050987400773e+19
        self.test(n, t, 0.0, known=t in (NEG0, 0.0))


if __name__ == '__main__':

    t = Tests(__file__, __version__, karney)
    if geographiclib:
        g = karney._wrapped.Geodesic_WGS84
        t.test('Geodesic', isinstance(g, karney._wrapped.Geodesic), True)
        t.testDirect(g)
        t.testInverse(g)
        t.testMask()
    else:
        t.skip('no geographiclib', n=103)
    t.testMath()
    t.results()
    t.exit()
