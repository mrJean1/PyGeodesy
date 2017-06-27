
# -*- coding: utf-8 -*-

# Test spherical earth model functions and methods.

__all__ = ('Tests',)
__version__ = '17.06.25'

from testLatLon import Tests as _TestsLL
from testVectorial import Tests as _TestsV

from pygeodesy import F_D, F_DMS, lonDMS


class Tests(_TestsLL, _TestsV):

    def testSpherical(self, module):

        self.subtitle(module, 'Spherical')

        LatLon = module.LatLon

        p = LatLon(51.8853, 0.2545)
        self.test('isSpherical', p.isSpherical, True)
        self.test('isEllipsoidal', p.isEllipsoidal, False)

        q = LatLon(49.0034, 2.5735)
        self.test('isSpherical', q.isSpherical, True)
        self.test('isEllipsoidal', q.isEllipsoidal, False)

        i = p.intersection(108.55, q, 32.44)
        self.test('intersection', i.toStr(F_D),  '50.907608°N, 004.508575°E')  # 50.9076°N, 004.5086°E  # Trig
        self.test('intersection', i.toStr(F_DMS), '50°54′27.39″N, 004°30′30.87″E')
        self.test('intersection', isinstance(i, LatLon), True)

        REO = LatLon(42.600, -117.866)
        BKE = LatLon(44.840, -117.806)
        i = REO.intersection(51, BKE, 137)
        self.test('intersection', i.toStr(F_D), '43.5719°N, 116.188757°W')  # 43.572°N, 116.189°W
        self.test('intersection', i.toStr(F_DMS), '43°34′18.84″N, 116°11′19.53″W')
        self.test('intersection', isinstance(i, LatLon), True)

        p = LatLon(0, 0)
        self.test('maxLat0',  p.maxLat( 0), '90.0')
        self.test('maxLat1',  p.maxLat( 1), '89.0')
        self.test('maxLat90', p.maxLat(90),  '0.0')

        if hasattr(LatLon, 'crossingParallels'):
            ps = p.crossingParallels(LatLon(60, 30), 30)
            t = ', '.join(map(lonDMS, ps))
            self.test('crossingParallels', t, '009°35′38.65″E, 170°24′21.35″E')

        p = LatLon(51.127, 1.338)
        q = LatLon(50.964, 1.853)
        b = p.rhumbBearingTo(q)
        self.test('rhumbBearingTo', b, 116.722, fmt='%.3f')  # 116.7

        d = p.rhumbDestination(40300, 116.7)
        self.test('rhumbDestination', d, '50.964155°N, 001.853°E')  # 50.9642°N, 001.8530°E
        self.test('rhumbDestination', isinstance(d, LatLon), True)

        d = p.rhumbDistanceTo(q)
        self.test('rhumbDistanceTo', d, 40307.8, fmt='%.1f')  # XXX 40310 ?

        m = p.rhumbMidpointTo(q)
        self.test('rhumbMidpointo', m, '51.0455°N, 001.595727°E')
        self.test('rhumbMidpointo', isinstance(m, LatLon), True)

        b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
        self.test('areaOf', module.areaOf(b), '8.6660587507e+09', fmt='%.10e')  # 8666058750.718977

        c = LatLon(0, 0), LatLon(1, 0), LatLon(0, 1)
        self.test('areaOf', module.areaOf(c), '6.18e+09', fmt='%.2e')

        if hasattr(module, 'isPoleEnclosedBy'):
            p = LatLon(85, 90), LatLon(85, 0), LatLon(85, -90), LatLon(85, -180)
            self.test('isPoleEnclosedBy', module.isPoleEnclosedBy(p), 'True')
            p = LatLon(85, 90), LatLon(85, 0), LatLon(85, -180)
            self.test('isPoleEnclosedBy', module.isPoleEnclosedBy(p), 'True', known=True)


if __name__ == '__main__':

    from pygeodesy import sphericalNvector as N, \
                          sphericalTrigonometry as T

    t = Tests(__file__, __version__)

    t.testLatLon(N, Sph=True)
    t.testSpherical(N)
    t.testVectorial(N)

    t.testLatLon(T, Sph=True)
    t.testSpherical(T)

    t.results()
    t.exit()
