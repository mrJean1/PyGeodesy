
# -*- coding: utf-8 -*-

# Test spherical earth model functions and methods.

__all__ = ('Tests',)
__version__ = '18.01.12'

from testLatLon import Tests as _TestsLL
from testVectorial import Tests as _TestsV

from pygeodesy import F_D, F_DMS, lonDMS


class Tests(_TestsLL, _TestsV):

    def testSpherical(self, module, Vct=False):  # MCCABE 13

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

        if hasattr(LatLon, 'isEnclosedBy'):
            p = LatLon(45.1, 1.1)

            b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
            for _ in self.testiter():
                self.test('isEnclosedBy', p.isEnclosedBy(b), True)

            b = LatLon(45, 1), LatLon(45, 3), LatLon(46, 2), LatLon(47, 3), LatLon(47, 1)
            for _ in self.testiter():
                try:
                    self.test('isEnclosedBy', p.isEnclosedBy(b), True)  # Nvector
                except ValueError as x:
                    t = ' '.join(str(x).split()[:3] + ['...)'])
                    self.test('isEnclosedBy', t, 'non-convex: (LatLon(45°00′00.0″N, 001°00′00.0″E), ...)')  # Trig

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

        if hasattr(module, 'nearestOn2'):
            c, d = module.nearestOn2(p, b)
            self.test('nearestOn2', c, '46.000996°N, 001.353049°E' if Vct else '46.0°N, 001.369324°E')
            self.test('nearestOn2', d, '569987.49' if Vct else '570101.83', fmt='%.2f')
            d = p.distanceTo(c)
            self.test('distanceTo', d, '569987.49' if Vct else '570101.82', fmt='%.2f')

            p = LatLon(47, 3)
            c, d = module.nearestOn2(p, b)
            self.test('nearestOn2', c, '46.0°N, 002.0°E' if Vct else '46.0°N, 002.0°E')
            self.test('nearestOn2', d, '134989.80' if Vct else '134992.48', fmt='%.2f')
            d = p.distanceTo(c)
            self.test('distanceTo', d, '134989.80' if Vct else '134989.80', fmt='%.2f')

            p = LatLon(45, 2)
            b = LatLon(45, 1), LatLon(47, 3)
            if Vct:
                c, d = module.nearestOn2(p, b)
                self.test('nearestOn2', c, '45.330691°N, 001.318551°E')
                self.test('distanceTo2', d, '64856.28', fmt='%.2f')
            else:
                c, d = module.nearestOn2(p, b, adjust=False)
                self.test('nearestOn2', c, '45.5°N, 001.5°E')
                self.test('distanceTo2', d, '78626.79', fmt='%.2f')
                c, d = module.nearestOn2(p, b, adjust=True)
                self.test('nearestOn2', c, '45.331319°N, 001.331319°E')
                self.test('distanceTo2', d, '64074.48', fmt='%.2f')
            d = p.distanceTo(c)
            self.test('distanceTo', d, '64856.28' if Vct else '64074.12', fmt='%.2f')
            # TrigTrue vs Nvector closests
            p = LatLon(45.330691, 001.318551)
            d = p.distanceTo(LatLon(45.331319, 001.331319))
            self.test('difference', d, '1000.53', fmt='%.2f')  # PYCHOK false?

            if not Vct:  # check nearestOn2 with closest on the segment
                b = LatLon(0, 1), LatLon(2, 3), LatLon(4, 5), LatLon(6, 7), LatLon(8, 9)
                for i in range(8):
                    p = LatLon(i + 2, i)
                    c, d = module.nearestOn2(p, b, adjust=False)
                    p.lat -= 1.5
                    p.lon += 1.5
                    self.test('nearestOn2', c.toStr(F_D, prec=6), p.toStr(F_D, prec=6))
                    self.test('neartesOn2', d, '235880.385', fmt='%.3f')

        if hasattr(module, 'isPoleEnclosedBy'):
            p = LatLon(85, 90), LatLon(85, 0), LatLon(85, -90), LatLon(85, -180)
            for _ in self.testiter():
                self.test('isPoleEnclosedBy', module.isPoleEnclosedBy(p), 'True')  # PYCHOK false?
            p = LatLon(85, 90), LatLon(85, 0), LatLon(85, -180)
            for _ in self.testiter():
                self.test('isPoleEnclosedBy', module.isPoleEnclosedBy(p), 'True', known=True)  # PYCHOK false?


if __name__ == '__main__':

    from pygeodesy import sphericalNvector as N, \
                          sphericalTrigonometry as T

    t = Tests(__file__, __version__)

    t.testLatLon(N, Sph=True)
    t.testVectorial(N)
    t.testSpherical(N, Vct=True)

    t.testLatLon(T, Sph=True)
    t.testSpherical(T, Vct=False)

    t.results()
    t.exit()
