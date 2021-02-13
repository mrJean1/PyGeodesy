
# -*- coding: utf-8 -*-

# Test datums, ellipsoids and transforms.

__all__ = ('Tests',)
__version__ = '21.02.11'

from base import TestsBase

from pygeodesy import EcefKarney, Ellipsoid, Ellipsoid2, Ellipsoids, \
                      a_b2f_, a_b2f2, a_b2n, a_f2Tuple, b_f2a, b_f_2a, \
                      circle4, ellipsoids, f_2f, fstr, hypot_, n2e2, \
                      n2f, PI_2, R_M, sincos2d


class Tests(TestsBase):

    def testEllipsoids(self):  # MCCABE 14
        # datum module tests
        E = Ellipsoid(1000, 1000, 0, name='TestEllipsoid')
        self.test('ellipsoid', E is Ellipsoids.TestEllipsoid, True)
#       print(Ellipsoid())

        e = Ellipsoids.unregister('TestEllipsoid')
        self.test(e.name, e, E)

        # to show WGS84 and NAD83 are identical up to 3 decimals
        for E in (Ellipsoids.WGS84, Ellipsoids.GRS80):  # NAD83
            self.subtitle(ellipsoids, E.name)
            self.test('R1', E.R1, R_M, fmt='%.4f')
            self.test('R2', E.R2, '6371007.2', fmt='%.1f')
            self.test('R3', E.R3, '6371000.8', fmt='%.1f')
            self.test('A',  E.A,  '6367449.1', fmt='%.1f')
            self.test('L',  E.L, '10001965.7', fmt='%.1f')

            self.test('Rrectifying', E.Rrectifying,     '6367449.1',   fmt='%.1f')  # 6367445.0
            self.test('Rgeometric',  E.Rgeometric,      '6367435.7',   fmt='%.1f')
            self.test('Rgeocentric', E.Rgeocentric(0),  '6378137.000', fmt='%.3f')
            self.test('Rgeocentric', E.Rgeocentric(45), '6367489.544', fmt='%.3f')
            self.test('Rgeocentric', E.Rgeocentric(90), '6356752.314', fmt='%.3f')

            self.test('Rlat', E.Rlat(0),  '6378137.000', fmt='%.3f')
            self.test('Rlat', E.Rlat(45), '6367444.657', fmt='%.3f')
            self.test('Rlat', E.Rlat(90), '6356752.314', fmt='%.3f')

            self.test('circle4.radius', E.circle4(0).radius,  '6378137.000', fmt='%.3f')
            self.test('circle4.radius', E.circle4(45).radius, '4517590.879', fmt='%.3f')
            self.test('circle4.radius', E.circle4(90).radius, '0.000',       fmt='%.3f')

            self.test('distance2', fstr(E.distance2( 0,  0,  1,  1), prec=3),  '156903.472, 45.192')
            self.test('distance2', fstr(E.distance2( 0,  0, 10, 10), prec=3), '1569034.719, 45.192')
            self.test('distance2', fstr(E.distance2(40, 40, 50, 50), prec=3), '1400742.676, 37.563')
            self.test('distance2', fstr(E.distance2(70, 70, 80, 80), prec=3), '1179164.848, 18.896')

            self.test('roc2', fstr(E.roc2(0),  prec=3), '6335439.327, 6378137.0')
            self.test('roc2', fstr(E.roc2(45), prec=3), '6367381.816, 6388838.29')
            self.test('roc2', fstr(E.roc2(90), prec=3), '6399593.626, 6399593.626')

            self.test('rocBearing', E.rocBearing( 0,  0), '6335439.327', fmt='%.3f')
            self.test('rocBearing', E.rocBearing(45, 45), '6378092.008', fmt='%.3f')
            self.test('rocBearing', E.rocBearing(90,  0), '6399593.626', fmt='%.3f')

            self.test('rocGauss', E.rocGauss(0),  '6356752.314', fmt='%.3f')
            self.test('rocGauss', E.rocGauss(45), '6378101.030', fmt='%.3f')
            self.test('rocGauss', E.rocGauss(90), '6399593.626', fmt='%.3f')

            self.test('rocMean', E.rocMean(0),  '6356716.465', fmt='%.3f')
            self.test('rocMean', E.rocMean(45), '6378092.008', fmt='%.3f')
            self.test('rocMean', E.rocMean(90), '6399593.626', fmt='%.3f')

            self.test('rocMeridional', fstr(E.rocMeridional(0),  prec=3), '6335439.327')
            self.test('rocMeridional', fstr(E.rocMeridional(45), prec=3), '6367381.816')
            self.test('rocMeridional', fstr(E.rocMeridional(90), prec=3), '6399593.626')

            self.test('rocPrimeVertical', fstr(E.rocPrimeVertical(0),  prec=3), '6378137.0')
            self.test('rocPrimeVertical', fstr(E.rocPrimeVertical(45), prec=3), '6388838.29')
            self.test('rocPrimeVertical', fstr(E.rocPrimeVertical(90), prec=3), '6399593.626')

        self.subtitle(ellipsoids, Ellipsoid.__init__)
        self.test('a, b, None',  Ellipsoid(1000, 500, None, name='_f_None').f_,  2.0)  # coverage
        self.test('a, None, f_', Ellipsoid(1000, None, 2,   name='_b_None').b, 500.0)  # coverage

        E = Ellipsoids.WGS84.copy()  # coverage
        self.subtitle(ellipsoids, E.name)
        self.test('WGS84.copy', E is not Ellipsoids.WGS84, True)
        self.test('WGS84.copy', E == Ellipsoids.WGS84, True)
        self.test('WGS84.find', Ellipsoids.find(E), None)

        self.test('WGS84.a2_b', E.a2_b, E.a2 / E.b,      prec=6)
        self.test('WGS84.b2_a', E.b2_a, E.b2 / E.a,      prec=6)
        self.test('WGS84.R2',   E.R2,   E.Rauthalic,     prec=6)  # obvious
        self.test('WGS84.c2',   E.c2,   E.R2**2,         fmt='%.0f')
        self.test('WGS84.es',   E.es,   E.e,             prec=6)
        self.test('WGS84.f2',   E.f2, (E.a - E.b) / E.b, prec=6)
        self.test('WGS84.m2degrees', int(E.m2degrees(E.a * PI_2)), 90)
        self.test('WGS84.area',   E.area,   '5.101e+14', fmt='%.3e')
        self.test('WGS84.volume', E.volume, '1.083e+21', fmt='%.3e')

        self.test('WGS84.ecef', E.ecef().__class__, EcefKarney)
        self.test('WGS84.ecef', E.ecef().name, E.name)

        t = E.toStr(prec=10)
        self.test('WGS84', t, "name='WGS84', a=6378137, b=6356752.3142451793, f_=298.257223563, f=0.0033528107, f2=0.0033640898, n=0.0016792204, e=0.0818191908, e2=0.00669438, e22=0.0067394967, e32=0.0033584313, A=6367449.1458234144, L=10001965.7293127235, R1=6371008.7714150595, R2=6371007.1809184747, R3=6371000.790009154")
        e = (E.a - E.b) / (E.a + E.b) - E.n
        t = 'A=%.10f, e=%.10f, f_=%.10f, n=%.10f(%.10e)' % (E.A, E.e, E.f_, E.n, e)
        self.test('WGS84.', t, 'A=6367449.1458234144, e=0.0818191908, f_=298.2572235630, n=0.0016792204(-2.1684043450e-19)')

        self.subtitle(ellipsoids, 'Kruegers')

        def _AB(E, K, A, B):
            E.KsOrder = K
            self.test('WGS84.AlphaKs', fstr(E.AlphaKs, prec=12, fmt='e', ints=True), A)
            self.test('WGS84.BetaKs ', fstr(E.BetaKs,  prec=12, fmt='e', ints=True), B)

        _AB(E, 8, '8.377318206245e-04, 7.608527773573e-07, 1.197645503243e-09, 2.429170680397e-12, 5.711818370428e-15, 1.47999793138e-17, 4.107624109371e-20, 1.210785038923e-22',
                  '8.37732164058e-04, 5.90587015222e-08, 1.673482665344e-10, 2.164798110491e-13, 3.787930968626e-16, 7.236769021816e-19, 1.493479824778e-21, 3.259522545838e-24')

        _AB(E, 6, '8.377318206245e-04, 7.608527773572e-07, 1.197645503329e-09, 2.429170607201e-12, 5.711757677866e-15, 1.491117731258e-17',
                  '8.37732164058e-04, 5.90587015222e-08, 1.673482665284e-10, 2.164798040063e-13, 3.787978046169e-16, 7.248748890695e-19')

        _AB(E, 4, '8.377318206304e-04, 7.608527714249e-07, 1.197638001561e-09, 2.443376194522e-12',
                  '8.377321640601e-04, 5.905869567934e-08, 1.673488880355e-10, 2.167737763022e-13')

        P = Ellipsoid(E.b, E.a, name='Prolate')
        _TOL = ellipsoids._TOL

        self.subtitle(ellipsoids, P.name)
        for p, e in ((P.a, E.b), (P.b, E.a), (P.n, -E.n),
                     (P.R1, 6363880.543), (P.R2, 6363878.941), (P.R3, 6363872.564),
                     (P.Rbiaxial, E.Rbiaxial), (P.Rgeometric, E.Rgeometric),
                     (P.c2, 40498955180263.188), (P.area, 508924880289508.500),
                     (P.volume, 1079575530747445379072.000)):
            t = '%s [%s]' % (p.name, p.units)
            self.test(t, p, e, fmt='%.3f', known=abs(p - e) < _TOL)

        for E, el, ob, pr in ((E, True, True, False),
                              (P, True, False, True),
                              (Ellipsoids.Sphere, False, False, False)):
            self.subtitle(ellipsoids, 'AuxiliaryLats ' + E.name)
            self.test('isEllipsoidal', E.isEllipsoidal, el)
            self.test('isOblate',      E.isOblate,      ob)
            self.test('isProlate',     E.isProlate,     pr)
            self.test('isSpherical',   E.isSpherical,   not el)
            for d in range(-90, 91, 30):
                x = 'lat (%d.0)' % (d,)
                for aux in (E.auxAuthalic,  E.auxConformal,  E.auxGeocentric,
                            E.auxIsometric, E.auxParametric, E.auxRectifying):
                    a = aux(d)
                    t = fstr(a, prec=9)
                    self.test('%s(%d)' % (aux.__name__, d), t, t)
                    self.test('name', a.name, aux.__name__)
                    self.test('inverse', aux(a, inverse=True).toRepr(prec=2), x)

        self.subtitle(ellipsoids, 'Flattenings')

        self.test('_TOL', _TOL, _TOL)
        for n, E in Ellipsoids.items():  # includes f_None, b_None, Prolate
            if E.f and E.f_:
                e = E.f_ - 1 / E.f
                self.test(n + '.f_ - 1 / .f', e, e if abs(e) < _TOL else _TOL)
                e = E.f - 1 / E.f_
                self.test(n + '.f - 1 / .f_', e, e if abs(e) < _TOL else _TOL)  # PYCHOK attr

        self.subtitle(ellipsoids, Ellipsoid2.__name__)
        for n, E in tuple(Ellipsoids.items()):  # includes f_None, b_None, Prolate
            E2 = Ellipsoid2(E.a, E.f, name='_2_' + n)
            n, E2 = str(E2).split(', ', 1)
            _, E  = str(E).split(', ', 1)
            self.test(n, E2, E)

        self.subtitle(ellipsoids, a_f2Tuple.__name__)
        for n, E in tuple(Ellipsoids.items()):  # includes f_None, b_None, Prolate
            a_f = a_f2Tuple(E.a, E.f)
            E2 = Ellipsoid(a_f.a, a_f.b, name='_a_b_' + n)  # PYCHOK a
            n, E2 = str(E2).split(', ', 1)
            _, E  = str(E).split(', ', 1)
            self.test(n, E2, E)

        self.subtitle(ellipsoids, 'Functions')
        for n, E in tuple(Ellipsoids.items()):  # includes f_None, b_None, Prolate
            f_ = a_b2f_(E.a, E.b)
            self.test('%s(%s)' % (a_b2f_.__name__, E.name), f_, E.f_, fmt='%.8f', known=abs(f_ - E.f_) < _TOL)
            f2 = a_b2f2(E.a, E.b)
            self.test('%s(%s)' % (a_b2f2.__name__, E.name), f2, E.f2, fmt='%.8f', known=abs(f2 - E.f2) < _TOL)
            n = a_b2n(E.a, E.b)
            self.test('%s(%s)' % (a_b2n.__name__, E.name), n, E.n, fmt='%.8f', known=abs(n - E.n) < _TOL)
            a = b_f2a(E.b, E.f)
            self.test('%s(%s)' % (b_f2a.__name__, E.name), a, E.a, fmt='%.3f', known=abs(a - E.a) < _TOL)  # millimeter
            a = b_f_2a(E.b, E.f_)
            self.test('%s(%s)' % (b_f_2a.__name__, E.name), a, E.a, fmt='%.3f', known=abs(a - E.a) < _TOL)  # millimeter
            f = f_2f(E.f_)
            self.test('%s(%s)' % (f_2f.__name__, E.name), f, E.f, fmt='%.8f', known=abs(f - E.f) < _TOL)
            e2 = n2e2(E.n)
            self.test('%s(%s)' % (n2e2.__name__, E.name), e2, E.e2, fmt='%.8f', known=abs(e2 - E.e2) < _TOL)
            f = n2f(E.n)
            self.test('%s(%s)' % (n2f.__name__, E.name), f, E.f, fmt='%.8f', known=abs(f - E.f) < _TOL)

        t = P.roc1_.__name__ + ' '
        for E, x in ((Ellipsoids.WGS84, 1.863e-9), (P, 1.863e-9),
                     (Ellipsoids.SphereAuthalic, '0.0')):
            self.subtitle(ellipsoids, E.name)
            for d in range(0, 91, 5):
                s, c = sincos2d(d)
                n = E.roc1_(s)
                h = E.roc1_(s, c)
                e = abs(n - h)  # delta in meter
                self.test(t + str(d), e, x, known=e < 1.864e-9)
                n = E.roc2(d).prime_vertical
                e = abs(n - h)  # delta in meter
                self.test(t + str(d), e, x, known=e < 1.864e-9)

        n = circle4.__name__ + ' '
        self.subtitle(ellipsoids, circle4.__name__)
        for E in (Ellipsoids.WGS84, Ellipsoids.Sphere):
            self.subtitle(ellipsoids, E.name)
            for d in range(0, 91, 10):
                r = E.Rgeocentric(d)
                t = E.circle4(d)
                self.test(n + str(d), hypot_(t.radius, t.height), r, prec=6)
                t = circle4(E, d)
                self.test(n + str(d), hypot_(t.radius, t.height), r, prec=6)


if __name__ == '__main__':

    t = Tests(__file__, __version__, ellipsoids)
    t.testEllipsoids()
    t.results()
    t.exit()
