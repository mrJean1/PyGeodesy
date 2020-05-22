
# -*- coding: utf-8 -*-

# Test datums, ellipsoids and transforms.

__all__ = ('Tests',)
__version__ = '20.05.15'

from base import TestsBase

from pygeodesy import Datum, Datums, EcefKarney, Ellipsoid, Ellipsoids, \
                      fstr, R_M, PI_2, Transform, Transforms


class Tests(TestsBase):

    def testDatum(self):
        # datum module tests
        E = Ellipsoid(1000, 1000, 0, name='TestEllipsiod')
        self.test('ellipsoid', E is Ellipsoids.TestEllipsiod, True)
#       print(Ellipsoid())

        T = Transform(name='TestTransform')
        self.test('transform', T is Transforms.TestTransform, True)
#       print(Transform())

        D = Datum(E, T, name='TestDatum')
        self.test('datum', D is Datums.TestDatum, True)
#       print(Datum())

        e = Ellipsoids.unregister('TestEllipsiod')
        self.test(e.name, e, E)
        t = Transforms.unregister('TestTransform')
        self.test(t.name, t, T)
        d = Datums.unregister('TestDatum')
        self.test(d.name, d, D)

        T = Transforms.ED50
        t = T.inverse().inverse("ED50_")
        self.test('ED50.inverse().inverse()', t == T, True)

        # to show WGS84 and NAD83 are identical up to 3 decimals
        for E in (Ellipsoids.WGS84, Ellipsoids.GRS80):  # NAD83
            self.test('R1', E.R1, R_M, fmt='%.4f')
            self.test('R2', E.R2, '6371007.2', fmt='%.1f')
            self.test('R3', E.R3, '6371000.8', fmt='%.1f')
            self.test('Rr', E.Rr, '6367449.1', fmt='%.1f')  # 6367445.0
            self.test('Rs', E.Rs, '6367435.7', fmt='%.1f')

            self.test('Rgeocentric', E.Rgeocentric(0),  '6378137.000', fmt='%.3f')
            self.test('Rgeocentric', E.Rgeocentric(45), '6367489.544', fmt='%.3f')
            self.test('Rgeocentric', E.Rgeocentric(90), '6356752.314', fmt='%.3f')

            self.test('Rlat', E.Rlat(0),  '6378137.000', fmt='%.3f')
            self.test('Rlat', E.Rlat(45), '6367444.657', fmt='%.3f')
            self.test('Rlat', E.Rlat(90), '6356752.314', fmt='%.3f')

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

        self.test('a, b, None',  Ellipsoid(1000, 500, None).f_, 2.0)  # coverage
        self.test('a, None, f_', Ellipsoid(1000, None, 2).b, 500.0)  # coverage

        E = Ellipsoids.WGS84.copy()  # coverage
        self.test('WGS84.copy', E is not Ellipsoids.WGS84, True)
        self.test('WGS84.copy', E == Ellipsoids.WGS84, True)
        self.test('WGS84.find', Ellipsoids.find(E), None)

        self.test('WGS84.a2_b', E.a2_b, E.a2 / E.b, fmt='%.6f')
        self.test('WGS84.b2_a', E.b2_a, E.b2 / E.a, fmt='%.6f')
        self.test('WGS84.c',    E.c,    E.R2,       fmt='%.6f')
        self.test('WGS84.es',   E.es,   E.e,        fmt='%.6f')
        self.test('WGS84.f2',   E.f2, (E.a - E.b) / E.b, fmt='%.6f')
        self.test('WGS84.m2degrees', int(E.m2degrees(E.a * PI_2)), 90)
        self.test('WGS84.area',   E.area,   '5.101e+14', fmt='%.3e')
        self.test('WGS84.volume', E.volume, '1.083e+21', fmt='%.3e')

        self.test('WGS84.ecef', E.ecef().__class__, EcefKarney)
        self.test('WGS84.ecef', E.ecef().name, E.name)

        t = E.toStr(prec=10)
        self.test('WGS84', t, "name='WGS84', a=6378137, b=6356752.3142499998, f_=298.257223563, f=0.0033528107, e=0.0818191908, e2=0.00669438, e12=0.99330562, e22=0.0067394967, n=0.0016792204, R1=6371008.7714166669, R2=6371007.180920884, R3=6371000.7900107643, Rr=6367449.1458250275, Rs=6367435.6797186071")
        e = (E.a - E.b) / (E.a + E.b) - E.n
        t = 'A=%.10f, e=%.10f, f_=%.10f, n=%.10f(%.10e)' % (E.A, E.e, E.f_, E.n, e)
        self.test('WGS84.', t, "A=6367449.1458234144, e=0.0818191908, f_=298.2572235630, n=0.0016792204(-3.7914875232e-13)")

        def _AB(E, K, A, B):
            E.KsOrder = K
            self.test('WGS84.AlphaKs', fstr(E.AlphaKs, prec=12, fmt='%.*e', ints=True), A)
            self.test('WGS84.BetaKs ', fstr(E.BetaKs,  prec=12, fmt='%.*e', ints=True), B)

        _AB(E, 8, '8.377318206245e-04, 7.608527773572e-07, 1.197645503242e-09, 2.429170680397e-12, 5.711818370428e-15, 1.47999793138e-17, 4.107624109371e-20, 1.210785038923e-22',
                  '8.377321640579e-04, 5.90587015222e-08, 1.673482665344e-10, 2.164798110491e-13, 3.787930968626e-16, 7.236769021816e-19, 1.493479824778e-21, 3.259522545838e-24')

        _AB(E, 6, '8.377318206245e-04, 7.608527773572e-07, 1.197645503329e-09, 2.429170607201e-12, 5.711757677866e-15, 1.491117731258e-17',
                  '8.377321640579e-04, 5.90587015222e-08, 1.673482665284e-10, 2.164798040063e-13, 3.787978046169e-16, 7.248748890694e-19')

        _AB(E, 4, '8.377318206304e-04, 7.608527714249e-07, 1.197638001561e-09, 2.443376194522e-12',
                  '8.377321640601e-04, 5.905869567934e-08, 1.673488880355e-10, 2.167737763022e-13')


if __name__ == '__main__':

    from pygeodesy import datum  # private

    t = Tests(__file__, __version__, datum)
    t.testDatum()
    t.results()
    t.exit()
