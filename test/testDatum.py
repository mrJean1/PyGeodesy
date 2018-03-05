
# -*- coding: utf-8 -*-

# Test datums, ellipsoids and transforms.

__all__ = ('Tests',)
__version__ = '18.03.04'

from base import TestsBase

from pygeodesy import R_M, Datum, Datums, Ellipsoid, Ellipsoids, \
                      fStr, Transform, Transforms


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

        E = Ellipsoids.WGS84
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

        self.test('distance2', fStr(E.distance2(0, 0,    1,  1), prec=3),  '156903.472, 45.192')
        self.test('distance2', fStr(E.distance2(0, 0,   10, 10), prec=3), '1569034.719, 45.192')
        self.test('distance2', fStr(E.distance2(40, 40, 50, 50), prec=3), '1400742.676, 37.563')
        self.test('distance2', fStr(E.distance2(70, 70, 80, 80), prec=3), '1179164.848, 18.896')

        self.test('roc2', fStr(E.roc2(0),  prec=3), '6335439.327, 6378137.0')
        self.test('roc2', fStr(E.roc2(45), prec=3), '6367381.816, 6388838.29')
        self.test('roc2', fStr(E.roc2(90), prec=3), '6399593.626, 6399593.626')

        self.test('rocBearing', E.rocBearing(0, 0),   '6335439.327', fmt='%.3f')
        self.test('rocBearing', E.rocBearing(45, 45), '6378092.008', fmt='%.3f')
        self.test('rocBearing', E.rocBearing(90, 0),  '6399593.626', fmt='%.3f')

        self.test('rocGauss', E.rocGauss(0),  '6356752.314', fmt='%.3f')
        self.test('rocGauss', E.rocGauss(45), '6378101.030', fmt='%.3f')
        self.test('rocGauss', E.rocGauss(90), '6399593.626', fmt='%.3f')

        self.test('rocMean', E.rocMean(0),  '6356716.465', fmt='%.3f')
        self.test('rocMean', E.rocMean(45), '6378092.008', fmt='%.3f')
        self.test('rocMean', E.rocMean(90), '6399593.626', fmt='%.3f')

        E = Datums.WGS84.ellipsoid
        e = (E.a - E.b) / (E.a + E.b) - E.n
        t = (E.toStr(prec=10),
            'A=%.10f, e=%.10f, f_=%.10f, n=%.10f(%.10e)' % (E.A, E.e, E.f_, E.n, e),
            'Alpha6=(%s)' % (fStr(E.Alpha6, prec=12, fmt='%.*e', ints=True),),
            'Beta6=(%s)'  % (fStr(E.Beta6,  prec=12, fmt='%.*e', ints=True),))
        self.test('WGS84', t[0], "name='WGS84', a=6378137, b=6356752.3142499998, f_=298.257223563, f=0.0033528107, e=0.0818191908, e2=0.00669438, e22=0.0067394967, n=0.0016792204, R1=6371008.7714166669, R2=6371007.180920884, R3=6371000.7900107643, Rr=6367449.1458250266, Rs=6367435.6797186071")
        self.test('WGS84', t[1], "A=6367449.1458234144, e=0.0818191908, f_=298.2572235630, n=0.0016792204(-3.7914875232e-13)")
        self.test('WGS84', t[2], "Alpha6=(0, 8.377318206245e-04, 7.608527773572e-07, 1.197645503329e-09, 2.429170607201e-12, 5.711757677866e-15, 1.491117731258e-17)")
        self.test('WGS84', t[3], "Beta6=(0, 8.377321640579e-04, 5.905870152220e-08, 1.673482665284e-1, 2.164798040063e-13, 3.787978046169e-16, 7.248748890694e-19)")


if __name__ == '__main__':

    from pygeodesy import datum  # private

    t = Tests(__file__, __version__, datum)
    t.testDatum()
    t.results(nl=0)
    t.exit()
