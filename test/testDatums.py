
# -*- coding: utf-8 -*-

# Test L{datums} and transforms.

__all__ = ('Tests',)
__version__ = '26.04.21'

from bases import TestsBase

from pygeodesy import Datum, Datums, Ellipsoid, Ellipsoids, map2, \
                      R_M, Similarity, Transform, Transforms
from pygeodesy.datums import _spherical_datum
from math import degrees


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
        t = T.inverse().inverse(name="ED50_")
        self.test('ED50.inverse().inverse()', t == T, True)

        S = Datums.Sphere
        self.test(S.name, _spherical_datum(R_M) is S, True)

    def testDatums(self):
        self.test(all.__name__, Datums.register(all), all.__name__, nl=1)
        n = 0
        for _, d in Datums.items(asorted=True):
            self.test(d.name, d, d, nl=1)
            e = d.ellipsoid
            self.test(e.name, e, e)
            t = d.transform
            self.test(t.name, t, t)
            n += 1
        self.test('total', n, 19, nl=1)

    def testSimilarity(self):
        for T, x in ((Similarity(tx=565.7381, ty=50.4018,  tz=465.2904, s=4.07244,
                                 rx=1.91514,  ry=-1.60363, rz=9.09546, name='_xRD2ETRS'),
                     '(3945490.847071, 325245.322152, 4984379.549295)'),
                     (Similarity(tx=-565.7346, ty=-50.4058, tz=-465.2895, s=-4.07242,
                                 rx=-1.91513,  ry=1.60365,  rz=-9.09546, name='_xETRS2RD'),
                     '(3944359.917728, 325145.721259, 4983449.176491)')):
            t = T.toStr(prec=8)
            self.test(T.name, t, t, nl=1)
            for r, s, t in zip((T.rx, T.ry, T.rz), (T.sx, T.sy, T.sz), 'xyz'):
                self.test(T.name + '.' + t, degrees(r), s / 3600, fmt='%.5e')
#           A = 3904046.1730176862, 368161.3118488123, 5013449.0509925615  # Amersfoort RD
            A = 3904046.17302,      368161.31185,      5013449.05099  # Amersfoort RD
            r = T.transform(3944925.380652793, 325195.5237056559, 4983914.362442657, A)
            self.test(T.name, r, x)
            r = T.transform(r.x, r.y, r.z, A, inverse=True)
            self.test(T.name, r, '(3944925.381, 325195.524, 4983914.362)', known=map2(int, r)==
                                  (3944925,     325195,     4983914))


if __name__ == '__main__':

    from pygeodesy import datums  # private

    t = Tests(__file__, __version__, datums)
    t.testDatum()
    t.testDatums()
    t.testSimilarity()
    t.results()
    t.exit()
