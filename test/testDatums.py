
# -*- coding: utf-8 -*-

# Test L{datums} and transforms.

__all__ = ('Tests',)
__version__ = '23.03.27'

from bases import TestsBase

from pygeodesy import Datum, Datums, Ellipsoid, Ellipsoids, \
                      R_M, Transform, Transforms
from pygeodesy.datums import _spherical_datum


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

        S = Datums.Sphere
        self.test(S.name, _spherical_datum(R_M) is S, True)

    def testDatums(self):
        n = 0
        for _, d in Datums.items(all=True, asorted=True):
            self.test(d.name, d, d, nl=1)
            e = d.ellipsoid
            self.test(e.name, e, e)
            t = d.transform
            self.test(t.name, t, t)
            n += 1
        self.test('total', n, 18, nl=1)


if __name__ == '__main__':

    from pygeodesy import datums  # private

    t = Tests(__file__, __version__, datums)
    t.testDatum()
    t.testDatums()
    t.results()
    t.exit()
