
# -*- coding: utf-8 -*-

# Test datums, ellipsoids and transforms.

__all__ = ('Tests',)
__version__ = '22.07.01'

from base import TestsBase

from pygeodesy import Datum, Datums, Ellipsoid, Ellipsoids, \
                      R_M, Transform, Transforms
from pygeodesy.datums import _spherical_datum


class Tests(TestsBase):

    def testDatums(self):
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


if __name__ == '__main__':

    from pygeodesy import datums  # private

    t = Tests(__file__, __version__, datums)
    t.testDatums()
    t.results()
    t.exit()
