
# -*- coding: utf-8 -*-

# Some basic C{geodesici} class C{Intersector} tests.

__all__ = ('Tests',)
__version__ = '24.06.27'

from bases import GeodSolve, geographiclib, TestsBase

from pygeodesy import geodesici, Intersector
from pygeodesy.interns import _DOT_

_xdict   = {}  # GeodesicExact results
_xpected = _xdict.setdefault


class Tests(TestsBase):

    def test_(self, name, v, **kwds):
        t = ('%9g' % (v,)).strip() if isinstance(v, float) else str(v)
        x = _xpected(name, t)
        TestsBase.test(self, name, t, x, known=True, **kwds)

    def testEnumerate(self, name, y):
        for i, t in enumerate(y):
            n = '%s[%s]' % (name, i)
            self.testItems(n, t)
            self.test_(_DOT_(n, 'iteration'), t.iteration)
            if i > 9:
                break

    def testItems(self, name, t):
        nl = 1
        for n, v in t.items():
            self.test_(_DOT_(name, n), v, nl=nl)
            nl = 0

    def testIntersector(self, G):
        self.subtitle(geodesici, G.__name__)

        I = Intersector(G(), name='Test')  # PYCHOK I
        self.test('Intersector', I, I)

        # <https://GeographicLib.sourceforge.io/C++/doc/classGeographicLib_1_1Intersect.html>
        a = I.Line( 0.0,  0.0,  45.0)
        b = I.Line(45.0, 10.0, 135.0)
        self.testItems('Closest.1',  I.Closest( a, b))
        self.testItems('Closest5.1', I.Closest5(a, b))

        self.testEnumerate('Next5s', I.Next5s(a, b))

        # <https://GeographicLib.sourceforge.io/C++/doc/IntersectTool.1.html>
        a = I.Line(50.0, -4.0, -147.7)
        b = I.Line( 0.0,  0.0,   90.0)
        self.testItems('Closest.2',  I.Closest( a, b))
        self.testItems('Closest5.2', I.Closest5(a, b))

        a = I.Line(50.0,  -4.0, -147.7)
        b = I.Line( 0.0, 180.0,    0.0)
        self.testEnumerate('All',  I.All( a, b))
        self.testEnumerate('All5', I.All5(a, b))


if __name__ == '__main__':

    from pygeodesy import GeodesicExact

    t = Tests(__file__, __version__, geodesici)

    if GeodSolve:  # first, to populate _xdict, ...
        from pygeodesy.geodsolve import GeodesicSolve
        t.testIntersector(GeodesicSolve)

    t.testIntersector(GeodesicExact)  # ... otherwise

    if geographiclib:
        from pygeodesy.geodesicw import Geodesic
        t.testIntersector(Geodesic)

    n = len(_xdict)  # not a test
    t.test('_xdict', n, n, nl=1)

    t.results()
    t.exit()
