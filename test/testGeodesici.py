
# -*- coding: utf-8 -*-

# Some C{geodesici} tests for classes C{Intersectool} and C{Intersector}.

__all__ = ('Tests',)
__version__ = '24.07.28'

from bases import GeodSolve, geographiclib, IntersectTool, TestsBase

from pygeodesy import euclid, geodesici, Intersectool, Intersector, LatLon_
from pygeodesy.interns import _DOT_


class _LL(LatLon_):

    def __abs__(self):
        return euclid(self.lat, self.lon)

    def __neg__(self):
        return _LL(-self.lat, -self.lon)

    def __sub__(self, other):
        return _LL(self.lat - other.lat,
                   self.lon - other.lon)


def _re(v, x):  # rel error
    e = 0 if v is None or x is None else abs(v - x)  # None iteration
    if e:
        m = max(abs(v), abs(x))
        if m:
            e = min(e, m, e / m)
    return e


class Tests(TestsBase):

    _Xdict2 = dict()
    _known1 = set(('a1M', 's1M', 's12', 'latM', 'lonX', 'azi1'))
    _known2 = set(('aMM', 'sMM', 'sX0', 'sB', 'iteration', 'latB', 'A', 'B'))

    def test_(self, name, v, **nl):
        t = ('%9g' % (v,)).strip() if isinstance(v, float) else str(v)
        X, x = self._Xdict2.setdefault(name, (v, t))
        if isinstance(X, dict):
            TestsBase.test(self, name, t, x, known=True, **nl)
            if t != x:
                if len(v) > len(X):
                    X, v = v, X
                for n, v in v.items():
                    x =  X.get(n, v)
                    e = _re(v, x)
                    if e > 1e-12:
                        m = _DOT_(name, n)
                        TestsBase.test(self, m, v, x, error=e, known=n in self._known1)
        else:
            e, n = _re(v,  X), name.split(_DOT_)[-1]  # _tailof
            TestsBase.test(self, name, t, x, error=e, known=(e < 9e-9  # max(e) 6e-16
                                                         or (n in self._known2)
                                                         or _re(v, -X) < 1e-12), **nl)

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

    def testIntersect(self, Intersect, G, **_C):
        self.subtitle(geodesici, G.__name__)

        I = Intersect(G())  # PYCHOK I
        self.test(I.named2, I, I)

        # <https://GeographicLib.sourceforge.io/C++/doc/classGeographicLib_1_1Intersect.html>
        a = I.Line( 0.0,  0.0,  45.0)
        b = I.Line(45.0, 10.0, 135.0)
        self.testItems('Closest.1',  I.Closest( a, b, **_C))
        self.testItems('Closest5.1', I.Closest5(a, b))

        # <https://GeographicLib.sourceforge.io/C++/doc/IntersectTool.1.html>
        a = I.Line(50.0, -4.0, -147.7)
        b = I.Line( 0.0,  0.0,   90.0)
        self.testItems('Closest.2',  I.Closest( a, b, **_C))
        self.testItems('Closest5.2', I.Closest5(a, b))

        a = I.Line(50.0,  -4.0, -147.7)
        b = I.Line( 0.0, 180.0,    0.0)
        self.testEnumerate('All',  I.All( a, b, **_C))
        self.testEnumerate('All5', I.All5(a, b))

        a = I.Line( 0.0,  0.0,  10.0, 10.0)
        b = I.Line(50.0, -4.0, -50.0, -4.0)
        self.testItems('Middle',  I.Middle(a, b, **_C))
        self.testItems('Middle5', I.Middle5(a, b))

        # % echo 0 0 10 10 50 -4 -50 -4 | IntersectTool -i -p 0 -C
        # -631414 5988887 0 -3
        # -4.05187 -4.00000 -4.05187 -4.00000 0
        self.testItems('Segment', I.Segment(a, b, **_C))
        self.testItems('Segment5', I.Segment5(a, b))

        self.testEnumerate('intersect7s', I.intersect7(LatLon_( 0,  0), LatLon_( 10, 10),
                                                       LatLon_(50, -4), LatLon_(-50, -4),
                                                       LatLon=_LL))

        # % echo 50N 4W 147.7W 0 0 90 | IntersectTool -e 6371e3 0 -c -p 0 -C
        # 6077191 -3318019 0
        # -0.00000 -29.83966 -0.00000 -29.83966 0
        S = Intersect(G(6371e3, 0.0), name='Test')  # PYCHOK I
        a = S.Line(50.0, -4.0, -147.7)
        b = S.Line( 0.0,  0.0,   90.0)
        self.testItems('Sphere.Closest', S.Closest(a, b, **_C))
        self.testItems('Sphere.Closest5', S.Closest5(a, b))

        try:
            n = I.invokation
            self.test('invokations', n, n, nl=1)
        except AttributeError:
            pass
        n = len(self._Xdict2)  # not a test
        self.test('_Xdict2', n, n, nl=1)


if __name__ == '__main__':

    from pygeodesy import GeodesicExact

    t = Tests(__file__, __version__, geodesici)

    for _C in (False, True):
        t._Xdict2 = dict()

        if IntersectTool:  # first, to populate _xdict, ...
            t.testIntersect(Intersectool, GeodesicExact, _C=_C)

        if GeodSolve:  # ... or ...
            from pygeodesy.geodsolve import GeodesicSolve
            t.testIntersect(Intersector, GeodesicSolve, _C=_C)

        t.testIntersect(Intersector, GeodesicExact)  # ... otherwise

        if geographiclib:
            from pygeodesy.geodesicw import Geodesic
            t.testIntersect(Intersector, Geodesic, _C=_C)

    t.results()
    t.exit()
