
# -*- coding: utf-8 -*-

# Test module attributes.

__all__ = ('Tests',)
__version__ = '21.06.01'

from base import coverage, GeodSolve, isPyPy, isPython2, \
                 numpy_version, TestsBase

from pygeodesy import EPS, F_D, fstr, IntersectionError, NEG0, \
                      sphericalNvector, trilaterate2d2, \
                      vector3d, VectorError
from pygeodesy.interns import _DOT_  # INTERNAL


class Tests(TestsBase):

    def testIntersection3d3(self):

        i3d = vector3d.intersection3d3
        V3d = vector3d.Vector3d
        self.subtitle(vector3d, i3d.__name__)

        # <https://www.MathOpenRef.com/coordintersection.html>
        s1, e1 = V3d(15, 10, 1), V3d(49, 25, 2)
        s2, e2 = V3d(29, 5, 3),  V3d(32, 32, 4)
        self.test('(30, 17)', i3d(s1, e1, s2, e2, useZ=False), '(Vector3d(30.30584, 16.75258, 0.0), 0, 0)')
        s2 = V3d(7, 10, 5)
        self.test('(-1,  3)', i3d(s1, e1, s2, e2, useZ=False), '(Vector3d(-1.0429, 2.92225, 0.0), -1, -2)')
        s2 = V3d(62, 32, 6)
        self.test('(65, 32)', i3d(s1, e1, s2, e2, useZ=False), '(Vector3d(64.86667, 32.0, 0.0), 1, -2)')
        try:
            s2 = V3d(32 - (49 - 15), 32 - (25 - 10), 7)
            self.test('(-2, 17)', i3d(s1, e1, s2, e2, useZ=False), IntersectionError.__name__)
        except Exception as x:
            self.test('(-2, 17)', x.__class__, IntersectionError)
        self.test('(49, 25)', i3d(s1, e1, e1, e1, useZ=False), '(Vector3d(49.0, 25.0, 0.0), 0, 0)')

    def testNvectorBase(self, module, **kwds):

        try:
            Nvector = module.Nvector
            c = Nvector.__name__
        except AttributeError:
            Nvector = module.NvectorBase
            c = 'Vector4Tuple'
        self.subtitle(module, Nvector.__name__)

        v = Nvector(0.500, 0.500, 0.707, **kwds)
        s = module.sumOf((v, v), h=0, name='sumOf')
        self.test('sumOf', s.__class__.__name__, c)

        p = v.toLatLon(LatLon=None)
        c = v.toCartesian(Cartesian=None)
        self.test('ecef.x, .y, .z', fstr(p[:3],  prec=5), fstr(c[:3],  prec=5))
        self.test('ecef.lat, .lon', fstr(p[3:5], prec=6), fstr(c[3:5], prec=6))
        self.test('ecef.height', fstr(p.height, prec=6), fstr(c.height, prec=6), known=True)
        if c.M is not None:
            self.test('ecef.M', fstr(p.M, prec=9), fstr(c.M, prec=9))

        if coverage:
            from pygeodesy.namedTuples import LatLon2Tuple, LatLon3Tuple, \
                                              PhiLam2Tuple, PhiLam3Tuple

            self.test('.isEllipsoidal', v.isEllipsoidal, not v.isSpherical)
            self.test('.isSpherical',   v.isSpherical,   not v.isEllipsoidal)

            self.test('.latlon', v.latlon, LatLon2Tuple(v.lat, v.lon))
            self.test('.philam', v.philam, PhiLam2Tuple(v.phi, v.lam))

            self.test('.latlonheight', v.latlonheight, LatLon3Tuple(v.lat, v.lon, v.h), known=v.h in (0, 0.0, NEG0))
            self.test('.philamheight', v.philamheight, PhiLam3Tuple(v.phi, v.lam, v.h), known=v.h in (0, 0.0, NEG0))

            t = v.parse('0.5, 0.5, 0.707')
            self.test('parse', t, v)
            self.test('cmp', t.cmp(v), 0)

            self.test('eq', t == v, True)
            self.test('ge', t >= v, True)
            self.test('gt', t >  v, False)
            self.test('le', t <= v, True)
            self.test('lt', t <  v, False)
            self.test('ne', t != v, False)

            m = t * 2
            self.test('*', m, '(1.0, 1.0, 1.414)')
            self.test('+', t + v, m)
            self.test('/', m / 2, t)
            self.test('-', m - t, t)

            m = v.__matmul__(t)
            self.test('@', m, '(0.0, 0.0, 0.0)')
            r = t.__rmatmul__(m)
            self.test('@', r, m)

            r = v.rotate(m, 45)
            self.test('rotate', r, '(0.26268, 0.26268, 0.37143)')

            r.crosserrors = True
            self.test('crosserrors', r.crosserrors, True)

            try:
                self.test('0', v.dividedBy(0), VectorError.__name__)
            except Exception as x:
                self.test('0', str(x), 'factor (0): float division by zero')

            t = vector3d.intersections2(Nvector(   0, 0, 0), 500,
                                        Nvector(1000, 0, 0), 500, sphere=False)
            self.test('intersections2', t[0], t[1])  # abutting

            p1, p2 = Nvector(-100, -100, -100), Nvector(100, 100, 100)
            t = vector3d.nearestOn(Nvector(0, 0, 0), p1, p2)
            self.test('nearestOn', t, Nvector(0, 0, 0))
            t = vector3d.nearestOn(Nvector(-200, -200, 0), p1, p2)
            self.test('nearestOn', t is p1, True)
            t = vector3d.nearestOn(Nvector(200, 200, 0), p1, p2)
            self.test('nearestOn', t, p2)
            self.test('nearestOn', t is p2, True)
            t = vector3d.iscolinearWith(Nvector(200, 200, 0), p1, p2)
            self.test('iscolinearWith', t, False)
            t = vector3d.iscolinearWith(Nvector(0, 0, 0), p1, p2)
            self.test('iscolinearWith', t, True)


    def testVectorial(self, module):  # MCCABE 14

        self.subtitle(module, 'Vectorial')

        LatLon, Nvector, meanOf, sumOf = module.LatLon, module.Nvector, module.meanOf, module.sumOf

        if hasattr(LatLon, 'crossTrackDistanceTo'):
            p = LatLon(53.2611, -0.7972)
            s = LatLon(53.3206, -1.7297)
            d = p.crossTrackDistanceTo(s, 96.0)
            self.test('crossTrackDistanceTo', d, '-305.67', prec=2)  # -305.7
            e = LatLon(53.1887, 0.1334)
            d = p.crossTrackDistanceTo(s, e)
            self.test('crossTrackDistanceTo', d, '-307.55', prec=2)  # -307.5

        if hasattr(LatLon, 'enclosedby'):
            r = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
            p = LatLon(45.1, 1.1)
            self.test('enclosedby', p.enclosedby(r), True)
            r = LatLon(45,1), LatLon(46,2), LatLon(45,2), LatLon(46,1)
            self.test('enclosedby', p.enclosedby(r), False)

        v = Nvector(0.500, 0.500, 0.707)
        p = v.toLatLon()
        self.test('toLatLon', p, '44.995674°N, 045.0°E')  # 45.0°N, 45.0°E
        c = p.toNvector()
        self.test('toNvector', c, '(0.50004, 0.50004, 0.70705)')  # 0.500, 0.500, 0.707
        self.test('isequalTo', c.isequalTo(v), False)
        self.test('isequalTo', c.isequalTo(v, units=True), True)
        self.test('length', v.length, '0.99992449715', prec=11)
        self.test('euclid', v.euclid, '0.99995577', prec=8)
        self.test('length', c.length, '1.00', prec=2)
        self.test('euclid', c.euclid, '1.0000', prec=4)

        s = meanOf((p, p, p, p), height=0, LatLon=LatLon)
        self.test('meanOf', s, '44.995674°N, 045.0°E')
        self.test('meanOf', s.__class__.__name__, LatLon.__name__)

        class Nv(Nvector):
            pass
        v = Nvector(52.205, 0.119, 0.0)
        s = sumOf((v, c), Vector=Nv, h=0, name='sumOf')
        self.test('sumOf', s, '(52.70504, 0.61904, 0.70705)')
        self.test('sumOf', s.__class__.__name__, 'Nv')
        self.test('sumOf', s._name, 'sumOf')
        self.test('length', s.length, '52.7134151513', prec=10)

        c = v.copy()
        self.test('copy', c.isequalTo(v), True)
        self.test('length', v.length, '52.2051356286', prec=10)
        self.test('length', c.length, '52.2051356286', prec=10)

        if module is sphericalNvector:  # coverage
            c = p.toCartesian()
            self.test('toCartesian', c, '[3185744.919, 3185744.919, 4504643.315]')
            self.test('toLatLon',  c.toLatLon(), p, known=True)  # '44.995674°N, 045.0°E, -0.00m'
            self.test('toNvector', c.toNvector(), '(0.50004, 0.50004, 0.70705, -0.00)', known=True)

        if hasattr(LatLon, 'intersection'):
            # <https://GitHub.com/ChrisVeness/geodesy/blob/master/test/latlon-vectors-tests.js>
            p = LatLon(1, 1)
            i = p.intersection(LatLon(2, 2), LatLon(1, 4), LatLon(2, 3))
            self.test('intersection', i, '02.499372°N, 002.5°E')

        if hasattr(LatLon, 'isenclosedBy'):
            p = LatLon(45.1, 1.1)

            b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
            for _ in self.testiter():
                self.test('isenclosedBy', p.isenclosedBy(b), True)

            b = LatLon(45, 1), LatLon(45, 3), LatLon(46, 2), LatLon(47, 3), LatLon(47, 1)
            for _ in self.testiter():
                try:
                    self.test('isenclosedBy', p.isenclosedBy(b), True)  # Nvector
                except ValueError as x:
                    t = ' '.join(str(x).split()[:3] + ['...)'])
                    self.test('isenclosedBy', t, 'non-convex: (LatLon(45°00′00.0″N, 001°00′00.0″E), ...)')  # Trig

        if hasattr(LatLon, 'iswithin'):
            # courtesy of Paulius Šarka  psarka  Aug 30, 2017
            p = LatLon(1, 1).iswithin(LatLon(2, 2), LatLon(2, 2))
            self.test('iswithin', p, False)
            p = LatLon(2, 2).iswithin(LatLon(2, 2), LatLon(2, 2))
            self.test('iswithin', p, True)

        if hasattr(LatLon, 'nearestOn'):
            s1 = LatLon(51.0, 1.0)
            s2 = LatLon(51.0, 2.0)

            s = LatLon(41.0, 0.0)
            p = s.nearestOn(s1, s2, within=True)
            self.test('nearestOn', p.toStr(F_D, prec=3), '51.0°N, 001.0°E')
            p = s.nearestOn(s1, s2, within=False)
            self.test('nearestOn', p.toStr(F_D, prec=3), '50.987°N, 000.298°W')

            s = LatLon(61.0, 3.0)
            p = s.nearestOn(s1, s2, within=True)
            self.test('nearestOn', p.toStr(F_D, prec=3), '51.0°N, 002.0°E')
            p = s.nearestOn(s1, s2, within=False)
            self.test('nearestOn', p.toStr(F_D, prec=3), '50.995°N, 002.655°E')

            s = LatLon(51.0, 1.9)
            p = s.nearestOn(s1, s2)  # 51.0004°N, 001.9000°E
            self.test('nearestOn', p.toStr(F_D, prec=3), '51.0°N, 001.9°E')
            self.test('nearestOn', isinstance(p, LatLon), True)

            d = p.distanceTo(s)  # 42.71 m
            self.test('distanceTo', d, 42.712 if module is sphericalNvector else 42.826, prec=3, known=int(d) == 42)
            s = LatLon(51.0, 2.1)
            p = s.nearestOn(s1, s2)  # 51.0000°N, 002.0000°E
            self.test('nearestOn', p.toStr(F_D), '51.0°N, 002.0°E')
            self.test('nearestOn', isinstance(p, LatLon), True)

            # courtesy AkimboEG on GitHub
            s1 = LatLon(0, 0)
            s2 = LatLon(0, 1)
            s = LatLon(1, 0)
            p = s.nearestOn(s1, s2)  # 0.0°N, 0.0°E
            self.test('nearestOn', p, '00.0°N, 000.0°E')
            self.test('nearestOn', isinstance(p, LatLon), True)

            p = LatLon(10, -140).nearestOn(LatLon(0, 20), LatLon(0, 40))
            self.test('nearestOn', p, '00.0°N, 020.0°E')
            self.test('nearestOn', isinstance(p, LatLon), True)

            # courtesy of Paulius Šarka  psarka  Aug 30, 2017
            p = LatLon(1, 1).nearestOn(LatLon(2, 2), LatLon(2, 2))
            self.test('nearestOn', p, '02.0°N, 002.0°E')
            p = LatLon(2, 2).nearestOn(LatLon(2, 2), LatLon(2, 2))
            self.test('nearestOn', p, '02.0°N, 002.0°E')  # PYCHOK test attr?

        if hasattr(LatLon, 'triangulate'):
            # courtesy of pvezid  Feb 10, 2017
            p = LatLon("47°18.228'N","002°34.326'W")  # Basse Castouillet
            self.test('BasseC', p, '47.3038°N, 002.5721°W')
            s = LatLon("47°18.664'N","002°31.717'W")  # Basse Hergo
            self.test('BasseH', s, '47.311067°N, 002.528617°W')
            t = p.triangulate(7, s, 295)
            self.test('triangulate', t, '47.323667°N, 002.568501°W')
            self.test('triangulate', isinstance(t, LatLon), True)  # PYCHOK test attr?

        if hasattr(LatLon, 'trilaterate'):
            # <https://GitHub.com/chrisveness/geodesy/blob/master/test/latlon-nvector-spherical-tests.js>
            p = LatLon(37.418436, -121.963477)
            t = p.trilaterate(265.710701754, LatLon(37.417243, -121.961889), 234.592423446,
                                             LatLon(37.418692, -121.960194), 54.8954278262)
            self.test('trilaterate', t, '37.419078°N, 121.960579°W')
            self.test('trilaterate', isinstance(t, LatLon), True)

            # courtesy Carlos Freitas <https://GitHub.com/mrJean1/PyGeodesy/issues/33>
            b1 = LatLon(-8.068361, -34.892722)
            b2 = LatLon(-8.075917, -34.894611)
            b3 = LatLon(-8.076361, -34.908000)
            p  = LatLon(-8.068912, -34.888699)
            d1 = b1.distanceTo(p)
            d2 = b2.distanceTo(p)
            d3 = b3.distanceTo(p)
            t = b1.trilaterate(d1, b2, d2, b3, d3)
            self.test('trilaterate', t, p)
            self.test('trilaterate', isinstance(t, LatLon), True)
            t = b1.trilaterate(d1, b2, d2, b3, d3, useZ=True)
            self.test('trilaterate', t, p, known=abs(t.lon - p.lon) < 1e-5)
            self.test('trilaterate', isinstance(t, LatLon), True)

            # courtesy AleixDev <https://GitHub.com/mrJean1/PyGeodesy/issues/43>
            d  = 5110  # meter
            p1 = LatLon(42.688839, 2.438857)
            p2 = LatLon(42.635421, 2.522570)
            t = p1.trilaterate(d, p2, d, LatLon(42.630788,2.500713), d)
            self.test('trilaterate', t.toStr(F_D, prec=8), '42.67456065°N, 002.49539502°E')
            try:
                t = p1.trilaterate(d, p2, d, LatLon(42.64540, 2.504811), d)
                self.test('trilaterate', t.toStr(F_D, prec=8), IntersectionError.__name__)
            except IntersectionError as x:
                self.test('trilaterate', str(x), str(x))  # PYCHOK test attr?

        self.testNvectorBase(module)

    def testTrilaterate3d(self, module, Vector):

        try:  # need numpy
            self.subtitle(module, Vector.__name__)

            # <https://GISPoint.de/artikelarchiv/avn/2008/avn-ausgabe-12008/
            # 2303-direkte-loesung-des-raeumlichen-bogenschnitts-mit-methoden-der-linearen-algebra.html>
            c1, r1 = Vector( 50.0, 150.0, 13.0), 74.760
            c2, r2 = Vector(100.0,  50.0, 11.0), 84.623
            c3, r3 = Vector(200.0, 120.0, 12.0), 82.608
            t = c1.trilaterate3d2(r1, c2, r2, c3, r3)
            k = numpy_version < 1.12  # numpy too old?

            n = _DOT_(c1.named3, Vector.trilaterate3d2.__name__)
            self.test(n, len(t), 2)
            self.test(t[0].named3, fstr(t[0].xyz, prec=4), '119.8958, 130.6508, -5.1451', known=k)
            self.test(t[1].named3, fstr(t[1].xyz, prec=4), '119.9999, 129.9999, 30.0019', known=k)

            # <https://www.ResearchGate.net/publication/
            # 275027725_An_Algebraic_Solution_to_the_Multilateration_Problem>
            c1, r1 = Vector(27.297, -4.953, 1.470), 3.851  # 3.857
            c2, r2 = Vector(25.475, -6.124, 2.360), 3.875  # 3.988
            c3, r3 = Vector(22.590,  0.524, 1.200), 3.514  # 3.497
            for C, V in ((vector3d, dict(Vector=Vector)),
                         (Vector, {})):  # method
                f = C.trilaterate3d2
                t = f(c1, r1, c2, r2, c3, r3, **V)
                n = _DOT_(C.__name__, f.__name__)

                self.test(n, len(t), 2)
                self.test(t[0].named3, fstr(t[0].xyz, prec=5), '24.31229, -2.52045, 1.53649', known=k)
                self.test(t[1].named3, fstr(t[1].xyz, prec=5), '24.35062, -2.48109, 1.66673', known=k)

                try:  # concentric
                    t = f(c3, r1, c2, r2, c3, r3, **V)
                    self.test(n, t, IntersectionError.__name__)
                except IntersectionError as x:
                    t = str(x)
                    self.test(n, t, t)

                try:  # too distant
                    t = f(c1, r1 * 0.1, c2, r2 * 0.1, c3, r3, **V)
                    self.test(n, t, IntersectionError.__name__)
                except IntersectionError as x:
                    t = str(x)
                    self.test(n, t, t)

                try:  # colinear
                    t = f(Vector(0, 0, 0), 10, Vector(0, 9, 0), 10, Vector(0, -9, 0), 10)
                    self.test(n, t, IntersectionError.__name__)
                except IntersectionError as x:
                    t = str(x)
                    self.test(n, t, t)

        except ImportError as x:
            self.skip(str(x), n=14)

    def testTrilaterate2d2(self):

        self.subtitle(vector3d, trilaterate2d2.__name__.capitalize())
        # courtesy MartinManganiello <https://GitHub.com/mrJean1/PyGeodesy/issues/49>
        t = trilaterate2d2(2, 2, 1, 3, 3, 1, 1, 4, 1.4142)
        self.test(trilaterate2d2.__name__, t.toStr(prec=1), '(2.0, 3.0)')  # XXX ?
        try:
            t = trilaterate2d2(2, 2, 1, 3, 3, 1, 1, 4, 1.4142, eps=EPS)
            self.test(trilaterate2d2.__name__, t.toStr(prec=1), IntersectionError.__name__)
        except ValueError as x:
            self.test(trilaterate2d2.__name__, str(x), str(x))

        t = trilaterate2d2(-500, -200, 450,
                            100, -100, 694.6221994724903,
                            500,  100, 1011.1874208078342)
        self.test(trilaterate2d2.__name__, t.toStr(prec=1), '(-500.0, 250.0)')

    def testTrilaterate3d2(self, Vector):

        try:  # need numpy
            trilaterate3d2 = vector3d.trilaterate3d2
            self.subtitle(vector3d, trilaterate3d2.__name__.capitalize())
            n = _DOT_(vector3d.__name__, trilaterate3d2.__name__)

            # coutesy Martin Manganiello <https://GitHub.com/mrJean1/PyGeodesy/issues/49>
            c1, r1 = Vector(-500, -200, 0),  450.0
            c2, r2 = Vector( 100, -100, 0),  694.6221994724903
            c3, r3 = Vector( 500,  100, 0), 1011.1874208078342

            try:  # no intersection, no perturbation for eps=0
                t = trilaterate3d2(c1, r1, c2, r2, c3, r3, eps=0)
                self.test(n, t, IntersectionError.__name__)
            except IntersectionError as x:
                t = str(x)
                self.test(n, t, t)

            # by default, perturbation at EPS and qsrt(EPS)
            t = trilaterate3d2(c1, r1, c2, r2, c3, r3)  # eps=EPS
            self.test(n, len(t), 2)
            self.test(t[0].named3, fstr(t[0].xyz, prec=5), '-500.0, 250.0, -0.00535', known=isPyPy and isPython2)  # -0.00531
            self.test(t[1].named3, fstr(t[1].xyz, prec=5), '-500.0, 250.0, 0.00535',  known=isPyPy and isPython2)  # 0.00531

        except ImportError as x:
            self.skip(str(x), n=4)


if __name__ == '__main__':

    from pygeodesy import Datums, cartesianBase, ellipsoidalExact, ellipsoidalKarney, \
                          ellipsoidalNvector, ellipsoidalVincenty, nvectorBase, \
                          sphericalTrigonometry

    t = Tests(__file__, __version__)

    t.testVectorial(ellipsoidalNvector)
    t.testVectorial(sphericalNvector)

    t.testNvectorBase(nvectorBase, datum=Datums.Sphere)
    t.testNvectorBase(nvectorBase, datum=Datums.WGS84)

    t.testTrilaterate3d(sphericalNvector,      sphericalNvector.Cartesian)
    t.testTrilaterate3d(sphericalTrigonometry, sphericalTrigonometry.Cartesian)

    t.testTrilaterate3d(ellipsoidalNvector,  ellipsoidalNvector.Cartesian)
    t.testTrilaterate3d(ellipsoidalVincenty, ellipsoidalVincenty.Cartesian)
    t.testTrilaterate3d(ellipsoidalKarney,   ellipsoidalKarney.Cartesian)
    t.testTrilaterate3d(ellipsoidalExact,    ellipsoidalExact.Cartesian)

    if GeodSolve:
        from pygeodesy import ellipsoidalGeodSolve
        t.testTrilaterate3d(ellipsoidalGeodSolve, ellipsoidalGeodSolve.Cartesian)

    t.testTrilaterate3d(cartesianBase,         cartesianBase.CartesianBase)
    t.testTrilaterate3d(nvectorBase,           nvectorBase.NvectorBase)
    t.testTrilaterate3d(vector3d,              vector3d.Vector3d)

    t.testTrilaterate2d2()
    t.testTrilaterate3d2(vector3d.Vector3d)

    t.testIntersection3d3()

    t.results()
    t.exit()
