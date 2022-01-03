
# -*- coding: utf-8 -*-

# Test azimuthal projections and intersections2 functions.

__all__ = ('Tests',)
__version__ = '22.01.03'

from base import GeodSolve, geographiclib, TestsBase, RandomLatLon

from pygeodesy import Equidistant, EquidistantExact, EquidistantGeodSolve, \
                      EquidistantKarney, Gnomonic, GnomonicExact, \
                      GnomonicGeodSolve, GnomonicKarney, LambertEqualArea, \
                      Orthographic, Stereographic, ellipsoidalExact, \
                      ellipsoidalKarney, ellipsoidalNvector, \
                      ellipsoidalVincenty, F_D, fstr, hypot, \
                      IntersectionError, latlonDMS, latlonDMS_, vincentys


class Tests(TestsBase):

    def testAzimuthal(self, Azimuthal, *xs):

        P = Azimuthal(48 + 50/60.0, 2 + 20/60.0, name='Paris')
        self.test(repr(P), P, P)

        f = P.forward(50.9, 1.8, name='Calais')
        self.test('forward', fstr(f[:6], prec=6), xs[0])
        r = P.reverse(f.x, f.y, name='Calais')
        self.test('reverse', fstr(r[:6], prec=6), xs[1])
        self.test('iteration', P.iteration, P.iteration)

        self.testCopy(P)

        r = P.reverse(-38e3, 230e3, name='Calais')
        self.test('reverse', fstr(r[:6], prec=6), xs[2])
        f = P.forward(r.lat, r.lon, name='Calais')
        self.test('forward', fstr(f[:6], prec=6), xs[3])

        for m in (ellipsoidalExact, ellipsoidalKarney,
                  ellipsoidalNvector, ellipsoidalVincenty):
            r = P.reverse(-38e3, 230e3, LatLon=m.LatLon)
            self.test('reverse', repr(r), xs[4])

        G = Azimuthal(51.4934, 0.0098, name='Greenwich')
        self.test(repr(G), G, G)
        f = G.forward(P.lat0, P.lon0, P.name)
        self.test('forward', fstr(f[:6], prec=6), xs[5])
        r = G.reverse(f.x, f.y)
        self.test('reverse', fstr(r[:6], prec=6), xs[6])
        self.test('iteration', G.iteration, G.iteration)

        h = hypot(f.x, f.y)  # easting + norting ~= distance
        d = vincentys(G.lat0, G.lon0, r.lat, r.lon)  # haversine
        self.test('hypot', h, d, fmt='%.3f', known=abs(d - h) < 1000, nt=1)

    def testDiscrepancies(self):
        # Compare ellipsoidal intersections2 for EquidistantKarney
        # and Equidistant showing the differences in degrees
        # and as percentage of the reference RandomLatLon.  Also
        # show the first, spherical intersections2 gu-/estimates.

        # Equidistant implements Snyder's formulas for the sphere
        # for ellipsoidal projections.  That plus the (high) accuracy
        # of EquidistantKarney likely cause the discrepancies.  An
        # other factor may be the innate distortions of azimuthal
        # equidistance projections for distances beyond 10,000 Km
        # (about one quarter of the earth circumference), see
        # <https://WikiPedia.org/wiki/Azimuthal_equidistant_projection>
        from pygeodesy.ellipsoidalBaseDI import _intersects2 as _ei2
        from pygeodesy.sphericalTrigonometry import _intersects2 as _si2

        def _100p(p, q, w):
            r = abs(100 * p / q) if q else 0
            return '%0*.3f%%' % (w + 4,r)

        def _max(i, r, t=''):
            s = latlonDMS(i, form=F_D, prec=-6, sep=', ')
            return '%s  %s, %s of Random%s' % (s, _100p(i.lat, r.lat, 2),
                                                  _100p(i.lon, r.lon, 3), t)

        for m in (ellipsoidalExact, ellipsoidalVincenty):  # ellipsoidalKarney
            LL = m.LatLon
            e = LL(0, 0)
            n = m.__name__
            # courtesy of U{Samuel Čavoj<https://GitHub.com/mrJean1/PyGeodesy/issues/41>}
            R = RandomLatLon(LL, 90, 90)  # +/- 45
            r = R()
            s = latlonDMS(r, form=F_D, prec=-6) + ' Random +/- 45'
            self.test(n, s, s)
            for _ in range(12):  # 100+
                p, q = R(), R()
                r1 = r.distanceTo(p)
                r2 = r.distanceTo(q)
                t = []
                for E in (None, EquidistantExact, Equidistant):
                    a = getattr(E, '__name__', 'Spherical')
                    try:
                        i1, i2 = _si2(p, r1, q, r2, LatLon=LL) if E is None else \
                                 _ei2(p, r1, q, r2, equidistant=E, LatLon=LL)
                        d, d2 = r.distanceTo(i1), r.distanceTo(i2)
                        if d2 < d:
                            d, i1, i2 = d2, i2, i1
                        s = latlonDMS_(i1, i2, form=F_D, prec=-6, sep=', ')
                        s = '%s  d %g meter  %s' % (s, d, a)
                        self.test(n, s, s)
                        if E is not None:
                            t.append(i1)
                    except (IntersectionError, TypeError, ValueError) as x:
                        self.test(n, str(x), a, known=True)
                if len(t) == 2:
                    i1, i2 = t
                    i = LL(i1.lat - i2.lat, i1.lon - i2.lon)
                    s = _max(i, r)
                    self.test(n, s, s)
                    e = LL(max(abs(i.lat), e.lat),
                           max(abs(i.lon), e.lon))
            s = _max(e, r, ', max')
            self.test(n, s, s, nt=1)

    def testSnyder(self, lat, lon, As, *xys):
        ll = lat, lon
        t = str(ll)
        for A, xy in zip(As, xys):
            n = A.__class__.__name__ + t
            f = A.forward(lat, lon)
            self.test(n, fstr(f[:2], prec=5), fstr(xy, prec=5))
            r = A.reverse(f.x, f.y)
            self.test(n, fstr(r[2:4], prec=5), fstr(ll, prec=5))


if __name__ == '__main__':

    from pygeodesy import azimuthal, equidistant, gnomonic, namedTuples

    t = Tests(__file__, __version__, azimuthal)

    t.testAzimuthal(Equidistant,
                   '-37467.812512, 230294.518853, 50.9, 1.8, 350.759218, 1.000223',
                   '-37467.812512, 230294.518853, 50.9, 1.8, 350.759218, 1.000223',
                   '-38000.0, 230000.0, 50.897321, 1.792455, 350.61849, 1.000222',
                   '-38000.0, 230000.0, 50.897321, 1.792455, 350.61849, 1.000222',
                   'LatLon(50°53′50.36″N, 001°47′32.84″E)',
                   '170420.92566, -293667.828613, 48.833333, 2.333333, 149.872606, 1.000472',
                   '170420.92566, -293667.828613, 48.833333, 2.333333, 149.872606, 1.000472')

    t.testAzimuthal(Gnomonic,
                   '-37484.520018, 230397.210923, 50.9, 1.8, 350.759218, 1.000669',
                   '-37484.520018, 230397.210923, 50.9, 1.8, 350.759218, 1.000223',
                   '-38000.0, 230000.0, 50.896405, 1.792706, 350.61849, 1.000222',
                   '-38000.0, 230000.0, 50.896405, 1.792706, 350.61849, 1.000667',
                   'LatLon(50°53′47.06″N, 001°47′33.74″E)',
                   '170581.851218, -293945.134107, 48.833333, 2.333333, 149.872606, 1.001416',
                   '170581.851218, -293945.134107, 48.833333, 2.333333, 149.872606, 1.000472')

    t.testAzimuthal(EquidistantExact,
                   '-37526.978232, 230000.911579, 50.9, 1.8, 350.325442, 0.999778',
                   '-37526.978232, 230000.911579, 50.9, 1.8, 350.325442, 0.999778',
                   '-38000.0, 230000.0, 50.899962, 1.793278, 350.205524, 0.999778',
                   '-38000.0, 230000.0, 50.899962, 1.793278, 350.205524, 0.999778',
                   'LatLon(50°53′59.86″N, 001°47′35.8″E)',
                   '170617.186469, -293210.754313, 48.833333, 2.333333, 151.589952, 0.999529',
                   '170617.186469, -293210.754313, 48.833333, 2.333333, 151.589952, 0.999529')

    t.testAzimuthal(GnomonicExact,
                   '-37543.665895, 230103.189403, 50.9, 1.8, 350.325442, 0.999333',
                   '-37543.665895, 230103.189403, 50.9, 1.8, 350.325442, 0.999333',
                   '-38000.0, 230000.0, 50.899044, 1.793528, 350.205718, 0.999333',
                   '-37999.995965, 229999.975581, 50.899044, 1.793528, 350.205718, 0.999334',
                   'LatLon(50°53′56.56″N, 001°47′36.7″E)',
                   '170778.089295, -293487.270649, 48.833333, 2.333333, 151.589952, 0.998587',
                   '170778.089295, -293487.270649, 48.833333, 2.333334, 151.589953, 0.998588')

    if GeodSolve:
        t.testAzimuthal(EquidistantGeodSolve,
                       '-37526.978232, 230000.911579, 50.9, 1.8, 350.325442, 0.999778',
                       '-37526.978232, 230000.911579, 50.9, 1.8, 350.325442, 0.999778',
                       '-38000.0, 230000.0, 50.899962, 1.793278, 350.205524, 0.999778',
                       '-38000.0, 230000.0, 50.899962, 1.793278, 350.205524, 0.999778',
                       'LatLon(50°53′59.86″N, 001°47′35.8″E)',
                       '170617.186469, -293210.754313, 48.833333, 2.333333, 151.589952, 0.999529',
                       '170617.186469, -293210.754313, 48.833333, 2.333333, 151.589952, 0.999529')

        t.testAzimuthal(GnomonicGeodSolve,
                       '-37543.665895, 230103.189403, 50.9, 1.8, 350.325442, 0.999333',
                       '-37543.665895, 230103.189403, 50.9, 1.8, 350.325442, 0.999333',
                       '-38000.0, 230000.0, 50.899044, 1.793528, 350.205718, 0.999333',
                       '-37999.995965, 229999.975581, 50.899044, 1.793528, 350.205718, 0.999334',
                       'LatLon(50°53′56.56″N, 001°47′36.7″E)',
                       '170778.089295, -293487.270649, 48.833333, 2.333333, 151.589952, 0.998587',
                       '170778.089295, -293487.270649, 48.833333, 2.333334, 151.589953, 0.998588')

    if geographiclib:
        t.testAzimuthal(EquidistantKarney,
                       '-37526.978232, 230000.911579, 50.9, 1.8, 350.325442, 0.999778',
                       '-37526.978232, 230000.911579, 50.9, 1.8, 350.325442, 0.999778',
                       '-38000.0, 230000.0, 50.899962, 1.793278, 350.205524, 0.999778',
                       '-38000.0, 230000.0, 50.899962, 1.793278, 350.205524, 0.999778',
                       'LatLon(50°53′59.86″N, 001°47′35.8″E)',
                       '170617.186469, -293210.754313, 48.833333, 2.333333, 151.589952, 0.999529',
                       '170617.186469, -293210.754313, 48.833333, 2.333333, 151.589952, 0.999529')

        t.testAzimuthal(GnomonicKarney,
                       '-37543.665895, 230103.189403, 50.9, 1.8, 350.325442, 0.999333',
                       '-37543.665895, 230103.189403, 50.9, 1.8, 350.325442, 0.999333',
                       '-38000.0, 230000.0, 50.899044, 1.793528, 350.205718, 0.999333',
                       '-37999.995965, 229999.975581, 50.899044, 1.793528, 350.205718, 0.999334',
                       'LatLon(50°53′56.56″N, 001°47′36.7″E)',
                       '170778.089295, -293487.270649, 48.833333, 2.333333, 151.589952, 0.998587',
                       '170778.089295, -293487.270649, 48.833333, 2.333334, 151.589953, 0.998588')
    else:
        t.skip('no geographiclib', n=14)

    t.testDiscrepancies()

    As = tuple(Azimuthal(0, 0, datum=1) for  # spherical datum of radius=1 for Snyder's Table ...
                        # ... 30, pp 196-197      26, pp 168           28, pp 188-189      22, pp 151          24, pp 158-159
               Azimuthal in  (Equidistant,        Gnomonic,            LambertEqualArea,   Orthographic,       Stereographic))
    t.testSnyder(10, 80, As, (1.37704, 0.24656), (5.67128,  1.01543), (1.26747, 0.22694), (0.96985, 0.17365), (1.65643, 0.29658))
    t.testSnyder(20, 20, As, (0.33454, 0.35601), (0.36397,  0.38733), (0.33123, 0.35248), (0.32139, 0.34202), (0.34136, 0.36327))
    t.testSnyder(40, 40, As, (0.57386, 0.74912), (0.8391,   1.09537), (0.55281, 0.72164), (0.4924,  0.64279), (0.62062, 0.81016))
    t.testSnyder(60, 60, As, (0.58948, 1.17896), (1.73205,  3.4641),  (0.54772, 1.09545), (0.43301, 0.86603), (0.69282, 1.38564))
    t.testSnyder(70, 80, As, (0.50997, 1.42273), (5.67128, 15.82209), (0.46280, 1.29114), (0.33682, 0.93969), (0.63588, 1.77402))
    t.testSnyder(80, 80, As, (0.26358, 1.51792), (5.67128, 32.65961), (0.23828, 1.37219), (0.17101, 0.98481), (0.33201, 1.91196))  # XXX 1.96962?
    t.testSnyder(80, 10, As, (0.04281, 1.39829), (0.17633,  5.75877), (0.03941, 1.28702), (0.03015, 0.98481), (0.05150, 1.68198))

    for A in (equidistant, gnomonic):
        t.test('function', A.__name__, A.__name__, nl=1)
        A = A(0, 0, datum=1, name='coverage')
        t.test('equatoradius', A.equatoradius, 1.0)
        t.test('flattening', A.flattening, 0.0)
        t.test('iteration', A.iteration, A.iteration)
        t.test('latlon0', A.latlon0, (0.0, 0.0))
        A.latlon0 = namedTuples.LatLon2Tuple(1, 2)
        t.test('latlon0', A.latlon0, (1.0, 2.0))
        t.test('name', A.name, 'coverage')
        t.test('radius', A.radius, 1.0)

    A = Stereographic(0, 0)  # coverage
    A.k0 = A.k0 + 1

    t.results()
    t.exit()
