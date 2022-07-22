
# -*- coding: utf-8 -*-

# Test spherical earth model functions and methods.

__all__ = ('Tests',)
__version__ = '22.07.08'

from base import isPython2, isWindows, RandomLatLon
from testLatLon import Tests as _TestsLL
from testVectorial import Tests as _TestsV

from pygeodesy import F_D, F_DEG_, F_DMS, PI_4, R_M, \
                      classname, Ellipsoids, IntersectionError, \
                      latlonDMS, lonDMS
from math import radians

# <https://GeographicLib.SourceForge.io/html/python/examples.html>
Antarctica = ((-63.1, -58),
              (-72.9, -74),
              (-71.9,-102),
              (-74.9,-102),
              (-74.3,-131),
              (-77.5,-163),
              (-77.4, 163),
              (-71.7, 172),
              (-65.9, 140),
              (-65.7, 113),
              (-66.6,  88),
              (-66.9,  59),
              (-69.8,  25),
              (-70.0,  -4),
              (-71.0, -14),
              (-77.3, -33),
              (-77.9, -46),
              (-74.7, -61))  # open


class Tests(_TestsLL, _TestsV):

    def testSpherical(self, module, Sph=True, Nv=False):  # MCCABE 13

        self.subtitle(module, 'Spherical')

        LatLon = module.LatLon

        p = LatLon(51.8853, 0.2545)
        self.test('isSpherical', p.isSpherical, True)
        self.test('isEllipsoidal', p.isEllipsoidal, False)

        q = LatLon(49.0034, 2.5735)
        self.test('datum', q.datum, p.datum)  # coverage
        q.datum = Ellipsoids.SphereAuthalic
        self.test('ellipsoid', q.datum.ellipsoid, Ellipsoids.SphereAuthalic)
        q.datum = p.datum
        self.test('datum', q.datum, p.datum)

        self.test('isSpherical', q.isSpherical, True)
        self.test('isEllipsoidal', q.isEllipsoidal, False)

        i = p.intersection(108.55, q, 32.44)
        self.test('intersection1', i.toStr(F_D), '50.907608°N, 004.508575°E')  # 50.9076°N, 004.5086°E  # Trig
        self.test('intersection1', i.toStr(F_DMS), '50°54′27.39″N, 004°30′30.87″E')
        self.test('intersection1', isinstance(i, LatLon), True)

        REO = LatLon(42.600, -117.866)
        BKE = LatLon(44.840, -117.806)
        i = REO.intersection(51, BKE, 137)
        self.test('intersection2', isinstance(i, LatLon), True)
        self.test('intersection2', i.toStr(F_D), '43.5719°N, 116.188757°W')  # 43.572°N, 116.189°W
        self.test('intersection2', i.toStr(F_DMS), '43°34′18.84″N, 116°11′19.53″W')

        # <https://GitHub.com/ChrisVeness/geodesy/issues/46>
        p = LatLon(51.8853, 0.2545)
        q = LatLon(51.8763, 0.2545)  # identical lon
        i = p.intersection(110.8878, q, 54.4525)
        self.test('intersection3', i, '51.882166°N, 000.267801°E')  # 51°52′55.8″N, 000°16′04.08″E?

        p = LatLon(+30, 0)
        q = LatLon(-30, 0)  # identical, zero lon
        i = p.intersection(135, q, 45)
        self.test('intersection4', i, '00.0°N, 026.565051°E', known=abs(i.lat) < 1e-6)

        p = LatLon(0, -30)
        q = LatLon(0, +30)  # identical, zero lat
        i = p.intersection(45, q, 315)
        self.test('intersection5', i, '26.565051°N, 000.0°W', known=abs(i.lon) < 1e-6)

        # <https://GitHub.com/ChrisVeness/geodesy/blob/master/test/latlon-vectors-tests.js>
        STN = LatLon(51.8853, 0.2545)
        CDG = LatLon(49.0034, 2.5735)
        i = STN.intersection(108.547, CDG, 32.435)
        self.test('intersection6', i, '50.907809°N, 004.50841°E')  # 50.9078°N, 004.5084°E

        # courtesy of U{sbonaime<https://GitHub.com/mrJean1/PyGeodesy/issues/58>}
        s1, s2 = LatLon(8, 0), LatLon(0, 8.4)
        self.test('intersection7', s1.intersection(150.06, s2, 55.61).toStr(), '''01°52′46.41″S, 005°39′06.87″E''')  # 01°52′46.41″N, 174°20′53.13″W
        s1, s2 = LatLon(80, 0), LatLon(0, 84)
        self.test('intersection8', s1.intersection(150.06, s2, 55.61).toStr(), '''28°15′40.35″S, 032°14′33.2″E''')  # 28°15′40.35″N, 147°45′26.8″W

        # <https://GitHub.com/ChrisVeness/geodesy/blob/master/test/latlon-vectors-tests.js>
        # <https://GitHub.com/ChrisVeness/geodesy/blob/master/test/latlon-spherical-tests.js>
        N, E, S, W, p, q = 0, 90, 180, 270, LatLon(0, 1), LatLon(1, 0)
        self.test('toward 1,1 N,E nearest',        p.intersection(N, q, E), '00.999848°N, 001.0°E')
        self.test('toward 1,1 E,N nearest',        q.intersection(E, p, N), '00.999848°N, 001.0°E')
        self.test('toward 1,1 N,E antipodal',      LatLon(2, 1).intersection(N, q, E), '00.999848°S, 179.0°W')
        self.test('toward/away 1,1 N,W antipodal', p.intersection(N, q, W), '00.999848°N, 001.0°E')
        self.test('toward/away 1,1 W,N antipodal', q.intersection(W, p, N), '00.999848°S, 179.0°W')
        self.test('toward/away 1,1 S,E antipodal', p.intersection(S, q, E), '00.999848°S, 179.0°W')
        self.test('toward/away 1,1 E,S antipodal', q.intersection(E, p, S), '00.999848°N, 001.0°E')
        self.test('away 1,1 S,W antipodal',        p.intersection(S, q, W), '00.999848°S, 179.0°W')
        self.test('away 1,1 W,S antipodal',        q.intersection(W, p, S), '00.999848°S, 179.0°W')
        self.test('1E/90E N,E antipodal',          p.intersection(N, LatLon(1, 90), E), '00.017454°N, 001.0°E')
        self.test('1E/90E N,E nearest',            p.intersection(N, LatLon(1, 92), E), '00.017454°N, 179.0°W')

        # <https://GitHub.com/ChrisVeness/geodesy/blob/master/test/latlon-vectors-tests.js>
        p, r = LatLon(1, 3), LatLon(2, 2)
        self.test('brng+end 1a', q.intersection(p, r, S), '01.000305°N, 002.0°E')
        self.test('brng+end 1b', r.intersection(S, q, p), '01.000305°N, 002.0°E')
        self.test('brng+end 2a', q.intersection(p, r, N), '01.000305°S, 178.0°W')
        self.test('brng+end 2b', r.intersection(N, q, p), '01.000305°S, 178.0°W')

        i = LatLon(1, 1).intersection(LatLon(2, 2), LatLon(1, 4), LatLon(2, 3))
        self.test('intersection7', i, '02.499372°N, 002.5°E')  # 02.4994°N, 002.5°E'

        p = LatLon(0, 0)
        self.test('maxLat0',  p.maxLat( 0), '90.0')
        self.test('maxLat1',  p.maxLat( 1), '89.0')
        self.test('maxLat90', p.maxLat(90),  '0.0')
        self.test('minLat0',  p.minLat( 0), '-90.0')
        self.test('minLat1',  p.minLat( 1), '-89.0')
        self.test('minLat90', p.minLat(90),  '-0.0', known=True)

        self.test('parse', p.parse('0, 0', name='parse'), p)  # coverage

        if hasattr(LatLon, 'crossingParallels'):
            ps = p.crossingParallels(LatLon(60, 30), 30)
            t = ', '.join(map(lonDMS, ps))
            self.test('crossingParallels', t, '009°35′38.65″E, 170°24′21.35″E')

        if hasattr(LatLon, 'intersections2') and Sph:

            n = 'intersections2 (%s)' % (LatLon.__module__,)

            def _100p2(t, r, *s):
                e = max(abs(a.distanceTo(b) - r) for a in s
                                                 for b in t) / r
                return e, '%g (%% of radius)' % (e,)  # percentages

            # <https://GIS.StackExchange.com/questions/48937/calculating-intersection-of-two-circles>
            p = LatLon(37.673442, -90.234036)  # (-0.00323306, -0.7915,   0.61116)
            q = LatLon(36.109997, -90.953669)  # (-0.0134464,  -0.807775, 0.589337)
            t = p.intersections2(0.0312705, q, 0.0421788, radius=None, height=0)  # radii in radians
            s = ', '.join(latlonDMS(t, form=F_D, prec=6))
            self.test(n, s, '36.98931°N, 088.151425°W, 38.23838°N, 092.390487°W')
            s = ', '.join(latlonDMS(t, form=F_DEG_, prec=9))
            # 36.989311051533505, -88.15142628069133, 38.2383796094578, -92.39048549120287
            self.test(n, s, '36.989310429, -088.151425243, 38.238379679, -092.390486808', known=True)

            t = LatLon(30, 0).intersections2(PI_4, LatLon(-30, 0), PI_4, radius=None)  # radii in radians
            s = ', '.join(latlonDMS(t, form=F_D))
            self.test(n, s, '00.0°N, 035.26439°W, 00.0°N, 035.26439°E', known='S, ' in s)

            t = LatLon(0, 40).intersections2(PI_4, LatLon(0, -40), PI_4, radius=None)  # radii in radians
            s = ', '.join(latlonDMS(t, form=F_D))
            self.test(n, s, '22.622036°N, 000.0°E, 22.622036°S, 000.0°E', known='W' in s)

            t = LatLon(30, 20).intersections2(PI_4, LatLon(-30, -20), PI_4, radius=None)  # radii in radians
            s = ', '.join(latlonDMS(t, form=F_D))
            self.test(n, s, '14.612841°N, 026.110934°W, 14.612841°S, 026.110934°E')

            t = LatLon(0, 0).intersections2(PI_4, LatLon(0, 22.5), PI_4 / 2, radius=None)  # abutting
            s = ', '.join(latlonDMS(t, form=F_D))
            self.test(n, s, '00.000001°S, 045.0°E, 00.000001°N, 045.0°E', known=True)  # N-S

            # centers at 2 opposite corners of a "square" and
            # radius equal to length of square side, expecting
            # the other 2 as the intersections ... but the
            # longitudes are farther and farther out
            for d in range(5, 66, 5):
                p = LatLon(d, -d)
                q = LatLon(-d, d)
                r = radians(2 * d) * R_M
                t = p.intersections2(r, q, r, radius=R_M)
                if t[0] is t[1]:
                    s = ', '.join(latlonDMS(t[:1], form=F_D)) + ' abutting'
                else:
                    s = ', '.join(latlonDMS(t, form=F_D))
                d = '%s %d' % (n, d)
                self.test(d, s, s)
                _, s = _100p2(t, r, q, p)
                self.test(d, s, s)

            d_m = 5e-5  # 50 micrometer
            # courtesy of U{Samuel Čavoj<https://GitHub.com/mrJean1/PyGeodesy/issues/41>}
            R = RandomLatLon(LatLon, 178, 178)  # +/- 89
            r = R()
            s = latlonDMS(r, form=F_D) + ' Random +/- 89'
            self.test(n, s, s)
            for _ in range(12):
                p, q = R(), R()
                try:  # see .testEllipsoidal
                    i1, i2 = p.intersections2(r.distanceTo(p),
                                           q, r.distanceTo(q), radius=R_M)
                    d, d2 = r.distanceTo(i1), r.distanceTo(i2)
                    if d2 < d:
                        d, i1, i2 = d2, i2, i1
                    s = '%s  d %g meter' % (', '.join(latlonDMS((i1, i2), form=F_D)), d)
                    self.test(n, s, s)
                    if d > d_m:
                        raise IntersectionError(d=d, fmt_name_value='%s (%g)', txt='over')
                except IntersectionError as x:
                    self.test(n, str(x), 'd < %g m' % (d_m), known=True)  # too distant, near concetric, etc.

            # courtesy of U{AleixDev<https://GitHub.com/mrJean1/PyGeodesy/issues/43>}
            def _known(p):
                return int(p.lat) == 42 and int(p.lon) == 2

            n = 'trilaterate5 (%s) .' % (LatLon.__module__,)
            d  = 5110  # meter
            p1 = LatLon(42.688839, 2.438857)
            p2 = LatLon(42.635421, 2.522570)

            p3 = LatLon(42.630788, 2.500713)
            t = p1.trilaterate5(d, p2, d, p3, d, area=True)  # overlap, default eps=EPS1
            self.test(n + 'min',   t.min, '313.671', prec=3, known= 300 < t.min < 320)
            p = t.minPoint
            self.test(n + 'point', p.toStr(F_D, prec=8), '42.66937229°N, 002.48639477°E', known=_known(p))
            self.test(n + 'max',   t.max, '1591.044', prec=3, known=1570 < t.max < 1610)
            p = t.maxPoint
            self.test(n + 'point', p.toStr(F_D, prec=8), '42.65153054°N, 002.46822157°E', known=_known(p))
            self.test(n + 'n', t.n, t.n)

            t = p1.trilaterate5(d, p2, d, p3, d, area=False, eps=200)  # intersection
            self.test(n + 'min',   t.min, '133.815', prec=3, known= 100 < t.min < 200)
            p = t.minPoint
            self.test(n + 'inter', p.toStr(F_D, prec=8), '42.6767291°N, 002.49916157°E', known=_known(p))
            self.test(n + 'n', t.n, t.n)

            p3 = LatLon(42.64540, 2.504811)
            t = p1.trilaterate5(d, p2, d, p3, d, area=True)  # overlap, default eps=EPS1
            self.test(n + 'min',   t.min, '2403.293', prec=3, known=2390 < t.min < 2420)
            self.test(n + 'max',   t.max, '2403.293', prec=3, known=2390 < t.max < 2420)
            p = t.maxPoint
            self.test(n + 'point', p.toStr(F_D, prec=8), '42.66135649°N, 002.47981645°E', known=_known(p))
            self.test(n + 'min- is .maxPoint', t.minPoint is t.maxPoint, True)
            self.test(n + 'n', t.n, t.n)
            try:
                t = p1.trilaterate5(d, p2, d, p3, d, area=False, eps=1000)  # no intersection in 1000 meter
                self.test(n + 'inter', t.minPoint, IntersectionError.__name__)
            except IntersectionError as x:
                self.test(n + 'inter', str(x), str(x))

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
                    t = str(x).replace(',)', ')')
                    self.test('isenclosedBy', t, 'points[3] (%s(47°00′00.0″N, 003°00′00.0″E)): not convex' % (classname(p),))

        p = LatLon(51.127, 1.338)
        q = LatLon(50.964, 1.853)
        b = p.rhumbBearingTo(q)  # rhumbAzimuthTo
        self.test('rhumbBearingTo', b, 116.722, fmt='%.3f')  # 116.7

        d = p.rhumbDestination(40300, 116.7)
        self.test('rhumbDestination', d, '50.964155°N, 001.853°E')  # 50.9642°N, 001.8530°E
        self.test('rhumbDestination', isinstance(d, LatLon), True)

        d = p.rhumbDistanceTo(q)
        self.test('rhumbDistanceTo', d, 40307.8, fmt='%.1f')  # XXX 40310 ?

        m = p.rhumbMidpointTo(q)
        self.test('rhumbMidpointo-0.5', m, '51.0455°N, 001.595727°E')
        self.test('rhumbMidpointo', isinstance(m, LatLon), True)
        m = p.rhumbMidpointTo(q, fraction=0)
        self.test('rhumbMidpointo-0.0', m, '51.127°N, 001.338°E')
        m = p.rhumbMidpointTo(q, fraction=0.25)
        self.test('rhumbMidpointo-0.25', m, '51.08625°N, 001.46692°E')
        m = p.rhumbMidpointTo(q, fraction=0.75)
        self.test('rhumbMidpointo-0.75', m, '51.00475°N, 001.72442°E')
        m = p.rhumbMidpointTo(q, fraction=1)
        self.test('rhumbMidpointo-1.0', m, '50.964°N, 001.853°E')
        m = p.rhumbMidpointTo(q, fraction=2)
        self.test('rhumbMidpointo-2.0', m, '50.801°N, 002.366196°E')

        b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
        self.test('areaOf', module.areaOf(b), '8.66605875e+09', fmt='%.8e')  # 8666058750.718977
        self.test('perimeterOf', module.perimeterOf(b, closed=True),  '3.78258541e+05', fmt='%.8e')
        self.test('perimeterOf', module.perimeterOf(b, closed=False), '2.67063461e+05', fmt='%.8e')

        c = LatLon(0, 0), LatLon(1, 0), LatLon(0, 1)
        self.test('areaOf', module.areaOf(c), '6.18e+09', fmt='%.2e')
        self.test('perimeterOf', module.perimeterOf(c, closed=True),  '3.79639757e+05', fmt='%.8e')
        self.test('perimeterOf', module.perimeterOf(c, closed=False), '2.68444678e+05', fmt='%.8e')

        if hasattr(module, 'nearestOn3'):
            c, d, _ = module.nearestOn3(p, b)
            self.test('nearestOn3', c, '46.000996°N, 001.353049°E' if Nv else '46.0°N, 001.369324°E')
            self.test('nearestOn3', d, '569987.49' if Nv else '570101.83', fmt='%.2f')
            d = p.distanceTo(c)
            self.test('distanceTo', d, '569987.49' if Nv else '570101.82', fmt='%.2f')

            p = LatLon(47, 3)
            c, d, _ = module.nearestOn3(p, b)
            self.test('nearestOn3', c, '46.0°N, 002.0°E' if Nv else '46.0°N, 002.0°E')
            self.test('nearestOn3', d, '134989.80' if Nv else '134992.48', fmt='%.2f')
            d = p.distanceTo(c)
            self.test('distanceTo', d, '134989.80' if Nv else '134989.80', fmt='%.2f')

            p = LatLon(45, 2)
            b = LatLon(45, 1), LatLon(47, 3)
            if Nv:
                c, d, _ = module.nearestOn3(p, b)
                self.test('nearestOn3', c, '45.330691°N, 001.318551°E')
                self.test('distance', d, '64856.28', fmt='%.2f')
            else:
                c, d, a = p.nearestOn3(b, adjust=False)
                self.test('nearestOn3', c, '45.5°N, 001.5°E')
                self.test('distance', d, '78626.79', fmt='%.2f')
                self.test('angle', a, '315.00', fmt='%.2f')
                a = p.compassAngleTo(c, adjust=False)
                self.test('compassAngleTo', a, '315.00', fmt='%.2f')
                c, d, a = p.nearestOn3(b, adjust=True)
                self.test('nearestOn3', c, '45.331319°N, 001.331319°E')
                self.test('distance', d, '64074.48', fmt='%.2f')
                self.test('angle', a, '305.10', fmt='%.2f')
            d = p.distanceTo(c)
            self.test('distanceTo', d, '64856.28' if Nv else '64074.12', fmt='%.2f')
            a = p.compassAngleTo(c)  # adjust=True
            self.test('compassAngleTo', a, '304.54' if Nv else '305.10', fmt='%.2f')
            # TrigTrue vs Nvector closests
            p = LatLon(45.330691, 001.318551)
            d = p.distanceTo(LatLon(45.331319, 001.331319))
            self.test('difference', d, '1000.53', fmt='%.2f')

            if Sph:  # check nearestOn2/3 with closest on the segment
                b = LatLon(0, 1), LatLon(2, 3), LatLon(4, 5), LatLon(6, 7), LatLon(8, 9)
                for i in range(8):
                    p = LatLon(i + 2, i)
                    c, d, a = p.nearestOn3(b, adjust=False)
                    t = LatLon(p.lat - 1.5, p.lon + 1.5).toStr(F_D, prec=6)
                    self.test('nearestOn3', c.toStr(F_D, prec=6), t)
                    self.test('distance', d, '235880.385', fmt='%.3f')
                    self.test('angle', a, '135.00', fmt='%.2f')

                n = module.meanOf(b)  # coverage
                self.test('meanOf', n.toStr(F_D, prec=6), '04.004858°N, 004.990226°E')
                n = module.nearestOn3(p, b, LatLon=LatLon, adjust=False)[0]  # coverage
                self.test('nearestOn3', n, '07.5°N, 008.5°E')
                c = p.toCartesian()  # coverage
                self.test('toCartesian', c, '[6245667.211, 766871.506, 996645.349]')

        if hasattr(module, 'ispolar'):
            p = LatLon(85, 90), LatLon(85, 0), LatLon(85, -90), LatLon(85, -180)
            for _ in self.testiter():
                self.test('ispolar', module.ispolar(p), 'True')
            p = LatLon(85, 90), LatLon(85, 0), LatLon(85, -180)
            for _ in self.testiter():
                self.test('ispolar', module.ispolar(p), 'True', known=True)
            p = [LatLon(*ll) for ll in Antarctica]
            for _ in self.testiter():
                self.test('ispolar', module.ispolar(p), 'True', known=Nv)

        if hasattr(LatLon, 'triangle7'):
            t = LatLon(10, 10).triangle7(LatLon(70, -20), LatLon(70, 40))  # radius=R_M
            self.test('triangle7', t.toStr(prec=6), '(22.142114, 2189776.256122, 84.91384, 6998970.168483, 84.91384, 6998970.168483, 8479698853790.72168)', known=isPython2 and isWindows)
            t = LatLon(10, 10).triangle7(LatLon(70, -20), LatLon(70, 40), radius=None)
            self.test('triangle7', t.toRepr(prec=5), 'Triangle8Tuple(A=0.38645, a=0.34371, B=1.48203, b=1.09857, C=1.48203, c=1.09857, D=3.74234, E=0.20891)')

        if hasattr(LatLon, 'nearestOn'):
            # <https://GitHub.com/mrJean1/PyGeodesy/issues/25>
            a = LatLon(1, 1, height=100)
            b = LatLon(2, 2, height=200)
            t = LatLon(1, 2).nearestOn(a, b).toStr(form=F_D, prec=1)
            self.test('nearestOn', t, '01.5°N, 001.5°E, +149.99m')
            t = LatLon(1, 2).nearestOn3([a, b])[0].toStr(form=F_D, prec=1)
            self.test('nearestOn3', t, '01.5°N, 001.5°E, +149.99m')

        if hasattr(LatLon, 'midpointTo'):
            a = LatLon(52.205, 0.119, height=100)
            b = LatLon(48.857, 2.351, height=200)
            t = a.midpointTo(b).toStr(form=F_D, prec=4)
            self.test('midpointTo', t, '50.5363°N, 001.2746°E, +150.00m')
            t = a.midpointTo(b, fraction=0.5).toStr(form=F_D, prec=4)
            self.test('intermediateTo', t, '50.5363°N, 001.2746°E, +150.00m')


if __name__ == '__main__':

    from pygeodesy import sphericalNvector as Nv, \
                          sphericalTrigonometry as Trig

    t = Tests(__file__, __version__)

    t.testLatLon(Nv, Sph=True, Nv=True)
    t.testVectorial(Nv, Sph=True)  # special case
    t.testSpherical(Nv, Sph=False, Nv=True)

    t.testLatLon(Trig, Sph=True)
    t.testSpherical(Trig)

    t.results()
    t.exit()
