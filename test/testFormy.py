
# -*- coding: utf-8 -*-

# Test L{formy} module.

__all__ = ('Tests',)
__version__ = '23.08.15'

from bases import TestsBase

from pygeodesy import PI, PI_2, R_M, antipode, bearing, cosineAndoyerLambert, \
                      cosineForsytheAndoyerLambert as _cosineForsythe_, \
                      cosineLaw, Datums, equirectangular, euclidean, \
                      excessAbc_, excessCagnoli_, excessGirard_, excessLHuilier_, \
                      excessKarney, excessQuad, flatLocal, flatPolar, formy, \
                      hartzell, haversine, heightOf, heightOrthometric, horizon, hubeny, \
                      IntersectionError, intersection2, intersections2, \
                      isantipode, isantipode_, isnormal, isnormal_, \
                      LatLon_, latlonDMS, LimitError, limiterrors, map1, normal, \
                      parseDMS, radical2, thomas, Vector3d as V3, vincentys

from math import degrees, radians


class Tests(TestsBase):

    def testDistance(self, t, f, lat1, lon1, lat2, lon2, x, m, **kwds):
        d = f(lat1, lon1, lat2, lon2, wrap=True, **kwds)
        e = 100.0 * abs(d - x) / x
        t = '%s%s (%.2f%%)' % (f.__name__, t, e)
        self.test(t, d, x, fmt='%.3f', known=e < m)

    def testDistances(self, t, lat1, lon1, lat2, lon2, x=0):
        # allow 0.1% margin
        self.testDistance(t, haversine,            lat1, lon1, lat2, lon2, x, 0.1)
        # assumed as reference
        self.testDistance(t, vincentys,            lat1, lon1, lat2, lon2, x, 0.0)
        self.testDistance(t, vincentys,            lat1, lon1, lat2, lon2, x, 0.0,
                                                   radius=Datums.Sphere)  # datums._mean_radius
        # allow 0.4% margin
        self.testDistance(t, cosineAndoyerLambert, lat1, lon1, lat2, lon2, x, 0.4)
        # allow 0.4% margin
        self.testDistance(t, _cosineForsythe_,     lat1, lon1, lat2, lon2, x, 0.4)
        # allow 0.1% margin
        self.testDistance(t, cosineLaw,            lat1, lon1, lat2, lon2, x, 0.1)
        # allow 3% margin, 10% for Auckland-Santiago
        self.testDistance(t, equirectangular,      lat1, lon1, lat2, lon2, x, 10, limit=None)  # =90)
        # allow 13% margin
        self.testDistance(t, euclidean,            lat1, lon1, lat2, lon2, x, 13)
        # allow .3% margin, 8% for long distances
        self.testDistance(t, flatLocal,            lat1, lon1, lat2, lon2, x, 8.1)
        # allow 21% margin, 57% for Auckland-Santiago
        self.testDistance(t, flatPolar,            lat1, lon1, lat2, lon2, x, 57)
        # allow 0.4% margin
        self.testDistance(t, thomas,               lat1, lon1, lat2, lon2, x, 0.4)
        # same as flatLocal
        self.test('hubeny' + t, hubeny, flatLocal, nt=1)

    def testDistances2(self, t, lat1, lon1, lat2, lon2, x=0, **datum):
        # allow 0.1% margin
        self.testDistance(t, haversine,            lat1, lon1, lat2, lon2, x, 0.1)
        # assumed as reference
        self.testDistance(t, vincentys,            lat1, lon1, lat2, lon2, x, 0.1)
        # allow 0.4% margin
        self.testDistance(t, cosineAndoyerLambert, lat1, lon1, lat2, lon2, x, 0.4, **datum)
        # allow 0.4% margin
        self.testDistance(t, _cosineForsythe_,     lat1, lon1, lat2, lon2, x, 0.4, **datum)
        # allow 0.1% margin
        self.testDistance(t, cosineLaw,            lat1, lon1, lat2, lon2, x, 0.1)
        # allow .3% margin, 16% for long distances
        self.testDistance(t, flatLocal,            lat1, lon1, lat2, lon2, x, 16, **datum)
        # allow 0.4% margin
        self.testDistance(t, thomas,               lat1, lon1, lat2, lon2, x, 0.4, **datum)
        # same as flatLocal
        self.test('hubeny' + t, hubeny, flatLocal, nt=1)

    def testFormy(self):

        self.test('antipode1', antipode( 89,  179), (-89.0, -1.0))
        self.test('antipode2', antipode(-89, -179), ( 89.0,  1.0))

        # roughly, Newport to New York
        self.test('bearing1', bearing(41.49, -71.31, 40.78, -73.97),              251.364, fmt='%.3f')
        self.test('bearing2', bearing(41.49, -71.31, 40.78, -73.97, final=False), 251.364, fmt='%.3f')
        self.test('bearing3', bearing(41.49, -71.31, 40.78, -73.97, final=True),  249.614, fmt='%.3f')

        # <https://codegolf.StackExchange.com/questions/63870/spherical-excess-of-a-triangle>
        # t = sphericalTrigonometry.triangle7(10, 10, 70, -20, 70, 40, radius=None) ... -8Tuple
        # (A=0.386453, a=0.34371, B=1.482026, b=1.098565, C=1.482026, c=1.098565, D=3.742345, E=0.208912)')
        # XXX with the full t.A, t.a, ... values, all 5 excess' produce the exact same results
        for f, t in ((excessAbc_,      (0.386453, 1.098565, 1.098565)),
                     (excessAbc_,      (1.482026, 0.34371,  1.098565)),
                     (excessCagnoli_,  (0.34371,  1.098565, 1.098565)),
                     (excessGirard_,   (0.386453, 1.482026, 1.482026)),
                     (excessLHuilier_, (0.34371,  1.098565, 1.098565))):
            self.test(f.__name__, degrees(f(*t)), '11.9698', prec=4)
        for f, t, x in ((excessKarney, (70, -20, 70, 40),  '56.9625'),
                        (excessQuad,   (70, -20, 70, 40),  '56.9625'),
                        (excessKarney, ( 0, -20, 70, 40),  '44.0235'),
                        (excessQuad,   ( 0, -20, 70, 40),  '44.0235'),
                        (excessKarney, (70, 40,  0, -20), '-44.0235'),
                        (excessQuad,   (70, 40,  0, -20), '-44.0235')):
            self.test(f.__name__, degrees(f(*t, radius=None)), x, prec=4)

        self.test('isantipode1', isantipode( 89,  179, -89,  -1), True, nl=1)
        self.test('isantipode2', isantipode(-89, -179,  89,   1), True)
        self.test('isantipode3', isantipode(-89, -179, -89,  -1), False)
        self.test('isantipode4', isantipode(  0,  -45, -0., 135), True)

        self.test('isantipode5', isantipode_(*map1(radians,  89,  179, -89,  -1)), True)
        self.test('isantipode6', isantipode_(*map1(radians, -89, -179,  89,   1)), True)
        self.test('isantipode7', isantipode_(*map1(radians, -89, -179, -89,  -1)), False)
        self.test('isantipode8', isantipode_(*map1(radians,   0,  -45, -0., 135)), True)

        self.test('isnormal1', isnormal( 90, -180), True, nl=1)
        self.test('isnormal2', isnormal(-99,  180), False)
        self.test('isnormal3', isnormal(-0.,  180), True)
        self.test('isnormal4', isnormal(-89,  179, eps=1), True)

        self.test('isnormal5', isnormal_( PI_2, -PI), True)
        self.test('isnormal6', isnormal_(-PI,    PI), False)
        self.test('isnormal7', isnormal_(-0.,    PI), True)
        self.test('isnormal8', isnormal_(-PI_2+0.1, PI-0.1, eps=0.1), True)

        pov = V3(10.1e6, 10.2e6, 10.3e6)  # 10+ km
        los = V3(-0.7274, -0.3637, -0.5819)
        self.test('hartzell', hartzell(pov, los).toStr(prec=6), '(884080.396945, 5592040.198472, 2927517.711001)', nl=1)
        ll = hartzell(pov, los, LatLon=LatLon_)
        self.test('hartzell', ll.toStr(prec=6), "27.500482°N, 081.016111°E, +0.00, 'hartzell'", known=abs(ll.height) < 1e-6)
        self.test('hartzell', hartzell(pov).toStr(prec=6), '(3642031.283571, 3678090.99925, 3714150.714929)')
        ll = hartzell(pov, LatLon=LatLon_)
        self.test('hartzell', ll.toStr(prec=6), "35.843738°N, 045.282243°E, +0.00, 'hartzell'", known=abs(ll.height) < 1e-6)

        self.test('heightOf0',   heightOf(0,   R_M), 2638958.23912, fmt='%.5f')
        self.test('heightOf45',  heightOf(45,  R_M), 5401080.43931, fmt='%.5f')
        self.test('heightOf90',  heightOf(90,  R_M), R_M,           fmt='%.5f')  # Height
        self.test('heightOf135', heightOf(135, R_M), 5401080.43931, fmt='%.5f')

        self.test('horizon0',     horizon(0), 0.0)
        self.test('horizon10Km',  horizon(10000), '357099.672', fmt='%.3f')
        self.test('horizon30Kft', horizon(10000, refraction=True), '392310.704', fmt='%.3f')
        self.test('horizon10Kft', horizon( 3000, refraction=True), '214877.422', fmt='%.3f')

        self.test('normal1', normal(-89,  179), (-89.0,  179.0))
        self.test('normal2', normal( 99,    0), ( 81.0,  180.0))
        self.test('normal3', normal( 99, -199), ( 81.0,  -19.0))
        self.test('normal4', normal(-99,  180), (-81.0,    0.0), nt=1)

        #              lat          lon
        Auckland   = -36.8485,    174.7633
        Boston     =  42.3541165, -71.0693514
        Cleveland  =  41.499498,  -81.695391
        LosAngeles =  34.0522,   -118.2437
        MtDiablo   =  37.8816,   -121.9142
        Newport    =  41.49008,   -71.312796
        NewYork    =  40.7791472, -73.9680804
        Santiago   = -33.4489,    -70.6693
        X          =  25.2522,     55.28
        Y          =  14.6042,    120.982
        # <https://GeographicLib.SourceForge.io/cgi-bin/GeodSolve>, <https://www.Distance.to>
        for i,  (ll1,        ll2,        expected) in enumerate((
                (Boston,     NewYork,      298009.404),    # ..328,722.580        370 km
                (Boston,     Newport,       98164.988),    # ..106,147.318         99 km
                (Cleveland,  NewYork,      651816.987),    # ..736,534.840        536 km
                (NewYork,    MtDiablo,    4084985.780),    # ..4,587,896.452    3,952 km
                (Auckland,   Santiago,    9670051.606),    # ..15,045,906.074   9,665 km
                (Auckland,   LosAngeles, 10496496.577),    # ..13,002,288.857  10,940 km
                (LosAngeles, Santiago,    8998396.669),    # ..10,578,638.162   8,993 km
                (X,          Y,           6906867.946))):  # ..6916,085.326     6,907 km
            self.testDistances(str(i + 1), *(ll1 + ll2), x=expected)

        # Thomas' S2, page 153 <https://apps.DTIC.mil/dtic/tr/fulltext/u2/703541.pdf>
        self.testDistances2('9', parseDMS('75 57 42.053'), 0,
                                 parseDMS('17  5 21.296'), parseDMS('85 31 54.631'),
                                 x=8044806.076, datum=Datums.NAD27)  # Clarke1866 ellipsoid

        f = intersection2
        self.test(f.__name__, f.__module__, formy.__name__)
        for datum in (None, R_M, -R_M, Datums.WGS84):
            t = str(intersection2(0, 0, 30, 0, 30, -30, datum))
            self.test('%s(%s)' % (f.__name__, datum), t, '(24.284126, 15.0)', known=t.endswith(' 15.0)'))

        f = intersections2
        self.test(f.__name__, f.__module__, formy.__name__, nl=1)
        for datum in (None, R_M, Datums.WGS84):
            self.testIntersections2(datum)

        n = radical2.__name__
        self.test(n, radical2(10, 4, 8), '(0.26, 2.6)', nl=1)
        self.test(n, radical2(10, 8, 4), '(0.74, 7.4)')
        self.test(n, radical2(10, 5, 5), '(0.5, 5.0)')
        self.test(n, radical2( 0, 5, 5), '(0.5, 0.0)')

        try:  # coverage
            self.test(n, radical2(10, 5, 4), IntersectionError.__name__, nt=1)
        except Exception as x:
            self.test(IntersectionError.__name__, str(x), 'distance (10.0), ...', known=True, nt=1)

        n = LimitError.__name__
        t = limiterrors(True)
        try:  # coverage
            self.test(n, equirectangular(0, 0, 60, 120), n)
        except Exception as x:
            self.test(n, str(x), 'delta exceeds ...', known=True)
        limiterrors(t)

        p = LatLon_(1, 2, height=3)
        self.test(heightOrthometric.__name__, heightOrthometric(p, 4), -1.0, nl=1)
        self.test(heightOrthometric.__name__, heightOrthometric(5, 4),  1.0)

    def testIntersections2(self, datum):
        # centers at 2 opposite corners of a "square" and
        # radius equal to length of square side, expecting
        # the other 2 as the intersections ... but the
        # longitudes are farther and farther out
        for d in (1, 2, 5, 10, 20, 30, 40):
            r = radians(2 * d) * R_M
            n = '%s(%s) %d' % (intersections2.__name__, getattr(datum, 'name', datum), d)
            try:
                t = intersections2(d, -d, r, -d, d, r, datum=datum)
                if t[0] is t[1]:
                    s = ', '.join(latlonDMS(t[:1], prec=4)) + ' abutting'
                else:
                    s = ', '.join(latlonDMS(t, prec=4))
                self.test(n, s, s)
            except IntersectionError as x:
                self.test(n, str(x), '2-tuple', known=True)


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testFormy()
    t.results()
    t.exit()
