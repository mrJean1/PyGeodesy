
# -*- coding: utf-8 -*-

# Test formulary functions.

__all__ = ('Tests',)
__version__ = '20.07.29'

from base import TestsBase

from pygeodesy import R_M, antipode, bearing, cosineAndoyerLambert, \
                      cosineForsytheAndoyerLambert as _cosineForsythe_, \
                      cosineLaw, equirectangular, Datums, euclidean, \
                      flatLocal, flatPolar, formy, haversine, heightOf, \
                      horizon, hubeny, IntersectionError, intersections2, \
                      isantipode, isantipode_, latlonDMS, map1, parseDMS, \
                      thomas, vincentys

from math import radians


class Tests(TestsBase):

    def testDistance(self, t, f, lat1, lon1, lat2, lon2, x, m, **kwds):
        d = f(lat1, lon1, lat2, lon2, wrap=True, **kwds)
        e = 100.0 * abs(d - x) / x
        t = '%s%s (%.1f%%)' % (f.__name__, t, e)
        self.test(t, d, x, fmt='%.3f', known=e < m)

    def testDistances(self, t, lat1, lon1, lat2, lon2, x=0):
        # allow 0.1% margin
        self.testDistance(t, haversine,            lat1, lon1, lat2, lon2, x, 0.1)
        # assumed as reference
        self.testDistance(t, vincentys,            lat1, lon1, lat2, lon2, x, 0.0)
        # allow 0.4% margin
        self.testDistance(t, cosineAndoyerLambert, lat1, lon1, lat2, lon2, x, 0.4)
        # allow 0.4% margin
        self.testDistance(t, _cosineForsythe_,     lat1, lon1, lat2, lon2, x, 0.4)
        # allow 0.1% margin
        self.testDistance(t, cosineLaw,            lat1, lon1, lat2, lon2, x, 0.1)
        # allow 3% margin, 10% for Auckland-Santiago
        self.testDistance(t, equirectangular,      lat1, lon1, lat2, lon2, x, 10, limit=None)  # =90)
        # allow 13% margin
        self.testDistance(t, euclidean,            lat1, lon1, lat2, lon2, x,  13)
        # allow .3% margin, 8% for long distances
        self.testDistance(t, flatLocal,            lat1, lon1, lat2, lon2, x, 8.1)
        # allow 21% margin, 57% for Auckland-Santiago
        self.testDistance(t, flatPolar,            lat1, lon1, lat2, lon2, x,  57)
        # allow 0.4% margin
        self.testDistance(t, thomas,               lat1, lon1, lat2, lon2, x, 0.4)
        # same as flatLocal
        self.test('hubeny', hubeny, flatLocal)

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
        self.testDistance(t, flatLocal,            lat1, lon1, lat2, lon2, x,  16, **datum)
        # allow 0.4% margin
        self.testDistance(t, thomas,               lat1, lon1, lat2, lon2, x, 0.4, **datum)
        # same as flatLocal
        self.test('hubeny', hubeny, flatLocal)

    def testFormy(self):

        self.test('antipode1', antipode( 89,  179), (-89, -1))
        self.test('antipode2', antipode(-89, -179), (89, 1))

        # roughly, Newport to New York
        self.test('bearing1', bearing(41.49, -71.31, 40.78, -73.97),              251.364, fmt='%.3f')
        self.test('bearing2', bearing(41.49, -71.31, 40.78, -73.97, final=False), 251.364, fmt='%.3f')
        self.test('bearing3', bearing(41.49, -71.31, 40.78, -73.97, final=True),  249.614, fmt='%.3f')

        self.test('isantipode1', isantipode( 89,  179, -89, -1), True)
        self.test('isantipode2', isantipode(-89, -179,  89,  1), True)
        self.test('isantipode3', isantipode(-89, -179, -89, -1), False)

        self.test('isantipode4', isantipode_(*map1(radians,  89,  179, -89, -1)), True)
        self.test('isantipode5', isantipode_(*map1(radians, -89, -179,  89,  1)), True)
        self.test('isantipode6', isantipode_(*map1(radians, -89, -179, -89, -1)), False)

        self.test('heightOf0',   heightOf(0,   R_M), 2638958.23912, fmt='%.5f')
        self.test('heightOf45',  heightOf(45,  R_M), 5401080.43931, fmt='%.5f')
        self.test('heightOf90',  heightOf(90,  R_M), R_M,           fmt='%.5f')  # Height
        self.test('heightOf135', heightOf(135, R_M), 5401080.43931, fmt='%.5f')

        self.test('horizon0',     horizon(0), 0.0)
        self.test('horizon10Km',  horizon(10000), '357099.672', fmt='%.3f')
        self.test('horizon30Kft', horizon(10000, refraction=True), '392310.704', fmt='%.3f')
        self.test('horizon10Kft', horizon( 3000, refraction=True), '214877.422', fmt='%.3f')

        Auckland   = -36.8485, 174.7633
        Boston     = 42.3541165, -71.0693514
        Cleveland  = 41.499498, -81.695391
        LosAngeles = 34.0522, -118.2437
        MtDiablo   = 37.8816, -121.9142
        Newport    = 41.49008, -71.312796
        NewYork    = 40.7791472, -73.9680804
        Santiago   = -33.4489, -70.6693
        X          = 25.2522, 55.28
        Y          = 14.6042, 120.982
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

        self.test(intersections2.__name__, intersections2.__module__, formy.__name__, nl=1)
        for datum in (None, R_M, Datums.WGS84):
            self.testIntersections2(datum)

    def testIntersections2(self, datum):
        # centers at 2 opposite corners of a "square" and
        # radius equal to length of square side, expecting
        # the other 2 as the intersections ... but the
        # longitudes are farther and farther out
        for d in (1, 2, 5, 10, 20, 30, 40):
            r = radians(2 * d) * R_M
            n = 'intersection2 (%s) %d' % (getattr(datum, 'name', datum), d)
            try:
                t = intersections2(d, -d, r, -d, d, r, datum=datum)
                if t[0] is t[1]:
                    s = latlonDMS(t[:1], prec=4, sep=', ') + ' abutting'
                else:
                    s = latlonDMS(t, prec=4, sep=', ')
                self.test(n, s, s)
            except IntersectionError as x:
                self.test(n, str(x), '2-tuple', known=True)


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testFormy()
    t.results()
    t.exit()
