
# -*- coding: utf-8 -*-

# Test LCC functions and methods.

__all__ = ('Tests',)
__version__ = '20.07.19'

from base import TestsBase, geographiclib

from pygeodesy import Equidistant, EquidistantKarney, Gnomonic, LambertEqualArea, \
                      Orthographic, Stereographic, \
                      ellipsoidalKarney, ellipsoidalNvector, ellipsoidalVincenty, \
                      fstr, haversine, hypot


class Tests(TestsBase):

    def testAzimuthal(self, Azimuthal, *x):

        P = Azimuthal(48 + 50/60.0, 2 + 20/60.0, name='Paris')
        self.test(repr(P), P, P)

        f = P.forward(50.9, 1.8, name='Calais')
        self.test('forward', fstr(f[:6], prec=6), x[0])
        r = P.reverse(f.x, f.y, name='Calais')
        self.test('reverse', fstr(r[:6], prec=6), x[0])

        self.testCopy(P)

        r = P.reverse(-38e3, 230e3, name='Calais')
        self.test('reverse', fstr(r[:6], prec=6), x[1])
        f = P.forward(r.lat, r.lon, name='Calais')
        self.test('forward', fstr(f[:6], prec=6), x[1])

        for m in (ellipsoidalKarney, ellipsoidalNvector, ellipsoidalVincenty):
            r = P.reverse(-38e3, 230e3, LatLon=m.LatLon)
            self.test('reverse', repr(r), x[2])

        G = Azimuthal(51.4934, 0.0098, name='Greenwich')
        self.test(repr(G), G, G)
        f = G.forward(P.lat0, P.lon0, P.name)
        self.test('forward', fstr(f[:6], prec=6), x[3])
        r = G.reverse(f.x, f.y)
        self.test('reverse', fstr(r[:6], prec=6), x[3])

        h = hypot(f.x, f.y)  # easting + norting ~= distance
        d = haversine(G.lat0, G.lon0, r.lat, r.lon)
        self.test('hypot', h, d, fmt='%.3f', known=abs(d - h) < 1000)

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

    from pygeodesy import azimuthal, equidistant, named

    t = Tests(__file__, __version__, azimuthal)

    t.testAzimuthal(Equidistant,
                   '-37467.812512, 230294.518853, 50.9, 1.8, 350.759218, 1.000223',
                   '-38000.0, 230000.0, 50.897321, 1.792455, 350.61849, 1.000222',
                   'LatLon(50°53′50.36″N, 001°47′32.84″E)',
                   '170420.92566, -293667.828613, 48.833333, 2.333333, 149.872606, 1.000472')
    if geographiclib:
        t.testAzimuthal(EquidistantKarney,
                       '-37526.978232, 230000.911579, 50.9, 1.8, 350.325442, 0.999778',
                       '-38000.0, 230000.0, 50.899962, 1.793278, 350.205524, 0.999778',
                       'LatLon(50°53′59.86″N, 001°47′35.8″E)',
                       '170617.186469, -293210.754313, 48.833333, 2.333333, 151.589952, 0.999529')
    else:
        t.skip('no geographiclib', n=14)

    As = tuple(Azimuthal(0, 0, datum=1)  # spherical datum of radius=1 for Snyder's Table ...
                            # 30, pp 196-197      26, pp 168           28, pp 188-189      22, pp 151          24, pp 158-159
           for Azimuthal in  (Equidistant,        Gnomonic,            LambertEqualArea,   Orthographic,       Stereographic))
    t.testSnyder(10, 80, As, (1.37704, 0.24656), (5.67128,  1.01543), (1.26747, 0.22694), (0.96985, 0.17365), (1.65643, 0.29658))
    t.testSnyder(20, 20, As, (0.33454, 0.35601), (0.36397,  0.38733), (0.33123, 0.35248), (0.32139, 0.34202), (0.34136, 0.36327))
    t.testSnyder(40, 40, As, (0.57386, 0.74912), (0.8391,   1.09537), (0.55281, 0.72164), (0.4924,  0.64279), (0.62062, 0.81016))
    t.testSnyder(60, 60, As, (0.58948, 1.17896), (1.73205,  3.4641),  (0.54772, 1.09545), (0.43301, 0.86603), (0.69282, 1.38564))
    t.testSnyder(70, 80, As, (0.50997, 1.42273), (5.67128, 15.82209), (0.46280, 1.29114), (0.33682, 0.93969), (0.63588, 1.77402))
    t.testSnyder(80, 80, As, (0.26358, 1.51792), (5.67128, 32.65961), (0.23828, 1.37219), (0.17101, 0.98481), (0.33201, 1.91196))  # XXX 1.96962?
    t.testSnyder(80, 10, As, (0.04281, 1.39829), (0.17633,  5.75877), (0.03941, 1.28702), (0.03015, 0.98481), (0.05150, 1.68198))

    A = equidistant(0, 0, datum=1, name='coverage')
    t.test('equatoradius', A.equatoradius, 1.0)
    t.test('flattening', A.flattening, 0)
    t.test('latlon0', A.latlon0, (0.0, 0.0))
    A.latlon0 = named.LatLon2Tuple(1, 2)
    t.test('latlon0', A.latlon0, (1.0, 2.0))
    t.test('name', A.name, 'coverage')
    t.test('radius', A.radius, 1.0)

    t.results()
    t.exit()
