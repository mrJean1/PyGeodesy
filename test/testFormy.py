
# -*- coding: utf-8 -*-

# Test formulary functions.

__all__ = ('Tests',)
__version__ = '20.03.23'

from base import TestsBase

from pygeodesy import R_M, antipode, bearing, equirectangular, euclidean, \
                      haversine, heightOf, horizon, isantipode, isantipode_, \
                      map1, vincentys

from math import radians


class Tests(TestsBase):

    def testDistances(self, a1, b1, a2, b2, x):
        k = x * 0.003  # allow 0.3% margin
        d = haversine(a1, b1, a2, b2)
        self.test('haversine', d, x, fmt='%.3f', known=abs(d - x) < k)

        d = vincentys(a1, b1, a2, b2)
        self.test('vincentys', d, x, fmt='%.3f', known=abs(d - x) < k)

        k = x * 0.02  # allow 2% margin
        d = equirectangular(a1, b1, a2, b2, limit=90)
        self.test('equirectangular', d, x, fmt='%.3f', known=abs(d - x) < k)

        k = x * 0.11  # allow 11% margin
        d = euclidean(a1, b1, a2, b2)
        self.test('euclidean', d, x, fmt='%.3f', known=abs(d - x) < k)

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
        self.test('heightOf90',  heightOf(90,  R_M), R_M)
        self.test('heightOf135', heightOf(135, R_M), 5401080.43931, fmt='%.5f')

        self.test('horizon0',     horizon(0), 0.0)
        self.test('horizon10Km',  horizon(10000), '357099.672', fmt='%.3f')
        self.test('horizon30Kft', horizon(10000, refraction=True), '392310.704', fmt='%.3f')
        self.test('horizon10Kft', horizon( 3000, refraction=True), '214877.422', fmt='%.3f')

        Boston    = 42.3541165, -71.0693514
        Cleveland = 41.499498, -81.695391
        MtDiablo  = 37.8816, -121.9142
        Newport   = 41.49008, -71.312796
        NewYork   = 40.7791472, -73.9680804
        # <https://GeographicLib.SourceForge.io/cgi-bin/GeodSolve>
        for ll1, ll2, x in ((Boston, NewYork,    298396.057),
                            (Boston, Newport,     98071.559),
                            (Cleveland, NewYork, 653456.173),
                            (NewYork, MtDiablo, 4094953.628)):
            abx = ll1 + ll2 + (x,)
            self.testDistances(*abx)


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testFormy()
    t.results()
    t.exit()
