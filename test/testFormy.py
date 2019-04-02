
# -*- coding: utf-8 -*-

# Test formulary functions.

__all__ = ('Tests',)
__version__ = '19.04.02'

from base import TestsBase

from pygeodesy import equirectangular, euclidean, haversine, vincentys


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
        Boston    = 42.3541165, -71.0693514
        Cleveland = 41.499498, -81.695391
        MtDiablo  = 37.8816, -121.9142
        Newport   = 41.49008, -71.312796
        NewYork   = 40.7791472, -73.9680804
        # <http://GeographicLib.SourceForge.io/cgi-bin/GeodSolve>
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
