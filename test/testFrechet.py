
# -*- coding: utf-8 -*-

# Test the Frechet distances.

__all__ = ('Tests',)
__version__ = '19.08.19'

from base import isPython3, isWindows, TestsBase

from pygeodesy import fStr, LatLon_, randomrangenerator

_rr = randomrangenerator('R')
_ms = [LatLon_(*_ll) for _ll in zip(_rr( -90,  90, 2),  # 90
                                    _rr(-180, 180, 4))]

_ps = [LatLon_(*_ll) for _ll in zip(_rr( -89,  90, 3),  # 60
                                    _rr(-179, 180, 6))]


class Tests(TestsBase):

    def test2(self, Frechet, x, y, **kwds):

        def _tstr(t):
            s = list(t[:5])
            s[0] = fStr(t.fd, prec=5)
            return '(%s)' % (', '.join(map(str, s)),)

        f = Frechet(_ms, fraction=None, **kwds)

        t = _tstr(f.discrete(_ps))
        self.test(f.named, t, x)  # + (t.units,)

        t = _tstr(f.discrete(_ps, fraction=0.5))
        self.test(f.named, t, y)  # + (t.units,)


if __name__ == '__main__':

    from pygeodesy import FrechetDegrees, FrechetRadians, \
                          FrechetEquirectangular, FrechetEuclidean, \
                          FrechetHaversine, FrechetVincentys

    class FrechetDegrees_(FrechetDegrees):
        '''Custom Frechet.'''
        def distance(self, p1, p2):
            dy, dx = abs(p1.lat - p2.lat), abs(p1.lon - p2.lon)
            if dx < dy:
                dx, dy = dy, dx
            return dx + dy * 0.5

    class FrechetRadians_(FrechetRadians):
        '''Custom Frechet.'''
        def distance(self, p1, p2):
            dy, dx = abs(p1.a - p2.a), abs(p1.b - p2.b)
            if dx < dy:
                dx, dy = dy, dx
            return dx + dy * 0.5

    t = Tests(__file__, __version__)

    if isPython3:  # XXX different Random?
        t.test2(FrechetDegrees_, (178.5, 74, 56,   19,  5400),
                                 (175.5, 74, 52.5, 29, 10710))

        t.test2(FrechetRadians_, (3.11541, 74, 56,   19,  5400),
                                 (3.06305, 74, 52.5, 29, 10710))

        t.test2(FrechetEquirectangular, (7.1331,  8, 3,   138,  5400),
                                        (7.01295, 0, 0.0, 208, 10710))

        t.test2(FrechetEuclidean, (2.84717, 8, 3,   138,  5400),
                                  (2.76523, 0, 0.0, 208, 10710))

        t.test2(FrechetHaversine, (2.63867, 0, 0,   149,  5400),
                                  (2.63867, 0, 0.0, 208, 10710))

        t.test2(FrechetVincentys, (2.63867, 0, 0,   149,  5400),
                                  (2.63867, 0, 0.0, 208, 10710))

    elif isWindows:  # Python 2
        t.test2(FrechetDegrees_, (288.0, 1, 1,   147,  5400),
                                 (288.0, 1, 1.0, 205, 10710))

        t.test2(FrechetRadians_, (5.02655, 1, 1,   147,  5400),
                                 (5.02655, 1, 1.0, 205, 10710))

        t.test2(FrechetEquirectangular, ( 7.53702, 1, 3,   145,  5400),
                                        (12.58507, 0, 2.5, 203, 10710))

        t.test2(FrechetEuclidean, (2.81941, 1, 3,   145,  5400),
                                  (3.95734, 0, 2.5, 203, 10710))

        t.test2(FrechetHaversine, (1.81341, 18, 14,   117,  5400),
                                  (1.83289,  3,  4.5, 196, 10710))

        t.test2(FrechetVincentys, (1.81341, 18, 14,   117,  5400),
                                  (1.83289,  3,  4.5, 196, 10710))

    else:  # Python 2, elsewhere
        t.test2(FrechetDegrees_, (288.0, 1, 1,   147,  5400),
                                 (288.0, 1, 1.0, 205, 10710))

        t.test2(FrechetRadians_, (5.02655, 1, 1,   147,  5400),
                                 (5.02655, 1, 1.0, 205, 10710))

        t.test2(FrechetEquirectangular, ( 7.53702, 1, 3,   145,  5400),
                                        (12.58507, 0, 2.5, 203, 10710))

        t.test2(FrechetEuclidean, (2.81941, 1, 3,   145,  5400),
                                  (3.95734, 0, 2.5, 203, 10710))

        t.test2(FrechetHaversine, (1.81341, 18, 14,   117,  5400),
                                  (1.83289,  3,  4.5, 196, 10710))

        t.test2(FrechetVincentys, (1.81341, 18, 14,   117,  5400),
                                  (1.83289,  3,  4.5, 196, 10710))
    t.results()
    t.exit()
