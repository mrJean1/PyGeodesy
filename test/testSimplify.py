
# -*- coding: utf-8 -*-

# Test the simplify functions.

__all__ = ('Tests',)
__version__ = '20.04.06'

from base import numpy, TestsBase, secs2str

from pygeodesy import EPS, R_M, LatLon_, Numpy2LatLon, pairs, \
                      simplify1, simplifyRW, \
                      simplifyRDP, simplifyRDPm, \
                      simplifyVW, simplifyVWm, \
                      simplify

from math import cos, degrees, radians
from time import time

_Simplifys = ()  # simplifyXYZ functions to run


class Tests(TestsBase):

    def test2(self, function, points, ms, typ=None, **kwds):

        if _Simplifys and function.__name__[8:] not in _Simplifys:
            return  # skip this simplify function

        n = len(points)
        t = ', '.join(pairs(kwds.items()))
        f = '%s(%s, %s)' % (function.__name__, n, t)

        S = s = 0
        for m in reversed(sorted(ms.keys())):
            s = time()
            r = function(points, m, **kwds)
            n = len(r)
            s = time() - s
            t = '%s %dm (%s)' % (f, m, secs2str(s))
            self.test(t, n, ms[m])
            S += s
            if typ:
                self.test(f + ' result', type(r), typ)

        if S > s:  # sub-total time
            self.test__('%s %s', f, secs2str(S), nt=1)


# for comparison, following are 2 other RDP implementations,
# modified to match signatures of simplifyRDP and -RDPm.

def simplifyRDPfw(points, epsilon, radius=R_M, adjust=False, shortest=False,  # MCCABE 20
                                   modified=False, indices=False):  # PYCHOK expected
    '''Iterative Ramer-Douglas-Peucker algorithm.

       <https://GitHub.com/FlorianWilhelm/gps_data_with_python>

       points[] -- Input coordinates as LatLon's in degrees
       epsilon -- distance tolerance in meter
       radius -- mean earth radius in meter
       adjust -- adjust lon's delta by cos(average lat's)
       shortest -- use shortest or perpendicular distance
       modified -- stop search at the first deviant point
       indices -- return points indices inlieu of points
    '''
    def _d2xy(p1, p2):
        # get deltas and hypot squared
        dx = p2.lon - p1.lon
        if dx > 180:
            dx -= 360
        elif dx < -180:
            dx += 360
        if adjust:
            dx *= cos(radians((p1.lat + p2.lat) * 0.5))
        dy = p2.lat - p1.lat
        d2 = dx**2 + dy**2
        return d2, dx, dy  # PYCHOK returns

    if shortest:
        raise NotImplementedError('shortest=%s' % (shortest,))

    if len(points) < 3:
        return points

    # convert epsilon meters to degrees
    ed2 = degrees(epsilon / radius)**2  # use squared distances

    use = [True] * len(points)
    stk = [(0, len(points) - 1)]

    while stk:
        s, e = stk.pop()
        if e > (s + 1):
            t2, t = ed2, 0

            p2, px, py = _d2xy(points[s], points[e])
            if p2 > EPS:
                p2 = 1.0 / p2
            else:  # distance to point
                p2 = 0

            for i in range(s + 1, e):
                if use[i]:
                    d2, dx, dy = _d2xy(points[i], points[s])
                    if p2 and d2 > EPS:  # distance to line
                        d2 = p2 * (px * dy - py * dx)**2
                    if d2 > t2:
                        t2, t = d2, i
                        if modified:
                            break

            if t2 > ed2:  # and t > 0:
                stk.append((t, e))
                if modified:
                    for i in range(s + 1, t):
                        use[i] = False
                else:
                    stk.append((s, t))
            else:
                for i in range(s + 1, e):
                    use[i] = False

    if indices:
        return [i for i, b in enumerate(use) if b]
    else:
        return [points[i] for i, b in enumerate(use) if b]


def simplifyRDPgr(source, kink, radius=R_M, adjust=True, shortest=True,  # MCCABE 17
                                modified=False, indices=False):  # PYCHOK expected
    '''Stack-based Douglas Peucker line simplification.

       Transcribed from JavaScript original after code by U{Dr. Gary J.
       Robinson<https://www.BDCC.co.UK/Gmaps/GDouglasPeuker.js>},
       Environmental Systems Science Centre, University of Reading,
       Reading, UK.

       source[] -- Input coordinates as LatLon's in degrees
       kink -- distance in metres, kinks above this depth are kept
       radius -- mean earth radius in metres
       adjust -- adjust lon's delta by cos(average lat's)
       shortest -- use shortest or perpendicular distance
       modified -- stop search at the first deviant point
       indices -- return source indices inlieu of points

       The kink depth is the height of the triangle abc where
       a-b and b-c are two consecutive line segments.
    '''

    def _d2xy(p1, p2):
        # get deltas and hypot squared
        dx = p2.lon - p1.lon
        if abs(dx) > 180:
            dx = 360 - abs(dx)
        if adjust:
            dx *= cos(radians((p1.lat + p2.lat) * 0.5))
        dy = p2.lat - p1.lat
        d2 = dx**2 + dy**2
        return d2, dx, dy  # PYCHOK returns

    if not shortest:
        raise NotImplementedError('shortest=%s' % (shortest,))

    n1 = len(source) - 1
    if n1 < 2:
        return source

    k2 = degrees(kink / radius) ** 2  # kink in degrees squared

    ixs = {0: True, n1: True}
    stk = [(0, n1)]

    while stk:
        s, e = stk.pop()
        if e > (s + 1):
            d12, x12, y12 = _d2xy(source[s], source[e])

            max2, ix = k2, s
            for i in range(s + 1, e):
                d13, x13, y13 = _d2xy(source[s], source[i])
                d23,   _,   _ = _d2xy(source[e], source[i])

                if d13 >= (d12 + d23):
                    d2 = d23
                elif d23 >= (d12 + d13):
                    d2 = d13
                elif d12 > EPS:  # perpendicular distance squared
                    d2 = (x13 * y12 - y13 * x12)**2 / d12
                else:  # null edge
                    d2 = (d13 + d23) * 0.5

                if d2 > max2:
                    max2, ix = d2, i
                    if modified:
                        break

            if max2 > k2:
                stk.append((ix, e))
                if not modified:
                    stk.append((s, ix))
                    s = 0
        ixs[s] = True

    if indices:
        return [i for i in sorted(ixs.keys())]
    else:
        return [source[i] for i in sorted(ixs.keys())]


if __name__ == '__main__':  # PYCHOK internal error?

    # usage: python testSimplify [[1-9] [RDP RDPm VW VWm ...]]

    import sys
    from testRoutes import Pts, PtsFFI  # PtsJS, RdpFFI

    # simplifyXYZ functions to run, all otherwise
    _Simplifys = sys.argv[2:]
    # number of meter values for each test
    m = 1 if len(sys.argv) < 2 else int(sys.argv[1])

    def _ms(ms):  # reduce the number of tests
        return dict(t for t in list(sorted(ms.items()))[:m])

    t = Tests(__file__, __version__, simplify)

    t.test2(simplify1, Pts, _ms({160: 6638, 80: 9362, 40: 12079, 20: 14245, 10: 15621, 1: 16597}), adjust=True)

    t.test2(simplifyRW, Pts, _ms({160: 1179, 80: 1737, 40: 2322, 20: 3121, 10: 4041, 1: 7095}), adjust=True, shortest=False, indices=True)
    t.test2(simplifyRW, Pts, _ms({160: 1179, 80: 1737, 40: 2322, 20: 3121, 10: 4041, 1: 7095}), adjust=True, shortest=False)
    t.test2(simplifyRW, Pts, _ms({160: 4350, 80: 5646, 40: 6744, 20: 7535, 10: 7995, 1: 8302}), adjust=True, shortest=True)

    t.test2(simplifyVWm, Pts, _ms({160: 1425, 80: 2648, 40: 4701, 20: 7678, 10: 11166, 1: 16328}), adjust=True, indices=True)
    t.test2(simplifyVWm, Pts, _ms({160: 1425, 80: 2648, 40: 4701, 20: 7678, 10: 11166, 1: 16328}), adjust=True)
    t.test2(simplifyVWm, Pts, _ms({160: 1755, 80: 3253, 40: 5590, 20: 8911, 10: 12331, 1: 16373}), adjust=False)

    t.test2(simplifyRDPm,  Pts, _ms({160: 3099, 80: 4924, 40: 7189, 20: 9663, 10: 11876, 1: 15864}), adjust=True, shortest=False)
    t.test2(simplifyRDPm,  Pts, _ms({160: 3099, 80: 4926, 40: 7195, 20: 9670, 10: 11884, 1: 15867}), adjust=True, shortest=True)
    # for comparison, the RDPgr results should be less than a few points different
    t.test2(simplifyRDPgr, Pts, _ms({160: 3099, 80: 4925, 40: 7196, 20: 9670, 10: 11884, 1: 15867}), adjust=True, shortest=True, modified=True)

    t.test2(simplifyRDPm,  Pts, _ms({160: 3166, 80: 5002, 40: 7259, 20:  9720, 10: 11939, 1: 15869}), adjust=False, shortest=False)
    # for comparison, the RDPfw results should be less than a few points different
    t.test2(simplifyRDPfw, Pts, _ms({160: 3166, 80: 5002, 40: 7259, 20:  9720, 10: 11939, 1: 15869}), adjust=False, shortest=False, modified=True)

    # run time may be too long
    t.test2(simplifyRDP,   Pts, _ms({100: 1185, 10: 4209, 1: 10960}), adjust=True, shortest=True, indices=True)
    t.test2(simplifyRDP,   Pts, _ms({100: 1185, 10: 4209, 1: 10960}), adjust=True, shortest=True)
    # for comparison, the RDPgr results should be less than a few points different
    t.test2(simplifyRDPgr, Pts, _ms({100: 1185, 10: 4209, 1: 10960}), adjust=True, shortest=True, modified=False)

    t.test2(simplifyRDP,   Pts, _ms({100: 1251, 10: 4461, 1: 11248}), adjust=False, shortest=False, indices=True)
    t.test2(simplifyRDP,   Pts, _ms({100: 1251, 10: 4461, 1: 11248}), adjust=False, shortest=False)
    # for comparison, the RDPfw results should be less than a few points different
    t.test2(simplifyRDPfw, Pts, _ms({100: 1251, 10: 4461, 1: 11248}), adjust=False, shortest=False, modified=False)

    # cut number of points (to shorten run time)
    n = len(Pts) // 10
    Ptsn = Pts[:n]

    t.test2(simplifyVW, Ptsn, _ms({160: 277, 80: 430, 40: 694, 20:  981, 10: 1266, 1: 1641}), adjust=True)
    t.test2(simplifyVW, Ptsn, _ms({160: 314, 80: 495, 40: 793, 20: 1081, 10: 1351, 1: 1646}), adjust=False)

    t.test2(simplifyRDP,   Ptsn, _ms({160: 98, 80: 147, 40: 226, 20: 354, 10: 501, 1: 1231}), adjust=True, shortest=False)
    t.test2(simplifyRDP,   Ptsn, _ms({160: 98, 80: 147, 40: 226, 20: 354, 10: 501, 1: 1231}), adjust=True, shortest=True)
    t.test2(simplifyRDPgr, Ptsn, _ms({160: 98, 80: 147, 40: 226, 20: 354, 10: 501, 1: 1231}), adjust=True, shortest=True)

    t.test2(simplifyRDP,   Ptsn, _ms({160: 111, 80: 161, 40: 256, 20: 387, 10: 542, 1: 1267}), adjust=False, shortest=False)
    # for comparison, the RDPfw results should be less than a few points different
    t.test2(simplifyRDPfw, Ptsn, _ms({160: 111, 80: 161, 40: 256, 20: 387, 10: 542, 1: 1267}), adjust=False, shortest=False)
    # for comparison, the RDPgr results should be less than a few points different
    t.test2(simplifyRDPgr, Ptsn, _ms({160: 111, 80: 161, 40: 256, 20: 387, 10: 542, 1: 1267}), adjust=False, shortest=True)

    # different points
    t.test2(simplifyVW,    PtsFFI, _ms({1000: 3, 100: 12, 10: 48, 1: 69}), adjust=False)
    t.test2(simplifyRDP,   PtsFFI, _ms({1000: 3, 100:  7, 10: 18, 1: 50}), adjust=False, shortest=False)  # XXX len(RdpFFI) = 7
    # for comparison, the RDPgr results should be less than a few points different
    t.test2(simplifyRDPfw, PtsFFI, _ms({1000: 3, 100:  7, 10: 18, 1: 50}), adjust=False, shortest=False)
    # for comparison, the RDPgr results should be less than a few points different
    t.test2(simplifyRDPgr, PtsFFI, _ms({1000: 3, 100:  7, 10: 18, 1: 50}), adjust=False, shortest=True)
    t.test2(simplifyRDPgr, PtsFFI, _ms({1000: 2, 100:  7, 10: 15, 1: 45}), adjust=True, shortest=True)

    # <https://docs.RS/geo/0.8.3/geo/algorithm/simplify/trait.Simplify.html>
#   t.test2(simplifyRDP, [LatLon_(*ll) for ll in ((0.0, 0.0), (5.0, 4.0), (11.0, 5.5), (17.3, 3.2), (27.8, 0.1))],
#                        _ms({1: 5}), adjust=False, shortest=True)  # (0.0, 0.0), (5.0, 4.0), (11.0, 5.5), (27.8, 0.1)
    t.test2(simplifyRDP, [LatLon_(*ll) for ll in ((0.0, 0.0), (5.0, 4.0), (11.0, 5.5), (17.3, 3.2), (EPS, EPS))],
                         _ms({1: 5}), adjust=False, shortest=True)  # coverage d2ih

    # <https://docs.RS/geo/0.8.3/geo/algorithm/simplifyvw/trait.SimplifyVW.html>
#   t.test2(simplifyVW, [LatLon_(*ll) for ll in ((5.0, 2.0), (3.0, 8.0), (6.0, 20.0), (7.0, 25.0), (10.0, 10.0))],
#                       _ms({30: 3}), adjust=False)  # (5.0, 2.0), (7.0, 25.0), (10.0, 10.0)
    t.test2(simplifyVW, [LatLon_(*ll) for ll in ((5.0, 2.0), (3.0, 8.0), (6.0, 20.0), (7.0, 25.0), (10.0, 10.0))],
                        _ms({30: 5}), adjust=False, attr='name')  # coverage attr, 3 set to 5

    if numpy:
        if 'numpy' in _Simplifys:
            _Simplifys.remove('numpy')

        t.test('numpy.__version__', numpy.__version__, numpy.__version__)

        npy = numpy.array([(ll.lon, 0, ll.lat, 0) for ll in PtsFFI], dtype=float)
        pts = Numpy2LatLon(npy, ilat=2, ilon=0)

        t.test2(simplify1,     pts, _ms({1000: 5, 100: 25, 10: 67, 1: 69}), adjust=False, typ=type(npy))
        t.test2(simplifyRW,    pts, _ms({1000: 4, 100:  9, 10: 22, 1: 33}), adjust=False, typ=type(npy))
        t.test2(simplifyRDP,   pts, _ms({1000: 3, 100:  7, 10: 18, 1: 50}), adjust=False, typ=type(npy))
        t.test2(simplifyRDPm,  pts, _ms({1000: 3, 100: 16, 10: 48, 1: 67}), adjust=False, typ=type(npy))
        t.test2(simplifyRDPfw, pts, _ms({1000: 3, 100:  7, 10: 18, 1: 50}), adjust=False)
        t.test2(simplifyRDPgr, pts, _ms({1000: 3, 100:  7, 10: 18, 1: 50}), adjust=False)
        t.test2(simplifyVW,    pts, _ms({1000: 3, 100: 12, 10: 48, 1: 69}), adjust=False, typ=type(npy))
        t.test2(simplifyVWm,   pts, _ms({1000: 2, 100:  7, 10: 45, 1: 69}), adjust=False, typ=type(npy))

        t.test2(simplify1,     pts, _ms({1000: 4, 100: 23, 10: 65, 1: 69}), adjust=True, typ=type(npy))
        t.test2(simplifyRW,    pts, _ms({1000: 3, 100:  8, 10: 21, 1: 31}), adjust=True, typ=type(npy))
        t.test2(simplifyRDP,   pts, _ms({1000: 2, 100:  7, 10: 15, 1: 45}), adjust=True, typ=type(npy))
        t.test2(simplifyRDPm,  pts, _ms({1000: 2, 100: 13, 10: 46, 1: 67}), adjust=True, typ=type(npy))
        t.test2(simplifyRDPfw, pts, _ms({1000: 2, 100:  7, 10: 15, 1: 45}), adjust=True)
        t.test2(simplifyRDPgr, pts, _ms({1000: 2, 100:  7, 10: 15, 1: 45}), adjust=True)
        t.test2(simplifyVW,    pts, _ms({1000: 3, 100: 11, 10: 43, 1: 69}), adjust=True, typ=type(npy))
        t.test2(simplifyVWm,   pts, _ms({1000: 2, 100:  6, 10: 39, 1: 69}), adjust=True, typ=type(npy))

        t.test2(simplify1,     pts, _ms({1000: 5, 100: 25, 10: 67, 1: 69}), adjust=False, indices=True, typ=list)
        t.test2(simplifyRW,    pts, _ms({1000: 4, 100:  9, 10: 22, 1: 33}), adjust=False, indices=True, typ=list)
        t.test2(simplifyRDP,   pts, _ms({1000: 3, 100:  7, 10: 18, 1: 50}), adjust=False, indices=True, typ=list)
        t.test2(simplifyRDPm,  pts, _ms({1000: 3, 100: 16, 10: 48, 1: 67}), adjust=False, indices=True, typ=list)
        t.test2(simplifyRDPfw, pts, _ms({1000: 2, 100:  7, 10: 15, 1: 45}), adjust=True, indices=True)
        t.test2(simplifyRDPgr, pts, _ms({1000: 2, 100:  7, 10: 15, 1: 45}), adjust=True, indices=True)
        t.test2(simplifyVW,    pts, _ms({1000: 3, 100: 12, 10: 48, 1: 69}), adjust=False, indices=True, typ=list)
        t.test2(simplifyVWm,   pts, _ms({1000: 2, 100:  7, 10: 45, 1: 69}), adjust=False, indices=True, typ=list)

    else:
        t.test('no module', 'numpy', 'numpy')

    t.results()
    t.exit()

#   Compare Vw routes
#   for i, p in enumerate(simplifyVW(Pts, 890, attr='vw2')):
#       print('%3d: %10.6f, %10.6f, %.15f' % (i, p.lon, p.lat, p.vw2))
