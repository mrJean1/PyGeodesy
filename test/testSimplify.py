
# -*- coding: utf-8 -*-

# Test the simplify functions.

__all__ = ('Tests',)
__version__ = '17.06.23'

from base import TestsBase, secs2str

from pygeodesy import simplify1, simplify2, \
                      simplifyRDP, simplifyRDPm, \
                      simplifyVW, simplifyVWm, \
                      simplify  # module

from time import time

_Simplifys = ()  # simplifyXYZ functions to run


class Tests(TestsBase):

    def test2(self, function, points, ms, **kwds):

        if _Simplifys and function.__name__[8:] not in _Simplifys:
            return  # skip this simplify function

        n = len(points)
        t = ', '.join('%s=%s' % t for t in sorted(kwds.items()))
        f = '%s(%s, %s)' % (function.__name__, n, t)

        for m in reversed(sorted(ms.keys())):
            s = time()
            r = function(points, m, **kwds)
            n = len(r)
            s = time() - s
            t = '%s %dm (%s)' % (f, m, secs2str(s))
            self.test(t, n, ms[m])

        self.printf('')


if __name__ == '__main__':  # PYCHOK internal error?

    # usage: python testSimplify [[1-9] [RDP RDPm VW VWm ...]]

    import sys
    from testRoutes import Pts, PtsFFI  # RdpFFI

    # simplifyXYZ functions to run, all otherwise
    _Simplifys = sys.argv[2:]
    # number of meter values for each test
    m = 1 if len(sys.argv) < 2 else int(sys.argv[1])

    def _ms(ms):  # reduce the number of tests
        return dict(t for t in list(sorted(ms.items()))[:m])


    t = Tests(__file__, __version__, simplify)  # PYCHOK expected

    t.test2(simplify1, Pts, _ms({320: 4423, 160: 6638, 80: 9362, 40: 12079, 20: 14245, 10: 15621, 1: 16597}), adjust=True)

    t.test2(simplify2, Pts, _ms({320: 2327, 160: 3565, 80: 4895, 40: 6130, 20: 7022, 10: 7598, 1: 8239}), adjust=True, shortest=False)
    t.test2(simplify2, Pts, _ms({320: 2440, 160: 3753, 80: 5116, 40: 6347, 20: 7188, 10: 7709, 1: 8247}), adjust=True, shortest=True)

    t.test2(simplifyVWm, Pts, _ms({320: 2575, 160: 4762, 80: 7924, 40: 11334, 20: 13968, 10: 15452, 1: 16488}), adjust=True)
    t.test2(simplifyVWm, Pts, _ms({320: 3225, 160: 5787, 80: 9139, 40: 12349, 20: 14559, 10: 15781, 1: 16490}), adjust=False)

    t.test2(simplifyRDPm, Pts, _ms({320: 2512, 160: 4106, 80: 6150, 40: 8620, 20: 11138, 10: 13239, 1: 16196}), adjust=True, shortest=False)
    t.test2(simplifyRDPm, Pts, _ms({320: 2526, 160: 4127, 80: 6179, 40: 8654, 20: 11174, 10: 13266, 1: 16201}), adjust=True, shortest=True)

    # cut number of points (to shorten run time)
    n = len(Pts) // 10
    Ptsn = Pts[:n]

    t.test2(simplifyVW, Ptsn, _ms({320: 447, 160: 728, 80: 1032, 40: 1290, 20: 1478, 10: 1575, 1: 1657}), adjust=True)
    t.test2(simplifyVW, Ptsn, _ms({320: 518, 160: 837, 80: 1127, 40: 1353, 20: 1510, 10: 1591, 1: 1657}), adjust=False)

    t.test2(simplifyRDP, Ptsn, _ms({320: 1605, 160: 1616, 80: 1630, 40: 1638, 20: 1647, 10: 1654, 1: 1660}), adjust=True, shortest=False)
    t.test2(simplifyRDP, Ptsn, _ms({320: 1605, 160: 1616, 80: 1631, 40: 1639, 20: 1649, 10: 1655, 1: 1661}), adjust=True, shortest=True)

    # different points
    t.test2(simplifyVW,  PtsFFI, _ms({1678:  3, 1000:  5, 100: 23, 10: 65, 1: 69}), adjust=False)
    t.test2(simplifyRDP, PtsFFI, _ms({1678: 11, 1000: 31, 100: 61, 10: 67, 1: 68}), adjust=False, shortest=False)  # XXX len(RdpFFI) = 7

    # <http://georust.github.io/rust-geo/geo/algorithm/simplify/trait.Simplify.html>
#   t.test2(simplifyRDP, [_LatLon(*ll) for ll in ((0.0, 0.0), (5.0, 4.0), (11.0, 5.5), (17.3, 3.2), (27.8, 0.1))],
#                         _ms({1: 4}), adjust=False, shortest=True)  # (0.0, 0.0), (5.0, 4.0), (11.0, 5.5), (27.8, 0.1)

    # <http://georust.github.io/rust-geo/geo/algorithm/simplifyvw/trait.SimplifyVW.html>
#   t.test2(simplifyVW, [_LatLon(*ll) for ll in ((5.0, 2.0), (3.0, 8.0), (6.0, 20.0), (7.0, 25.0), (10.0, 10.0))],
#                        _ms({30: 3}), adjust=False)  # (5.0, 2.0), (7.0, 25.0), (10.0, 10.0)
    t.results(nl=0)
    t.exit()

#   Compare Vw routes
#   for i, p in enumerate(simplifyVW(Pts, 890, attr='vw2')):
#       print('%3d: %10.6f, %10.6f, %.15f' % (i, p.lon, p.lat, p.vw2))
