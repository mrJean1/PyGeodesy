
# -*- coding: utf-8 -*-

# Test the simplify functions.

__all__ = ('Tests',)
__version__ = '17.04.06'

from tests import Tests as _Tests

from geodesy import simplify1, simplify2, \
                    simplifyRDP, simplifyRDPm, \
                    simplifyVW, simplifyVWm

from time import time


class Tests(_Tests):

    def test2(self, simplify, points, ms, **kwds):

        n = len(points)
        t = ', '.join('%s=%s' % t for t in sorted(kwds.items()))
        s = '%s(%s, %s)' % (simplify.__name__, n, t)

        for m in reversed(sorted(ms.keys())):
            t = time()
            n = len(simplify(points, m, **kwds))
            t = time() - t
            t = '%s %dm (%.3f sec)' % (s, m, t)
            self.test(t, n, str(ms[m]))

        self.printf('')


if __name__ == '__main__':  # PYCHOK internal error?

    import sys
    from testRoutes import Pts

    # number of meter values for each test
    m = 1 if len(sys.argv) < 2 else int(sys.argv[1])

    def _ms(ms):  # reduce the number of tests
        return dict(t for t in list(sorted(ms.items()))[:m])

    t = Tests(__file__, __version__)

    t.test2(simplify1, Pts, _ms({320: 4423, 160: 6638, 80: 9362, 40: 12079, 20: 14245, 10: 15621, 1: 16597}), adjust=True)

    t.test2(simplify2, Pts, _ms({320: 2327, 160: 3565, 80: 4895, 40: 6130, 20: 7022, 10: 7598, 1: 8239}), adjust=True, shortest=False)
    t.test2(simplify2, Pts, _ms({320: 2206, 160: 3272, 80: 4581, 40: 5838, 20: 6850, 10: 7548, 1: 8247}), adjust=True, shortest=True)

    t.test2(simplifyVWm, Pts, _ms({320: 2005, 160: 3833, 80: 6596, 40: 9931, 20: 12904, 10: 14868, 1: 16482}), adjust=True)

    t.test2(simplifyRDPm, Pts, _ms({320: 2512, 160: 4106, 80: 6150, 40: 8620, 20: 11138, 10: 13239, 1: 16196}), adjust=True, shortest=False)
    t.test2(simplifyRDPm, Pts, _ms({320: 2526, 160: 4127, 80: 6179, 40: 8654, 20: 11174, 10: 13266, 1: 16201}), adjust=True, shortest=True)

    # cut number of points (to shorten run time)
    n = len(Pts) // 10
    Ptsn = Pts[:n]

    t.test2(simplifyVW, Ptsn, _ms({320: 347, 160: 590, 80: 883, 40: 1172, 20: 1392, 10: 1528, 1: 1657}), adjust=True)

    t.test2(simplifyRDP, Ptsn, _ms({320: 1605, 160: 1616, 80: 1630, 40: 1638, 20: 1647, 10: 1654, 1: 1660}), adjust=True, shortest=False)
    t.test2(simplifyRDP, Ptsn, _ms({320: 1605, 160: 1616, 80: 1631, 40: 1639, 20: 1649, 10: 1655, 1: 1661}), adjust=True, shortest=True)

    t.results()
    t.exit()

#   Compare Vw routes
#   for i, p in enumerate(simplifyVW(Pts, 890, attr='vw2')):
#       print('%3d: %10.6f, %10.6f, %.15f' % (i, p.lon, p.lat, p.vw2))
