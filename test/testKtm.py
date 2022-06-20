# -*- coding: utf-8 -*-

u'''Some basic C{ktm} L{KTransverseMercator} tests.
'''
__all__ = ('Tests',)
__version__ = '22.06.15'

from base import TestsBase

from pygeodesy import fstr, hypot, NN, sincos2d
from pygeodesy.interns import _SLASH_


class Tests(TestsBase):

    def testKtm(self, xTM, known, fw=NN, rv=NN):

        t = xTM.toRepr()
        self.test('_TM', t, t, nl=1)
        e = 0
        for i in range(0, 361, 3):
            lon, lat = sincos2d(i)
            lat *= 80
            lon *= 160
            try:
                w, x, y = 'forward ', lat, lon
                n = '%s%dN (%.3f, %.3f) ' % (w, i, x, y)
                x, y, g, k = xTM.forward(lat, lon)
                z = getattr(xTM, '_zetaC', '-')
                z = _SLASH_(z, xTM.iteration or '-')
                t = fstr((x, y, g, k), prec=6)
                self.test(n + str(z), t, t)

                w = 'reverse '
                x, y, _, _ = xTM.reverse(x, y)
                s = getattr(xTM, '_sigmaC', '-')
                s = _SLASH_(s, xTM.iteration or '-')
                d = hypot(lat - x, lon - y)
                n = '%s%dN (%.3f, %.3f) %s' % (w, i, x, y, s)
                self.test(n, d, min(d, known), known=True)
                e = max(e, d)
            except Exception as X:
                t = repr(X)
                n = '%s%dN (%.3f, %.3f) ERROR' % (w, i, x, y)
                self.test(n, t, known, known=True)
        self.test('max', e, known, known=True)

        t = xTM.forward(30, 60)
        self.test('forward', fstr(t, prec=6), fw, nl=1)
        t = xTM.reverse(62e5, 55e5)
        self.test('reverse', fstr(t, prec=6), rv)


if __name__ == '__main__':

    from pygeodesy import Ellipsoids, ExactTransverseMercator as _ETM, \
                          ktm, KTransverseMercator as _KTM
    from sys import argv

    _t = len(argv) > 1 and argv[1] == '-t'

    t = Tests(__file__, __version__, ktm)

    # echo 30 60 | TransverseMercatorProj -s
    # 6208422.537400 5452954.287187 41.077484300765 1.511911171199
    # echo 62e5 55e5 | TransverseMercatorProj -s -r
    # 30.24422805624 60.16966413219 41.479185132982 1.510345924158
    t.testKtm(_KTM(Ellipsoids.WGS84, TMorder=7),  6e-14,
               fw='6208422.5374, 5452954.287187, 41.077484, 1.511911',
               rv='30.244228, 60.169664, 41.479185, 1.510346')

    # echo 30 60 | TransverseMercatorProj -s -e 6371008.771415 0
    # 6196225.831883 5458228.732328 40.893394649131 1.511253148880
    # echo 62e5 55e5 | TransverseMercatorProj -s -r -e 6371008.771415 0
    # 30.17255567482 60.20831594572 41.279729402162 1.511925129658
    t.testKtm(_KTM(Ellipsoids.Sphere, TMorder=4), 5e-14,
               fw='6196225.831883, 5458228.732328, 40.893395, 1.511253',
               rv='30.172556, 60.208316, 41.279729, 1.511925')

    # echo 30 60 | TransverseMercatorProj
    # 6208422.537400 5452954.287187 41.077484300760 1.511911171199
    # echo 62e5 55e5 | TransverseMercatorProj -r
    # 30.24422805624 60.16966413219 41.479185132982 1.510345924158
    t.testKtm(_ETM(extendp=_t), 9e-14,
               fw='6208422.5374, 5452954.287187, 41.077484, 1.511911',
               rv='30.244228, 60.169664, 41.479185, 1.510346')

    t.results()
    t.exit()
