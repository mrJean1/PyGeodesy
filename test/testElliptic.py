
# -*- coding: utf-8 -*-

u'''Test Elliptic Python implementation.
'''

__all__ = ('Tests',)
__version__ = '20.01.17'

from base import TestsBase

from pygeodesy import elliptic, EPS, fStr, PI_2, radians, sincos2


class Tests(TestsBase):

    def testElliptic(self):

        _RC = elliptic._RC
        _RD = elliptic._RD
        _RF = elliptic._RF
        _RG = elliptic._RG
        _RJ = elliptic._RJ

        eps = EPS * 4
        self.test('eps', eps, eps, fmt='%.12e')

        y = 0.1
        for x in range(2, 100):
            x *= 0.01
            rc = _RC(x, y)
            rf = _RF(x, y, y)
            k = abs(rc - rf) < eps
            t = '_RC, _RF(%.3f, ...)' % (x,)
            self.test(t, rc, rf, fmt='%.12f', known=k)
            y = x

        for kp2 in range(1, 100):
            kp2 *= 0.01
            rd = _RD(0, kp2, 1)
            rj = _RJ(0, kp2, 1, 1)
            k = abs(rd - rj) < eps
            t = '_RD, _RJ(%.3f, ...)' % (kp2,)
            self.test(t, rd, rj, fmt='%.12f', known=k)

        self.test('eps', eps, eps, fmt='%.12e')

        # <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1EllipticFunction.html>
        e = elliptic.Elliptic(0.1)
        self.test('cK', e.cK, '1.612441349', fmt='%.9f')
        self.test('cE', e.cE, '1.530757637', fmt='%.9f')
        phi = radians(20)
        sn, cn = sincos2(phi)
        self.test('fE(phi)', e.fE(phi), '0.348372822', fmt='%.9f')
        dn = e.fDelta(sn, cn)
        self.test('fDelta', dn, '0.994133906', fmt='%.9f')
        self.test('fE(sn, cn, Delta)', e.fE(sn, cn, dn), '0.348372822', fmt='%.9f')

        self.test('fEd(PI_2)',       e.fEd(PI_2),       '0.027415224', fmt='%.9f')
        self.test('fEinv(PI_2)',     e.fEinv(PI_2),     '1.612999420', fmt='%.9f')
        self.test('fF(sn, cn, dn)',  e.fF(sn, cn, dn),  '0.349761345', fmt='%.9f')
        self.test('fG(sn, cn, dn)',  e.fG(sn, cn, dn),  '0.348372822', fmt='%.9f')
        self.test('fH(sn, cn, dn)',  e.fH(sn, cn, dn),  '0.363646580', fmt='%.9f')
        self.test('fPi(sn, cn, dn)', e.fPi(sn, cn, dn), '0.349761345', fmt='%.9f')

        self.test('_RF(1, 2, 0)', _RF(1, 2, 0), '1.311028777', fmt='%.9f')
        self.test('_RF(2, 3, 4)', _RF(2, 3, 4), '0.584082842', fmt='%.9f')

        self.test('_RC(0, 1/4)', _RC(0, 0.25), '3.141592654', fmt='%.9f')
        self.test('_RC(9/4, 2)', _RC(2.25, 2), '0.693147181', fmt='%.9f')
        self.test('_RC(1/4,-2)', _RC(0.25,-2), '0.231049060', fmt='%.9f')

        self.test('_RJ(0, 1, 2, 3)', _RJ(0, 1, 2, 3), '0.776886238', fmt='%.9f')
        self.test('_RJ(2, 3, 4, 5)', _RJ(2, 3, 4, 5), '0.142975797', fmt='%.9f')

        self.test('_RD(0, 2, 1)', _RD(0, 2, 1), '1.797210352', fmt='%.9f')
        self.test('_RD(2, 3, 4)', _RD(2, 3, 4), '0.165105273', fmt='%.9f')

        self.test('_RG(0, 16, 16)',     _RG(0, 16, 16), '3.141592654', fmt='%.9f')
        self.test('_RG(2,  3,  4)',     _RG(2,  3,  4), '1.725503028', fmt='%.9f')
        self.test('_RG(0,  0.0796, 4)', _RG(0,  0.0796, 4), '1.028475809', fmt='%.9f')

        e.reset(0, 0)
        self.test('sncndn(x)', fStr(e.sncndn(0), prec=9), '0.0, 1.0, 1.0')
        self.test('sncndn(x)', fStr(e.sncndn(PI_2), prec=9), '1.0, -0.0, 1.0', known=True)
        e.reset(1, 1)
        self.test('sncndn(x)', fStr(e.sncndn(0), prec=9), '0.0, 1.0, 1.0')
        self.test('sncndn(x)', fStr(e.sncndn(PI_2), prec=9), '0.917152336, 0.398536815, 0.398536815')

        self.testCopy(e)


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testElliptic()
    t.results()
    t.exit()
