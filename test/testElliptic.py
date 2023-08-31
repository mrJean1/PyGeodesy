
# -*- coding: utf-8 -*-

# Test L{elliptic} Python implementation.

__all__ = ('Tests',)
__version__ = '23.08.31'

from bases import endswith, TestsBase

from pygeodesy import elliptic, Elliptic, EllipticError, Elliptic3Tuple, \
                      EPS, fstr, PI_2, PI_4, radians, Scalar, sincos2


class Tests(TestsBase):

    def testElliptic(self):

        RC = Elliptic.fRC
        RD = Elliptic.fRD
        RF = Elliptic.fRF
        RG = Elliptic.fRG
        RJ = Elliptic.fRJ

        eps4 = EPS * 4
        self.test('eps4', eps4, eps4, fmt='%.9e')

        y = 0.1
        for x in range(2, 100):
            x *= 0.01
            rc = RC(x, y)
            rf = RF(x, y, y)
            k  = abs(rc - rf) < eps4
            t  = 'RC, RF(%.3f, ...)' % (x,)
            self.test(t, rc, rf, fmt='%.12f', known=k)
            y = x

        for kp2 in range(1, 100):
            kp2 *= 0.01
            rd = RD(0, kp2, 1)
            rj = RJ(0, kp2, 1, 1)
            k  = abs(rd - rj) < eps4
            t  = 'RD, RJ(%.3f, ...)' % (kp2,)
            self.test(t, rd, rj, fmt='%.12f', known=k)

        self.test('eps4', eps4, eps4, fmt='%.9e', nl=1)

        # <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1EllipticFunction.html>
        e = Elliptic(0.1)
        self.test('k2',  e.k2,  '0.1000000',   fmt='%.7f')
        self.test('kp2', e.kp2, '0.9000000',   fmt='%.7f')
        self.test('eps', e.eps, '0.0263340',   fmt='%.7f')
        self.test('cD',  e.cD,  '0.816837118', fmt='%.9f')  # ?
        self.test('cE',  e.cE,  '1.530757637', fmt='%.9f')
        self.test('cG',  e.cG,  '1.530757637', fmt='%.9f')
        self.test('cH',  e.cH,  '0.795604230', fmt='%.9f')  # ?
        self.test('cK',  e.cK,  '1.612441349', fmt='%.9f')
        self.test('cKE', e.cKE, '0.081683712', fmt='%.9f')  # ?
        self.test('cPi', e.cPi, '1.612441349', fmt='%.9f')
        phi = radians(20)
        sn, cn = sincos2(phi)
        self.test('fE(phi)', e.fE(phi), '0.348372822', fmt='%.9f')
        dn = e.fDelta(sn, cn)
        self.test('fDelta(sn, cn)', dn, '0.994133906', fmt='%.9f')

        self.test('fD(sn, cn, dn)',  e.fD(sn, cn, dn),  '0.013885234', fmt='%.9f')
        self.test('fE(sn, cn, dn)',  e.fE(sn, cn, dn),  '0.348372822', fmt='%.9f')
        self.test('fEd(PI_2)',       e.fEd(PI_2),       '0.027415224', fmt='%.9f')
        self.test('fEinv(PI_2)',     e.fEinv(PI_2),     '1.612999420', fmt='%.9f')
        self.test('fF(sn, cn, dn)',  e.fF(sn, cn, dn),  '0.349761345', fmt='%.9f')
        self.test('fG(sn, cn, dn)',  e.fG(sn, cn, dn),  '0.348372822', fmt='%.9f')
        self.test('fH(sn, cn, dn)',  e.fH(sn, cn, dn),  '0.335876111', fmt='%.9f')
        self.test('fPi(sn, cn, dn)', e.fPi(sn, cn, dn), '0.349761345', fmt='%.9f')
        try:
            t = e.fPi(0, None, 1)
        except EllipticError as x:
            t = str(x)
        self.test('fPi(sn, None, dn)', t, 'invokation Elliptic.fPi(0, None, 1): invalid')
        try:
            t = e.fH(0, 1, None)
        except EllipticError as x:
            t = str(x)
        self.test('fPi(sn, dn, None)', t, 'invokation Elliptic.fH(0, 1, None): invalid')

        self.test('deltaD(sn, cn, dn)',  e.deltaD(sn, cn, dn),  '-0.3223642', fmt='%.7f', nl=1)
        self.test('deltaE(sn, cn, dn)',  e.deltaE(sn, cn, dn),   '0.0084191', fmt='%.7f')
        self.test('deltaEinv(sn, cn)',   e.deltaEinv(sn, cn),   '-0.0082518', fmt='%.7f')
        self.test('deltaF(sn, cn, dn)',  e.deltaF(sn, cn, dn),  '-0.0083379', fmt='%.7f')
        self.test('deltaG(sn, cn, dn)',  e.deltaG(sn, cn, dn),   '0.0084191', fmt='%.7f')
        self.test('deltaH(sn, cn, dn)',  e.deltaH(sn, cn, dn),   '0.3140691', fmt='%.7f')
        self.test('deltaPi(sn, cn, dn)', e.deltaPi(sn, cn, dn), '-0.0083379', fmt='%.7f')
        try:
            t = e.deltaPi(0, None, 1)
        except EllipticError as x:
            t = str(x)
        self.test('deltaPi(sn, None, dn)', t, 'invokation Elliptic.deltaPi(0, None, 1): invalid')
        try:
            t = e.deltaPi(0, 1, None)
        except EllipticError as x:
            t = str(x)
        self.test('deltaPi(sn, dn, None)', t, 'invokation Elliptic.deltaPi(0, 1, None): invalid')

        # U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>} 3 Numerical checks
        self.test('RF(1,   2, 0)', RF(1,   2, 0), '1.3110287771461',  fmt='%.13f', nl=1)
        self.test('RF(0.5, 1, 0)', RF(0.5, 1, 0), '1.8540746773014',  fmt='%.13f')
        self.test('RF(2,   3, 4)', RF(2,   3, 4), '0.58408284167715', fmt='%.14f')

        self.test('RC(0,  1/4)', RC(0,  0.25), '3.1415926535898',  fmt='%.13f', nl=1)  # PI
        self.test('RC(9/4,  2)', RC(2.25,  2), '0.69314718055995', fmt='%.14f')  # ln(2)
        self.test('RC(1/4, -2)', RC(0.25, -2), '0.23104906018665', fmt='%.14f')  # ln(2) / 3

        self.test('RJ(0, 1, 2, 3)',    RJ(0, 1, 2, 3),    '0.77688623778582', fmt='%.14f', nl=1)
        self.test('RJ(2, 3, 4, 5)',    RJ(2, 3, 4, 5),    '0.14297579667157', fmt='%.14f')
#       self.test('RJ(2, 3, 4, -0.5)', RJ(2, 3, 4, -0.5), '0.24723819703052', fmt='%.14f')
#       self.test('RJ(2, 3, 4, -5)',   RJ(2, 3, 4, -5),  '−0.12711230042964', fmt='%.14f')

        self.test('RD(0, 2, 1)', RD(0, 2, 1), '1.7972103521034',  fmt='%.13f', nl=1)
        self.test('RD(2, 3, 4)', RD(2, 3, 4), '0.16510527294261', fmt='%.14f')

        self.test('RG(0, 16, 16)',     RG(0, 16, 16),     '3.1415926535898', fmt='%.13f', nl=1)  # PI
        self.test('RG(2,  3,  4)',     RG(2,  3,  4),     '1.7255030280692', fmt='%.13f')
        self.test('RG(0,  0.0796, 4)', RG(0,  0.0796, 4), '1.0284758090288', fmt='%.13f')  # E(0.99)?

        e.reset(0, 0)
        self.test('reset', len(e.__dict__), 4, nl=1)  # len(_k2, _kp2, _alpha2, _alphap2)

        self.test('sncndn(x)', fstr(e.sncndn(0), prec=9),     '0.0, 1.0, 1.0')
        t = fstr(e.sncndn(PI_2), prec=9)
        self.test('sncndn(x)', t, '1.0, 0.0, 1.0', known=t == '1.0, -0.0, 1.0')
        e.reset(1, 1)
        self.test('sncndn(x)', fstr(e.sncndn(0), prec=9),     '0.0, 1.0, 1.0')
        self.test('sncndn(x)', fstr(e.sncndn(PI_2), prec=9),  '0.917152336, 0.398536815, 0.398536815')
        self.test('sncndn(x)', type(e.sncndn(PI_4)), Elliptic3Tuple, nt=1)

        self.testCopy(e)

        for f in range(4):  # coverage
            t = (f / 4.0,) * 4
            e = Elliptic(*t)
            s = 'k2 alpha2 kp2 alphap2'
            e = tuple(getattr(e, a) for a in s.split())
            t = tuple(Scalar(f, name=n) for n, f in zip(s.split(), t))
            self.test(s, e, t)

        f = elliptic._convergenceError
        x = f(1, 2)
        self.test(f.__name__, str(x), "no convergence (1), tolerance (2)", nl=1)

        f = elliptic._ellipticError
        x = f(self.testElliptic, None, txt='test')
        self.test(f.__name__, str(x), "invokation Elliptic.testElliptic(None): test", nl=1)
        x = f(f.__name__, None, cause=x)
        self.test(f.__name__, str(x), "invokation Elliptic.testElliptic(None): test", known=endswith)


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testElliptic()
    t.results()
    t.exit()
