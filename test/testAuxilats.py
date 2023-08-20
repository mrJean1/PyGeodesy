
# -*- coding: utf-8 -*-

# Some basic L{auxilats} tests.

__all__ = ('Tests',)
__version__ = '23.08.16'

from bases import numpy, TestsBase

from pygeodesy import NN, PI_2, PI_4, Fsum, fsum, printf, sincos2
from pygeodesy.auxilats import Aux, AuxAngle, AuxDST, AuxLat, \
                               AuxBeta, AuxChi, AuxMu, AuxPhi, AuxTheta, AuxXi
from math import fabs


class Tests(TestsBase):

    def testAngles(self, deltas):
        de = {}
        # <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AuxLatitude.html>
        aL = AuxLat()
        for A in (AuxBeta, AuxChi, AuxMu, AuxPhi, AuxTheta, AuxXi):
            AUX = A._AUX
            for d in range(0, 360, 3):
                a = A.fromDegrees(d)
                for exact in (False, True):
                    for auxout in range(len(Aux)):
                        r = aL.convert(auxout, a, exact=exact)
                        assert r._AUX == auxout
                        b = aL.convert(AUX, r, exact=exact)
                        assert b._AUX == AUX
                        i =  b.iteration
                        i =  NN if i is None else (', iteration=' + str(i))
                        n = '%2d %.12f %s%s' % (d, (r.toDegrees % 360.0), r, i)
                        self.test(n, b.tan, a.tan, prec=12, nl=int(not auxout))
                        if deltas and b != a:  # or i:
                            e = fabs(b.tan - a.tan)
                            printf('%s tan=%s, e=%.3e', b, b.tan, e)
                            printf('%s tan=%s, iteration=%s', a, a.tan, b.iteration)
                            if e > de.get(d, 0):
                                de[d] = e
        n = 1
        for d, e in sorted(de.items()):
            printf('%2d error %.3e ', d, e, nl=n)
            n = 0

        A = AuxAngle(2)
        self.test('abs', abs(A), 'AuxAngle(tan=2.0, x=1.0, y=2.0)', nl=1)
        self.test('add', A + A,  'AuxAngle(tan=-1.33333, x=-3, y=4.0)')
        self.test('eq ', A == A, True)
        self.test('float', float(A), 2.0)
        self.test('sub', A - A,  'AuxAngle(tan=0.0, x=1.0, y=0.0)')
        self.test('neg', -A,     'AuxAngle(tan=-2, x=1.0, y=-2)')
        self.test('ne ', A != A, False)
        self.test('pos', +A, A)  # PYCHOK no effect
        A += A
        self.test('iadd', A,     'AuxAngle(tan=-1.33333, x=-3, y=4.0)')
        A -= A
        self.test('isub', A,     'AuxAngle(tan=0.0, x=1.0, y=0.0)')
        self.test('radd', A.__radd__(A), 'AuxAngle(tan=0.0, x=1.0, y=0.0)')
        self.test('rsub', A.__rsub__(A), 'AuxAngle(tan=0.0, x=1.0, y=0.0)')

    def testCoeffs(self):
        a = AuxLat()
        for al in (4, 6, 8):
            a.ALorder = al
            self.test('Aux', al, al, nl=1)
            for ain in range(Aux.N):
                for aout in range(Aux.N):
                    t = a._coeffs(ain, aout)
                    self.test('Aux', len(t), al)

    def testDST(self, enums):
        # <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1DST.html>

        def _f(a):  # sawtooth
            return a + PI_4

        def _j(i):
            return (2 * i + 1)**2 * (1 - ((i & 1) << 1))

        N = 5
        dst = AuxDST(N)
        self.test('N', dst.N, N, prec=0, nl=1)

        tx = dst.transform(_f)
        if enums:
            for i, t in enumerate(tx):
                printf('N /%s: %g %g', i, t, t * _j(i), nl=not i)
        self.test('N /sum', fsum(tx), '2.748844788926', prec=12)

        tx = dst.refine(_f, tx)
        if enums:
            for i, t in enumerate(tx):
                printf('+N/%s: %g %g', i, t, t * _j(i), nl=not i)
        self.test('+N/sum', fsum(tx), '3.071245975238', prec=12)

        K  = dst.reset(N*2)
        tx = dst.transform(_f)
        if enums:
            for i, t in enumerate(tx):
                printf('2N/%s: %g %g', i, t, t * _j(i), nl=not i)
        self.test('2N/sum', fsum(tx), '3.071245975238', prec=12)

        Te, Tg = Fsum(), Fsum()
        for i in range(K):
            x = PI_2 * i / K
            s, c = sincos2(x)
            e = AuxDST.evaluate(s, c, tx)
            g = AuxDST.integral(s, c, tx)
            if enums:
                printf('T/%s: f(%g)=%g, e=%g, g=%g', i, x, _f(x), e, g, nl=not i)
            Te += e
            Tg += g
        self.test('Te/sum', float(Te),   '57.582664067074', prec=12)
        self.test('Tg/sum', float(Tg), '-182.807444594653', prec=12)


if __name__ == '__main__':

    from pygeodesy import auxilats

    t = Tests(__file__, __version__, auxilats)

    t.testCoeffs()
    t.testAngles(False)
    if numpy:
        t.testDST(False)
    else:
        t.skip('no numpy', 6)

    t.results()
    t.exit()
