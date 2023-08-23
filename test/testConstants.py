
# -*- coding: utf-8 -*-

# Test L{constants} module.

__all__ = ('Tests',)
__version__ = '23.08.23'

from bases import isPyPy, TestsBase

from pygeodesy import Float, Int, Radius, basics, constants, interns, \
                      float_, isinf, isint0, isnan, \
                      EPS, EPS0, EPS02, EPS1, EPS2, EPS_2, EPS4, \
                      INF, INT0, NAN, NEG0, NINF


class Tests(TestsBase):

    def testConstants(self):

        _all    = set(constants.__all__)
        _DOT_   = interns._DOT_
        _UNDER_ = interns._UNDER_
        _off90  = constants._off90
        _0_0    = constants._0_0
        _0_0s   = constants._0_0s

        for n in sorted(dir(constants), key=str.lower):
            v = getattr(constants, n, None)
            if isinstance(v, (Float, Int, Radius)):
                r = v.toRepr(std=False)
                self.test(n, r, r, nl=1)
                self.test(n, v.name, n)
                if not isPyPy:
                    self.test(n, v.name is n, True)
                if not n.startswith(_UNDER_):
                    self.test(n, n in _all, True)
            elif isinstance(v, float):
                a = n.strip(_UNDER_).replace(_UNDER_, _DOT_)
                self.test(n, v, a, known=n.startswith(_UNDER_))

        self.test('EPS',        EPS  > 0, True, nl=1)
        self.test('EPS+1', (1 + EPS) > 1, True)
        self.test('EPS-1', (1 - EPS) < 1, True)

        self.test('EPS0',  0 < EPS0  < EPS,  True, nl=1)
        self.test('EPS02', 0 < EPS02 < EPS0, True)

        self.test('EPS_2', 0 < EPS_2 < EPS, True, nl=1)
        self.test('EPS_2',     EPS_2,  EPS / 2)

        self.test('EPS1',    0 < EPS1  < 1, True, nl=1)
        self.test('EPS1+1', (1 - EPS1) > 0, True)
        self.test('EPS1-1', (EPS1 - 1) < 0, True)

        self.test('EPS2', EPS2 > EPS, True, nl=1)
        self.test('EPS2', EPS2,  EPS  * 2)

        self.test('EPS4', EPS4 > EPS2, True, nl=1)
        self.test('EPS4', EPS4,  EPS2 * 2)

        self.test('INF',  isinf(INF),  True, nl=1)
        self.test('INF',  INF == NINF, False)
        self.test('NINF', isinf(NINF), True)
        self.test('NINF',       NINF, -INF)

        self.test('INT0',        INT0,  '0', nl=1)
        self.test('INT0', isint0(INT0),  True)
        self.test('INT0', isint0(0),     True)
        self.test('INT0', isint0(0.0),   False)
        self.test('INT0', isint0(True),  False)
        self.test('INT0', isint0(False), False)

        self.test('NAN',  isnan(NAN),  True, nl=1)
        self.test('NAN',  NAN == INF,  False)
        self.test('NAN',  NAN == NINF, False)

        self.test('NEG0',       NEG0,  '-0.0', nl=1)
        self.test('NEG0',  0 == NEG0,   True)
        self.test('NEG0',  bool(NEG0),  False)
        self.test('NEG0',   abs(NEG0), '0.0')

        self.test(_off90.__name__, _off90(90) < 90, True, nl=1)
        self.test(_off90.__name__, _off90(90) > 89.999999, True)

        t = float_(1, 2, 3)
        self.test(float_.__name__, t, '(1.0, 2.0, 3.0)')
        t = float_(3.14, sets=True)
        self.test(float_.__name__, t is float(3.14), True)

        self.test('_0_0', _0_0 is basics._0_0, True, nl=1, nt=1)

        n = _0_0s.__name__
        for i in (0, 1, 2, 3, 5, 8, 9, 10, 12, 25, 49, 129, 257):
            t = _0_0s(i)
            self.test(n, len(t), i)
            self.test(n, t.count(_0_0), i)


if __name__ == '__main__':

    t = Tests(__file__, __version__, constants)
    t.testConstants()
    t.results()
    t.exit()
