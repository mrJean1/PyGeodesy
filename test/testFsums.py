
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '22.04.18'

from base import isPython2, startswith, TestsBase

from pygeodesy import Fsum, fsum, fsum_, fsums, NN, ResidualError

from math import ceil, floor
from random import random, gauss, shuffle
from sys import getsizeof

_dot0 = '.0' if isPython2 else NN


class Tests(TestsBase):

    def testFsums(self):  # MCCABE 26

        # U{Neumaier<https://WikiPedia.org/wiki/Kahan_summation_algorithm>}
        t = 1, 1e101, 1, -1e101
        for i in range(1, 11):
            s = float(len(t) / 2)  # number of ones
            self.test('sum' + str(i), sum(t), s, known=True)
            self.test('fsum', fsum(t), s)
            self.test('Fsum', Fsum().fsum(t), s)
            t += t

        # <https://GitHub.com/ActiveState/code/blob/master/recipes/Python/
        #        393090_Binary_floating_point_summatiaccurate_full/recipe-393090.py>
        t = 1.0, 0.0050, 0.0000000000010
        s = 1.0050000000010
        self.test('sum', sum(t), s, known=True)
        self.test('fsum', fsum(t), s)
        self.test('Fsum', Fsum().fsum(t), s)

        # <https://GitHub.com/python/cpython/blob/master/Modules/mathmodule.c>
        t = 1e-16, 1, 1e16
        s = '1.0000000000000002e+16'  # half-even rounded
        self.test('fsum', fsum(t), s, prec=-16)
        self.test('Fsum', Fsum().fsum(t), s, prec=-16, nt=1)

        # <https://GitHub.com/ActiveState/code/blob/master/recipes/Python/
        #        393090_Binary_floating_point_summatiaccurate_full/recipe-393090.py>
        for i in range(1, 101):
            t = [7, 1e100, -9e-20, -7, -1e100, 8e-20] * 10
            s = 0
            for _ in range(20):
                v = gauss(0, random())**7 - s
                s += v
                t.append(v)
            shuffle(t)
            f = Fsum().fadd(t)
            s = f.fsum()
            self.test('sum' + str(i), sum(t), s, known=True)
            self.test('fsum', s, s)

            p = f * f * f * f
            self.test('pow(4)', f.pow(4), p, known=True)  # abs(p - f**4) < 1e-3
            p = f.pow(1)
            self.test('pow(1)', p, f, known=p == f)
            self.test('pow(0)', f.pow(0), "fsums.Fsum[1] (1.0, 0)")

            i = f.ceil  # Python 2, ceil(f) != f.__ceil__
            self.test('ceil', (i - 1) < f <= i, True)
            i = f.floor  # Python 2, floor(f) != f.__floor__
            self.test('floor', i <= f < (i + 1), True)
            n = len(t) // 2
            i, m = divmod(f, n)  # == f.__divmod__ in Python 2
            f -= (m + i * n)
            self.test('divmod', f, '0.0', known=f == 0)

            r = f.residual
            self.test('residual', r, '0')
            self.test('is_exact', f.is_exact(), (not r))

            f = Fsum()
            f.fsum(t[:n])  # test ps
            self.test('Fsum', f.fsum(t[n:]), s)
            f = Fsum()
            f.fsum(t[n:])  # test ps
            self.test('Fsum', f.fsum(t[:n]), s)
        p = f * (f * 1e10)  # coverage Fsum.__imul__
        f *= f * 1e10
        self.test('fmul', p.fsum(), f.fsum(), prec=8, nt=1)

        t = 1, 1e101, 1, -1e101
        for i in range(1, 9):
            a = Fsum(*t)
            self.test('len' + str(i), len(a), len(a))
            b = a.fcopy()
            c = a + b
            self.test('FSum+', c.fsum(), a.fsum() + b.fsum())
            c -= a
            self.test('FSum-', c.fsum(), b.fsum())
            c -= b
            self.test('FSum-', c.fsum(), 0.0)
            b = a * 2
            a += a
            self.test('FSum*', a.fsum(), b.fsum())
            t += t
            self.testCopy(a, '_fint2', '_fprs', '_fprs2', '_n', '_ps', deep=True)
        self.test('len', len(a), 513)
        t = a.partials
        self.test('partials', t, t, nt=1)

        c = a - b
        self.test('FSum0', c.fsum(), 0.0)
        c -= 0
        self.test('FSum0', c.fsum(), 0.0)
        c -= c
        self.test('FSum0', c.fsum(), 0.0)
        c *= Fsum(1.0)
        self.test('FSum0', c.fsum(), 0.0)
        t = a + a.fcopy().fmul(-1)
        self.test('FSum0', t.fsum(), 0.0)
        a.fsub_(*a._ps)
        self.test('FSum0', a.fsum(), 0.0)
        self.test('Fsum#', len(a), 514)
        self.test('Fsum#', len(a._ps), 1)
        self.test('FSum.', a, 'fsums.Fsum[4097] (0.0, 0)', known=not a)

        self.test('FsumI', c.imag, 0.0)
        self.test('FsumR', c.real, float(c))

        self.test('radd', float(2 + b),  '514.0')
        self.test('rdiv', float(2 / b),  '3.91e-03', prec=-2)
        self.test('rmul', float(2 * b),  '1024.0')
        self.test('rpow', float(2**a),      '1.0')
        self.test('rsub', float(2 - b), '-510.0')
        z = getsizeof(t)
        self.test('sizeof', z, 372, known=100 < z < 600)  # or 456 Python 2.
        try:
            self.test('_2sum', fsums._2sum(1e308, 1e803), OverflowError.__name__)
        except OverflowError as x:
            self.test('_2sum', repr(x), repr(x))
        try:
            self.test('F(None)', Fsum(None).fsum(), TypeError.__name__)
        except Exception as X:
            self.test('F(None)', repr(X), TypeError.__name__, known=startswith, nt=1)

        # Knuth/Kulisch, TAOCP, vol 2, p 245, sec 4.2.2, item 31, see also .testKarney.py
        # <https://SeriousComputerist.Atariverse.com/media/pdf/book/
        #  Art%20of%20Computer%20Programming%20-%20Volume%202%20(Seminumerical%20Algorithms).pdf>
        x = 408855776
        y = 708158977
        self.test('ints', 2*y**2 +  9*x**4 - y**4, 1)
        self.test('ints', 2*y**2 + (3*x**2 - y**2) * (3*x**2 + y**2), 1)
        t = 2*float(y)**2, 9*float(x)**4, -(float(y)**4)
        self.test('fsum ', fsum(t),          '1.0', prec=-8, known=True)  # -3.589050987400773e+19
        self.test('fsum_', fsum_(*t),        '1.0', prec=-8, known=True)
        self.test('Fsum ', Fsum().fsum_(*t), '1.0', prec=-8, known=True)
        self.test('sum  ',  sum(t),          '1.0', prec=-8, known=True)

        t = Fsum(1, 1e101, 1, -1e101)
        d = t * 2
        self.test('F * 2', float(d), 4.0, prec=4, known=True)
        s = d / 2  # /= chockes PyChecker
        self.test('F / 2', float(s), float(t), prec=4, known=True)
        self.test('F / F', s / s == Fsum(1.0), True, known=True)  # PYCHOK "s / s is always 1 or ZeroDivisionError"
        self.test('F / F', float(s / t), '1.0', known=True)
        self.test('F / F', float(d / t), '2.0', known=True)
        m = abs(t)
        self.test('abs  ', m, +t, known=m == t)  # PYCHOK "Unary positive (+) usually has no effect"
        self.test('int  ', int(t), 2, known=True)

        self.test('eq F', t == s,  True)
        self.test('ge F', t >= s,  True)
        self.test('gt F', t > s,  False)
        self.test('le F', t <= s,  True)
        self.test('lt F', t <  s, False)
        self.test('ne F', t != s, False)
        self.test('if F', bool(s), True)

        self.test('gt 0', t > 0, True)
        self.test('lt 0', t < 0, False)
        self.test('eq 0', t == 0, False)
        m = -t  # t is positive
        self.test('lt 0', m < 0, True)
        self.test('gt 0', t > 0, True)
        self.test('gt 0', m > 0, False)

        self.test('signOf', t.signOf(),  1)
        self.test('signOf', m.signOf(), -1)

        self.test('ceil ', ceil(t.fcopy().fadd_(1e-15)), '3', known=startswith)
        self.test('floor', floor(t), '2', known=startswith, nt=1)

        self.test('divmod ', divmod(d, 2),  ("(2%s, <fsums.Fsum '__divmod__'[2] (0.0, 0)" % _dot0),  known=startswith)
        self.test('divmod ', d.fcopy().divmod(2), ("(2%s, <fsums.Fsum 'divmod'[2] (0.0, 0)" % _dot0), known=startswith)
        self.test('rdivmod ', divmod(2, d), ("(0%s, <fsums.Fsum '__rdivmod__'[1] (2.0, 0)" % _dot0), known=startswith)
        self.test('divmod ', divmod(Fsum(-3), 2), ("(-2%s, <fsums.Fsum '__divmod__'[2] (1.0, 0)" % _dot0), known=startswith)
        m  = d.fcopy(name='__imod__')
        m %= 2
        self.test('imod', m, "fsums.Fsum '__imod__'[2] (0.0, 0)")
        self.test('mod ', d % 2, "fsums.Fsum '__mod__'[2] (0.0, 0)")
        self.test('rmod', 2 % d, "fsums.Fsum '__rmod__'[1] (2.0, 0)")
        m = -t
        self.test('neg ', m, -t, known=m == -t)
        m = +t  # PYCHOK "Unary positive (+) usually has no effect"
        self.test('pos ', m, t, known=m == t)
        self.test('is_int', t.is_integer(), True, known=True)

        x = Fsum(1, 1e-101, -1, -1e-102)
        self.test('float', float(x), '9e-102', known=True)
        self.test('is_int', x.is_integer(), False, known=True)

        m = Fsum(1, 1e-101, -4, -1e-102)  # about -3
        self.test('F //', m // 3,  "fsums.Fsum '__floordiv__'[1] (-1, 0)")
        try:
            self.test('// F', 5 // m, "fsums.Fsum '__rfloordiv__'[1] (-2, 0)")
        except Exception as X:
            self.test('// F', repr(X), ResidualError.__name__, known=startswith)
        m.__ifloordiv__(2)  # //= chockes PyChecker
        self.test('F //=', m, "fsums.Fsum[1] (-2, 0)")
        try:
            self.test('F / 0', x / 0, ZeroDivisionError.__name__)
        except Exception as X:
            self.test('F / 0', repr(X), ZeroDivisionError.__name__, known=startswith)

        try:
            self.test('pow(F, +)', pow(x, 2.1), ResidualError.__name__)
        except Exception as X:
            self.test('pow(F, +)', repr(X), ResidualError.__name__, known=startswith)
        try:
            self.test('pow(F, -)', pow(x, -1), "fsums.Fsum '__pow__'[1] (1.11111e+101, 0)")
        except Exception as X:
            self.test('pow(F, -)', repr(X), ResidualError.__name__, known=startswith)
        try:
            self.test('pow(-F, F)', pow(m, x), ValueError.__name__)
        except Exception as X:
            self.test('pow(-F, F)', repr(X), ValueError.__name__, known=startswith)
        try:
            self.test('pow(F, F)', pow(-m, x), Fsum(1))  # -m = 2, x = 0.+
        except Exception as X:
            self.test('pow(F, F)', repr(X), ResidualError.__name__, known=startswith)
        try:
            self.test('pow(F, f, i)', pow(x, 2.1, 2), ResidualError.__name__)
        except Exception as X:
            self.test('pow(F, f, i)', repr(X), ResidualError.__name__, known=startswith)
        try:
            self.test('pow(F, F, i)', m.pow(Fsum(2.1), 2), TypeError.__name__)
        except Exception as X:
            self.test('pow(F, F, i)', repr(X), TypeError.__name__, known=startswith)
        self.test('pow(F, F, None)', pow(-m, 2, None), ("fsums.Fsum[1] (4%s, 0)" % _dot0))
        try:
            self.test('Z**-2', Fsum(0.0)**-2, ZeroDivisionError.__name__)
        except Exception as X:
            self.test('Z**-2', repr(X), ZeroDivisionError.__name__, known=startswith)

        x = Fsum(1, 1e-101, -4, -1e-102)  # about -3
        self.test('pow(0)',  x**0,             '1.000', prec=3)
        self.test('pow(1)',  x**1,            '-3.000', prec=3)
        self.test('pow(2)',  x**2,             '9.000', prec=3)
        self.test('pow(21)', x**21, '-10460353203.000', prec=3)
        x **= 2
        self.test('**= 2',  x,                 '9.000', prec=3)

        x = Fsum()
        self.test('F0**0',  x**0,            "fsums.Fsum '__pow__'[1] (1.0, 0)")
        self.test('F0**0.', x**0.,           "fsums.Fsum '__pow__'[1] (1.0, 0)")
        self.test('0**F0',  0**x,            "fsums.Fsum '__rpow__'[1] (1.0, 0)")
        self.test('0.**F0', 0.**x,           "fsums.Fsum '__rpow__'[1] (1.0, 0)")
        self.test('F0**0', x.pow(0),         "fsums.Fsum 'pow'[1] (1.0, 0)")
        self.test('F0**2', x.pow(2),         "fsums.Fsum 'pow'[1] (0.0, 0)")
        self.test('F0**0.', x.pow(0.),       "fsums.Fsum 'pow'[1] (1.0, 0)")
        self.test('F0**3.', x.pow(3.),       "fsums.Fsum 'pow'[1] (0.0, 0)")
        self.test('F0**0.', x.pow(0., None), "fsums.Fsum 'pow'[1] (1, 0)")

        # preserve C{type(base)}
        t = Fsum()._fset(2, asis=True)
        self.test('2**F0',  (2**Fsum()).toStr(),   "fsums.Fsum '__rpow__'[1] (1.0, 0)")
        self.test('2.**F0', (2.**Fsum()).toStr(),  "fsums.Fsum '__rpow__'[1] (1.0, 0)")
        self.test('F2**0',  (t**0).toStr(),        "fsums.Fsum '__pow__'[1] (1.0, 0)")
        self.test('F2.**0', (Fsum(2.)**0).toStr(), "fsums.Fsum '__pow__'[1] (1.0, 0)")
        self.test('F2**F2',  t**t,                 "fsums.Fsum '__pow__'[1] (4, 0)")
        self.test('F2**F2',  t.__rpow__(t),        "fsums.Fsum '__rpow__'[1] (4, 0)")

        x = Fsum(2, 3)
        self.test('F**2',     x**x,     '3125.000', prec=3)
        self.test('F**-1',    x**-1,       '0.200', prec=3)
        self.test('F**-2',    x**-2,       '0.040', prec=3)
        self.test('F**-2.5',  x**-2.5,     '0.018', prec=3)
        self.test('F** 2.5',  x** 2.5,    '55.902', prec=3)
        self.test('pow(2)',   x.pow(2),   '25.000', prec=3)
        self.test('pow(2.5)', x.pow(2.5), '55.902', prec=3)
        self.test('pow(F)',   x.pow(x), '3125.000', prec=3)

        self.test('3pow(2, None)',   x.pow(2, None).toStr(),         "fsums.Fsum 'pow'[1] (25, 0)")
        self.test('3pow(2.5, None)', x.pow(2.5, None).toStr(prec=5), "fsums.Fsum 'pow'[1] (55.902, 0)")
        self.test('3pow(2, 20)',     x.pow(2, 20).toStr(),           "fsums.Fsum 'pow'[1] (5, 0)")

        t = x.fsum()
        self.test('fsum()',  t, '5.0')
        self.test('fsum()',  t is x.fsum(), True)  # Property_RO
        t = x.fsum2()
        self.test('fsum2()', t, '(5.0, 0)')
        self.test('fsum2()', t is x.fsum2(), True)  # Property_RO
        self.test('fsum2()', t.toRepr(), 'Fsum2Tuple(fsum=5.0, residual=0)')

        s, r = t.toUnits()
        t = t.__class__.__name__
        self.test(t, (s.name, s, s.__class__), "('fsum', 5.0, <class 'pygeodesy.units.Float'>)")
        self.test(t, (r.name, r, r.__class__), "('residual', 0, <class 'pygeodesy.units.Int'>)")

        self.test('fmul(x)', x.fmul(x),          '25.0', prec=1)
        self.test('fmul(F)', x.fmul(Fsum(2.5)),  '62.5', prec=1)
        self.test('fadd(F)', x.fadd(Fsum(2.5)),  '65.0', prec=1)
        self.test('fsub(F)', x.fsub(Fsum(2.5)),  '62.5', prec=1)  # iter(Fsum)
        self.test('Fsum(F)', Fsum(x, x).fsum(), '125.0', prec=1, nt=1)

        f = Fsum(1, 1e-11, -4, -1e-12)  # about -3
        r = f.as_integer_ratio()
        self.test('ratio', str(r).replace('L', NN), '(-27021597764141911, 9007199254740992)')  # L on Windows
        t = Fsum(r[0] / r[1])  # C{int} in Python 2
        self.test('ratio', t, f, known=abs(t.fsum() - f.fsum()) < 1e-9)  # python special
        self.test('int_float', t.int_float(), '-3.000', prec=3)
        self.test('fint',  t.fint(raiser=False), ("fsums.Fsum 'fint'[1] (-3, 0)" if isPython2
                                             else "fsums.Fsum 'fint'[1] (-2, 0)"))
        self.test('fint2', t.fint2(), ('(-3, 0' if isPython2 else '(-2, -1.0)'), known=startswith)

        f = Fsum(3) // 1
        try:
            t = pow(f, 3, 4)
            self.test('pow3', t, "fsums.Fsum '__pow__'[1] (3", known=startswith)
        except Exception as X:
            self.test('pow3', repr(X), TypeError.__name__, known=startswith)

        self.test('math_fsum', f.is_math_fsum(), True, known=isPython2)

        self.test('RESIDUAL', f.RESIDUAL(), False)
        self.test('RESIDUAL', f.RESIDUAL(True), False)
        self.test('RESIDUAL', f.RESIDUAL(None), True)


if __name__ == '__main__':

    t = Tests(__file__, __version__, fsums)
    t.testFsums()
    t.results()
    t.exit()
