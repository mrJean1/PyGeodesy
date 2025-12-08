
# -*- coding: utf-8 -*-

# Test L{angles} module.

__all__ = ('Tests',)
__version__ = '25.12.01'

from bases import startswith, TestsBase

from pygeodesy import EPS, INF, PI, Ang, Deg, Fsum, \
                      Lambertian, Rad, angles

from math import ceil, degrees, floor, radians


class Tests(TestsBase):

    def testAngles(self):  # MCCABE 26

        for d in range(-716, 716, 3):
            D = Deg(d)
            self.test('D', D.degrees, d, prec=6)
            self.test('r', D.radians, radians(d), prec=6)
            self.test('n', D.n, D.n)

            r = d * 0.01
            R = Rad(r)
            self.test('R', R.radians, r, prec=6)
            self.test('d', R.degrees, degrees(r), prec=6)
            self.test('n', R.n, R.n)

        D = Deg(30)
        self.test('radd', float(2 + D), '32.0', prec=1, nl=1)
        self.test('rdiv', float(2 / D),  '6.67e-02', prec=-2)
        self.test('rmul', float(2 * D), '60.0', prec=1)
        self.test('rpow', float(2**D),  '1073741824.0', prec=1)
        self.test('rsub', float(2 - D), '-28.0', prec=1)

        R = Rad(2)
        r = R * 2
        self.test('R * 2', float(r), 4.0, prec=4, known=True, nl=1)
        s = r / 2  # /= chockes PyChecker
        self.test('R / 2', float(s), float(R), prec=4, known=True)
        self.test('R / R', s / s == Ang(1.0), True, known=True)  # PYCHOK "s / s is always 1 or ZeroDivisionError"
        self.test('R / R', float(s / R), '1.0', known=True)
        self.test('R / R', float(r / R), '2.0', known=True)
        m = abs(R)
        self.test('abs  ', m, +R, known=m == R)  # PYCHOK "Unary positive (+) usually has no effect"
        self.test('int  ', int(R), 2, known=True)

        self.test('eq R', R == s,  True)
        self.test('ge R', R >= s,  True)
        self.test('gt R', R > s,  False)
        self.test('le R', R <= s,  True)
        self.test('lt R', R <  s, False)
        self.test('ne R', R != s, False)
        self.test('if R', bool(s), True)

        self.test('gt 0', R > 0, True)
        self.test('lt 0', R < 0, False)
        self.test('eq 0', R == 0, False)
        m = -R  # R is positive
        self.test('lt 0', m < 0, True)
        self.test('gt 0', R > 0, True)
        self.test('gt 0', m > 0, False)

        self.test('signOf', R.signOf(),  1)
        self.test('signOf', m.signOf(), -1)

        self.test('ceil ', ceil(R.copy() + 1e-15), '3', known=startswith)
        self.test('floor', floor(R), '2', known=startswith, nt=1)

        self.test('divmod ',  divmod(r, 2),        '(2.0, Radians(0.0))')
        self.test('divmod ',  r.copy().divmod(2),  '(2.0, Radians(0.0))')
        self.test('rdivmod ', divmod(2, r),        '(0.0, Radians(2.0))')
        self.test('divmod ',  divmod(Ang(-3), 2), '(-2.0, Radians(1.0))')
        m  = r.copy(name='__imod__')
        m %= 2
        self.test('imod', m,     '0.0')
        self.test('mod ', r % 2, '0.0')
        self.test('rmod', 2 % r, '2.0')
        m = -R
        self.test('neg ', m, -R, known=m == -R)
        m = +R  # PYCHOK "Unary positive (+) usually has no effect"
        self.test('pos ', m, R, known=m == R)
        self.test('is_int', R.is_integer(), False, known=True)

        x = Rad(float(Fsum(1, 1e-101, -1, -1e-102)))
        self.test('float', float(x), '9e-102', known=True)
        self.test('is_int', x.is_integer(), False, known=True)
        self.test('round1', float(round(x, 1)), '0.0')  # Python2

        m = Deg(float(Fsum(1, 1e-101, -4, -1e-102, name='m')))  # about -3
        self.test('R //', m // 3,  '-2.0')
        self.test('// R', 5 // m, '-2.0')
        m.__ifloordiv__(2)  # //= chockes PyChecker
        self.test('R //=', m, '-2.0')  # -3.0 // 2 = -2.0, 3.0 // 2 = 1.0
        try:
            self.test('R / 0', x / 0, ZeroDivisionError.__name__)
        except Exception as X:
            self.test('R / 0', repr(X), ZeroDivisionError.__name__, known=startswith)

        try:
            self.test('pow(R, +)', pow(x, 2.1), '0.0', nl=1)
        except Exception as X:
            self.test('pow(R, +)', repr(X), ZeroDivisionError.__name__, known=startswith)
        try:
            self.test('pow(R, -)', pow(x, -1), '0.4767037')
        except Exception as X:
            self.test('pow(R, -)', repr(X), ZeroDivisionError.__name__, known=startswith)
        try:
            self.test('pow(-R, R)', pow(m, float(x)), ZeroDivisionError.__name__)
        except Exception as X:  # TypeError('fromDegrees((1+2.8274333882308138e-101j))')
            self.test('pow(-R, R)', repr(X), TypeError.__name__, known=True)  # =startswith)
        try:
            self.test('pow(R, R)', pow(-m, x), '1.0')  # -m = 2, x = 0.+
        except Exception as X:
            self.test('pow(R, R)', repr(X), ZeroDivisionError.__name__, known=startswith)
        try:
            self.test('pow(R, f, i)', pow(x, 2.1, 2), ZeroDivisionError.__name__)
        except Exception as X:
            self.test('pow(R, f, i)', repr(X), TypeError.__name__, known=startswith)
        try:
            self.test('pow(R, R, i)', m.pow(Fsum(2.1), 2), TypeError.__name__)
        except Exception as X:
            self.test('pow(R, R, i)', repr(X), TypeError.__name__, known=startswith)
        self.test('pow(R, i, None)', pow(-m, 2, None), '4.0')
        try:
            self.test('Z**-2', Fsum(0.0)**-2, ZeroDivisionError.__name__)
        except Exception as X:
            self.test('Z**-2', repr(X), ZeroDivisionError.__name__, known=startswith)

        x = Rad(float(Fsum(1, 1e-101, -4, -1e-102)))  # about -3
        self.test('pow(0)',  x**0,             '1.000', prec=3)
        self.test('pow(1)',  x**1,            '-3.000', prec=3)
        self.test('pow(2)',  x**2,             '9.000', prec=3)
        self.test('pow(21)', x**21, '-10460353203.000', prec=3)
        self.test('pow(-5)', x**-5,           '-0.004', prec=3)
        x **= 2
        self.test('**= 2',  x,                 '9.000', prec=3)

        x = Deg(0)
        self.test('F0**0',  x**0,            '1.0', nl=1)
        self.test('F0**0.', x**0.,           '1.0')
        self.test('0**F0',  0**x,            '1.0')
        self.test('0.**F0', 0.**x,           '1.0')
        self.test('F0**0', x.pow(0),         '1.0')
        self.test('F0**2', x.pow(2),         '0.0')
        self.test('F0**0.', x.pow(0.),       '1.0')
        self.test('F0**3.', x.pow(3.),       '0.0')
        self.test('F0**0.', x.pow(0., None), '1.0')

        x = Rad(float(Fsum(2, 3)))
        self.test('R**2',     x**x,     '3125.000', prec=3, nl=1)
        self.test('R**-1',    x**-1,       '0.200', prec=3)
        self.test('R**-2',    x**-2,       '0.040', prec=3)
        self.test('R**-2.5',  x**-2.5,     '0.018', prec=3)
        self.test('R** 2.5',  x** 2.5,    '55.902', prec=3)
        self.test('pow(2)',   x.pow(2),   '25.000', prec=3)
        self.test('pow(2.5)', x.pow(2.5), '55.902', prec=3)
        self.test('pow(R)',   x.pow(x), '3125.000', prec=3)

        self.test('3pow(2, None)',   x.pow(2, None).toStr(), '25.0')
        self.test('3pow(2.5, None)', x.pow(2.5, None).toStr(prec=5), '55.9017')
        self.test('3pow(2, 20)', pow(int(x), 2, 20), '5')

        x = float(x)
        self.test('x * x', x * x,        '25.0', prec=1)
        self.test('x * R', x * Ang(2.5), '12.5', prec=1)
        self.test('x + R', x + Ang(2.5), '7.5', prec=1)
        self.test('x - R', x - Ang(2.5), '2.5', prec=1)

        t = Ang(EPS)
        self.test('abs(T)',   abs(t),   '0.', known=startswith, nl=1)
        self.test('bool(T)',  bool(t),   True)
        self.test('float(T)', float(t), '2.',  known=startswith)
        self.test('int(T)',   int(t),   '0')
        self.test('-T',       (-t),     '-0.', known=startswith)
        self.test('+T',       (+t),     '0.', known=startswith)  # PYCHOK no effect

        p = Rad(PI)
        self.test('R==T', p == t, False, nl=1)
        self.test('R>=T', p >= t, True)
        self.test('R> T', p >  t, True)
        self.test('R<=T', p <= t, False)
        self.test('R< T', p <  t, False)
        self.test('R!=T', p != t, True)

        if True:  # coverage:
            D = Deg(61)
            t = D.base(D)
            self.test('base', t, t, nl=1)
            t = D.flipsign()
            self.test('flipsign', t, -D)
            t = D.lambertian
            self.test('lambertian', t, t)
            t = Ang.fromScalar(t, unit=Lambertian)
            self.test('fromScalar', t, t)  # D?
            t = D.copy()
            t.n0 = D.n0 + 2
            self.test('n0', t, t)
            t = D.nearest(1)
            self.test('nearest', t, t)
            t = D.normalize(2).normalize(0)
            self.test('normalize', t, D)
            t = D.quadrant
            D.quadrant = t + 4
            self.test('quadrant', D.quadrant, t)
            t = D.copy().reflect(True, True, True).reflect(True, True, True)
            self.test('reflect', t, D)
            t = D.round()
            self.test('round', t, t)
            t = D.shift(INF)
            self.test('shift', t, t)
            self.test('t', D.t, D.t)
            t = D.toLambertian()
            self.test('toLambertian', t, t)
            t = D.toTuple()
            self.test('toTuple', t, t)


if __name__ == '__main__':

    t = Tests(__file__, __version__, angles)
    t.testAngles()
    t.results()
    t.exit()
