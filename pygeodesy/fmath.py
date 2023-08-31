
# -*- coding: utf-8 -*-

u'''Utilities using precision floating point summation.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import _copysign, copysign0, isint, isscalar, len2
from pygeodesy.constants import EPS0, EPS02, EPS1, NAN, PI, PI_2, PI_4, \
                               _0_0, _0_125, _0_25, _0_5, _1_0, _N_1_0, \
                               _1_3rd, _1_5, _1_6th, _2_0, _2_3rd, _3_0, \
                               _isfinite, isnear0, isnear1, _over, remainder
from pygeodesy.errors import _IsnotError, LenError, _TypeError, _ValueError, \
                             _xError, _xkwds_get, _xkwds_pop
from pygeodesy.fsums import _2float, _Powers, Fsum, _fsum, fsum, fsum1_, \
                            _pow_op_,  Fmt, unstr
from pygeodesy.interns import MISSING, _few_, _h_, _invokation_, _negative_, \
                             _not_scalar_, _SPACE_, _too_
from pygeodesy.lazily import _ALL_LAZY, _sys_version_info2
# from pygeodesy.streprs import Fmt, unstr  # from .fsums
from pygeodesy.units import Int_,  Float_  # PYCHOK for .heights

from math import fabs, sqrt  # pow
from operator import mul as _mul  # in .triaxials

__all__ = _ALL_LAZY.fmath
__version__ = '23.08.31'

# sqrt(2) <https://WikiPedia.org/wiki/Square_root_of_2>
_0_4142 = 0.414213562373095  # sqrt(_2_0) - _1_0


class Fdot(Fsum):
    '''Precision dot product.
    '''
    def __init__(self, a, *b, **name):
        '''New L{Fdot} precision dot product M{sum(a[i] * b[i]
           for i=0..len(a))}.

           @arg a: Iterable, list, tuple, etc. (C{scalar}s).
           @arg b: Other values (C{scalar}s), all positional.
           @kwarg name: Optional name (C{str}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise LenError: Unequal C{len(B{a})} and C{len(B{b})}.

           @see: Function L{fdot} and method L{Fsum.fadd}.
        '''
        Fsum.__init__(self, **name)
        self.fadd(_map_mul(a, b, Fdot))


class Fhorner(Fsum):
    '''Precision polynomial evaluation using the Horner form.
    '''
    def __init__(self, x, *cs, **name):
        '''New L{Fhorner} evaluation of the polynomial
           M{sum(cs[i] * x**i for i=0..len(cs))}.

           @arg x: Polynomial argument (C{scalar}).
           @arg cs: Polynomial coeffients (C{scalar} or C{Fsum}
                    instances), all positional.
           @kwarg name: Optional name (C{str}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{x}}.

           @raise ValueError: Non-finite B{C{x}}.

           @see: Function L{fhorner} and methods L{Fsum.fadd} and L{Fsum.fmul}.
        '''
        Fsum.__init__(self, *cs[-1:], **name)
        if len(cs) > 1:
            x  = _2float(x=x)
            _a =  self._fadd    # (other, op)
            _f =  self._finite  # (other, op)
            op =  Fhorner.__name__
            ps =  self._ps
            for c in reversed(cs[:-1]):  # multiply-accumulate
                ps[:] = [_f(p * x, op) for p in ps]
                _a(c, op)
            # assert self._ps is ps


class Fhypot(Fsum):
    '''Precision hypotenuse of summation.
    '''
    def __init__(self, *xs, **power_name_RESIDUAL):
        '''New L{Fhypot} hypotenuse of (the I{power} of) several
           C{scalar} or C{Fsum} values.

           @arg xs: One or more values to include (each C{scalar}
                    or an C{Fsum} instance).
           @kwarg power_name_RESIDUAL: Optional exponent and root
                        order C{B{power}=2}, C{B{name}=NN} and
                        C{B{RESIDUAL}=None}, see L{Fsum.__init__}.
        '''
        try:
            p = _xkwds_pop(power_name_RESIDUAL, power=2)
            Fsum.__init__(self, **power_name_RESIDUAL)
            if xs:
                self._facc(_Powers(p, xs), up=False)  # PYCHOK yield
            self._fset(self._fpow(_1_0 / p, _pow_op_), asis=True)
        except Exception as X:
            raise self._ErrorX(X, xs, power=p)


class Fpolynomial(Fsum):
    '''Precision polynomial evaluation.
    '''
    def __init__(self, x, *cs, **name):
        '''New L{Fpolynomial} evaluation of the polynomial
           M{sum(cs[i] * x**i for i=0..len(cs))}.

           @arg x: Polynomial argument (C{scalar}).
           @arg cs: Polynomial coeffients (C{scalar}s), all
                    positional.
           @kwarg name: Optional name (C{str}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{x}}.

           @raise ValueError: Non-finite B{C{x}}.

           @see: Function L{fpolynomial} and method L{Fsum.fadd}.
        '''
        Fsum.__init__(self, *cs[:1], **name)
        n = len(cs) - 1
        if n > 0:
            self.fadd(_map_mul(cs[1:], fpowers(x, n), Fpolynomial))


class Fpowers(Fsum):
    '''Precision summation or powers, optimized for C{power=2}.
    '''
    def __init__(self, power, *xs, **name_RESIDUAL):
        '''New L{Fpowers} sum of (the I{power} of) several C{scalar}
           or C{Fsum} values.

           @arg power: The exponent (C{scalar} or C{Fsum}).
           @arg xs: One or more values to include (each C{scalar}
                    or an C{Fsum} instance).
           @kwarg power_name_RESIDUAL: Optional exponent and root
                        order C{B{power}=2}, C{B{name}=NN} and
                        C{B{RESIDUAL}=None}, see L{Fsum.__init__}.
        '''
        try:
            Fsum.__init__(self, **name_RESIDUAL)
            if xs:
                self._facc(_Powers(power, xs), up=False)  # PYCHOK yield
        except Exception as X:
            raise self._ErrorX(X, xs, power=power)


class Fn_rt(Fsum):
    '''Precision n-th root of summation.
    '''
    def __init__(self, root, *xs, **name_RESIDUAL):
        '''New L{Fn_rt} root of the precision sum of several
           C{scalar} or C{Fsum} values.

           @arg root: The order (C{scalar} or C{Fsum}).
           @arg xs: Values to include (each C{scalar} or an
                    C{Fsum} instance).
           @kwarg name_RESIDUAL: See L{Fsum.__init__}.
        '''
        try:
            Fsum.__init__(self, *xs, **name_RESIDUAL)
            self._fset(self._fpow(_1_0 / root, _pow_op_), asis=True)
        except Exception as X:
            raise self._ErrorX(X, xs, root=root)


class Fcbrt(Fn_rt):
    '''Precision cubic root of summation.
    '''
    def __init__(self, *xs, **name_RESIDUAL):
        '''New L{Fcbrt} cubic root of the precision sum of
           several C{scalar} or C{Fsum} values.

           @arg xs: Values to include (each C{scalar} or an
                    C{Fsum} instance).
           @kwarg name_RESIDUAL: See L{Fsum.__init__}.
        '''
        Fn_rt.__init__(self, _3_0, *xs, **name_RESIDUAL)


class Fsqrt(Fn_rt):
    '''Precision square root of summation.
    '''
    def __init__(self, *xs, **name_RESIDUAL):
        '''New L{Fsqrt} square root of the precision sum of
           several C{scalar} or C{Fsum} values.

           @arg xs: Values to include (each C{scalar} or an
                    C{Fsum} instance).
           @kwarg name_RESIDUAL: See L{Fsum.__init__}.
        '''
        Fn_rt.__init__(self, _2_0, *xs, **name_RESIDUAL)


def bqrt(x):
    '''Return the 4-th, I{bi-quadratic} or I{quartic} root, M{x**(1 / 4)}.

       @arg x: Value (C{scalar}).

       @return: I{Quartic} root (C{float}).

       @raise ValueError: Negative B{C{x}}.

       @see: Functions L{zcrt} and L{zqrt}.
    '''
    return _root(x, _0_25, bqrt)


try:
    from math import cbrt  # Python 3.11+

    def cbrt2(x):
        '''Compute the cube root I{squared} M{x**(2/3)}.
        '''
        return cbrt(x)**2  # cbrt(-0.0*2) == -0.0

except ImportError:  # Python 3.10-

    def cbrt(x):
        '''Compute the cube root M{x**(1/3)}.

           @arg x: Value (C{scalar}).

           @return: Cubic root (C{float}).

           @see: Functions L{cbrt2} and L{sqrt3}.
        '''
        # <https://archive.lib.MSU.edu/crcmath/math/math/r/r021.htm>
        # simpler and more accurate than Ken Turkowski's CubeRoot, see
        # <https://People.FreeBSD.org/~lstewart/references/apple_tr_kt32_cuberoot.pdf>
        return _copysign(pow(fabs(x), _1_3rd), x)  # cbrt(-0.0) == -0.0

    def cbrt2(x):  # PYCHOK attr
        '''Compute the cube root I{squared} M{x**(2/3)}.

           @arg x: Value (C{scalar}).

           @return: Cube root I{squared} (C{float}).

           @see: Functions L{cbrt} and L{sqrt3}.
        '''
        return pow(fabs(x), _2_3rd)  # XXX pow(fabs(x), _1_3rd)**2


def euclid(x, y):
    '''I{Appoximate} the norm M{sqrt(x**2 + y**2)} by
       M{max(abs(x), abs(y)) + min(abs(x), abs(y)) * 0.4142...}.

       @arg x: X component (C{scalar}).
       @arg y: Y component (C{scalar}).

       @return: Appoximate norm (C{float}).

       @see: Function L{euclid_}.
    '''
    x, y = fabs(x), fabs(y)
    if x < y:
        x, y = y, x
    return x + y * _0_4142  # XXX * _0_5 before 20.10.02


def euclid_(*xs):
    '''I{Appoximate} the norm M{sqrt(sum(x**2 for x in xs))}
       by cascaded L{euclid}.

       @arg xs: X arguments, positional (C{scalar}s).

       @return: Appoximate norm (C{float}).

       @see: Function L{euclid}.
    '''
    e = _0_0
    for x in sorted(map(fabs, xs)):  # XXX not reverse=True
        # e = euclid(x, e)
        if e < x:
            e, x = x, e
        if x:
            e += x * _0_4142
    return e


def facos1(x):
    '''Fast approximation of L{pygeodesy.acos1}C{(B{x})}.

       @see: U{ShaderFastLibs.h<https://GitHub.com/michaldrobot/
             ShaderFastLibs/blob/master/ShaderFastMathLib.h>}.
    '''
    a = fabs(x)
    if a < EPS0:
        r = PI_2
    elif a < EPS1:
        H = Fhorner(-a, 1.5707288, 0.2121144, 0.0742610, 0.0187293)
        r = float(H * sqrt(_1_0 - a))
        if x < 0:
            r = PI - r
    else:
        r = PI if x < 0 else _0_0
    return r


def fasin1(x):  # PYCHOK no cover
    '''Fast approximation of L{pygeodesy.asin1}C{(B{x})}.

       @see: L{facos1}.
    '''
    return PI_2 - facos1(x)


def fatan(x):
    '''Fast approximation of C{atan(B{x})}.
    '''
    a = fabs(x)
    if a < _1_0:
        r = fatan1(a) if a else _0_0
    elif a > _1_0:
        r = PI_2 - fatan1(_1_0 / a)  # == fatan2(a, _1_0)
    else:
        r = PI_4
    if x < 0:  # copysign0(r, x)
        r = -r
    return r


def fatan1(x):
    '''Fast approximation of C{atan(B{x})} for C{0 <= B{x} <= 1}, I{unchecked}.

       @see: U{ShaderFastLibs.h<https://GitHub.com/michaldrobot/ShaderFastLibs/
             blob/master/ShaderFastMathLib.h>} and U{Efficient approximations
             for the arctangent function<http://www-Labs.IRO.UMontreal.CA/
             ~mignotte/IFT2425/Documents/EfficientApproximationArctgFunction.pdf>},
             IEEE Signal Processing Magazine, 111, May 2006.
    '''
    # Eq (9): PI_4 * x - x * (abs(x) - 1) * (0.2447 + 0.0663 * abs(x)), for -1 < x < 1
    #         PI_4 * x - (x**2 - x) * (0.2447 + 0.0663 * x), for 0 < x - 1
    #         x * (1.0300981633974482 + x * (-0.1784 - x * 0.0663))
    H = Fhorner(x, _0_0, 1.0300982, -0.1784, -0.0663)
    return float(H)


def fatan2(y, x):
    '''Fast approximation of C{atan2(B{y}, B{x})}.

       @see: U{fastApproximateAtan(x, y)<https://GitHub.com/CesiumGS/cesium/blob/
             master/Source/Shaders/Builtin/Functions/fastApproximateAtan.glsl>}
             and L{fatan1}.
    '''
    b, a = fabs(y), fabs(x)
    if a < b:
        r = (PI_2 - fatan1(a / b)) if a else PI_2
    elif b < a:
        r = fatan1(b / a) if b else _0_0
    elif a:  # == b != 0
        r = PI_4
    else:  # a == b == 0
        return _0_0
    if x < 0:
        r = PI - r
    if y < 0:  # copysign0(r, y)
        r = -r
    return r


def favg(v1, v2, f=_0_5):
    '''Return the average of two values.

       @arg v1: One value (C{scalar}).
       @arg v2: Other value (C{scalar}).
       @kwarg f: Optional fraction (C{float}).

       @return: M{v1 + f * (v2 - v1)} (C{float}).
    '''
#      @raise ValueError: Fraction out of range.
#   '''
#   if not 0 <= f <= 1:  # XXX restrict fraction?
#       raise _ValueError(fraction=f)
    # v1 + f * (v2 - v1) == v1 * (1 - f) + v2 * f
    return fsum1_(v1, -f * v1, f * v2)


def fdot(a, *b):
    '''Return the precision dot product M{sum(a[i] * b[i] for
       i=0..len(a))}.

       @arg a: Iterable, list, tuple, etc. (C{scalar}s).
       @arg b: All positional arguments (C{scalar}s).

       @return: Dot product (C{float}).

       @raise LenError: Unequal C{len(B{a})} and C{len(B{b})}.

       @see: Class L{Fdot} and U{Algorithm 5.10 B{DotK}
             <https://www.TUHH.De/ti3/paper/rump/OgRuOi05.pdf>}.
    '''
    return fsum(_map_mul(a, b, fdot))


def fdot3(a, b, c, start=0):
    '''Return the precision dot product M{start +
       sum(a[i] * b[i] * c[i] for i=0..len(a))}.

       @arg a: Iterable, list, tuple, etc. (C{scalar}s).
       @arg b: Iterable, list, tuple, etc. (C{scalar}s).
       @arg c: Iterable, list, tuple, etc. (C{scalar}s).
       @kwarg start: Optional bias (C{scalar}).

       @return: Dot product (C{float}).

       @raise LenError: Unequal C{len(B{a})}, C{len(B{b})}
                        and/or C{len(B{c})}.

       @raise OverflowError: Partial C{2sum} overflow.
    '''
    def _mul3(a, b, c):  # map function
        return a * b * c

    def _muly(a, b, c, start):
        yield start
        for abc in map(_mul3, a, b, c):
            yield abc

    if not len(a) == len(b) == len(c):
        raise LenError(fdot3, a=len(a), b=len(b), c=len(c))

    return fsum(_muly(a, b, c, start) if start else map(_mul3, a, b, c))


def fhorner(x, *cs):
    '''Evaluate the polynomial M{sum(cs[i] * x**i for
       i=0..len(cs))} using the Horner form.

       @arg x: Polynomial argument (C{scalar}).
       @arg cs: Polynomial coeffients (C{scalar}s).

       @return: Horner value (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{x}}.

       @raise ValueError: No B{C{cs}} coefficients or B{C{x}} is not finite.

       @see: Function L{fpolynomial} and class L{Fhorner}.
    '''
    H = Fhorner(x, *cs)
    return float(H)


def fidw(xs, ds, beta=2):
    '''Interpolate using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW).

       @arg xs: Known values (C{scalar}s).
       @arg ds: Non-negative distances (C{scalar}s).
       @kwarg beta: Inverse distance power (C{int}, 0, 1, 2, or 3).

       @return: Interpolated value C{x} (C{float}).

       @raise LenError: Unequal or zero C{len(B{ds})} and C{len(B{xs})}.

       @raise ValueError: Invalid B{C{beta}}, negative B{C{ds}} value,
                          weighted B{C{ds}} below L{EPS}.

       @note: Using C{B{beta}=0} returns the mean of B{C{xs}}.
    '''
    n, xs = len2(xs)
    d, ds = len2(ds)
    if n != d or n < 1:
        raise LenError(fidw, xs=n, ds=d)

    d, x = min(zip(ds, xs))
    if d > EPS0 and n > 1:
        b = -Int_(beta=beta, low=0, high=3)
        if b < 0:
            ws = tuple(float(d)**b for d in ds)
            t =  fsum(_map_mul1(xs, ws))  # fdot(xs, *ws)
            x = _over(t, fsum(ws, floats=True))
        else:  # b == 0
            x = fsum(xs) / n  # fmean(xs)
    elif d < 0:  # PYCHOK no cover
        n = Fmt.SQUARE(ds=ds.index(d))
        raise _ValueError(n, d, txt=_negative_)
    return x


def fmean(xs):
    '''Compute the accurate mean M{sum(xs[i] for
       i=0..len(xs)) / len(xs)}.

       @arg xs: Values (C{scalar} or L{Fsum} instances).

       @return: Mean value (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise ValueError: No B{C{xs}} values.
    '''
    n, xs = len2(xs)
    if n > 0:
        return fsum(xs) / n  # if n > 1 else _2float(index=0, xs=xs[0])
    raise _ValueError(xs=xs)


def fmean_(*xs):
    '''Compute the accurate mean M{sum(xs[i] for
       i=0..len(xs)) / len(xs)}.

       @arg xs: Values (C{scalar} or L{Fsum} instances).

       @return: Mean value (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise ValueError: No B{C{xs}} values.
    '''
    return fmean(xs)


def fpolynomial(x, *cs, **over):
    '''Evaluate the polynomial M{sum(cs[i] * x**i for
       i=0..len(cs)) [/ over]}.

       @arg x: Polynomial argument (C{scalar}).
       @arg cs: Polynomial coeffients (C{scalar}s), all
                positional.
       @kwarg over: Optional, final divisor (C{scalar}

       @return: Polynomial value (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{x}}.

       @raise ValueError: No B{C{cs}} coefficients or B{C{x}} is not finite.

       @see: Function L{fhorner} and class L{Fpolynomial}.
    '''
    P =  Fpolynomial(x, *cs)
    d = _xkwds_get(over, over=0) if over else 0
    return P.fover(d) if d else float(P)


def fpowers(x, n, alts=0):
    '''Return a series of powers M{[x**i for i=1..n]}.

       @arg x: Value (C{scalar}).
       @arg n: Highest exponent (C{int}).
       @kwarg alts: Only alternating powers, starting with
                    this exponent (C{int}).

       @return: Powers of B{C{x}} (C{float}s or C{int}s).

       @raise TypeError: Non-scalar B{C{x}} or B{C{n}} not C{int}.

       @raise ValueError: Non-finite B{C{x}} or non-positive B{C{n}}.
    '''
    if not isint(n):
        raise _IsnotError(int.__name__, n=n)
    elif n < 1:
        raise _ValueError(n=n)

    p  = t = x if isint(x) else _2float(x=x)
    ps = [p]
    _a = ps.append
    for _ in range(1, n):
        p *= t
        _a(p)

    if alts > 0:  # x**2, x**4, ...
        # ps[alts-1::2] chokes PyChecker
        ps = ps[slice(alts-1, None, 2)]

    return ps


try:
    from math import prod as fprod  # Python 3.8
except ImportError:

    def fprod(xs, start=1):
        '''Iterable product, like C{math.prod} or C{numpy.prod}.

           @arg xs: Terms to be multiplied, an iterable, list,
                    tuple, etc. (C{scalar}s).
           @kwarg start: Initial term, also the value returned
                         for an empty B{C{xs}} (C{scalar}).

           @return: The product (C{float}).

           @see: U{NumPy.prod<https://docs.SciPy.org/doc/
                 numpy/reference/generated/numpy.prod.html>}.
        '''
        return freduce(_mul, xs, start)


def frange(start, number, step=1):
    '''Generate a range of C{float}s.

       @arg start: First value (C{float}).
       @arg number: The number of C{float}s to generate (C{int}).
       @kwarg step: Increment value (C{float}).

       @return: A generator (C{float}s).

       @see: U{NumPy.prod<https://docs.SciPy.org/doc/
             numpy/reference/generated/numpy.arange.html>}.
    '''
    if not isint(number):
        raise _IsnotError(int.__name__, number=number)
    for i in range(number):
        yield start + i * step


try:
    from functools import reduce as freduce
except ImportError:
    try:
        freduce = reduce  # PYCHOK expected
    except NameError:  # Python 3+

        def freduce(f, xs, *start):
            '''For missing C{functools.reduce}.
            '''
            if start:
                r = v = start[0]
            else:
                r, v = 0, MISSING
            for v in xs:
                r = f(r, v)
            if v is MISSING:
                raise _TypeError(xs=(), start=MISSING)
            return r


def fremainder(x, y):
    '''Remainder in range C{[-B{y / 2}, B{y / 2}]}.

       @arg x: Numerator (C{scalar}).
       @arg y: Modulus, denominator (C{scalar}).

       @return: Remainder (C{scalar}, preserving signed
                0.0) or C{NAN} for any non-finite B{C{x}}.

       @raise ValueError: Infinite or near-zero B{C{y}}.

       @see: I{Karney}'s U{Math.remainder<https://PyPI.org/
             project/geographiclib/>} and Python 3.7+
             U{math.remainder<https://docs.Python.org/3/
             library/math.html#math.remainder>}.
    '''
    # with Python 2.7.16 and 3.7.3 on macOS 10.13.6 and
    # with Python 3.10.2 on macOS 12.2.1 M1 arm64 native
    #  fmod( 0,   360) ==  0.0
    #  fmod( 360, 360) ==  0.0
    #  fmod(-0,   360) ==  0.0
    #  fmod(-0.0, 360) == -0.0
    #  fmod(-360, 360) == -0.0
    # however, using the % operator ...
    #    0   % 360 == 0
    #  360   % 360 == 0
    #  360.0 % 360 == 0.0
    #   -0   % 360 == 0
    # -360   % 360 == 0   == (-360)   % 360
    #   -0.0 % 360 == 0.0 ==   (-0.0) % 360
    # -360.0 % 360 == 0.0 == (-360.0) % 360

    # On Windows 32-bit with python 2.7, math.fmod(-0.0, 360)
    # == +0.0.  This fixes this bug.  See also Math::AngNormalize
    # in the C++ library, Math.sincosd has a similar fix.
    if _isfinite(x):
        try:
            r = remainder(x, y) if x else x
        except Exception as e:
            raise _xError(e, unstr(fremainder, x, y))
    else:  # handle x INF and NINF as NAN
        r = NAN
    return r


if _sys_version_info2 < (3, 8):  # PYCHOK no cover
    from math import hypot  # OK in Python 3.7-

    def hypot_(*xs):
        '''Compute the norm M{sqrt(sum(x**2 for x in xs))}.

           Similar to Python 3.8+ n-dimension U{math.hypot
           <https://docs.Python.org/3.8/library/math.html#math.hypot>},
           but exceptions, C{nan} and C{infinite} values are
           handled differently.

           @arg xs: X arguments (C{scalar}s), all positional.

           @return: Norm (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise ValueError: Invalid or no B{C{xs}} values.

           @note: The Python 3.8+ Euclidian distance U{math.dist
                  <https://docs.Python.org/3.8/library/math.html#math.dist>}
                  between 2 I{n}-dimensional points I{p1} and I{p2} can be
                  computed as M{hypot_(*((c1 - c2) for c1, c2 in zip(p1, p2)))},
                  provided I{p1} and I{p2} have the same, non-zero length I{n}.
        '''
        h, x2 = _h_x2(xs)
        return (h * sqrt(x2)) if x2 else _0_0

elif _sys_version_info2 < (3, 10):
    # In Python 3.8 and 3.9 C{math.hypot} is inaccurate, see
    # U{agdhruv<https://GitHub.com/geopy/geopy/issues/466>},
    # U{cffk<https://Bugs.Python.org/issue43088>} and module
    # U{geomath.py<https://PyPI.org/project/geographiclib/1.52>}

    def hypot(x, y):
        '''Compute the norm M{sqrt(x**2 + y**2)}.

           @arg x: X argument (C{scalar}).
           @arg y: Y argument (C{scalar}).

           @return: C{sqrt(B{x}**2 + B{y}**2)} (C{float}).
        '''
        if x:
            h = sqrt(x**2 + y**2) if y else fabs(x)
        elif y:
            h = fabs(y)
        else:
            h = _0_0
        return h

    from math import hypot as hypot_  # PYCHOK in Python 3.8 and 3.9
else:
    from math import hypot  # PYCHOK in Python 3.10+
    hypot_ = hypot


def _h_x2(xs):
    '''(INTERNAL) Helper for L{hypot_} and L{hypot2_}.
    '''
    if xs:
        n, xs = len2(xs)
        if n > 0:
            h = float(max(map(fabs, xs)))
            if h < EPS0:
                x2 = _0_0
            else:  # math.fsum, see C{_hypot21_} below
                x2 = _fsum(_x2_h2(_1_0, xs, h, _N_1_0))
            return h, x2

    raise _ValueError(xs=xs, txt=_too_(_few_))


def hypot1(x):
    '''Compute the norm M{sqrt(1 + x**2)}.

       @arg x: Argument (C{scalar}).

       @return: Norm (C{float}).
    '''
    return hypot(_1_0, x) if x else _1_0


def hypot2(x, y):
    '''Compute the I{squared} norm M{x**2 + y**2}.

       @arg x: X argument (C{scalar}).
       @arg y: Y argument (C{scalar}).

       @return: C{B{x}**2 + B{y}**2} (C{float}).
    '''
    if x:
        if y:
            if fabs(x) < fabs(y):
                x, y = y, x
            h2 = x**2 * ((y / x)**2 + _1_0)
        else:
            h2 = x**2
    elif y:
        h2 = y**2
    else:
        h2 = _0_0
    return h2


def hypot2_(*xs):
    '''Compute the I{squared} norm C{sum(x**2 for x in B{xs})}.

       @arg xs: X arguments (C{scalar}s), all positional.

       @return: Squared norm (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise ValueError: Invalid or no B{C{xs}} value.

       @see: Function L{hypot_}.
    '''
    h, x2 = _h_x2(xs)
    return (h**2 * x2) if x2 else _0_0


def _map_mul(a, b, where):
    '''(INTERNAL) Yield each B{C{a * b}}.
    '''
    n = len(b)
    if len(a) != n:  # PYCHOK no cover
        raise LenError(where, a=len(a), b=n)
    return map(_mul, a, b) if n > 3 else _map_mul1(a, b)


def _map_mul1(a, b):
    '''(INTERNAL) Yield each B{C{a * b}}, 1-primed.
    '''
    yield _1_0
    for ab in map(_mul, a, b):
        if ab:
            yield ab
    yield _N_1_0


def norm2(x, y):
    '''Normalize a 2-dimensional vector.

       @arg x: X component (C{scalar}).
       @arg y: Y component (C{scalar}).

       @return: 2-Tuple C{(x, y)}, normalized.

       @raise ValueError: Invalid B{C{x}} or B{C{y}}
              or zero norm.
    '''
    h = hypot(x, y)
    if not h:
        x = y = _0_0  # pass?
    elif not isnear1(h):
        try:
            x, y = x / h, y / h
        except Exception as e:
            raise _xError(e, x=x, y=y, h=h)
    return x, y


def norm_(*xs):
    '''Normalize all n-dimensional vector components.

       @arg xs: Components (C{scalar}s), all positional.

       @return: Yield each component, normalized.

       @raise ValueError: Invalid or insufficent B{C{xs}}
              or zero norm.
    '''
    h = hypot_(*xs)
    if h:
        try:
            for i, x in enumerate(xs):
                yield x / h
        except Exception as e:
            raise _xError(e, Fmt.SQUARE(xs=i), x, _h_, h)
    else:
        for _ in xs:
            yield _0_0


def _root(x, p, where):
    '''(INTERNAL) Raise C{x} to power C{0 < p < 1}.
    '''
    if x < 0:
        t = _SPACE_(_invokation_, where.__name__)
        raise _ValueError(unstr(t, x), txt=_negative_)
    return pow(x, p) if x else _0_0


def sqrt0(x):
    '''Return the square root iff C{B{x} >} L{EPS02}.

       @arg x: Value (C{scalar}).

       @return: Square root (C{float}) or C{0.0}.

       @note: Any C{B{x} < }L{EPS02} I{including} C{B{x} < 0}
              returns C{0.0}.
    '''
    return sqrt(x) if x > EPS02 else (_0_0 if x < EPS02 else EPS0)


def sqrt3(x):
    '''Return the square root, I{cubed} M{sqrt(x)**3} or M{sqrt(x**3)}.

       @arg x: Value (C{scalar}).

       @return: Square root I{cubed} (C{float}).

       @raise ValueError: Negative B{C{x}}.

       @see: Functions L{cbrt} and L{cbrt2}.
    '''
    return _root(x, _1_5, sqrt3)


def sqrt_a(h, b):
    '''Compute C{I{a}} side of a right-angled triangle from
       C{sqrt(B{h}**2 - B{b}**2)}.

       @arg h: Hypotenuse or outer annulus radius (C{scalar}).
       @arg b: Triangle side or inner annulus radius (C{scalar}).

       @return: C{copysign(I{a}, B{h})} or C{unsigned 0.0} (C{float}).

       @raise TypeError: Non-scalar B{C{h}} or B{C{b}}.

       @raise ValueError: If C{abs(B{h}) < abs(B{b})}.

       @see: Inner tangent chord B{I{d}} of an U{annulus
             <https://WikiPedia.org/wiki/Annulus_(mathematics)>}
             and function U{annulus_area<https://People.SC.FSU.edu/
             ~jburkardt/py_src/geometry/geometry.py>}.
    '''
    try:
        if not (isscalar(h) and isscalar(b)):
            raise TypeError(_not_scalar_)
        elif isnear0(h):  # PYCHOK no cover
            c, b = fabs(h), fabs(b)
            d = c - b
            if d < 0:
                raise ValueError('abs(h) < abs(b)')
            a = copysign0(sqrt((c + b) * d), h) if d > 0 else _0_0
        else:
            c =  float(h)
            s = _1_0 - (b / c)**2
            if s < 0:
                raise ValueError('abs(h) < abs(b)')
            a = (sqrt(s) * c) if 0 < s < 1 else (c if s else _0_0)
    except Exception as x:
        raise _xError(x, h=h, b=b)
    return a


def _x2_h2(s, xs, h, e):
    '''(INTERNAL) Yield M{(x / h)**2 for x in xs}.
    '''
    yield s
    if h in (_0_0, _1_0):
        for x in xs:
            if x:
                yield x**2
    else:
        for x in xs:
            if x:
                yield (x / h)**2
    yield e


def zcrt(x):
    '''Return the 6-th, I{zenzi-cubic} root, M{x**(1 / 6)}.

       @arg x: Value (C{scalar}).

       @return: I{Zenzi-cubic} root (C{float}).

       @see: Functions L{bqrt} and L{zqrt}.

       @raise ValueError: Negative B{C{x}}.
    '''
    return _root(x, _1_6th, zcrt)


def zqrt(x):
    '''Return the 8-th, I{zenzi-quartic} or I{squared-quartic} root, M{x**(1 / 8)}.

       @arg x: Value (C{scalar}).

       @return: I{Zenzi-quartic} root (C{float}).

       @see: Functions L{bqrt} and L{zcrt}.

       @raise ValueError: Negative B{C{x}}.
    '''
    return _root(x, _0_125, zqrt)

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
