
# -*- coding: utf-8 -*-

u'''Utilities using precision floating point summation.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import _copysign, copysign0, isbool, isint, isscalar, \
                              len2, map1, _xiterable
from pygeodesy.constants import EPS0, EPS02, EPS1, NAN, PI, PI_2, PI_4, \
                               _0_0, _0_125, _1_6th, _0_25, _1_3rd, _0_5, _1_0, \
                               _N_1_0, _1_5, _copysign_0_0, _isfinite, remainder
from pygeodesy.errors import _IsnotError, LenError, _TypeError, _ValueError, \
                             _xError, _xkwds_get1, _xkwds_pop2
from pygeodesy.fsums import _2float, Fsum, fsum, fsum1_, _isFsumTuple, _1primed, \
                             Fmt, unstr
from pygeodesy.interns import MISSING, _negative_, _not_scalar_
from pygeodesy.lazily import _ALL_LAZY, _sys_version_info2
# from pygeodesy.streprs import Fmt, unstr  # from .fsums
from pygeodesy.units import Int_, _isHeight, _isRadius,  Float_  # PYCHOK for .heights

from math import fabs, sqrt  # pow
import operator as _operator  # in .datums, .trf, .utm

__all__ = _ALL_LAZY.fmath
__version__ = '24.05.29'

# sqrt(2) <https://WikiPedia.org/wiki/Square_root_of_2>
_0_4142  =  0.41421356237309504880  # ... sqrt(2) - 1
_2_3rd   = _1_3rd * 2
_h_lt_b_ = 'abs(h) < abs(b)'


class Fdot(Fsum):
    '''Precision dot product.
    '''
    def __init__(self, a, *b, **name_RESIDUAL):
        '''New L{Fdot} precision dot product M{sum(a[i] * b[i] for
           i=0..len(a)-1)}.

           @arg a: Iterable of values (each C{scalar} or an L{Fsum} or L{Fsum2Tuple}
                   instance).
           @arg b: Other values (each C{scalar} or an L{Fsum} or L{Fsum2Tuple} instance),
                   all positional.
           @kwarg name_RESIDUAL: Optional C{B{name}=NN} (C{str}) and the C{B{RESIDUAL}=0.0}
                       threshold (C{scalar}) for raising L{ResidualError}s, see class
                       L{Fsum<Fsum.__init__>}.

           @raise LenError: Unequal C{len(B{a})} and C{len(B{b})}.

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Invalid B{C{x}}.

           @raise ValueError: Non-finite B{C{x}}.

           @see: Function L{fdot} and method L{Fsum.fadd}.
        '''
        Fsum.__init__(self, **name_RESIDUAL)
        self.fadd(_map_mul(a, b, Fdot))


class Fhorner(Fsum):
    '''Precision polynomial evaluation using the Horner form.
    '''
    def __init__(self, x, *cs, **name_RESIDUAL):
        '''New L{Fhorner} evaluation of polynomial M{sum(cs[i] * x**i for
           i=0..len(cs)-1)}.

           @arg x: Polynomial argument (C{scalar} or an L{Fsum} or L{Fsum2Tuple}).
           @arg cs: Polynomial coeffients (each C{scalar} or an L{Fsum} or L{Fsum2Tuple}
                    instance), all positional.
           @kwarg name_RESIDUAL: Optional C{B{name}=NN} (C{str}) and the C{B{RESIDUAL}=0.0}
                       threshold (C{scalar}) for raising L{ResidualError}s, see class
                       L{Fsum<Fsum.__init__>}.

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Invalid B{C{x}}.

           @raise ValueError: Non-finite B{C{x}}.

           @see: Function L{fhorner} and methods L{Fsum.fadd} and L{Fsum.fmul}.
        '''
        Fsum.__init__(self, **name_RESIDUAL)
        if cs:
            self._fhorner(x, cs, Fhorner.__name__)
        else:
            self._fset_ps(_0_0)


class Fhypot(Fsum):
    '''Precision summation and hypotenuse, default C{root=2}.
    '''
    def __init__(self, *xs, **root_name_RESIDUAL_raiser):
        '''New L{Fhypot} hypotenuse of (the I{root} of) several components
           (raised to the power I{root}).

           @arg xs: Components (each C{scalar} or an L{Fsum} or L{Fsum2Tuple} instance),
                    all positional.
           @kwarg root_name_RESIDUAL_raiser: Optional, exponent and C{B{root}=2} order
                       (C{scalar}), C{B{name}=NN} (C{str}), the C{B{RESIDUAL}=0.0}
                       threshold (C{scalar}) and C{B{raiser}=True} (C{bool}) for
                       raising L{ResidualError}s, see class L{Fsum<Fsum.__init__>} and
                       method L{root<Fsum.root>}.
        '''
        r = None  # _xkwds_pop2 error
        try:
            r, kwds = _xkwds_pop2(root_name_RESIDUAL_raiser, root=2)
            r, kwds = _xkwds_pop2(kwds, power=r)  # for backward compatibility
            raiser  = _Fsum__init__(self, **kwds)
            if xs:
                self._facc_power(r, xs, Fhypot, **raiser)
            self._fset(self.root(r, **raiser))
        except Exception as X:
            raise self._ErrorXs(X, xs, root=r)


class Fpolynomial(Fsum):
    '''Precision polynomial evaluation.
    '''
    def __init__(self, x, *cs, **name_RESIDUAL):
        '''New L{Fpolynomial} evaluation of the polynomial
           M{sum(cs[i] * x**i for i=0..len(cs)-1)}.

           @arg x: Polynomial argument (C{scalar} or an L{Fsum} or L{Fsum2Tuple}).
           @arg cs: Polynomial coeffients (each C{scalar} or an L{Fsum} or L{Fsum2Tuple}
                    instance), all positional.
           @kwarg name_RESIDUAL: Optional C{B{name}=NN} (C{str}) and the C{B{RESIDUAL}=0.0}
                       threshold (C{scalar}) for raising L{ResidualError}s, see class
                       L{Fsum<Fsum.__init__>}.

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Invalid B{C{x}}.

           @raise ValueError: Non-finite B{C{x}}.

           @see: Class L{Fhorner}, function L{fpolynomial} and method L{Fsum.fadd}.
        '''
        Fsum.__init__(self, *cs[:1], **name_RESIDUAL)
        n = len(cs) - 1
        if n > 0:
            self.fadd(_1map_mul(cs[1:], _powers(x, n)))
        elif n < 0:
            self._fset_ps(_0_0)


class Fpowers(Fsum):
    '''Precision summation of powers, optimized for C{power=2, 3 and 4}.
    '''
    def __init__(self, power, *xs, **name_RESIDUAL_raiser):
        '''New L{Fpowers} sum of (the I{power} of) several bases.

           @arg power: The exponent (C{scalar} or an L{Fsum} or L{Fsum2Tuple}).
           @arg xs: One or more bases (each C{scalar} or an L{Fsum} or L{Fsum2Tuple} instance),
                    all positional.
           @kwarg name_RESIDUAL_raiser: Optional C{B{name}=NN} (C{str}), the C{B{RESIDUAL}=0.0}
                       threshold (C{scalar}) and C{B{raiser}=True} (C{bool}) for raising
                       L{ResidualError}s, see class L{Fsum<Fsum.__init__>} and method
                       L{fpow<Fsum.fpow>}.
        '''
        try:
            raiser = _Fsum__init__(self, **name_RESIDUAL_raiser)
            if xs:
                self._facc_power(power, xs, Fpowers, **raiser)  # x**0 == 1
        except Exception as X:
            raise self._ErrorXs(X, xs, power=power)


class Froot(Fsum):
    '''The root of a precision summation.
    '''
    def __init__(self, root, *xs, **name_RESIDUAL_raiser):
        '''New L{Froot} root of a precision sum.

           @arg root: The order (C{scalar} or an L{Fsum} or L{Fsum2Tuple}), non-zero.
           @arg xs: Items to summate (each a C{scalar} or an L{Fsum} or L{Fsum2Tuple} instance),
                    all positional.
           @kwarg name_RESIDUAL_raiser: Optional C{B{name}=NN} (C{str}), the C{B{RESIDUAL}=0.0}
                       threshold (C{scalar}) and C{B{raiser}=True} (C{bool}) for raising
                       L{ResidualError}s, see class L{Fsum<Fsum.__init__>} and method
                       L{fpow<Fsum.fpow>}.
        '''
        try:
            raiser = _Fsum__init__(self, **name_RESIDUAL_raiser)
            if xs:
                self.fadd(xs)
            self._fset(self.root(root, **raiser))
        except Exception as X:
            raise self._ErrorXs(X, xs, root=root)


class Fcbrt(Froot):
    '''Cubic root of a precision summation.
    '''
    def __init__(self, *xs, **name_RESIDUAL_raiser):
        '''New L{Fcbrt} cubic root of a precision sum.

           @see: Class L{Froot} for further details.
        '''
        Froot.__init__(self, 3, *xs, **name_RESIDUAL_raiser)


class Fsqrt(Froot):
    '''Square root of a precision summation.
    '''
    def __init__(self, *xs, **name_RESIDUAL_raiser):
        '''New L{Fsqrt} square root of a precision sum.

           @see: Class L{Froot} for further details.
        '''
        Froot.__init__(self, 2, *xs, **name_RESIDUAL_raiser)


def _Fsum__init__(inst, raiser=MISSING, **name_RESIDUAL):
    '''(INTERNAL) Init an C{F...} instance above.
    '''
    Fsum.__init__(inst, **name_RESIDUAL)  # PYCHOK self
    inst._fset_ps(_0_0)
    return {} if raiser is MISSING else dict(raiser=raiser)


def bqrt(x):
    '''Return the 4-th, I{bi-quadratic} or I{quartic} root, M{x**(1 / 4)},
       preserving C{type(B{x})}.

       @arg x: Value (C{scalar} or an L{Fsum} or L{Fsum2Tuple}).

       @return: I{Quartic} root (C{float} or an L{Fsum}).

       @raise TypeeError: Invalid B{C{x}}.

       @raise ValueError: Negative B{C{x}}.

       @see: Functions L{zcrt} and L{zqrt}.
    '''
    return _root(x, _0_25, bqrt)


try:
    from math import cbrt as _cbrt  # Python 3.11+

except ImportError:  # Python 3.10-

    def _cbrt(x):
        '''(INTERNAL) Compute the I{signed}, cube root M{x**(1/3)}.
        '''
        # <https://archive.lib.MSU.edu/crcmath/math/math/r/r021.htm>
        # simpler and more accurate than Ken Turkowski's CubeRoot, see
        # <https://People.FreeBSD.org/~lstewart/references/apple_tr_kt32_cuberoot.pdf>
        return _copysign(pow(fabs(x), _1_3rd), x)  # to avoid complex


def cbrt(x):
    '''Compute the cube root M{x**(1/3)}, preserving C{type(B{x})}.

       @arg x: Value (C{scalar} or an L{Fsum} or L{Fsum2Tuple}).

       @return: Cubic root (C{float} or L{Fsum}).

       @see: Functions L{cbrt2} and L{sqrt3}.
    '''
    if _isFsumTuple(x):
        r = abs(x).fpow(_1_3rd)
        if x.signOf() < 0:
            r = -r
    else:
        r = _cbrt(x)
    return r  # cbrt(-0.0) == -0.0


def cbrt2(x):  # PYCHOK attr
    '''Compute the cube root I{squared} M{x**(2/3)}, preserving C{type(B{x})}.

       @arg x: Value (C{scalar} or an L{Fsum} or L{Fsum2Tuple}).

       @return: Cube root I{squared} (C{float} or L{Fsum}).

       @see: Functions L{cbrt} and L{sqrt3}.
    '''
    return abs(x).fpow(_2_3rd) if _isFsumTuple(x) else _cbrt(x**2)


def euclid(x, y):
    '''I{Appoximate} the norm M{sqrt(x**2 + y**2)} by
       M{max(abs(x), abs(y)) + min(abs(x), abs(y)) * 0.4142...}.

       @arg x: X component (C{scalar} or L{Fsum} instance).
       @arg y: Y component (C{scalar} or L{Fsum} instance).

       @return: Appoximate norm (C{float} or L{Fsum}).

       @see: Function L{euclid_}.
    '''
    x, y = abs(x), abs(y)  # NOT fabs!
    if y > x:
        x, y = y, x
    return x + y * _0_4142  # XXX * _0_5 before 20.10.02


def euclid_(*xs):
    '''I{Appoximate} the norm M{sqrt(sum(x**2 for x in xs))} by
       cascaded L{euclid}.

       @arg xs: X arguments (each C{scalar} or an L{Fsum}
                instance), all positional.

       @return: Appoximate norm (C{float} or L{Fsum}).

       @see: Function L{euclid}.
    '''
    e = _0_0
    for x in sorted(map(abs, xs)):  # NOT fabs, reverse=True!
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
        H  = Fhorner(-a, 1.5707288, 0.2121144, 0.0742610, 0.0187293)
        H *= Fsqrt(_1_0, -a)
        r  = float(H)
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
    '''Fast approximation of C{atan(B{x})} for C{0 <= B{x} < 1}, I{unchecked}.

       @see: U{ShaderFastLibs.h<https://GitHub.com/michaldrobot/ShaderFastLibs/
             blob/master/ShaderFastMathLib.h>} and U{Efficient approximations
             for the arctangent function<http://www-Labs.IRO.UMontreal.CA/
             ~mignotte/IFT2425/Documents/EfficientApproximationArctgFunction.pdf>},
             IEEE Signal Processing Magazine, 111, May 2006.
    '''
    # Eq (9): PI_4 * x - x * (abs(x) - 1) * (0.2447 + 0.0663 * abs(x)), for -1 < x < 1
    #         PI_4 * x - (x**2 - x) * (0.2447 + 0.0663 * x), for 0 < x < 1
    #         x * (1.0300981633974482 + x * (-0.1784 - x * 0.0663))
    H = Fhorner(x, _0_0, 1.0300981634, -0.1784, -0.0663)
    return float(H)


def fatan2(y, x):
    '''Fast approximation of C{atan2(B{y}, B{x})}.

       @see: U{fastApproximateAtan(x, y)<https://GitHub.com/CesiumGS/cesium/blob/
             master/Source/Shaders/Builtin/Functions/fastApproximateAtan.glsl>}
             and L{fatan1}.
    '''
    a, b = fabs(x), fabs(y)
    if b > a:
        r = (PI_2 - fatan1(a / b)) if a else PI_2
    elif a > b:
        r = fatan1(b / a) if b else _0_0
    elif a:  # a == b != 0
        r = PI_4
    else:  # a == b == 0
        return _0_0
    if x < 0:
        r = PI - r
    if y < 0:  # copysign0(r, y)
        r = -r
    return r


def favg(a, b, f=_0_5):
    '''Return the precision average of two values.

       @arg a: One (C{scalar} or an L{Fsum} or L{Fsum2Tuple}).
       @arg b: Other (C{scalar} or an L{Fsum} or L{Fsum2Tuple}).
       @kwarg f: Optional fraction (C{float}).

       @return: M{a + f * (b - a)} (C{float}).
    '''
#      @raise ValueError: Fraction out of range.
#   '''
#   if not 0 <= f <= 1:  # XXX restrict fraction?
#       raise _ValueError(fraction=f)
    # a + f * (b - a) == a * (1 - f) + b * f
    return fsum1_(a, a * (-f), b * f)


def fdot(a, *b):
    '''Return the precision dot product M{sum(a[i] * b[i] for
       i=0..len(a))}.

       @arg a: Iterable of values (each C{scalar}).
       @arg b: Other values (each C{scalar}), all positional.

       @return: Dot product (C{float}).

       @raise LenError: Unequal C{len(B{a})} and C{len(B{b})}.

       @see: Class L{Fdot} and U{Algorithm 5.10 B{DotK}
             <https://www.TUHH.De/ti3/paper/rump/OgRuOi05.pdf>}.
    '''
    return fsum(_map_mul(a, b, fdot))


def fdot3(xs, ys, zs, start=0):
    '''Return the precision dot product M{start +
       sum(a[i] * b[i] * c[i] for i=0..len(a)-1)}.

       @arg xs: Iterable (each C{scalar} or an L{Fsum} or
                L{Fsum2Tuple} instance).
       @arg ys: Iterable (each C{scalar} or an L{Fsum} or
                L{Fsum2Tuple} instance).
       @arg zs: Iterable (each C{scalar} or an L{Fsum} or
                L{Fsum2Tuple} instance).
       @kwarg start: Optional bias (C{scalar} or an L{Fsum}
                     or L{Fsum2Tuple}).

       @return: Dot product (C{float}).

       @raise LenError: Unequal C{len(B{xs})}, C{len(B{ys})}
                        and/or C{len(B{zs})}.

       @raise OverflowError: Partial C{2sum} overflow.
    '''
    def _mul3(xs, ys, zs, s, p):
        if s:
            yield s
        if p:
            yield _1_0
        _F = Fsum
        for x, y, z in zip(xs, ys, zs):
            yield (_F(x) * y) * z
        if p:
            yield _N_1_0

    n = len(xs)
    if not n == len(ys) == len(zs):
        raise LenError(fdot3, xs=n, ys=len(ys), zs=len(zs))

    return fsum(_mul3(xs, ys, zs, start, n < 4))


def fhorner(x, *cs):
    '''Evaluate the polynomial M{sum(cs[i] * x**i for
       i=0..len(cs)-1)} using the Horner form.

       @return: Horner sum (C{float}).

       @see: Class L{Fhorner} for further details.
    '''
    H = Fhorner(x, *cs)
    return float(H)


def fidw(xs, ds, beta=2):
    '''Interpolate using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW).

       @arg xs: Known values (each C{scalar} or an L{Fsum} or
                L{Fsum2Tuple} instance).
       @arg ds: Non-negative distances (each C{scalar} or an L{Fsum}
                or L{Fsum2Tuple} instance).
       @kwarg beta: Inverse distance power (C{int}, 0, 1, 2, or 3).

       @return: Interpolated value C{x} (C{float}).

       @raise LenError: Unequal or zero C{len(B{ds})} and C{len(B{xs})}.

       @raise TypeError: An invalid B{C{ds}} or B{C{xs}}.

       @raise ValueError: Invalid B{C{beta}}, negative B{C{ds}} or
                          weighted B{C{ds}} below L{EPS}.

       @note: Using C{B{beta}=0} returns the mean of B{C{xs}}.
    '''
    n, xs = len2(xs)
    if n > 1:
        b = -Int_(beta=beta, low=0, high=3)
        if b < 0:
            try:  # weighted
                _F =  Fsum
                W  = _F()
                X  = _F()
                for i, d in enumerate(_xiterable(ds)):
                    x = xs[i]
                    D = _F(d)
                    if D < EPS0:
                        if D < 0:
                            raise ValueError(_negative_)
                        x = float(x)
                        i = n
                        break
                    if D.fpow(b):
                        W += D
                        X += D.fmul(x)
                else:
                    x  = X.fover(W, raiser=False)
                    i += 1  # len(xs) >= len(ds)
            except IndexError:
                i += 1  # len(xs) < i < len(ds)
            except Exception as X:
                _I = Fmt.INDEX
                raise _xError(X, _I(xs=i), x, _I(ds=i), d)
        else:  # b == 0
            x = fsum(xs) / n  # fmean(xs)
            i = n
    elif n:
        x = float(xs[0])
        i = n
    else:
        x    = _0_0
        i, _ =  len2(ds)
    if i != n:
        raise LenError(fidw, xs=n, ds=i)
    return x


def fmean(xs):
    '''Compute the accurate mean M{sum(xs) / len(xs)}.

       @arg xs: Values (C{scalar} or L{Fsum} instances).

       @return: Mean value (C{float}).

       @raise LenError: No B{C{xs}} values.

       @raise OverflowError: Partial C{2sum} overflow.
    '''
    n, xs = len2(xs)
    if n < 1:
        raise LenError(fmean, xs=xs)
    return Fsum(*xs).fover(n) if n > 1 else _2float(index=0, xs=xs[0])


def fmean_(*xs):
    '''Compute the accurate mean M{sum(xs) / len(xs)}.

       @see: Function L{fmean} for further details.
    '''
    return fmean(xs)


def fpolynomial(x, *cs, **over):
    '''Evaluate the polynomial M{sum(cs[i] * x**i for
       i=0..len(cs)) [/ over]}.

       @kwarg over: Optional final, I{non-zero} divisor (C{scalar}).

       @return: Polynomial value (C{float}).

       @see: Class L{Fpolynomial} for further details.
    '''
    P =  Fpolynomial(x, *cs)
    d = _xkwds_get1(over, over=0) if over else 0
    return P.fover(d) if d else float(P)


def fpowers(x, n, alts=0):
    '''Return a series of powers M{[x**i for i=1..n]}.

       @arg x: Value (C{scalar} or an L{Fsum} or L{Fsum2Tuple}).
       @arg n: Highest exponent (C{int}).
       @kwarg alts: Only alternating powers, starting with this
                    exponent (C{int}).

       @return: Tuple of powers of B{C{x}} (each C{type(B{x})}).

       @raise TypeError: Invalid B{C{x}} or B{C{n}} not C{int}.

       @raise ValueError: Non-finite B{C{x}} or invalid B{C{n}}.
    '''
    if not isint(n):
        raise _IsnotError(int.__name__, n=n)
    elif n < 1:
        raise _ValueError(n=n)

    p  = x if isint(x) or _isFsumTuple(x) else _2float(x=x)
    ps = tuple(_powers(p, n))

    if alts > 0:  # x**2, x**4, ...
        # ps[alts-1::2] chokes PyChecker
        ps = ps[slice(alts-1, None, 2)]

    return ps


try:
    from math import prod as fprod  # Python 3.8
except ImportError:

    def fprod(xs, start=1):
        '''Iterable product, like C{math.prod} or C{numpy.prod}.

           @arg xs: Iterable of values to be multiplied (each
                    C{scalar} or an L{Fsum}).
           @kwarg start: Initial value, also the value returned
                         for an empty B{C{xs}} (C{scalar}).

           @return: The product (C{float} or an L{Fsum}).

           @see: U{NumPy.prod<https://docs.SciPy.org/doc/
                 numpy/reference/generated/numpy.prod.html>}.
        '''
        return freduce(_operator.mul, xs, start)


def frandoms(n, seeded=None):
    '''Generate C{n} (long) lists of random C{floats}.

       @arg n: Number of lists to generate (C{int}, non-negative).
       @kwarg seeded: If C{scalar}, use C{random.seed(B{seeded})} or
                      if C{True}, seed using today's C{year-day}.

       @see: U{Hettinger<https://GitHub.com/ActiveState/code/tree/master/recipes/
             Python/393090_Binary_floating_point_summatiaccurate_full/recipe-393090.py>}.
    '''
    from random import gauss, random, seed, shuffle

    if seeded is None:
        pass
    elif seeded and isbool(seeded):
        from time import localtime
        seed(localtime().tm_yday)
    elif isscalar(seeded):
        seed(seeded)

    c = (7, 1e100, -7, -1e100, -9e-20, 8e-20) * 7
    for _ in range(n):
        s  = 0
        t  = list(c)
        _a = t.append
        for _ in range(n * 8):
            v = gauss(0, random())**7 - s
            _a(v)
            s += v
        shuffle(t)
        yield t


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
        yield start + (step * i)


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
        return float(Fhypot(*xs, raiser=False))

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
        return float(Fhypot(x, y, raiser=False))

    from math import hypot as hypot_  # PYCHOK in Python 3.8 and 3.9
else:
    from math import hypot  # PYCHOK in Python 3.10+
    hypot_ = hypot


def hypot1(x):
    '''Compute the norm M{sqrt(1 + x**2)}.

       @arg x: Argument (C{scalar} or L{Fsum} or L{Fsum2Tuple}).

       @return: Norm (C{float}).
    '''
    if _isFsumTuple(x):
        h = float(Fhypot(_1_0, x)) if x else _1_0
    else:
        h = hypot(_1_0, x) if x else _1_0
    return h


def hypot2(x, y):
    '''Compute the I{squared} norm M{x**2 + y**2}.

       @arg x: X (C{scalar} or L{Fsum} or L{Fsum2Tuple}).
       @arg y: Y (C{scalar} or L{Fsum} or L{Fsum2Tuple}).

       @return: C{B{x}**2 + B{y}**2} (C{float}).
    '''
    x, y = map1(abs, x, y)  # NOT fabs!
    if y > x:
        x, y = y, x
    if x:
        h2 = x**2
        if y:
            h2 *= (y / x)**2 + _1_0
        h2 =  float(h2)
    else:
        h2 = _0_0
    return h2


def hypot2_(*xs):
    '''Compute the I{squared} norm C{fsum(x**2 for x in B{xs})}.

       @arg xs: Components (each C{scalar} or an L{Fsum} or
                L{Fsum2Tuple} instance), all positional.

       @return: Squared norm (C{float}).

       @see: Class L{Fpowers} for further details.
    '''
    h2 = float(max(map(abs, xs))) if xs else _0_0
    if h2:
        _h = _1_0 / h2
        h2 =  Fpowers(2, *((x * _h) for x in xs))
        h2 =  h2.fover(_h**2)
    return h2


def _map_mul(xs, ys, where):
    '''(INTERNAL) Yield each B{C{x * y}}.
    '''
    n = len(ys)
    if len(xs) != n:  # PYCHOK no cover
        raise LenError(where, xs=len(xs), ys=n)
    return _1map_mul(xs, ys) if n < 4 else map(
           _operator.mul, map(Fsum, xs), ys)


def _1map_mul(xs, ys):
    '''(INTERNAL) Yield each B{C{x * y}}, 1-primed.
    '''
    return _1primed(map(_operator.mul, map(Fsum, xs), ys))


def norm2(x, y):
    '''Normalize a 2-dimensional vector.

       @arg x: X component (C{scalar}).
       @arg y: Y component (C{scalar}).

       @return: 2-Tuple C{(x, y)}, normalized.

       @raise ValueError: Invalid B{C{x}} or B{C{y}}
              or zero norm.
    '''
    try:
        h = hypot(x, y)
        if h:
            x, y = (x / h), (y / h)
        else:
            x = _copysign_0_0(x)  # pass?
            y = _copysign_0_0(y)
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
    try:
        i  = x = h = None
        h  = hypot_(*xs)
        _h = (_1_0 / h) if h else _0_0
        for i, x in enumerate(xs):
            yield x * _h
    except Exception as X:
        raise _xError(X, Fmt.SQUARE(xs=i), x, h=h)


def _powers(x, n):
    '''(INTERNAL) Yield C{x**i for i=1..n}.
    '''
    p = 1  # type(p) == type(x)
    for _ in range(n):
        p *= x
        yield p


def _root(x, p, where):
    '''(INTERNAL) Raise C{x} to power C{0 < p < 1}.
    '''
    try:
        if x > 0:
            return Fsum(x).fpow(p).as_iscalar
        elif x < 0:
            raise ValueError(_negative_)
    except Exception as X:
        raise _xError(X, unstr(where, x))
    return _0_0


def sqrt0(x, Error=None):
    '''Return the square root C{sqrt(B{x})} iff C{B{x} > }L{EPS02},
       preserving C{type(B{x})}.

       @arg x: Value (C{scalar} or an L{Fsum} or L{Fsum2Tuple}).
       @kwarg Error: Error to raise for negative B{C{x}}.

       @return: Square root (C{float} or L{Fsum}) or C{0.0}.

       @raise TypeeError: Invalid B{C{x}}.

       @note: Any C{B{x} < }L{EPS02} I{including} C{B{x} < 0}
              returns C{0.0}.
    '''
    if Error and x < 0:
        raise Error(unstr(sqrt0, x))
    return _root(x, _0_5, sqrt0) if x > EPS02 else (_0_0 if x < EPS02 else EPS0)


def sqrt3(x):
    '''Return the square root, I{cubed} M{sqrt(x)**3} or M{sqrt(x**3)},
       preserving C{type(B{x})}.

       @arg x: Value (C{scalar} or an L{Fsum} or L{Fsum2Tuple}).

       @return: Square root I{cubed} (C{float} or L{Fsum}).

       @raise TypeeError: Invalid B{C{x}}.

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
        if not (_isHeight(h) and _isRadius(b)):
            raise TypeError(_not_scalar_)
        c = fabs(h)
        if c > EPS0:
            s = _1_0 - (b / c)**2
            if s < 0:
                raise ValueError(_h_lt_b_)
            a = (sqrt(s) * c) if 0 < s < 1 else (c if s else _0_0)
        else:  # PYCHOK no cover
            b = fabs(b)
            d = c - b
            if d < 0:
                raise ValueError(_h_lt_b_)
            d *= c + b
            a  = sqrt(d) if d else _0_0
    except Exception as x:
        raise _xError(x, h=h, b=b)
    return copysign0(a, h)


def zcrt(x):
    '''Return the 6-th, I{zenzi-cubic} root, M{x**(1 / 6)},
       preserving C{type(B{x})}.

       @arg x: Value (C{scalar} or an L{Fsum} or L{Fsum2Tuple}).

       @return: I{Zenzi-cubic} root (C{float} or L{Fsum}).

       @see: Functions L{bqrt} and L{zqrt}.

       @raise TypeeError: Invalid B{C{x}}.

       @raise ValueError: Negative B{C{x}}.
    '''
    return _root(x, _1_6th, zcrt)


def zqrt(x):
    '''Return the 8-th, I{zenzi-quartic} or I{squared-quartic} root,
       M{x**(1 / 8)}, preserving C{type(B{x})}.

       @arg x: Value (C{scalar} or an L{Fsum} or L{Fsum2Tuple}).

       @return: I{Zenzi-quartic} root (C{float} or L{Fsum}).

       @see: Functions L{bqrt} and L{zcrt}.

       @raise TypeeError: Invalid B{C{x}}.

       @raise ValueError: Negative B{C{x}}.
    '''
    return _root(x, _0_125, zqrt)

# **) MIT License
#
# Copyright (C) 2016-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
