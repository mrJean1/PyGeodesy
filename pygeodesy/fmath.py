
# -*- coding: utf-8 -*-

u'''Precision floating point functions, utilities and constants.

@newfield example: Example, Examples
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division

from pygeodesy.basics import copysign, isfinite, isint, isscalar, \
                             len2, _xcopy
from pygeodesy.errors import _IsnotError, LenError, _OverflowError, \
                             _TypeError, _ValueError
from pygeodesy.interns import EPS0, EPS1, MISSING, NN, PI, PI_2, PI_4, \
                             _EPS0__2, _finite_, _few_, _negative_,\
                             _not_, _singular_, _SPACE_, _too_, \
                             _0_0, _1_0, _1_5 as _3_2nd, _2_0, _3_0
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.streprs import Fmt, unstr
from pygeodesy.units import Int_

from math import hypot, sqrt  # pow
from operator import mul as _mul

__all__ = _ALL_LAZY.fmath
__version__ = '21.04.14'

# sqrt(2) <https://WikiPedia.org/wiki/Square_root_of_2>
_0_4142 =  0.414213562373095  # sqrt(_2_0) - _1_0
_1_3rd  = _1_0 / _3_0
_2_3rd  = _2_0 / _3_0


def _2even(s, r, p):
    '''(INTERNAL) Half-even rounding.
    '''
    if (r > 0 and p > 0) or \
       (r < 0 and p < 0):  # signs match
        t, p = _2sum(s, p * 2)
        if not p:
            s = t
    return s


def _2sum(a, b):  # by .testFmath
    '''(INTERNAL) Precision C{2sum} of M{a + b}.
    '''
    s = a + b
    if not isfinite(s):
        raise _OverflowError(unstr(_2sum.__name__, a, b), txt=str(s))
    if abs(a) < abs(b):
        a, b = b, a
    return s, b - (s - a)


class Fsum(object):
    '''Precision summation similar to standard Python function C{math.fsum}.

       Unlike C{math.fsum}, this class accumulates the values and provides
       intermediate, precision running sums.  Accumulation may continue
       after intermediate summations.

       @note: Handling of exceptions, C{inf}, C{INF}, C{nan} and C{NAN}
              values is different from C{math.fsum}.

       @see: U{Hettinger<https://GitHub.com/ActiveState/code/blob/master/recipes/Python/
             393090_Binary_floating_point_summatiaccurate_full/recipe-393090.py>},
             U{Kahan<https://WikiPedia.org/wiki/Kahan_summation_algorithm>},
             U{Klein<https://Link.Springer.com/article/10.1007/s00607-005-0139-x>},
             Python 2.6+ file I{Modules/mathmodule.c} and the issue log
             U{Full precision summation<https://Bugs.Python.org/issue2819>}.
    '''
    _fsum2_ = None
    _n      = 0
    _ps     = []

    def __init__(self, *starts):
        '''Initialize a new accumulator with one or more start values.

           @arg starts: No, one or more start values (C{scalar}s).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{starts}} value.

           @raise ValueError: Invalid or non-finite B{C{starts}} value.
        '''
        self._n = 0
        self._ps = []
        if starts:
            self.fadd(starts)

    def __add__(self, other):
        '''Sum of this and an other instance or a scalar.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: The sum, a new instance (L{Fsum}).

           @see: Method L{Fsum.__iadd__}.
        '''
        f = self.fcopy()
        f += other
        return f  # self.fcopy().__iadd__(other)

    def __iadd__(self, other):
        '''Add a scalar or an other instance to this instance.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.fadd}.
        '''
        if isscalar(other):
            self.fadd_(other)
        elif other is self:
            self.fmul(2)
        elif isinstance(other, Fsum):
            self.fadd(other._ps)
        else:
            raise _TypeError(_SPACE_(self, '+=', repr(other)))
        return self

    def __imul__(self, other):
        '''Multiply this instance by a scalar or an other instance.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.fmul}.
        '''
        if isscalar(other):
            self.fmul(other)
        elif isinstance(other, Fsum):
            ps = list(other._ps)  # copy
            if ps:
                s = self.fcopy()
                self.fmul(ps.pop())
                while ps:  # self += s * ps.pop()
                    p = s.fcopy()
                    p.fmul(ps.pop())
                    self.fadd(p._ps)
            else:
                self._ps = []  # zero
                self._fsum2_ = None
        else:
            raise _TypeError(_SPACE_(self, '*=', repr(other)))
        return self

    def __isub__(self, other):
        '''Subtract a scalar or an other instance from this instance.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.fadd}.
        '''
        if isscalar(other):
            self.fadd_(-other)
        elif other is self:
            self._ps = []  # zero
            self._fsum2_ = None
        elif isinstance(other, Fsum):
            self.fadd(-p for p in other._ps)
        else:
            raise _TypeError(_SPACE_(self, '-=', repr(other)))
        return self

    def __len__(self):
        '''Return the number of accumulated values.
        '''
        return self._n

    def __mul__(self, other):
        '''Product of this and an other instance or a scalar.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: The product, a new instance (L{Fsum}).

           @see: Method L{Fsum.__imul__}.
        '''
        f = self.fcopy()
        f *= other
        return f

    def __str__(self):
        from pygeodesy.named import classname
        return Fmt.PAREN(classname(self, prefixed=True), NN)

    def __sub__(self, other):
        '''Difference of this and an other instance or a scalar.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: The difference, a new instance (L{Fsum}).

           @see: Method L{Fsum.__isub__}.
        '''
        f = self.fcopy()
        f -= other
        return f

    def fadd(self, xs):
        '''Accumulate more values from an iterable.

           @arg xs: Iterable, list, tuple, etc. (C{scalar}s).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.
        '''
        if isscalar(xs):  # for backward compatibility
            xs = (xs,)

        ps = self._ps
        for n, x in enumerate(map(float, xs)):  # _iter()
            if not isfinite(x):
                n = Fmt.SQUARE(xs=n)
                raise _ValueError(n, x, txt=_not_(_finite_))
            i = 0
            for p in ps:
                x, p = _2sum(x, p)
                if p:
                    ps[i] = p
                    i += 1
            ps[i:] = [x]
        # assert self._ps is ps
        self._n += n + 1
        self._fsum2_ = None

    def fadd_(self, *xs):
        '''Accumulate more values from positional arguments.

           @arg xs: Values to add (C{scalar}s), all positional.

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.
        '''
        self.fadd(xs)

    def fcopy(self, deep=False):
        '''Copy this instance, C{shallow} or B{C{deep}}.

           @return: The copy, a new instance (L{Fsum}).
         '''
        f = _xcopy(self, deep=deep)
        # f._fsum2_ = self._fsum2_
        # f._n = self._n
        f._ps = list(self._ps)  # separate copy
        return f

    copy = fcopy

    def fmul(self, factor):
        '''Multiple the current, partial sum by a factor.

           @arg factor: The multiplier (C{scalar}).

           @raise TypeError: Non-scalar B{C{factor}}.

           @raise ValueError: Invalid or non-finite B{C{factor}}.

           @see: Method L{Fsum.fadd}.
        '''
        if not isfinite(factor):
            raise _ValueError(factor=factor, txt=_not_(_finite_))

        f, ps = float(factor), self._ps
        if ps:  # multiply and adjust partial sums
            ps[:] = [p * f for p in ps]
            self.fadd_(ps.pop())
            self._n -= 1
        # assert self._ps is ps

    def fsub(self, iterable):
        '''Accumulate more values from an iterable.

           @arg iterable: Sequence, list, tuple, etc. (C{scalar}s).

           @see: Method L{Fsum.fadd}.
        '''
        if iterable:
            self.fadd(-s for s in iterable)

    def fsub_(self, *xs):
        '''Accumulate more values from positional arguments.

           @arg xs: Values to subtract (C{scalar}s), all positional.

           @see: Method L{Fsum.fadd}.
        '''
        self.fsub(xs)

    def fsum(self, iterable=()):
        '''Accumulate more values from an iterable and sum all.

           @kwarg iterable: Sequence, list, tuple, etc. (C{scalar}s), optional.

           @return: Accurate, running sum (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{iterable}} value.

           @raise ValueError: Invalid or non-finite B{C{iterable}} value.

           @note: Accumulation can continue after summation.
        '''
        if iterable:
            self.fadd(iterable)

        ps = self._ps
        i = len(ps) - 1
        if i < 0:
            s = _0_0
        else:
            s = ps[i]
            while i > 0:
                i -= 1
                s, p = _2sum(s, ps[i])
                ps[i:] = [s]
                if p:  # sum(ps) became inexact
                    ps.append(p)
                    if i > 0:  # half-even round if signs match
                        s = _2even(s, ps[i-1], p)
                    break
            # assert self._ps is ps
        self._fsum2_ = s
        return s

    def fsum_(self, *xs):
        '''Accumulate more values from positional arguments and sum all.

           @arg xs: Values to add (C{scalar}s), all positional.

           @return: Accurate, running sum (C{float}).

           @see: Method L{Fsum.fsum}.

           @note: Accumulation can continue after summation.
        '''
        return self.fsum(xs)

    def fsum2_(self, *xs):
        '''Accumulate more values from positional arguments, sum all
           and provide the sum and delta.

           @arg xs: Values to add (C{scalar}s), all positional.

           @return: 2-Tuple C{(sum, delta)} with the accurate,
                    running C{sum} and the C{delta} with the
                    previous running C{sum}, both (C{float}).

           @see: Method L{Fsum.fsum_}.

           @note: Accumulation can continue after summation.
        '''
        p = self._fsum2_
        if p is None:
            p = self.fsum()
        s = self.fsum(xs)  # if xs else self._fsum2_
        return s, s - p


class Fdot(Fsum):
    '''Precision dot product.
    '''
    def __init__(self, a, *b):
        '''New L{Fdot} precision dot product M{sum(a[i] * b[i]
           for i=0..len(a))}.

           @arg a: List, sequence, tuple, etc. (C{scalar}s).
           @arg b: All positional arguments (C{scalar}s).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise LenError: Unequal C{len(B{a})} and C{len(B{b})}.

           @see: Function L{fdot} and method L{Fsum.fadd}.
        '''
        if len(a) != len(b):
            raise LenError(Fdot, a=len(a), b=len(b))

        Fsum.__init__(self)
        self.fadd(map(_mul, a, b))


class Fhorner(Fsum):
    '''Precision polynomial evaluation using the Horner form.
    '''
    def __init__(self, x, *cs):
        '''New L{Fhorner} evaluation of the polynomial
           M{sum(cs[i] * x**i for i=0..len(cs))}.

           @arg x: Polynomial argument (C{scalar}).
           @arg cs: Polynomial coeffients (C{scalar}[]).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{x}}.

           @raise ValueError: No B{C{cs}} coefficients or B{C{x}} is not finite.

           @see: Function L{fhorner} and methods L{Fsum.fadd} and L{Fsum.fmul}.
        '''
        if not isfinite(x):
            raise _ValueError(x=x, txt=_not_(_finite_))
        if not cs:
            raise _ValueError(cs=cs, txt=MISSING)

        x, cs = float(x), list(cs)

        Fsum.__init__(self, cs.pop())
        while cs:
            self.fmul(x)
            c = cs.pop()
            if c:
                self.fadd_(c)


class Fpolynomial(Fsum):
    '''Precision polynomial evaluation.
    '''
    def __init__(self, x, *cs):
        '''New L{Fpolynomial} evaluation of the polynomial
           M{sum(cs[i] * x**i for i=0..len(cs))}.

           @arg x: Polynomial argument (C{scalar}).
           @arg cs: Polynomial coeffients (C{scalar}[]).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{x}}.

           @raise ValueError: No B{C{cs}} coefficients or B{C{x}} is not finite.

           @see: Function L{fpolynomial} and method L{Fsum.fadd}.
        '''
        if not isfinite(x):
            raise _ValueError(x=x, txt=_not_(_finite_))
        if not cs:
            raise _ValueError(cs=cs, txt=MISSING)

        x, cs, xp = float(x), list(cs), _1_0

        Fsum.__init__(self, cs.pop(0))
        while cs:
            xp *= x
            c = cs.pop(0)
            if c:
                self.fadd_(xp * c)


def cbrt(x):
    '''Compute the cubic root M{x**(1/3)}.

       @arg x: Value (C{scalar}).

       @return: Cubic root (C{float}).

       @see: Functions L{cbrt2} and L{sqrt3}.
    '''
    # simpler and more accurate than Ken Turkowski's CubeRoot, see
    # <https://People.FreeBSD.org/~lstewart/references/apple_tr_kt32_cuberoot.pdf>
    return copysign(pow(abs(x), _1_3rd), x)


def cbrt2(x):
    '''Compute the cubic root squared M{x**(2/3)}.

       @arg x: Value (C{scalar}).

       @return: Cubic root squared (C{float}).

       @see: Functions L{cbrt} and L{sqrt3}.
    '''
    return pow(abs(x), _2_3rd)  # XXX pow(abs(x), _1_3rd)**2


def euclid(x, y):
    '''I{Appoximate} the norm M{sqrt(x**2 + y**2)} by
       M{max(abs(x), abs(y)) + min(abs(x), abs(y)) * 0.4142...}.

       @arg x: X component (C{scalar}).
       @arg y: Y component (C{scalar}).

       @return: Appoximate norm (C{float}).

       @see: Function L{euclid_}.
    '''
    x, y = abs(x), abs(y)
    if y > x:
        x, y = y, x
    return x + y * _0_4142  # XXX _0_5 before 20.10.02


def euclid_(*xs):
    '''I{Appoximate} the norm M{sqrt(sum(x**2 for x in xs))}
       by cascaded L{euclid}.

       @arg xs: X arguments, positional (C{scalar}[]).

       @return: Appoximate norm (C{float}).

       @see: Function L{euclid}.
    '''
    e = _0_0
    for x in sorted(map(abs, xs)):
        # e = euclid(x, e)
        if x > e:
            e, x = x, e
        e += x * _0_4142
    return e


def facos1(x):
    '''Fast approximation of L{acos1}C{(B{x})}.

       @see: U{ShaderFastLibs.h<https://GitHub.com/michaldrobot/
             ShaderFastLibs/blob/master/ShaderFastMathLib.h>}.
    '''
    a = abs(x)
    if a < EPS0:
        r = PI_2
    elif a < EPS1:
        r  = 1.5707288 - a * (0.2121144 - a * (0.0742610 - a * 0.0187293))
        r *= sqrt(_1_0 - a)
    else:
        r = _0_0
    return (PI - r) if x < 0 else r


def fasin1(x):  # PYCHOK no cover
    '''Fast approximation of L{asin1}C{(B{x})}.

       @see: L{facos1}.
    '''
    return PI_2 - facos1(x)


def fatan(x):
    '''Fast approximation of C{atan(B{x})}.
    '''
    a = abs(x)
    if a < _1_0:
        r = fatan1(a) if a else _0_0
    elif a > _1_0:
        r = PI_2 - fatan1(_1_0 / a)  # == fatan2(a, _1_0)
    else:
        r = PI_4
    return -r if x < 0 else r


def fatan1(x):
    '''Fast approximation of C{atan(B{x})} for C{0 <= B{x} <= 1}, I{unchecked}.

       @see: U{ShaderFastLibs.h<https://GitHub.com/michaldrobot/
             ShaderFastLibs/blob/master/ShaderFastMathLib.h>} and
             U{Efficient approximations for the arctangent function
             <http://www-Labs.IRO.UMontreal.CA/~mignotte/IFT2425/Documents/
             EfficientApproximationArctgFunction.pdf>}, IEEE Signal
             Processing Magazine, 111, May 2006.
    '''
    # Eq (9): PI_4 * x - x * (x - 1) * (0.2447 + 0.0663 * x**2)
    return x * (1.0300982 - x * (0.1784 + 0.0663 * x))  # w/o x**4


def fatan2(y, x):
    '''Fast approximation of C{atan2(B{y}, B{x})}.

       @see: U{fastApproximateAtan(x, y)<https://GitHub.com/CesiumGS/cesium/blob/
             master/Source/Shaders/Builtin/Functions/fastApproximateAtan.glsl>}
             and L{fatan1}.
    '''
    b, a = abs(y), abs(x)
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
    return -r if y < 0 else r


def favg(v1, v2, f=0.5):
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
    return v1 + f * (v2 - v1)  # v1 * (1 - f) + v2 * f


def fdot(a, *b):
    '''Return the precision dot product M{sum(a[i] * b[i] for
       i=0..len(a))}.

       @arg a: List, sequence, tuple, etc. (C{scalar}s).
       @arg b: All positional arguments (C{scalar}s).

       @return: Dot product (C{float}).

       @raise LenError: Unequal C{len(B{a})} and C{len(B{b})}.

       @see: Class L{Fdot}.
    '''
    if len(a) != len(b):
        raise LenError(fdot, a=len(a), b=len(b))

    return fsum(map(_mul, a, b))


def fdot3(a, b, c, start=0):
    '''Return the precision dot product M{start +
       sum(a[i] * b[i] * c[i] for i=0..len(a))}.

       @arg a: List, sequence, tuple, etc. (C{scalar}[]).
       @arg b: List, sequence, tuple, etc. (C{scalar}[]).
       @arg c: List, sequence, tuple, etc. (C{scalar}[]).
       @kwarg start: Optional bias (C{scalar}).

       @return: Dot product (C{float}).

       @raise LenError: Unequal C{len(B{a})}, C{len(B{b})}
                        and/or C{len(B{c})}.

       @raise OverflowError: Partial C{2sum} overflow.
    '''
    def _mul3(a, b, c):  # map function
        return a * b * c  # PYCHOK returns

    if not len(a) == len(b) == len(c):
        raise LenError(fdot3, a=len(a), b=len(b), c=len(c))

    if start:
        f = Fsum(start)
        return f.fsum(map(_mul3, a, b, c))
    else:
        return fsum(map(_mul3, a, b, c))


def fhorner(x, *cs):
    '''Evaluate the polynomial M{sum(cs[i] * x**i for
       i=0..len(cs))} using the Horner form.

       @arg x: Polynomial argument (C{scalar}).
       @arg cs: Polynomial coeffients (C{scalar}[]).

       @return: Horner value (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{x}}.

       @raise ValueError: No B{C{cs}} coefficients or B{C{x}} is not finite.

       @see: Function L{fpolynomial} and class L{Fhorner}.
    '''
    h = Fhorner(x, *cs)
    return h.fsum()


def fidw(xs, ds, beta=2):
    '''Interpolate using using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW).

       @arg xs: Known values (C{scalar}[]).
       @arg ds: Non-negative distances (C{scalar}[]).
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
            ds = tuple(d**b for d in ds)
            d  = fsum(ds)
            if abs(d) < EPS0:
                n = Fmt.PAREN(fsum='ds')
                raise _ValueError(n, d, txt=_singular_)
            x = fdot(xs, *ds) / d
        else:  # b == 0
            x = fsum(xs) / n  # fmean(xs)
    elif d < 0:
        n = Fmt.SQUARE(ds=ds.index(d))
        raise _ValueError(n, d, txt=_negative_)
    return x


def fmean(xs):
    '''Compute the accurate mean M{sum(xs[i] for
       i=0..len(xs)) / len(xs)}.

       @arg xs: Values (C{scalar}s).

       @return: Mean value (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise ValueError: No B{C{xs}} values.
    '''
    n, xs = len2(xs)
    if n > 0:
        return fsum(xs) / n
    raise _ValueError(xs=xs)


def fmean_(*xs):
    '''Compute the accurate mean M{sum(xs[i] for
       i=0..len(xs)) / len(xs)}.

       @arg xs: Values (C{scalar}s).

       @return: Mean value (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise ValueError: No B{C{xs}} values.
    '''
    return fmean(xs)


def fpolynomial(x, *cs):
    '''Evaluate the polynomial M{sum(cs[i] * x**i for
       i=0..len(cs))}.

       @arg x: Polynomial argument (C{scalar}).
       @arg cs: Polynomial coeffients (C{scalar}[]).

       @return: Polynomial value (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{x}}.

       @raise ValueError: No B{C{cs}} coefficients or B{C{x}} is not finite.

       @see: Function L{fhorner} and class L{Fpolynomial}.
    '''
    p = Fpolynomial(x, *cs)
    return p.fsum()


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
    if not isfinite(x):
        raise _ValueError(x=x, txt=_not_(_finite_))
    if not isint(n):
        raise _IsnotError(int.__name__, n=n)
    elif n < 1:
        raise _ValueError(n=n)

    xs = [x]
    for _ in range(1, n):
        xs.append(xs[-1] * x)

    if alts > 0:  # x**2, x**4, ...
        # xs[alts-1::2] chokes PyChecker
        xs = xs[slice(alts-1, None, 2)]

    return xs


try:
    from math import prod as fprod  # Python 3.8
except ImportError:

    def fprod(iterable, start=_1_0):
        '''Iterable product, like C{math.prod} or C{numpy.prod}.

           @arg iterable: Terms to be multiplied (C{scalar}[]).
           @kwarg start: Initial term, also the value returned
                         for an empty iterable (C{scalar}).

           @return: The product (C{float}).

           @see: U{NumPy.prod<https://docs.SciPy.org/doc/
                 numpy/reference/generated/numpy.prod.html>}.
        '''
        return freduce(_mul, iterable, start)


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
except ImportError:  # PYCHOK no cover
    try:
        freduce = reduce  # PYCHOK expected
    except NameError:  # Python 3+

        def freduce(f, iterable, *start):
            '''For missing C{functools.reduce}.
            '''
            if start:
                r = v = start[0]
            else:
                r, v = 0, MISSING
            for v in iterable:
                r = f(r, v)
            if v is MISSING:
                raise _TypeError(iterable=(), start=MISSING)
            return r


def fsum_(*xs):
    '''Precision summation of the positional argument vulues.

       @arg xs: Values to be added (C{scalar}[]).

       @return: Accurate L{fsum} (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{xs}} value.

       @raise ValueError: Invalid or non-finite B{C{xs}} value.
    '''
    return fsum(map(float, xs))


try:
    from math import fsum  # precision IEEE-754 sum, Python 2.6+

    # make sure fsum works as expected (XXX check
    # float.__getformat__('float')[:4] == 'IEEE'?)
    if fsum_(1, 1e101, 1, -1e101) != 2:
        del fsum  # nope, remove fsum ...
        raise ImportError  # ... use fsum below

except ImportError:  # PYCHOK no cover

    def fsum(iterable):
        '''Precision summation similar to standard Python function C{math.fsum}.

           Exception and I{non-finite} handling differs from C{math.fsum}.

           @arg iterable: Values to be added (C{scalar}[]).

           @return: Accurate C{sum} (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{iterable}} value.

           @raise ValueError: Invalid or non-finite B{C{iterable}} value.

           @see: Class L{Fsum}.
        '''
        f = Fsum()
        return f.fsum(iterable)


try:
    _ = hypot(1, 2, 3)  # new in Python 3.8+
    hypot_ = hypot
    del _
except TypeError:  # Python 3.7-

    def hypot_(*xs):
        '''Compute the norm M{sqrt(sum(x**2 for x in xs))}.

           @arg xs: X arguments, positional (C{scalar}[]).

           @return: Norm (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise ValueError: Invalid or no B{C{xs}} value.

           @see: Similar to Python 3.8+ U{math.hypot
                 <https://docs.Python.org/3.8/library/math.html#math.hypot>},
                 but handling of exceptions, C{nan} and C{infinite} values
                 is different.

           @note: The Python 3.8+ U{math.dist
                  <https://docs.Python.org/3.8/library/math.html#math.dist>}
                  Euclidian distance between 2 I{n}-dimensional points I{p1}
                  and I{p2} can be computed as M{hypot_(*((c1 - c2) for c1,
                  c2 in zip(p1, p2)))}, provided I{p1} and I{p2} have the
                  same, non-zero length I{n}.
        '''
        h, x2 = _h_x2(xs)
        return (h * sqrt(x2)) if x2 else _0_0


def _h_x2(xs):
    '''(INTERNAL) Helper for L{hypot_} and L{hypot2_}.
    '''
    if xs:
        n, xs = len2(xs)
        if n > 0:
            h = float(max(abs(x) for x in xs))
            if h > 0:
                if n > 1:
                    X = Fsum(_1_0)
                    X.fadd((x / h)**2 for x in xs)
                    x2 = X.fsum_(-_1_0)
                else:
                    x2 = _1_0
            else:
                h = x2 = _0_0
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

       @arg x: Argument (C{scalar}).
       @arg y: Argument (C{scalar}).

       @return: C{B{x}**2 + B{y}**2} (C{float}).
    '''
    x, y = x**2, y**2
    if x < y:
        x, y = y, x
    if y:  # and x
        x *= _1_0 + y / x
    return x


def hypot2_(*xs):
    '''Compute the I{squared} norm C{sum(x**2 for x in B{xs})}.

       @arg xs: X arguments, positional (C{scalar}[]).

       @return: Squared norm (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise ValueError: Invalid or no B{C{xs}} value.

       @see: Function L{hypot_}.
    '''
    h, x2 = _h_x2(xs)
    return (h**2 * x2) if x2 else _0_0


def sqrt0(x):
    '''Compute the square root iff C{B{x} > EPS0**2}.

       @arg x: Value (C{scalar}).

       @return: Square root (C{float}) or C{0.0}.
    '''
    return sqrt(x) if x > _EPS0__2 else _0_0


def sqrt3(x):
    '''Compute the square root, cubed M{sqrt(x)**3} or M{sqrt(x**3)}.

       @arg x: Value (C{scalar}).

       @return: Cubed square root (C{float}).

       @raise ValueError: Negative B{C{x}}.

       @see: Functions L{cbrt} and L{cbrt2}.
    '''
    if x < 0:
        raise _ValueError(x=x)
    return pow(x, _3_2nd) if x else _0_0

# **) MIT License
#
# Copyright (C) 2016-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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
