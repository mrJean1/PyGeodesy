
# -*- coding: utf-8 -*-

u'''Precision floating point functions, utilities and constants.

@newfield example: Example, Examples
'''
# make sure int/int division yields float quotient
from __future__ import division
division = 1 / 2  # double check int division, see .datum.py, .utily.py
if not division:
    raise ImportError('%s 1/2 == %d' % ('division', division))
del division

from pygeodesy.lazily import _ALL_LAZY

from math import acos, copysign, hypot, isinf, isnan, sqrt  # pow
from operator import mul
from sys import float_info as _float_info

# all public contants, classes and functions
__all__ = _ALL_LAZY.fmath
__version__ = '19.09.17'

try:  # Luciano Ramalho, "Fluent Python", page 395, O'Reilly, 2016
    from numbers import Integral as _Ints  #: (INTERNAL) Int objects
except ImportError:
    try:  # _Ints imported by .utily
        _Ints = int, long  #: (INTERNAL) Int objects (C{tuple})
    except NameError:  # Python 3+
        _Ints = int,  #: (INTERNAL) Int objects (C{tuple})

try:  # similarly ...
    from numbers import Real as _Scalars  #: (INTERNAL) Scalar objects
except ImportError:
    try:
        _Scalars = int, long, float  #: (INTERNAL) Scalar objects (C{tuple})
    except NameError:
        _Scalars = int, float  #: (INTERNAL) Scalar objects (C{tuple})

try:  # _Seqs imported by .utily
    from collections import Sequence as _Seqs  #: (INTERNAL) incl MutableSequence
except ImportError:
    _Seqs = list, tuple, range  # XXX also set?

try:
    EPS    = _float_info.epsilon   #: System's epsilon (C{float})
    MANTIS = _float_info.mant_dig  #: System's mantissa bits (C{int})
    MAX    = _float_info.max       #: System's float max (C{float})
    MIN    = _float_info.min       #: System's float min (C{float})
except AttributeError:
    EPS    = 2.220446049250313e-16  #: Epsilon (C{float}) 2**-52?
    MANTIS = 53  #: Mantissa bits ≈53 (C{int})
    MAX    = pow(2.0,  1023) * (2 - EPS)  #: Float max (C{float}) ≈10**308, 2**1024?
    MIN    = pow(2.0, -1021)  # Float min (C{float}) ≈10**-308, 2**-1021?
EPS_2  = EPS / 2      #: M{EPS / 2}   ≈1.110223024625e-16 (C{float})
EPS1   = 1.0 - EPS    #: M{1 - EPS}   ≈0.9999999999999998 (C{float})
EPS1_2 = 1.0 - EPS_2  #: M{1 - EPS_2} ≈0.9999999999999999 (C{float})
# _1EPS  = 1.0 + EPS  #: M{1 + EPS}   ≈1.0000000000000002 (C{float})

INF  = float('inf')  #: Infinity (C{float}), see C{isinf}, C{isfinite}
NAN  = float('nan')  #: Not-A-Number (C{float}), see C{isnan}
NEG0 = -0.0          #: Negative 0.0 (C{float}), see C{isneg0}

_1_3rd = 1.0 / 3.0  #: (INTERNAL) One third (C{float})
_2_3rd = 2.0 / 3.0  #: (INTERNAL) Two thirds (C{float})
_3_2nd = 3.0 / 2.0  #: (INTERNAL) Three halfs (C{float})


def _IsNotError(*names, **pair):  # Error=TypeError
    '''(INTERNAL) Format a C{TypeError} for a C{name=value} pair.
    '''
    Error = pair.pop('Error', TypeError)
    for n, v in pair.items():
        break
    else:
        n, v = 'pair', 'N/A'
    t = ' or ' .join(names)
    return Error('%s not %s: %r' % (n, t, v))


def _2even(s, r, p):
    '''(INTERNAL) Half-even rounding.
    '''
    if (r > 0 and p > 0) or \
       (r < 0 and p < 0):  # signs match
        t, p = _2sum(s, p * 2)
        if not p:
            s = t
    return s


def _2sum(a, b):
    '''(INTERNAL) Precision C{2sum} of M{a + b}.
    '''
    s = a + b
    if not isfinite(s):
        raise OverflowError('%s: %r' % ('2sum', s))
    if abs(a) < abs(b):
        a, b = b, a
    return s, b - (s - a)


class Fsum(object):
    '''Precision summation similar to standard Python function C{math.fsum}.

       Unlike C{math.fsum}, this class accumulates the values repeatedly
       and provides intermediate, precision running sums.  Accumulation
       may continue after intermediate summations.

       @note: Handling of exceptions, C{nan} and C{finite} values is
              different from C{math.fsum}.

       @see: U{Hettinger<https://code.ActiveState.com/recipes/393090>},
             U{Kahan<https://WikiPedia.org/wiki/Kahan_summation_algorithm>},
             U{Klein<https://Link.Springer.com/article/10.1007/s00607-005-0139-x>},
             Python 2.6+ file I{Modules/mathmodule.c} and the issue log
             U{Full precision summation<https://Bugs.Python.org/issue2819>}.
    '''
    _fsum2_ = None

    def __init__(self, *starts):
        '''Initialize a new accumulator with one or more start values.

           @param starts: No, one or more start values (C{scalar}s).

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

           @param other: L{Fsum} instance or C{scalar}.

           @return: The sum, a new instance (L{Fsum}).

           @see: Method L{Fsum.__iadd__}.
        '''
        f = self.fcopy()
        f += other
        return f  # self.fcopy().__iadd__(other)

    def __iadd__(self, other):
        '''Add a scalar or an other instance to this instance.

           @param other: L{Fsum} instance or C{scalar}.

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
            raise TypeError('%s += %r' % (self, other))
        return self

    def __imul__(self, other):
        '''Multiply this instance by a scalar or an other instance.

           @param other: L{Fsum} instance or C{scalar}.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.fmul}.
        '''
        if isscalar(other):
            self.fmul(other)
        elif isinstance(other, Fsum):
            ps = list(other._ps)
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
            raise TypeError('%s *= %r' % (self, other))
        return self

    def __isub__(self, other):
        '''Subtract a scalar or an other instance from this instance.

           @param other: L{Fsum} instance or C{scalar}.

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
            raise TypeError('%s -= %r' % (self, other))
        return self

    def __len__(self):
        '''Return the number of accumulated values.
        '''
        return self._n

    def __mul__(self, other):
        '''Product of this and an other instance or a scalar.

           @param other: L{Fsum} instance or C{scalar}.

           @return: The product, a new instance (L{Fsum}).

           @see: Method L{Fsum.__imul__}.
        '''
        f = self.fcopy()
        f *= other
        return f

    def __str__(self):
        # m = self.__module__.split('.')[-1]
        return '%s.%s()' % (self.__module__, self.__class__.__name__)

    def __sub__(self, other):
        '''Difference of this and an other instance or a scalar.

           @param other: L{Fsum} instance or C{scalar}.

           @return: The difference, a new instance (L{Fsum}).

           @see: Method L{Fsum.__isub__}.
        '''
        f = self.fcopy()
        f -= other
        return f

    def fadd(self, iterable):
        '''Accumulate more values from an iterable.

           @param iterable: Sequence, list, tuple, etc. (C{scalar}s).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{iterable}} value.

           @raise ValueError: Invalid or non-finite B{C{iterable}} value.
        '''
        if isscalar(iterable):  # for backward compatibility
            iterable = tuple(iterable)

#       def _iter():
#           for a in iterable:
#               if isinstance(a, Fsum):
#                   if a is self:
#                       self.fmul(2)
#                   else:
#                       for a in a._ps:
#                           yield a
#               else:
#                   yield a

        ps = self._ps
        for a in iterable:  # _iter()
            if not isfinite(a):
                raise ValueError('not %s: %r' %('finite', a))
            i = 0
            for p in ps:
                a, p = _2sum(a, p)
                if p:
                    ps[i] = p
                    i += 1
            ps[i:] = [a]
            self._n += 1
        # assert self._ps is ps
        self._fsum2_ = None

    def fadd_(self, *xs):
        '''Accumulate more values from positional arguments.

           @param xs: Values to add (C{scalar}s), all positional.

           @see: Method L{Fsum.fadd}.
        '''
        self.fadd(xs)

    def fcopy(self):
        '''Copy this instance.

           @return: The copy, a new instance (L{Fsum}).
         '''
        f = Fsum()
        f._n  = self._n
        f._ps = list(self._ps)
        return f

    __copy__ = fcopy

    def fmul(self, factor):
        '''Multiple the current, partial sum by a factor.

           @param factor: The multiplier (C{scalar}).

           @raise TypeError: Non-scalar B{C{factor}}.

           @raise ValueError: Invalid or non-finite B{C{factor}}.

           @see: Method L{Fsum.fadd}.
        '''
        if not isfinite(factor):
            raise ValueError('not %s: %r' %('finite', factor))

        ps = self._ps
        if ps:  # multiply and adjust partial sums
            ps[:] = [p * factor for p in ps]
            self.fadd_(ps.pop())
            self._n -= 1
        # assert self._ps is ps

    def fsub(self, iterable):
        '''Accumulate more values from an iterable.

           @param iterable: Sequence, list, tuple, etc. (C{scalar}s).

           @see: Method L{Fsum.fadd}.
        '''
        if iterable:
            self.fadd(-s for s in iterable)

    def fsub_(self, *xs):
        '''Accumulate more values from positional arguments.

           @param xs: Values to subtract (C{scalar}s), all positional.

           @see: Method L{Fsum.fadd}.
        '''
        self.fsub(xs)

    def fsum(self, iterable=()):
        '''Accumulate more values from an iterable and sum all.

           @keyword iterable: Sequence, list, tuple, etc. (C{scalar}s), optional.

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
            s = 0.0
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

           @param xs: Values to add (C{scalar}s), all positional.

           @return: Accurate, running sum (C{float}).

           @see: Method L{Fsum.fsum}.

           @note: Accumulation can continue after summation.
        '''
        return self.fsum(xs)

    def fsum2_(self, *xs):
        '''Accumulate more values from positional arguments, sum all
           and provide the sum and delta.

           @param xs: Values to add (C{scalar}s), all positional.

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

           @param a: List, sequence, tuple, etc. (C{scalar}s).
           @param b: All positional arguments (C{scalar}s).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise ValueError: Unequal C{len}(B{a}) and C{len}(B{b}).

           @see: Function L{fdot} and method L{Fsum.fadd}.
        '''
        if not len(a) == len(b):
            raise ValueError('%s: %s vs %s' % ('len', len(a), len(b)))

        Fsum.__init__(self)
        self.fadd(map(mul, a, b))


class Fhorner(Fsum):
    '''Precision polynomial evaluation using the Horner form.
    '''
    def __init__(self, x, *cs):
        '''New L{Fhorner} evaluation of the polynomial
           M{sum(cs[i] * x**i for i=0..len(cs))}.

           @param x: Polynomial argument (C{scalar}).
           @param cs: Polynomial coeffients (C{scalar}[]).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{x}}.

           @raise ValueError: No B{C{cs}} coefficients or B{C{x}} is not finite.

           @see: Function L{fhorner} and methods L{Fsum.fadd} and L{Fsum.fmul}.
        '''
        if not isfinite(x):
            raise ValueError('not %s: %r' %('finite', x))
        if not cs:
            raise ValueError('no %s: %r' % ('coefficents', cs))

        x, cs = float(x), list(cs)

        Fsum.__init__(self, cs.pop())
        while cs:
            self.fmul(x)
            self.fadd_(cs.pop())


class Fpolynomial(Fsum):
    '''Precision polynomial evaluation.
    '''
    def __init__(self, x, *cs):
        '''New L{Fpolynomial} evaluation of the polynomial
           M{sum(cs[i] * x**i for i=0..len(cs))}.

           @param x: Polynomial argument (C{scalar}).
           @param cs: Polynomial coeffients (C{scalar}[]).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{x}}.

           @raise ValueError: No B{C{cs}} coefficients or B{C{x}} is not finite.

           @see: Function L{fpolynomial} and method L{Fsum.fadd}.
        '''
        if not isfinite(x):
            raise ValueError('not %s: %r' %('finite', x))
        if not cs:
            raise ValueError('no %s: %r' % ('coefficents', cs))

        x, cs, xp = float(x), list(cs), 1

        Fsum.__init__(self, cs.pop(0))
        while cs:
            xp *= x
            self.fadd_(xp * cs.pop(0))


def acos1(x):
    '''Return M{math.acos(max(-1, min(1, x)))}.
    '''
    return acos(max(-1.0, min(1.0, x)))


def cbrt(x):
    '''Compute the cubic root M{x**(1/3)}.

       @param x: Value (C{scalar}).

       @return: Cubic root (C{float}).

       @see: Functions L{cbrt2} and L{sqrt3}.
    '''
    # simpler and more accurate than Ken Turkowski's CubeRoot, see
    # <https://People.FreeBSD.org/~lstewart/references/apple_tr_kt32_cuberoot.pdf>
    return copysign(pow(abs(x), _1_3rd), x)


def cbrt2(x):
    '''Compute the cubic root squared M{x**(2/3)}.

       @param x: Value (C{scalar}).

       @return: Cubic root squared (C{float}).

       @see: Functions L{cbrt} and L{sqrt3}.
    '''
    return pow(abs(x), _2_3rd)


def favg(v1, v2, f=0.5):
    '''Return the weighted average of two values.

       @param v1: One value (C{scalar}).
       @param v2: Other value (C{scalar}).
       @keyword f: Optional fraction (C{float}).

       @return: M{v1 + f * (v2 - v1)} (C{float}).
    '''
#      @raise ValueError: Fraction out of range.
#   '''
#   if not 0 <= f <= 1:  # XXX restrict fraction?
#       raise ValueError('%s invalid: %r' % ('fraction', f))
    return v1 + f * (v2 - v1)  # v1 * (1 - f) + v2 * f


def fdot(a, *b):
    '''Return the precision dot product M{sum(a[i] * b[i] for
       i=0..len(a))}.

       @param a: List, sequence, tuple, etc. (C{scalar}s).
       @param b: All positional arguments (C{scalar}s).

       @return: Dot product (C{float}).

       @raise ValueError: Unequal C{len(B{a})} and C{len(B{b})}.

       @see: Class L{Fdot}.
    '''
    if not len(a) == len(b):
        raise ValueError('%s: %s vs %s' % ('len', len(a), len(b)))

    return fsum(map(mul, a, b))


def fdot3(a, b, c, start=0):
    '''Return the precision dot product M{start +
       sum(a[i] * b[i] * c[i] for i=0..len(a))}.

       @param a: List, sequence, tuple, etc. (C{scalar}[]).
       @param b: List, sequence, tuple, etc. (C{scalar}[]).
       @param c: List, sequence, tuple, etc. (C{scalar}[]).
       @keyword start: Optional bias (C{scalar}).

       @return: Dot product (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise ValueError: Unequal C{len(B{a})}, C{len(B{b})}
                          and/or C{len(B{c})}.
    '''
    def _mul3(a, b, c):  # map function
        return a * b * c  # PYCHOK returns

    if not len(a) == len(b) == len(c):
        raise ValueError('%s: %s vs %s vs %s' % ('len', len(a), len(b), len(c)))

    if start:
        f = Fsum(start)
        return f.fsum(map(_mul3, a, b, c))
    else:
        return fsum(map(_mul3, a, b, c))


def fhorner(x, *cs):
    '''Evaluate the polynomial M{sum(cs[i] * x**i for
       i=0..len(cs))} using the Horner form.

       @param x: Polynomial argument (C{scalar}).
       @param cs: Polynomial coeffients (C{scalar}[]).

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

       @param xs: Known values (C{scalar}[]).
       @param ds: Non-negative distances (C{scalar}[]).
       @keyword beta: Inverse distance power (C{int}, 0, 1, 2, or 3).

       @return: Interpolated value C{x} (C{float}).

       @raise ValueError: Invalid B{C{beta}}, negative B{C{ds}} value,
                          weighted B{C{ds}} below L{EPS} or unequal
                          C{len(B{ds})} and C{len(B{xs})}.

       @note: Using C{B{beta}=0} returns the mean of B{C{xs}}.
    '''
    n, xs = len2(xs)
    d, ds = len2(ds)
    if n != d or n < 1:
        raise ValueError('%s: %s vs %s' % ('len', n, d))

    d, x = min(zip(ds, xs))
    if d > EPS and n > 1:
        b = -int(beta)
        if -4 < b < 0:  # and b == -beta
            ds = tuple(d**b for d in ds)
            d = fsum(ds)
            if d < EPS:
                raise ValueError('%s[%s] invalid: %r' % ('ds', '', d))
            x = fdot(xs, *ds) / d
        elif b == 0:
            x = fmean(xs)
        else:
            raise ValueError('%s=%r invalid' % ('beta', beta))
    elif d < 0:
        i = ds.index(d)
        raise ValueError('%s[%s] invalid: %r' % ('ds', i, d))
    return x


def fmean(xs):
    '''Compute the accurate mean M{sum(xs[i] for
       i=0..len(xs)) / len(xs)}.

       @param xs: Values (C{scalar}s).

       @return: Mean value (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise ValueError: No B{C{xs}} values.
    '''
    n, xs = len2(xs)
    if n > 0:
        return fsum(xs) / n
    raise ValueError('no %s: %r' % ('xs', xs))


def fpolynomial(x, *cs):
    '''Evaluate the polynomial M{sum(cs[i] * x**i for
       i=0..len(cs))}.

       @param x: Polynomial argument (C{scalar}).
       @param cs: Polynomial coeffients (C{scalar}[]).

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

       @param x: Value (C{scalar}).
       @param n: Highest exponent (C{int}).
       @keyword alts: Only alternating powers, starting
                      with this exponent (C{int}).

       @return: Powers of B{C{x}} (C{float}[]).

       @raise TypeError: Non-scalar B{C{x}} or B{C{n}} not C{int}.

       @raise ValueError: Non-finite B{C{x}} or non-positive B{C{n}}.
    '''
    if not isfinite(x):
        raise ValueError('not %s: %r' %('finite', x))
    if not isinstance(n, _Ints):
        raise _IsNotError(int.__name_, n=n)
    elif n < 1:
        raise ValueError('%s invalid: %r' % ('n', n))

    xs = [x]
    for _ in range(1, n):
        xs.append(xs[-1] * x)

    if alts > 0:  # x**2, x**4, ...
        # XXX PyChecker chokes on xs[alts-1::2]
        xs = xs[slice(alts-1, None, 2)]

    # XXX PyChecker claims result is None
    return xs


try:
    from math import prod as fprod  # Python 3.8
except ImportError:

    def fprod(iterable, start=1.0):
        '''Iterable product, like C{math.prod} or C{numpy.prod}.

           @param iterable: Values to be multiplied (C{scalar}[]).
           @keyword start: Initial product, also the value returned
                           for an empty iterable (C{scalar}).

           @return: The product (C{float}).

           @see: U{NumPy.prod<https://docs.SciPy.org/doc/
                 numpy/reference/generated/numpy.prod.html>}.
        '''
        return freduce(mul, iterable, start)


def frange(start, number, step=1):
    '''Generate a range of C{float}s.

       @param start: First value (C{float}).
       @param number: The number of C{float}s to generate (C{int}).
       @keyword step: Increment value (C{float}).

       @return: A generator (C{float}s).

       @see: U{NumPy.prod<https://docs.SciPy.org/doc/
             numpy/reference/generated/numpy.arange.html>}.
    '''
    if not isinstance(number, _Ints):
        raise _IsNotError(int.__name_, number=number)
    for i in range(number):
        yield start + i * step


try:
    from functools import reduce as freduce
except ImportError:
    try:
        freduce = reduce  # PYCHOK expected
    except NameError:  # Python 3+
        _EMPTY = object()

        def freduce(f, iterable, *start):
            '''For missing C{functools.reduce}.
            '''
            if start:
                r = v = start[0]
            else:
                r, v = 0, _EMPTY
            for v in iterable:
                r = f(r, v)
            if v is _EMPTY:
                raise TypeError('%s() empty, no start' % ('freduce',))
            return r


def fStr(floats, prec=6, sep=', ', fmt='%.*f', ints=False):
    '''Convert floats to string, optionally with trailing zero
       decimals stripped.

       @param floats: List, sequence, tuple, etc. (C{scalar}s).
       @keyword prec: Optional precision, number of decimal digits (0..9).
                      Trailing zero decimals are stripped for B{C{prec}} values
                      of 1 and above, but kept for negative B{C{prec}} values.
       @keyword sep: Optional, separator to join (string).
       @keyword fmt: Optional, float format (string).
       @keyword ints: Optionally, remove decimal dot (C{bool}).

       @return: The floats as 'f, f, ... f' (string).
    '''
    def _fstr(p, f):
        t = fmt % (abs(p), float(f))
        if ints and (isint(f, both=True) or  # for ...
                     # corner case testLcc lon0=-96.0
                     t.rstrip('0').endswith('.')):
            t = t.split('.')[0]
        elif p > 1:
            t = fStrzs(t)
        return t  # PYCHOK returns

    if isscalar(floats):
        return _fstr(prec, floats)
    else:
        return sep.join(_fstr(prec, f) for f in floats)


def fStrzs(fstr):
    '''Strip trailing zero decimals from a float string.

       @param fstr: Float (string).

       @return: Float (string).
    '''
    if fstr.endswith('0'):
        z = fstr.find('.') + 2  # keep 1st zero decimal
        if z > 1 and fstr[z:].isdigit():  # don't strip 'e+0..'
            fstr = fstr[:z] + fstr[z:].rstrip('0')
    return fstr


def fsum_(*xs):
    '''Precision summation of the positional argument vulues.

       @param xs: Values to be added (C{scalar}[]).

       @return: Accurate L{fsum} (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{xs}} value.

       @raise ValueError: Invalid or non-finite B{C{xs}} value.
    '''
    return fsum(xs)


try:
    from math import fsum  # precision IEEE-754 sum, Python 2.6+

    # make sure fsum works as expected (XXX check
    # float.__getformat__('float')[:4] == 'IEEE'?)
    if fsum_(1, 1e101, 1, -1e101) != 2:
        del fsum  # no, remove fsum ...
        raise ImportError  # ... use fsum below

except ImportError:

    def fsum(iterable):
        '''Precision summation similar to standard Python function C{math.fsum}.

           Exception and I{non-finite} handling differs from C{math.fsum}.

           @param iterable: Values to be added (C{scalar}[]).

           @return: Accurate C{sum} (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{iterable}} value.

           @raise ValueError: Invalid or non-finite B{C{iterable}} value.

           @see: Class L{Fsum}.
        '''
        f = Fsum()
        return f.fsum(iterable)


def hypot1(x):
    '''Compute the norm M{sqrt(1 + x**2)}.

       @param x: Argument (C{scalar}).

       @return: Norm (C{float}).
    '''
    return hypot(1.0, x)


try:
    _ = hypot(1, 2, 3)  # new in Python 3.8+
    hypot_ = hypot
    del _
except TypeError:  # Python 3.7-

    def hypot_(*xs):
        '''Compute the norm M{sqrt(sum(xs[i]**2)) for i=0..len(xs)}.

           @param xs: X arguments, positional (C{scalar}[]).

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
        if xs:
            n, xs = len2(xs)
            if n > 0:
                h = float(max(abs(x) for x in xs))
                if h > 0 and n > 1:
                    X = Fsum(1.0)
                    X.fadd((x / h)**2 for x in xs)
                    h *= sqrt(X.fsum_(-1.0))
                return h
        else:
            n = ''
        raise ValueError('%s(): %r[%s]' % ('hypot_', xs, n))


try:
    from math import isfinite  # new in Python 3+
except ImportError:

    def isfinite(obj):
        '''Check for C{Inf} and C{NaN} values.

           @param obj: Value (C{scalar}).

           @return: C{False} if B{C{obj}} is C{INF} or C{NAN},
                    C{True} otherwise.

           @raise TypeError: Non-scalar B{C{obj}}.
        '''
        if not isscalar(obj):
            raise _IsNotError(isscalar.__name__, obj=obj)
        return not (isinf(obj) or isnan(obj))


def isint(obj, both=False):
    '''Check for integer type or integer value.

       @param obj: The object (any C{type}).
       @keyword both: Optionally, check both type and value (C{bool}).

       @return: C{True} if B{C{obj}} is C{int}, C{False} otherwise.
    '''
    if both and isinstance(obj, float):  # NOT _Scalars!
        try:
            return obj.is_integer()
        except AttributeError:
            return False  # XXX float(int(obj)) == obj?
    return isinstance(obj, _Ints)


def isneg0(obj):
    '''Check for NEG0, negative 0.0.

       @param obj: Value (C{scalar}).

       @return: C{True} if B{C{obj}} is C{NEG0} or -0.0,
                C{False} otherwise.
    '''
    return obj in (0.0, NEG0) and copysign(1, obj) < 0
#                             and str(obj).rstrip('0') == '-0.'


def isscalar(obj):
    '''Check for scalar types.

       @param obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is C{scalar}, C{False} otherwise.
    '''
    return isinstance(obj, _Scalars)


def len2(seq):
    '''Make built-in function L{len} work for generators, iterators,
       etc. since those can only be started exactly once.

       @param seq: Generator, iterator, list, range, tuple, etc.

       @return: 2-Tuple (number, ...) of items (C{int}, C{list} or
                C{range} or C{tuple}).
    '''
    if not isinstance(seq, _Seqs):  # not hasattr(seq, '__len__'):
        seq = list(seq)
    return len(seq), seq


def map1(func, *xs):
    '''Apply each argument to a single-argument function and
       return a tuple of results.

       @param func: Function to apply (C{callable}).
       @param xs: Arguments to apply (C{any positional}).

       @return: Function results (C{tuple}).
    '''
    return tuple(map(func, xs))


def map2(func, *xs):
    '''Apply arguments to a function and return a tuple of results.

       Unlike Python 2's built-in L{map}, Python 3+ L{map} returns a
       L{map} object, an iterator-like object which generates the
       results only once.  Converting the L{map} object to a tuple
       maintains Python 2 behavior.

       @param func: Function to apply (C{callable}).
       @param xs: Arguments to apply (C{list, tuple, ...}).

       @return: Function results (C{tuple}).
    '''
    return tuple(map(func, *xs))


def scalar(value, low=EPS, high=1.0, name='scalar', Error=ValueError):
    '''Validate a scalar.

       @param value: The value (C{scalar}).
       @keyword low: Optional lower bound (C{scalar}).
       @keyword high: Optional upper bound (C{scalar}).
       @keyword name: Optional name of value (C{str}).
       @keyword Error: Exception to raise (C{ValueError}).

       @return: New value (C{type} of B{C{low}}).

       @raise TypeError: Non-scalar B{C{value}}.

       @raise ValueError: Out-of-bounds B{C{value}}.
    '''
    if not isscalar(value):
        raise _IsNotError('scalar', **{name: value})
    try:
        if low is None:
            v = float(value)
        else:
            v = type(low)(value)
            if v < low or v > high:
                raise ValueError
    except (TypeError, ValueError):
        raise Error(str(_IsNotError('valid', **{name: value})))
    return v


def sqrt3(x):
    '''Compute the square root cubed M{sqrt(x)**3} or M{sqrt(x**3)}.

       @param x: Value (C{scalar}).

       @return: Cubed square root (C{float}).

       @raise ValueError: Negative B{C{x}}.

       @see: Functions L{cbrt} and L{cbrt2}.
    '''
    if x < 0:
        raise ValueError('%s(%r)' % ('sqrt3', x))
    return pow(x, _3_2nd)

# **) MIT License
#
# Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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
