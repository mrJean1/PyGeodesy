
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

from pygeodesy.basics import EPS, isfinite, isint, isscalar, \
                             len2, _xcopy
from pygeodesy.errors import _IsnotError, LenError, _OverflowError, \
                             _TypeError, _ValueError
from pygeodesy.interns import _beta_, _item_ps, _item_sq, _Missing, \
                               NN, _too_few_
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.streprs import unstr
from pygeodesy.units import Int_

from math import copysign, hypot, sqrt  # pow
from operator import mul as _mul_

__all__ = _ALL_LAZY.fmath
__version__ = '20.07.19'

_not_finite_ = 'not finite'

_1_3rd = 1 / 3.0  #: (INTERNAL) One third (C{float})
_2_3rd = 2 / 3.0  # PYCHOK exported to .datum
_3_2nd = 3 / 2.0  #: (INTERNAL) Three halfs (C{float})


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

       @see: U{Hettinger<https://code.ActiveState.com/recipes/393090>},
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
            raise _TypeError('%s += %r' % (self, other))
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
            raise _TypeError('%s *= %r' % (self, other))
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
            raise _TypeError('%s -= %r' % (self, other))
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
        return _item_ps(classname(self, prefixed=True), NN)

    def __sub__(self, other):
        '''Difference of this and an other instance or a scalar.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: The difference, a new instance (L{Fsum}).

           @see: Method L{Fsum.__isub__}.
        '''
        f = self.fcopy()
        f -= other
        return f

    def fadd(self, iterable):
        '''Accumulate more values from an iterable.

           @arg iterable: Sequence, list, tuple, etc. (C{scalar}s).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{iterable}} value.

           @raise ValueError: Invalid or non-finite B{C{iterable}} value.
        '''
        if isscalar(iterable):  # for backward compatibility
            iterable = tuple(iterable)

        ps = self._ps
        for a in iterable:  # _iter()
            if not isfinite(a):
                raise _ValueError(iterable=a, txt=_not_finite_)
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

           @arg xs: Values to add (C{scalar}s), all positional.

           @see: Method L{Fsum.fadd}.
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
            raise _ValueError(factor=factor, txt=_not_finite_)

        ps = self._ps
        if ps:  # multiply and adjust partial sums
            ps[:] = [p * factor for p in ps]
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

           @raise ValueError: Unequal C{len}(B{C{a}}) and C{len}(B{C{b}}).

           @see: Function L{fdot} and method L{Fsum.fadd}.
        '''
        if len(a) != len(b):
            raise LenError(Fdot, a=len(a), b=len(b))

        Fsum.__init__(self)
        self.fadd(map(_mul_, a, b))


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
            raise _ValueError(x=x, txt=_not_finite_)
        if not cs:
            raise _ValueError(cs=cs, txt=_Missing)

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

           @arg x: Polynomial argument (C{scalar}).
           @arg cs: Polynomial coeffients (C{scalar}[]).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{x}}.

           @raise ValueError: No B{C{cs}} coefficients or B{C{x}} is not finite.

           @see: Function L{fpolynomial} and method L{Fsum.fadd}.
        '''
        if not isfinite(x):
            raise _ValueError(x=x, txt=_not_finite_)
        if not cs:
            raise _ValueError(cs=cs, txt=_Missing)

        x, cs, xp = float(x), list(cs), 1

        Fsum.__init__(self, cs.pop(0))
        while cs:
            xp *= x
            self.fadd_(xp * cs.pop(0))


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
    return pow(abs(x), _1_3rd)**2  # XXX pow(abs(x), _2_3rd)


def favg(v1, v2, f=0.5):
    '''Return the weighted average of two values.

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

       @raise ValueError: Unequal C{len}(B{C{a}}) and C{len}(B{C{b}}).

       @see: Class L{Fdot}.
    '''
    if len(a) != len(b):
        raise LenError(fdot, a=len(a), b=len(b))

    return fsum(map(_mul_, a, b))


def fdot3(a, b, c, start=0):
    '''Return the precision dot product M{start +
       sum(a[i] * b[i] * c[i] for i=0..len(a))}.

       @arg a: List, sequence, tuple, etc. (C{scalar}[]).
       @arg b: List, sequence, tuple, etc. (C{scalar}[]).
       @arg c: List, sequence, tuple, etc. (C{scalar}[]).
       @kwarg start: Optional bias (C{scalar}).

       @return: Dot product (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise ValueError: Unequal C{len}C{(}B{C{a}}C{)},
                          C{len}C{(}B{C{b}}C{)} and/or
                          C{len}C{(}B{C{c}}C{)}.
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

       @raise ValueError: Invalid B{C{beta}}, negative B{C{ds}} value,
                          weighted B{C{ds}} below L{EPS} or unequal
                          C{len}C{(}B{C{ds}}C{)} and C{len}C{(}B{C{xs}}C{)}.

       @note: Using B{C{beta}}C{=0} returns the mean of B{C{xs}}.
    '''
    n, xs = len2(xs)
    d, ds = len2(ds)
    if n != d or n < 1:
        raise LenError(fidw, xs=n, ds=d)

    d, x = min(zip(ds, xs))
    if d > EPS and n > 1:
        b = -Int_(beta, name=_beta_, low=0, high=3)
        if b < 0:
            ds = tuple(d**b for d in ds)
            d = fsum(ds)
            if d < EPS:
                raise _ValueError(ds=d)
            x = fdot(xs, *ds) / d
        else:
            x = fmean(xs)
    elif d < 0:
        raise _ValueError(_item_sq('ds', ds.index(d)), d)
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

       @return: Powers of B{C{x}} (C{float}[]).

       @raise TypeError: Non-scalar B{C{x}} or B{C{n}} not C{int}.

       @raise ValueError: Non-finite B{C{x}} or non-positive B{C{n}}.
    '''
    if not isfinite(x):
        raise _ValueError(x=x, txt=_not_finite_)
    if not isint(n):
        raise _IsnotError(int.__name__, n=n)
    elif n < 1:
        raise _ValueError(n=n)

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

           @arg iterable: Values to be multiplied (C{scalar}[]).
           @kwarg start: Initial product, also the value returned
                         for an empty iterable (C{scalar}).

           @return: The product (C{float}).

           @see: U{NumPy.prod<https://docs.SciPy.org/doc/
                 numpy/reference/generated/numpy.prod.html>}.
        '''
        return freduce(_mul_, iterable, start)


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
                r, v = 0, _Missing
            for v in iterable:
                r = f(r, v)
            if v is _Missing:
                raise _TypeError(iterable=(), start=_Missing)
            return r


def fsum_(*xs):
    '''Precision summation of the positional argument vulues.

       @arg xs: Values to be added (C{scalar}[]).

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
        '''Compute the norm M{sqrt(sum(xs[i]**2)) for i=0..len(xs)}.

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
        if xs:
            n, xs = len2(xs)
            if n > 0:
                h = float(max(abs(x) for x in xs))
                if h > 0 and n > 1:
                    X = Fsum(1.0)
                    X.fadd((x / h)**2 for x in xs)
                    h *= sqrt(X.fsum_(-1.0))
                return h
        raise _ValueError(xs=xs, txt=_too_few_)


def hypot1(x):
    '''Compute the norm M{sqrt(1 + x**2)}.

       @arg x: Argument (C{scalar}).

       @return: Norm (C{float}).
    '''
    return hypot(1.0, x)


def hypot2(x, y):
    '''Compute the norm, I{squared} M{x**2 + y**2}.

       @arg x: Argument (C{scalar}).
       @arg y: Argument (C{scalar}).

       @return: B{C{x}}C{**2 + }B{C{y}}C{**2} (C{float}).
    '''
    x, y = x**2, y**2
    if x < y:
        x, y = y, x
    if y:  # and x
        x *= 1 + y / x
    return x


def sqrt3(x):
    '''Compute the square root, cubed M{sqrt(x)**3} or M{sqrt(x**3)}.

       @arg x: Value (C{scalar}).

       @return: Cubed square root (C{float}).

       @raise ValueError: Negative B{C{x}}.

       @see: Functions L{cbrt} and L{cbrt2}.
    '''
    if x < 0:
        raise _ValueError(x=x)
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
