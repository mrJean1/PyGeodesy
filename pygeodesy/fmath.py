
# -*- coding: utf-8 -*-

u'''Precision mathematical functions, utilities and constants.

@newfield example: Example, Examples
'''

from math import acos, copysign, hypot, sqrt  # pow
from operator import mul
import sys

# all public contants, classes and functions
__all__ = ('EPS', 'EPS1',  # constants
           'Fsum',  # classes
           'acos1',
           'cbrt', 'cbrt2',  # functions
           'favg', 'fdot', 'fdot3', 'fmean',
           'fhorner', 'fpolynomial', 'fpowers',
           'fStr', 'fStrzs', 'fsum', 'fsum_',
           'hypot', 'hypot1', 'hypot3',
           'isfinite', 'isint', 'isscalar',
           'len2',
           'map1', 'map2',
           'scalar', 'sqrt3')
__version__ = '18.10.26'

try:  # Luciano Ramalho, "Fluent Python", page 395, O'Reilly, 2016
    from numbers import Integral as _Ints  #: (INTERNAL) Int objects
except ImportError:
    try:
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

try:
    from collections import Sequence as _Seqs  #: (INTERNAL) incl MutableSequence
except ImportError:
    _Seqs = list, tuple, range  # XXX also set?

try:
    EPS = sys.float_info.epsilon  #: System's epsilon (C{float})
except AttributeError:
    EPS = 2.220446049250313e-16  #: Approximate epsilon (C{float})
EPS1  = 1.0 - EPS  #: M{1 - EPS} (C{float}), about 0.9999999999999998
_1EPS = 1.0 + EPS  #: M{1 + EPS} (C{float})

_1_3rd = 1.0 / 3.0  #: (INTERNAL) One third (C{float})
_2_3rd = 2.0 / 3.0  #: (INTERNAL) Two third (C{float})
_3_2nd = 3.0 / 2.0  #: (INTERNAL) Three halfs (C{float})


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
        return s, a - (s - b)
    else:
        return s, b - (s - a)


class Fsum(object):
    '''Precision floating point summation similar to standard
       Python function C{math.fsum}.

       Unlike C{math.fsum}, this class accumulates the values
       incrementally and provides intermediate, precision, running
       sums.  Accumulation may continue after intermediate
       summations.

       @note: Exception and I{non-finite} handling differ from C{math.fsum}.

       @see: U{Hettinger<http://code.ActiveState.com/recipes/393090>},
             U{Kahan<http://WikiPedia.org/wiki/Kahan_summation_algorithm>},
             U{Klein<http://link.Springer.com/article/10.1007/s00607-005-0139-x>},
             Python 2.6+ file I{Modules/mathmodule.c} and the issue log
             U{Full precision summation<http://bugs.Python.org/issue2819>}.
    '''
    def __init__(self, *starts):
        '''Initialize a new accumulator with one or more start values.

           @param starts: No, one or more start values (scalars).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar I{starts} value.

           @raise ValueError: Invalid or infinite I{starts} value.
        '''
        self._n = 0
        self._ps = []
        if starts:
            self._add(starts)

    def __len__(self):
        '''Return the number of accumulated values.
        '''
        return self._n

    def _add(self, iterable):
        '''(INTERNAL) Accumulate more values.

           @param iterable: Sequence, list, tuple, etc. (scalars).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar  I{iterable} value.

           @raise ValueError: Invalid or infinite I{iterable} value.
        '''
        ps = self._ps
        for a in iterable:
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

    def fadd(self, *args):
        '''Accumulate more values from positional arguments.

           @param args: Values to add (scalars), all positional.

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar I{args} value.

           @raise ValueError: Invalid or infinite I{arg}.
        '''
        self._add(args)

    def fmul(self, factor):
        '''Multiple the current, partial sum by a factor.

           @param factor: The multiplier (C{scalar}).

           @raise TypeError: Non-scalar I{factor}.

           @raise ValueError: Invalid or infinite I{factor}.
        '''
        if not isfinite(factor):
            raise ValueError('not %s: %r' %('finite', factor))

        ps = self._ps
        if ps:  # multiply and adjust partial sums
            ps[:] = [p * factor for p in ps]
            self.fadd(ps.pop())
            self._n -= 1
        # assert self._ps is ps

    def fsum_(self, *args):
        '''Accumulated more values from positional arguments and sum all.

           @param args: Values to add (scalars), all positional.

           @return: Accurate, running sum (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar I{args} value.

           @raise ValueError: Invalid or infinite I{args} value.

           @note: Accumulation can continue after summation.
        '''
        return self.fsum(args)

    def fsum(self, iterable=()):
        '''Accumulated more values from the iterable and sum all.

           @param iterable: Sequence, list, tuple, etc. (scalars), optional.

           @return: Accurate, running sum (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar I{iterable} value.

           @raise ValueError: Invalid or infinite I{iterable} value.

           @note: Accumulation can continue after summation.
        '''
        if iterable:
            self._add(iterable)

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
                if p:  # sum(ps) became inexaxt
                    ps.append(p)
                    if i > 0:  # half-even round if signs match
                        s = _2even(s, ps[i-1], p)
                    break
            # assert self._ps is ps
        return s


def acos1(x):
    '''Return M{math.acos(max(-1, min(1, x)))}.
    '''
    return acos(max(-1.0, min(1.0, x)))


def cbrt(x):
    '''Compute the cubic root M{x**(1/3)}.

       @param x: Value (C{scalar}).

       @return: Cubic root (C{float}).
    '''
    # simpler and more accurate than Ken Turkowski's CubeRoot, see
    # <http://People.FreeBSD.org/~lstewart/references/apple_tr_kt32_cuberoot.pdf>
    return copysign(pow(abs(x), _1_3rd), x)


def cbrt2(x):
    '''Compute the cubic root squared M{x**(2/3)}.

       @param x: Value (C{scalar}).

       @return: Cubic root squared (C{float}).

       @see: Function L{sqrt3}.
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

       @param a: List, sequence, tuple, etc. (scalars).
       @param b: All positional arguments (scalars).

       @return: Dot product (C{float}).

       @raise ValueError: Unequal len(a) and len(b).
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

       @raise ValueError: Unequal C{len}(I{a}), C{len}(I{b})
                          and/or C{len}(I{c}).
    '''
    def _mul3(a, b, c):  # map function
        return a * b * c

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

       @raise TypeError: Non-scalar I{x}.

       @raise ValueError: No I{cs} coefficients or I{x} is not finite.

       @see: Function L{fpolynomial}.
    '''
    if not isfinite(x):
        raise ValueError('not %s: %r' %('finite', x))
    if not cs:
        raise ValueError('no %s: %r' % ('coefficents', cs))

    x, cs = float(x), list(cs)

    h = Fsum(cs.pop())
    while cs:
        h.fmul(x)
        h.fadd(cs.pop())
    return h.fsum()


def fmean(xs):
    '''Compute the accurate mean M{sum(xs[i] for
       i=0..len(xs)) / len(xs)}.

       @param xs: Values (scalars).

       @return: Mean value (C{float}).

       @raise ValueError: No I{xs} values.
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

       @raise TypeError: Non-scalar I{x}.

       @raise ValueError: No I{cs} coefficients or I{x} is not finite.

       @see: Function L{fhorner}.
    '''
    if not isfinite(x):
        raise ValueError('not %s: %r' %('finite', x))
    if not cs:
        raise ValueError('no %s: %r' % ('coefficents', cs))

    def _terms(x, c0, *cs):
        xp = 1
        yield float(c0)
        for c in cs:
            xp *= x
            yield xp * c

    return fsum(_terms(float(x), *cs))


def fpowers(x, n, alts=0):
    '''Return a series of powers M{[x**i for i=1..n]}.

       @param x: Value (C{scalar}).
       @param n: Highest exponent (C{int}).
       @keyword alts: Only alternating powers, starting
                      with this exponent (C{int}).

       @return: Powers of I{x} (C{float}[]).

       @raise TypeError: Non-scalar I{x} or I{n} not C{int}.

       @raise ValueError: Non-positive I{n} or I{x} is not finite.
    '''
    if not isfinite(x):
        raise ValueError('not %s: %r' %('finite', x))
    if not isinstance(n, _Ints):
        raise TypeError('%s invalid: %r' % ('n', n))
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


def fStr(floats, prec=6, sep=', ', fmt='%.*f', ints=False):
    '''Convert floats to string, optionally with trailing zero
       decimals stripped.

       @param floats: List, sequence, tuple, etc. (scalars).
       @keyword prec: Optional precision, number of decimal digits (0..9).
                      Trailing zero decimals are stripped for I{prec} values
                      of 1 and above, but kept for negative I{prec} values.
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
        return t

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


def fsum_(*args):
    '''Precision floating point sum of the positional argument vulues.

       @param args: Values to be added (C{scalar}[]).

       @return: Accurate L{fsum} (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar I{arg} value.

       @raise ValueError: Invalid or infinite I{arg} value.
    '''
    return fsum(args)


try:
    from math import fsum  # precision IEEE-754 sum, Python 2.6+

    # make sure fsum works as expected (XXX check
    # float.__getformat__('float')[:4] == 'IEEE'?)
    if fsum_(1, 1e101, 1, -1e101) != 2:
        del fsum  # no, remove fsum ...
        raise ImportError  # ... use fsum below

except ImportError:

    def fsum(iterable):
        '''Precision floating point sum similar to standard
           Python function C{math.fsum}.

           Exception and I{non-finite} handling differs from C{math.fsum}.

           @param iterable: Values to be added (C{scalar}[]).

           @return: Accurate C{sum} (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar I{iterable} value.

           @raise ValueError: Invalid or infinite I{iterable} value.

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


def hypot3(x, y, z):
    '''Compute the norm M{sqrt(x**2 + y**2 + z**2)}.

       @param x: X argument (C{scalar}).
       @param y: Y argument (C{scalar}).
       @param z: Z argument (C{scalar}).

       @return: Norm (C{float}).
    '''
    x, y, z = map1(abs, x, y, z)
    if x < z:
        x, z = z, x
    if y < z:
        y, z = z, y
    if z:
        if x < y:
            x, y = y, x
        h = float(x)
        if h > EPS:
            # XXX PyChecker chokes on /= and *=
            t = fsum_(1.0, (y / h)**2, (z / h)**2)
            if t > _1EPS:
                h *= sqrt(t)
    else:
        h = hypot(x, y)
    return h


try:
    from math import isfinite  # new in Python 3+
except ImportError:
    from math import isinf, isnan

    def isfinite(obj):
        '''Check for C{Inf} and C{NaN} values.

           @param obj: Value (C{scalar}).

           @return: C{False} if I{obj} is C{Inf} or C{NaN},
                    C{True} otherwise.

           @raise TypeError: Non-scalar I{obj}.
        '''
        if isscalar(obj):
            return not (isinf(obj) or isnan(obj))
        raise TypeError('not %s: %r' % ('scalar', obj))


def isint(obj, both=False):
    '''Check for integer type or integer value.

       @param obj: The object (any C{type}).
       @keyword both: Optionally, check both type and value (C{bool}).

       @return: C{True} if I{obj} is C{int}, C{False} otherwise.
    '''
    if both and isinstance(obj, float):  # NOT _Scalars!
        try:
            return obj.is_integer()
        except AttributeError:
            return False  # XXX float(int(obj)) == obj?
    return isinstance(obj, _Ints)


def isscalar(obj):
    '''Check for scalar types.

       @param obj: The object (any C{type}).

       @return: C{True} if I{obj} is C{scalar}, C{False} otherwise.
    '''
    return isinstance(obj, _Scalars)


def len2(seq):
    '''Make built-in function L{len} work for generators, iterators,
       etc. since those can only be started exactly once.

       @param seq: Generator, iterator, list, range, tuple, etc.

       @return: 2-Tuple (number, ...) of items (C{int}, C{list} or
                C{range} or C{tuple}).
    '''
    if not isinstance(seq, _Seqs):
        seq = list(seq)
    return len(seq), seq


def map1(func, *args):
    '''Apply each argument to a single-argument function and
       return a tuple of results.

       @param func: Function to apply (callable).
       @param args: Arguments to apply (any positional).

       @return: Function results (C{tuple}).
    '''
    return tuple(map(func, args))


def map2(func, *args):
    '''Apply arguments to a function and return a tuple of results.

       Unlike Python 2's built-in L{map}, Python 3+ L{map} returns a
       L{map} object, an iterator-like object which generates the
       results only once.  Converting the L{map} object to a tuple
       maintains Python 2 behavior.

       @param func: Function to apply (callable).
       @param args: Arguments to apply (list, tuple, ...).

       @return: N-Tuple of function results (C{tuple}).
    '''
    return tuple(map(func, *args))


def scalar(value, low=EPS, high=1.0, name='scalar'):
    '''Validate a scalar.

       @param value: The value (C{scalar}).
       @keyword low: Optional lower bound (C{scalar}).
       @keyword high: Optional upper bound (C{scalar}).
       @keyword name: Optional name of value (C{str}).

       @return: New value (C{type} of I{low}).

       @raise TypeError: Non-scalar I{value}.

       @raise ValueError: Out-of-bounds I{value}.
    '''
    if not isscalar(value):
        raise TypeError('%s invalid: %r' % (name, value))
    try:
        if low is None:
            v = float(value)
        else:
            v = type(low)(value)
            if v < low or v > high:
                raise ValueError
    except (TypeError, ValueError):
        raise ValueError('%s invalid: %r' % (name, value))
    return v


def sqrt3(x):
    '''Compute the square root cubed M{sqrt(x)**3} or M{sqrt(x**3)}.

       @param x: Value (C{scalar}).

       @return: Cubid square root (C{float}).

       @raise ValueError: Negative I{x}.

       @see: Function L{cbrt2}.
    '''
    if x < 0:
        raise ValueError('%s(%r)' % ('sqrt3', x))
    return pow(x, _3_2nd)

# **) MIT License
#
# Copyright (C) 2016-2018 -- mrJean1 at Gmail dot com
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
