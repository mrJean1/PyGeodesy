
# -*- coding: utf-8 -*-

u'''Precision mathematical functions, utilities and constants.

@newfield example: Example, Examples
'''

from math import hypot, sqrt  # pow
from operator import mul
import sys

# all public contants, classes and functions
__all__ = ('EPS', 'EPS1', 'EPS2',  # constants
           'Fsum',  # classes
           'cbrt', 'cbrt2',  # functions
           'favg', 'fdot', 'fdot3', 'fmean', 'fpolynomial', 'fpowers',
           'fStr', 'fStrzs', 'fsum', 'fsum_',
           'hypot', 'hypot1', 'hypot3',
           'isfinite', 'isint', 'isscalar',
           'len2',
           'map1', 'map2',
           'scalar')
__version__ = '18.02.09'

try:  # Luciano Ramalho, "Fluent Python", page 395, O'Reilly, 2016
    from numbers import Integral as _Ints  #: (INTERNAL) Int objects
except ImportError:
    try:
        _Ints = int, long  #: (INTERNAL) Int objects (tuple)
    except NameError:  # Python 3+
        _Ints = int  #: (INTERNAL) Int objects (tuple)

try:  # similarly ...
    from numbers import Real as _Scalars  #: (INTERNAL) Scalar objects
except ImportError:
    try:
        _Scalars = int, long, float  #: (INTERNAL) Scalar objects (tuple)
    except NameError:
        _Scalars = int, float  #: (INTERNAL) Scalar objects (tuple)

try:
    from collections import Sequence as _Seqs  #: (INTERNAL) incl MutableSequence
except ImportError:
    _Seqs = list, tuple, range  # XXX also set?

try:
    EPS = sys.float_info.epsilon  #: System's epsilon (float)
except AttributeError:
    EPS = 2.220446049250313e-16  #: Approximate epsilon (float)
EPS1 = 1.0 - EPS  #: M{1.0 - EPS} (float), about 0.9999999999999998
EPS2 = sqrt(EPS)  #: M{sqrt(EPS)} (float)

_1_3rd = 1.0 / 3.0  #: (INTERNAL) One third (float)
_2_3rd = 2.0 / 3.0  #: (INTERNAL) Two third (float)


def _2sum(a, b):
    '''(INTERNAL) Precision I{2sum} of M{a + b}.
    '''
    if abs(b) > abs(a):
        a, b = b, a
    s = a + b
    if not isfinite(s):
        raise OverflowError('partial %s: %r' % ('2sum', s))
    return s, b - (s - a)


class Fsum(object):
    '''Precision floating point summation similar to standard
       Python function I{math.fsum}.

       Unlike I{math.fsum}, this class accumulates the values to
       be added incrementally and provides intermediate, precision
       summations.  Accumulation can continue after intermediate
       summations.

       Exception and I{non-finite} handling differs from I{math.fsum}.

       @see: U{Hettinger<http://code.activestate.com/recipes/393090/>},
             U{Klein<http://link.springer.com/article/10.1007/s00607-005-0139-x>},
             U{Kahan<http://wikipedia.org/wiki/Kahan_summation_algorithm>},
             Python 2.6+ file I{Modules/mathmodule.c} and the issue log
             U{Full precision summation<http://bugs.python.org/issue2819>}.
    '''
    def __init__(self, start=0):
        '''Initialize a new accumulator.
        '''
        self._n = 0
        self._ps = [self._arg(start)] if start else []

    def __len__(self):
        '''Return the number of accumulated values.
        '''
        return self._n

    def _arg(self, arg):
        try:
            a = float(arg)
            if isfinite(a):
                self._n += 1
                return a
        except (TypeError, ValueError):
            pass
        raise ValueError('%s invalid: %r' % ('Fsum', arg))

    def fadd(self, arg):
        '''Accumulate another value.

           @param arg: Value to add (scalar).

           @raise OverflowError: Partial I{2sum} overflow.

           @raise ValueError: Invalid or infinite I{arg}.
        '''
        a = self._arg(arg)
        i = 0
        for p in self._ps:
            a, p = _2sum(a, p)
            if p:
                self._ps[i] = p
                i += 1
        self._ps[i:] = [a]

    def fadd2(self, iterable):
        '''Accumulate multiple values.

           @param iterable: Sequence, list, tuple, etc. (scalars).

           @raise OverflowError: Partial I{2sum} overflow.

           @raise ValueError: Invalid or infinite I{iterable} value.
        '''
        for a in iterable:
            self.fadd(a)

    def fsum(self):
        '''Summation of the values accumulated so far.

           @return: Precision sum (float).

           @raise OverflowError: Partial I{2sum} overflow.

           @note: Accumulation can continue after summation.
        '''
        s = 0.0
        i = len(self._ps)
        if i > 0:
            i -= 1
            s = self._ps[i]
            while i > 0:
                i -= 1
                s, p = _2sum(s, self._ps[i])
                if p:
                    if i > 0:  # half-even round if signs match
                        r = self._ps[i - 1]
                        if (r > 0 and p > 0) or \
                           (r < 0 and p < 0):
                            r, p = _2sum(s, p * 2)
                            if not p:
                                s = r
                    break
        return s


def cbrt(x):
    '''Compute the cubic root M{x**(1/3)}.

       @param x: Value (scalar).

       @return: Cubic root (float).
    '''
    # simpler and more accurate than Ken Turkowski's CubeRoot, see
    # <http://people.freebsd.org/~lstewart/references/apple_tr_kt32_cuberoot.pdf>
    if x < 0:
        return -pow(-x, _1_3rd)
    else:
        return  pow( x, _1_3rd)


def cbrt2(x):
    '''Compute the cubic root squared M{x**(2/3)}.

       @param x: Value (scalar).

       @return: Cubic root squared (float).
    '''
    return pow(abs(x), _2_3rd)


def favg(v1, v2, f=0.5):
    '''Return the weighted average of two values.

       @param v1: One value (scalar).
       @param v2: Other value (scalar).
       @keyword f: Optional fraction (scalar).

       @return: M{v1 + f * (v2 - v1)} (float).
    '''
#      @raise ValueError: Fraction out of range.
#   '''
#   if not 0 <= f <= 1:  # XXX restrict fraction?
#       raise ValueError('%s invalid: %r' % ('fraction', f))
    return v1 + f * (v2 - v1)  # v1 * (1 - f) + v2 * f


def fdot(a, *b):
    '''Return the precision dot product M{sum(a[i] * b[i]),
       i=0..len(a)}.

       @param a: List, sequence, tuple, etc. (scalars).
       @param b: List, sequence, tuple, etc. (scalars).

       @return: Dot product (float).

       @raise ValueError: Unequal len(a) and len(b).
    '''
    if not len(a) == len(b):
        raise ValueError('%s: %s vs %s' % ('len', len(a), len(b)))

    return fsum(map(mul, a, b))


def fdot3(a, b, c, start=0):
    '''Return the precision dot product M{start +
       sum(a[i] * b[i] * c[i]), i=0..len(a)}.

       @param a: List, sequence, tuple, etc. (scalars).
       @param b: List, sequence, tuple, etc. (scalars).
       @param c: List, sequence, tuple, etc. (scalars).
       @keyword start: Optional bias (scalar).

       @return: Dot product (float).

       @raise ValueError: Unequal len(a), len(b) and/or len(c).
    '''
    def _mul3(a, b, c):  # map function
        return a * b * c

    if not len(a) == len(b) == len(c):
        raise ValueError('%s: %s vs %s vs %s' % ('len', len(a), len(b), len(c)))

    if start:
        f = Fsum(start)
        f.fadd2(map(_mul3, a, b, c))
        return f.fsum()
    else:
        return fsum(map(_mul3, a, b, c))


def fmean(xs):
    '''Compute the accurate mean M{sum(xs[i]), i=0..len(xs) /
       len(xs)}.

       @param xs: Values (scalars).

       @return: Mean value (float).

       @raise ValueError: No I{xs} values.
    '''
    n, xs = len2(xs)
    if n > 0:
        return fsum(xs) / n
    raise ValueError('%s missing: %r' % ('xs', xs))


def fpolynomial(x, *cs):
    '''Accuately evaluate the polynomial M{sum(cs[i] * x**i),
       i=0..len(cs))}.

       @param x: Polynomial argument (scalar).
       @param cs: Polynomial coeffients (scalars).

       @return: Polynomial value (float).

       @raise TypeError: Argument I{x} not scalar.

       @raise ValueError: No I{cs} coefficients.
    '''
    if not isscalar(x):
        raise TypeError('%s invalid: %r' % ('argument', x))
    if not cs:
        raise ValueError('%s missing: %r' % ('coefficents', cs))

    def _terms(x, c0, *cs):
        xp = x
        yield float(c0)
        for c in cs:
            yield xp * c
            xp *= x

    return fsum(_terms(float(x), *cs))


def fpowers(x, n, alts=0):
    '''Return a series of powers M{x**i for i=1..n}.

       @param x: Value (scalar).
       @param n: Highest exponent (int).
       @keyword alts: Only alternating powers, starting
                      with this exponent (int).

       @return: Powers of I{x} (list of floats).
    '''
    xs = [x]
    for _ in range(1, n):
        xs.append(xs[-1] * x)

    if alts > 0:  # x**2, x**4, ...
        xs = xs[alts-1: :2]

    # XXX PyChecker falsely claims result is None
    return xs


def fStr(floats, prec=6, sep=', ', fmt='%.*f', ints=False):
    '''Convert floats to string, optionally with trailing zero
       decimals stripped.

       @param floats: List, sequence, tuple, etc. (scalars).
       @keyword prec: Optional precision, number of decimal digits (0..9).
                      Trailing zero decimals are stripped for prec values
                      of 1 and above, but kept for negative prec values.
       @keyword sep: Optional, separator to join (string).
       @keyword fmt: Optional, float format (string).
       @keyword ints: Optionally, remove decimal dot (bool).

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
        if z > 1:
            fstr = fstr[:z] + fstr[z:].rstrip('0')
    return fstr


def fsum_(*args):
    '''Precision floating point summation of the positional arguments.

       @param args: Values to be added (scalars).

       @return: Precision L{fsum} (float).
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
        '''Precision floating point summation similar to standard
           Python function I{math.fsum}.

           Exception and I{non-finite} handling differs from I{math.fsum}.

           @param iterable: Sequence, list, tuple, etc. (scalars).

           @return: Precision sum (float).

           @raise OverflowError: Partial I{2sum} overflow.

           @raise ValueError: Invalid or infinite I{iterable} value.

           @see: Class L{Fsum}.
        '''
        f = Fsum()
        f.fadd2(iterable)
        return f.fsum()


def hypot1(x):
    '''Compute the norm M{sqrt(1 + x**2)}.

       @param x: Argument (scalar).

       @return: Norm (float).
    '''
    return hypot(1.0, x)


def hypot3(x, y, z):
    '''Compute the norm M{sqrt(x**2 + y**2 + z**2)}.

       @param x: X argument (scalar).
       @param y: Y argument (scalar).
       @param z: Z argument (scalar).

       @return: Norm (float).
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
            # XXX PyChecker chokes on some /= and *=
            t = (y / h)**2 + (z / h)**2
            if t > EPS:
                h *= sqrt(1.0 + t)
    else:
        h = hypot(x, y)
    return h


try:
    from math import isfinite  # new in Python 3+
except ImportError:
    from math import isinf, isnan

    def isfinite(obj):
        '''Check for I{Inf} and I{NaN} values.

           @param obj: Value (scalar).

           @return: False if I{Inf} or I{NaN}, True otherwise (bool).

           @raise TypeError: If I{obj} value is not scalar.
        '''
        if isscalar(obj):
            return not (isinf(obj) or isnan(obj))
        raise TypeError('%s invalid: %r' % ('isfinite', obj))


def isint(obj, both=False):
    '''Check for integer type or integer value.

       @param obj: The object (any).
       @keyword both: Optionally, check both type and value (bool).

       @return: True if I{obj} is integer (bool).
    '''
    if both and isinstance(obj, float):  # NOT _Scalars!
        try:
            return obj.is_integer()
        except AttributeError:
            return False  # XXX float(int(obj)) == obj?
    return isinstance(obj, _Ints)


def isscalar(obj):
    '''Check for scalar types.

       @param obj: The object (any).

       @return: True if I{obj} is scalar (bool).
    '''
    return isinstance(obj, _Scalars)


def len2(seq):
    '''Make built-in function L{len} work for generators,
       iterators, etc. since those can only be started once.

       @param seq: Generator, iterator, list, range, tuple, etc.

       @return: 2-Tuple (number, list) of items (int, list).
    '''
    if not isinstance(seq, _Seqs):
        seq = list(seq)
    return len(seq), seq


def map1(func, *args):
    '''Apply each argument to a single-argument function and
       return a tuple of results.

       @param func: Function to apply (callable).
       @param args: Arguments to apply (any).

       @return: Function results (tuple).
    '''
    return tuple(map(func, args))


def map2(func, *args):
    '''Apply arguments to a function and return a tuple of results.

       Unlike Python 2 built-in L{map}, Python 3+ L{map} returns a
       L{map} object, an iterator-like object which generates the
       results only once.  Converting the L{map} object to a tuple
       maintains Python 2 behavior.

       @param func: Function to apply (callable).
       @param args: Arguments to apply (list, tuple, ...).

       @return: N-Tuple of function results (tuple).
    '''
    return tuple(map(func, *args))


def scalar(value, low=EPS, high=1.0, name='scalar'):
    '''Validate a scalar.

       @param value: The value (scalar).
       @keyword low: Optional lower bound (scalar).
       @keyword high: Optional upper bound (scalar).
       @keyword name: Optional name of value (string).

       @return: New value (type(low)).

       @raise TypeError: The I{value} is not scalar.

       @raise ValueError: The I{value} is out of bounds.
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
