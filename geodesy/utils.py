
# -*- coding: utf-8 -*-

'''Mathematical and utility functions and constants.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**,
see U{http://www.movable-type.co.uk/scripts/latlong.html}
and U{http://www.movable-type.co.uk/scripts/latlong-vectors.html}.

@newfield example: Example, Examples
'''

from math import degrees, pi as PI, radians, sin, sqrt, tan
try:
    from math import fsum  # precision sum, Python 2.6+
except ImportError:
    fsum = sum  # use standard, built-in sum (or Kahan's summation
    # <https://en.wikipedia.org/wiki/Kahan_summation_algorithm> or
    # Hettinger's <https://code.activestate.com/recipes/393090/>)
from operator import mul
import sys

# all public contants, classes and functions
__all__ = ('EPS', 'EPS1', 'EPS2', 'PI', 'PI2', 'PI_2',  # constants
           'cbrt',
           'degrees', 'degrees90', 'degrees180', 'degrees360',
           'false2f', 'favg', 'fdot', 'fdot3', 'fStr', 'fsum',
           'halfs', 'hsin', 'hypot1', 'hypot3',
           'isint', 'isscalar', 'len2', 'map1', 'map2',
           'radians', 'radiansPI', 'radiansPI2', 'radiansPI_2',
           'tanPI_2_2',
           'wrap90', 'wrap180', 'wrapPI', 'wrapPI2', 'wrapPI_2')
__version__ = '17.03.12'

try:
    _Ints = int, long  #: (INTERNAL) Int objects (tuple)
    _Scalars = int, long, float  #: (INTERNAL) Scalar objects (tuple)
except NameError:  # Python 3+
    _Ints = int  #: (INTERNAL) Int objects (tuple)
    _Scalars = int, float  #: (INTERNAL) Scalar objects (tuple)

try:
    EPS = sys.float_info.epsilon  #: System's epsilon (float)
except AttributeError:
    EPS = 2.2204460492503131e-16  #: Approximate epsilon (float)
EPS1 = 1.0 - EPS  #: 1.0 - EPS (float)
EPS2 = sqrt(EPS)  #: M{sqrt(EPS)} (float)

PI2  = PI * 2  #: Two PI, M{PI * 2} (float)  # PYCHOK expected
PI_2 = PI / 2  #: Half PI, M{PI / 2} (float)

_3rd = 1.0 / 3.0  #: (INTERNAL) One third (float)


def cbrt(x):
    '''Computes the cubic root M{x**(1/3)}.

       @param x: Scalar (float or int).

       @return: Cubic root (float).
    '''
    # simpler and more accurate than Ken Turkowski's CubeRoot, see
    # <http://people.freebsd.org/~lstewart/references/apple_tr_kt32_cuberoot.pdf>
    if x < 0:
        return -(float(-x) ** _3rd)
    else:
        return float(x) ** _3rd


def degrees90(rad):
    '''Converts and wraps radians to degrees M{-270..+90}.

       @param rad: Angle (radians).

       @return: Degrees, wrapped (degrees90).
    '''
    return _drap(degrees(rad), 90)


def degrees180(rad):
    '''Converts and wraps radians to degrees M{-180..+180}.

       @param rad: Angle (radians).

       @return: Degrees, wrapped (degrees180).
    '''
    return _drap(degrees(rad), 180)


def degrees360(rad):
    '''Converts and wraps radians to degrees M{0..+360}.

       @param rad: Angle (radians).

       @return: Degrees, wrapped (degrees360).
    '''
    return _drap(degrees(rad), 360)


def _drap(deg, wrap):
    '''(INTERNAL) Degree wrapper M{(wrap-360)..+wrap}.

       @param deg: Angle (degrees).
       @param wrap: Limit (degrees).

       @returns: Degrees, wrapped (degrees).
    '''
    d = deg % 360  # -1.5 % 360 == 358.5
    if d > wrap:
        d -= 360
    return d


def false2f(value, name='value', false=True):
    '''Converts false east-/northing to non-negative float.

       @param value: Value to convert (scalar).
       @keyword name: Name of the value (string).
       @keyword false: Value must include false origin (bool).

       @return: The value (float).

       @raise ValueError: Invalid or negative value.
    '''
    try:
        f = float(value)
        if f < 0 and false:
            raise ValueError
    except (TypeError, ValueError):
        raise ValueError('%s invalid: %r' % (name, value))
    return f


def favg(v1, v2, f=0.5):
    '''Weighted average of two values.

       @param v1: One value (scalar).
       @param v2: Other value (saclar).
       @keyword f: Fraction (scalar).

       @return: M{v1 + f * (v2 - v1)} (float).

       @raise ValueError: Fraction out of range.
    '''
    if not 0 <= f <= 1:
        raise ValueError('%s invalid: %r' % ('fraction', f))
    return v1 + f * (v2 - v1)


def fdot(a, *b):
    '''Precision dot product M{sum(a[i] * b[i]
       for i in range(len(a)))}.

       @param a: List, sequence, tuple, etc. (scalars).
       @param b: List, sequence, tuple, etc. (scalars).

       @return: Dot product (float).

       @raise ValueError: Unequal len(a) and len(b).
    '''
    if not len(a) == len(b):
        raise ValueError('unequal len: %s vs %s' % (len(a), len(b)))

    return fsum(map(mul, a, b))


def fdot3(a, b, c, start=0):
    '''Precision dot product M{sum(a[i] * b[i] * c[i]
       for i in range(len(a))) + start}.

       @param a: List, sequence, tuple, etc. (scalars).
       @param b: List, sequence, tuple, etc. (scalars).
       @param c: List, sequence, tuple, etc. (scalars).

       @return: Dot product (float).

       @raise ValueError: Unequal len(a), len(b) and/or len(c).
    '''
    def mul3(a, b, c):  # map function
        return a * b * c

    if not len(a) == len(b) == len(c):
        raise ValueError('unequal len: %s vs %s vs %s' % (len(a), len(b), len(c)))

    m3 = map(mul3, a, b, c)
    if start:
        m3 = (start,) + tuple(m3)
    return fsum(m3)


def fStr(floats, prec=6, sep=', ', fmt='%.*f', ints=False):
    '''Converts floats to string, optionally with trailing
       zero decimals stripped.

       @param floats: List, sequence, tuple, etc. (scalars).
       @keyword prec: Optional, number of decimals, unstripped.  Trailing
                      zero decimals are not stripped if prec is 2 or negative.
       @keyword sep: Optional, separator to join (string).
       @keyword fmt: Optional, float format (string).
       @keyword ints: Optionally, remove decimal dot (bool).

       @return: The floats as 'f, f, ... f' (string).
    '''
    def _fstr(p, f):
        t = fmt % (abs(p), float(f))
        if ints and isint(f):
            t = t.split('.')[0]
        elif p > 1 and t.endswith('0'):
            z = len(t) - p + 1
            t = t[:z] + t[z:].rstrip('0')
        return t

    if isscalar(floats):
        return _fstr(prec, floats)
    else:
        return sep.join(_fstr(prec, f) for f in floats)


def halfs(str2):
    '''Splits a string in 2 halfs.

       @param str2: String to split (string).

       @return: 2-Tuple (1st, 2nd) halfs (strings).

       @raise ValueError: Zero or odd len(str2).
    '''
    h, r = divmod(len(str2), 2)
    if r or not h:
        raise ValueError('%s invalid: %r' % ('str2', str2))
    return str2[:h], str2[h:]


def hsin(rad):
    '''Compute the Haversine value of an angle.

       @param rad: Angle (radians).

       @return: M{sin(rad / 2)**2} (float).

       @note: U{http://en.wikipedia.org/wiki/Haversine_formula}.
    '''
    h = sin(rad * 0.5)
    return h * h


def hypot1(x):
    '''Computes the norm M{sqrt(1 + x**2)}.

       @param x: Argument (scalar).

       @return: Norm (float).
    '''
    h = abs(x)
    if h > 1:
        x = 1.0 / h
        if x > EPS2:
            h *= sqrt(1 + x * x)
    elif h > EPS:
        h = sqrt(1 + x * x)
    else:
        h = 1.0
    return h


def hypot3(x, y, z):
    '''Computes the norm M{sqrt(x**2 + y**2 + z**2)}.

       @param x: X argument (scalar).
       @param y: Y argument (scalar).
       @param z: Z argument (scalar).

       @return: Norm (float).
    '''
    x, y, z = abs(x), abs(y), abs(z)
    if x < y:
        x, y = y, x
    if x < z:
        x, z = z, x
    x = float(x)
    if x > EPS:
        y /= x
        z /= x
        h = x * sqrt(1 + y * y + z * z)
    elif x:
        h = x  # EPS
    else:
        h = 0
    return h


def isint(obj, both=False):
    '''Checks for integer types and value.

       @param obj: The object (any type).
       @keyword both: Check both type and value (bool).

       @return: True if obj is integer (bool).
    '''
    if both and isinstance(obj, float):
        try:
            return obj.is_integer()
        except AttributeError:
            return False
    return isinstance(obj, _Ints)


def isscalar(obj):
    '''Checks for scalar types.

       @param obj: The object (any type).

       @return: True if obj is scalar (bool).
    '''
    return isinstance(obj, _Scalars)


def len2(xtor):
    '''Makes built-in L{len}() function work for generators,
       iterators, etc. since those can only be started once.

       @param xtor: Generator, iterator, list, sequence, tuple, etc.

       @return: 2-Tuple (number, list) of items (int, list).
    '''
    if not isinstance(xtor, (list, tuple)):
        xtor = list(xtor)
    return len(xtor), xtor


def map1(func, *args):
    '''Applies each argument to a single-argument function.

       @param func: Function to apply (callable).
       @param args: Arguments to apply (any).

       @return: Function results (tuple).
    '''
    return tuple(map(func, args))


def map2(func, *args):
    '''Applies arguments to a function, like built-in L{map}.

       Python 3+ L{map} returns a L{map} object, an iterator-like
       object which generates the results only once.  Converting
       the L{map} object into a tuple maintains Python 2 behavior.

       @param func: Function to apply (callable).
       @param args: Arguments to apply (list, tuple, ...).

       @return: Function results (tuple).
    '''
    return tuple(map(func, *args))


def radiansPI(deg):
    '''Converts and wraps degrees to radians M{-PI..+PI}.

       @param deg: Angle (degrees).

       @return: Radians, wrapped (radiansPI)
    '''
    return _wrap(radians(deg), PI)


def radiansPI2(deg):
    '''Converts and wraps degrees to radians M{0..+2PI}.

       @param deg: Angle (degrees).

       @return: Radians, wrapped (radiansPI2)
    '''
    return _wrap(radians(deg), PI2)


def radiansPI_2(deg):
    '''Converts and wraps degrees to radians M{-3PI/2..+PI/2}.

       @param deg: Angle (degrees).

       @return: Radians, wrapped (radiansPI_2)
    '''
    return _wrap(radians(deg), PI_2)


def tanPI_2_2(rad):
    '''Computes tan of half angle, rotated.

       @param rad: Angle (radians).

       @return: M{tan((rad + PI/2) / 2)} (float).
    '''
    return tan((rad + PI_2) * 0.5)


def wrap90(deg):
    '''Wraps degrees to M{-270..+90}.

       @param deg: Angle (degrees).

       @return: Degrees, wrapped (degrees90).
    '''
    return _drap(deg, 90)


def wrap180(deg):
    '''Wraps degrees to M{-180..+180}.

       @param deg: Angle (degrees).

       @return: Degrees, wrapped (degrees180).
    '''
    return _drap(deg, 180)


def _wrap(rad, wrap):
    '''(INTERNAL) Radians wrapper M{(wrap-2PI)..+wrap}.

       @param rad: Angle (radians).
       @param wrap: Limit (radians).

       @returns: Radians, wrapped (radians).
    '''
    r = rad % PI2  # -1.5 % 3.14 == 1.64
    if r > wrap:
        r -= PI2
    return r


def wrapPI(rad):
    '''Wraps radians to M{-PI..+PI}.

       @param rad: Angle (radians).

       @return: Radians, wrapped (radiansPI).
    '''
    return _wrap(rad, PI)


def wrapPI2(rad):
    '''Wraps radians to M{0..+2PI}.

       @param rad: Angle (radians).

       @return: Radians, wrapped (radiansPI2).
    '''
    return _wrap(rad, PI2)


def wrapPI_2(rad):
    '''Wraps radians to M{-3PI/2..+PI/2}.

       @param rad: Angle (radians).

       @return: Radians, wrapped (radiansPI_2).
    '''
    return _wrap(rad, PI_2)

# **) MIT License
#
# Copyright (C) 2016-2017 -- mrJean1@Gmail.com
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
