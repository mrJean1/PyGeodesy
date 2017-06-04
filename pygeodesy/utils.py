
# -*- coding: utf-8 -*-

u'''Utility and mathematical functions and constants.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**,
see U{http://www.movable-type.co.uk/scripts/latlong.html}
and U{http://www.movable-type.co.uk/scripts/latlong-vectors.html}.

@newfield example: Example, Examples
'''

# force int division to yield float quotient
from __future__ import division as _
if not 1/2:  # PYCHOK 1/2 == 0
    raise ImportError('1/2 == %d' % (1/2,))

try:
    from math import fsum  # precision sum, Python 2.6+
except ImportError:
    fsum = sum  # use standard, built-in sum (or Kahan's summation
    # <http://wikipedia.org/wiki/Kahan_summation_algorithm> or
    # Hettinger's <http://code.activestate.com/recipes/393090/>)
from math import atan2, cos, degrees, hypot, \
                 pi as PI, radians, sin, sqrt, tan  # pow
from operator import mul
import sys

# all public contants, classes and functions
__all__ = ('EPS', 'EPS1', 'EPS2', 'PI', 'PI2', 'PI_2', 'R_M',  # constants
           'cbrt', 'cbrt2',
           'degrees', 'degrees90', 'degrees180', 'degrees360',
           'false2f', 'favg', 'fdot', 'fdot3', 'fStr', 'fStrzs',
           'fsum', 'ft2m',
           'halfs', 'hsin', 'hsin3', 'hypot', 'hypot1', 'hypot3',
           'isint', 'isscalar', 'len2',
           'm2ft', 'm2km', 'm2NM', 'm2SM', 'map1', 'map2',
           'radians', 'radiansPI', 'radiansPI2', 'radiansPI_2',
           'tanPI_2_2',
           'wrap90', 'wrap180', 'wrap360',
           'wrapPI', 'wrapPI2', 'wrapPI_2')
__version__ = '17.06.04'

try:  # Luciano Ramalho, "Fluent Python", page 395, O'Reilly, 2016
    from numbers import Real as _Scalars  #: (INTERNAL) Scalar objects
except ImportError:
    try:
        _Scalars = int, long, float  #: (INTERNAL) Scalar objects (tuple)
    except NameError:
        _Scalars = int, float  #: (INTERNAL) Scalar objects (tuple)

try:  # similarly ...
    from numbers import Integral as _Ints  #: (INTERNAL) Int objects
except ImportError:
    try:
        _Ints = int, long  #: (INTERNAL) Int objects (tuple)
    except NameError:  # Python 3+
        _Ints = int  #: (INTERNAL) Int objects (tuple)

try:
    EPS = sys.float_info.epsilon  #: System's epsilon (float)
except AttributeError:
    EPS = 2.220446049250313e-16  #: Approximate epsilon (float)
EPS1 = 1.0 - EPS  #: M{1.0 - EPS} (float), about 0.9999999999999998
EPS2 = sqrt(EPS)  #: M{sqrt(EPS)} (float)

PI2  = PI * 2  #: Two PI, M{PI * 2} (float)  # PYCHOK expected
PI_2 = PI / 2  #: Half PI, M{PI / 2} (float)

# R_M moved here to avoid circular import for bases and datum
R_M = 6371008.771415  #: Mean, spherical earth radius (meter).

_1_3rd = 1.0 / 3.0  #: (INTERNAL) One third (float)
_2_3rd = 2.0 / 3.0  #: (INTERNAL) Two third (float)


def cbrt(x):
    '''Compute the cubic root M{x**(1/3)}.

       @param x: Argument (scalar).

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

       @param x: Argument (scalar).

       @return: Cubic root squared (float).
    '''
    return pow(abs(x), _2_3rd)


def degrees90(rad):
    '''Convert and wrap radians to degrees M{-270..+90}.

       @param rad: Angle (radians).

       @return: Degrees, wrapped (degrees90).
    '''
    return _drap(degrees(rad), 90)


def degrees180(rad):
    '''Convert and wrap radians to degrees M{-180..+180}.

       @param rad: Angle (radians).

       @return: Degrees, wrapped (degrees180).
    '''
    return _drap(degrees(rad), 180)


def degrees360(rad):
    '''Convert and wrap radians to degrees M{0..+360}.

       @param rad: Angle (radians).

       @return: Degrees, wrapped (degrees360).
    '''
    return _drap(degrees(rad), 360)


def _drap(deg, wrap):
    '''(INTERNAL) Degree wrapper M{(wrap-360)..+wrap}.

       @param deg: Angle (degrees).
       @param wrap: Limit (degrees).

       @return: Degrees, wrapped (degrees).
    '''
    d = deg % 360  # -1.5 % 360 == 358.5
    if d > wrap:
        d -= 360
    return d


def false2f(value, name='value', false=True):
    '''Convert a false east-/northing to non-negative float.

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
    '''Return the weighted average of two values.

       @param v1: One value (scalar).
       @param v2: Other value (scalar).
       @keyword f: Fraction (scalar).

       @return: M{v1 + f * (v2 - v1)} (float).
    '''
#      @raise ValueError: Fraction out of range.
#   '''
#   if not 0 <= f <= 1:  # XXX restrict fraction?
#       raise ValueError('%s invalid: %r' % ('fraction', f))
    return v1 + f * (v2 - v1)  # v1 * (1 - f) + v2 * f


def fdot(a, *b):
    '''Return the precision dot product M{sum(a[i] * b[i]
       for i in range(len(a)))}.

       @param a: List, sequence, tuple, etc. (scalars).
       @param b: List, sequence, tuple, etc. (scalars).

       @return: Dot product (float).

       @raise ValueError: Unequal len(a) and len(b).
    '''
    if not len(a) == len(b):
        raise ValueError('%s: %s vs %s' % ('len', len(a), len(b)))

    return fsum(map(mul, a, b))


def fdot3(a, b, c, start=0):
    '''Return the precision dot product M{sum(a[i] * b[i] * c[i]
       for i in range(len(a))) + start}.

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

    m3 = map(_mul3, a, b, c)
    if start:
        m3 = (start,) + tuple(m3)
    return fsum(m3)


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
        if ints and isint(f):
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

       @param fstr: Float (string)

       @return: Float (string).
    '''
    if fstr.endswith('0'):
        z = fstr.find('.') + 2  # keep 1st zero decimal
        if z > 1:
            fstr = fstr[:z] + fstr[z:].rstrip('0')
    return fstr


def ft2m(feet):
    '''Convert feet to meter (m).

       @param feet: Value in feet (scalar).

       @return: Value in m (float).
    '''
    return feet * 0.3048


def halfs(str2):
    '''Split a string in 2 halfs.

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

       @see: U{http://wikipedia.org/wiki/Haversine_formula}.
    '''
    h = sin(rad * 0.5)
    return h * h


def hsin3(a2, a1, b21):
    '''Compute the angular distance using the Haversine formula.

       @param a2: Latitude2 (radians).
       @param a1: Latitude1 (radians).
       @param b21: Longitude delta (radians).

       @return: 3-Tuple (angle, cos(a2), cos(a1))

       @see: U{http://www.movable-type.co.uk/scripts/latlong.html} and
             U{http://www.edwilliams.org/avform.htm#Dist}
    '''
    ca2, ca1 = map1(cos, a2, a1)
    h = hsin(a2 - a1) + ca1 * ca2 * hsin(b21)  # haversine
    try:
        r = atan2(sqrt(h), sqrt(1 - h)) * 2  # == asin(sqrt(h)) * 2
    except ValueError:
        r = 0 if h < 0.5 else PI
    return r, ca2, ca1


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
            # XXX PyChecker chokes on /= and *=
            y = y / h
            z = z / h
            t = y * y + z * z
            if t > EPS:
                h *= sqrt(1.0 + t)
    else:
        h = hypot(x, y)
    return h


def isint(obj, both=False):
    '''Check for integer types and value.

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
    '''Check for scalar types.

       @param obj: The object (any type).

       @return: True if obj is scalar (bool).
    '''
    return isinstance(obj, _Scalars)


def len2(xtor):
    '''Make built-in L{len}() function work for generators,
       iterators, etc. since those can only be started once.

       @param xtor: Generator, iterator, list, sequence, tuple, etc.

       @return: 2-Tuple (number, list) of items (int, list).
    '''
    if not isinstance(xtor, (list, tuple)):
        xtor = list(xtor)
    return len(xtor), xtor


def m2ft(meter):
    '''Convert meter to feet (ft).

       @param meter: Value in meter (scalar).

       @return: Value in ft (float).
    '''
    return meter * 3.2808399


def m2km(meter):
    '''Convert meter to kilo meter (km).

       @param meter: Value in meter (scalar).

       @return: Value in km (float).
    '''
    return meter * 1.0e-3


def m2NM(meter):
    '''Convert meter to nautical miles (NM).

       @param meter: Value in meter (scalar).

       @return: Value in NM (float).
    '''
    return meter * 5.39956804e-4  # == * 1.0 / 1852.0


def m2SM(meter):
    '''Convert meter to statute miles (SM).

       @param meter: Value in meter (scalar).

       @return: Value in SM (float).
    '''
    return meter * 6.21369949e-4  # XXX 6.213712e-4


def map1(func, *args):
    '''Apply each argument to a single-argument function and
       returns a tuple of results.

       @param func: Function to apply (callable).
       @param args: Arguments to apply (any).

       @return: Function results (tuple).
    '''
    return tuple(map(func, args))


def map2(func, *args):
    '''Apply arguments to a function and returns a tuple of results.

       Unlike Python 2 built-in L{map}, Python 3+ L{map} returns a L{map}
       object, an iterator-like object which generates the results only
       once.  Converting the L{map} object to a tuple maintains Python 2
       behavior.

       @param func: Function to apply (callable).
       @param args: Arguments to apply (list, tuple, ...).

       @return: Function results (tuple).
    '''
    return tuple(map(func, *args))


def radiansPI(deg):
    '''Convert and wrap degrees to radians M{-PI..+PI}.

       @param deg: Angle (degrees).

       @return: Radians, wrapped (radiansPI)
    '''
    return _wrap(radians(deg), PI)


def radiansPI2(deg):
    '''Convert and wrap degrees to radians M{0..+2PI}.

       @param deg: Angle (degrees).

       @return: Radians, wrapped (radiansPI2)
    '''
    return _wrap(radians(deg), PI2)


def radiansPI_2(deg):
    '''Convert and wrap degrees to radians M{-3PI/2..+PI/2}.

       @param deg: Angle (degrees).

       @return: Radians, wrapped (radiansPI_2)
    '''
    return _wrap(radians(deg), PI_2)


def tanPI_2_2(rad):
    '''Compute tan of half angle, rotated.

       @param rad: Angle (radians).

       @return: M{tan((rad + PI/2) / 2)} (float).
    '''
    return tan((rad + PI_2) * 0.5)


def wrap90(deg):
    '''Wrap degrees to M{-270..+90}.

       @param deg: Angle (degrees).

       @return: Degrees, wrapped (degrees90).
    '''
    return _drap(deg, 90)


def wrap180(deg):
    '''Wrap degrees to M{-180..+180}.

       @param deg: Angle (degrees).

       @return: Degrees, wrapped (degrees180).
    '''
    return _drap(deg, 180)


def wrap360(deg):
    '''Wrap degrees to M{0..+360}.

       @param deg: Angle (degrees).

       @return: Degrees, wrapped (degrees360).
    '''
    return _drap(deg, 360)


def _wrap(rad, wrap):
    '''(INTERNAL) Radians wrapper M{(wrap-2PI)..+wrap}.

       @param rad: Angle (radians).
       @param wrap: Limit (radians).

       @return: Radians, wrapped (radians).
    '''
    r = rad % PI2  # -1.5 % 3.14 == 1.64
    if r > wrap:
        r -= PI2
    return r


def wrapPI(rad):
    '''Wrap radians to M{-PI..+PI}.

       @param rad: Angle (radians).

       @return: Radians, wrapped (radiansPI).
    '''
    return _wrap(rad, PI)


def wrapPI2(rad):
    '''Wrap radians to M{0..+2PI}.

       @param rad: Angle (radians).

       @return: Radians, wrapped (radiansPI2).
    '''
    return _wrap(rad, PI2)


def wrapPI_2(rad):
    '''Wrap radians to M{-3PI/2..+PI/2}.

       @param rad: Angle (radians).

       @return: Radians, wrapped (radiansPI_2).
    '''
    return _wrap(rad, PI_2)

# **) MIT License
#
# Copyright (C) 2016-2017 -- mrJean1 at Gmail dot com
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
