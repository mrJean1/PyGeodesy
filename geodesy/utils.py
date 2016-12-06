
# -*- coding: utf-8 -*-

# Common base classes and functions.

# After (C) Chris Veness 2011-2015 published under the same MIT Licence**,
# see <http://www.movable-type.co.uk/scripts/latlong.html>
# and <http://www.movable-type.co.uk/scripts/latlong-vectors.html>

from math import degrees, pi as PI, radians, sin, sqrt, tan
from operator import mul
import sys

# all public contants, classes and functions
__all__ = ('EPS', 'EPS1', 'EPS2', 'PI', 'PI2', 'PI_2',  # constants
           'cbrt',
           'degrees', 'degrees90', 'degrees180', 'degrees360',
           'fdot', 'fdot3', 'fStr', 'fsum',
           'halfs', 'hypot1', 'hypot3',
           'isint', 'isscalar', 'len2', 'map2',
           'radians', 'radiansPI', 'radiansPI_2',
           'sin_2', 'tanPI_2_2',
           'wrap90', 'wrap180', 'wrapPI', 'wrapPI2', 'wrapPI_2')
__version__ = '16.12.04'

try:
    from math import fsum  # precision sum
except ImportError:
    fsum = sum  # use standard, built-in sum (or Kahan's sum
    # <https://en.wikipedia.org/wiki/Kahan_summation_algorithm> or
    # Hettinger's <https://code.activestate.com/recipes/393090/>)

try:
    _Ints = int, long
    _Scalars = int, long, float
except NameError:  # Python 3+
    _Ints = int
    _Scalars = int, float

try:
    EPS = sys.float_info.epsilon
except AttributeError:
    EPS = 2.2204460492503131e-16
EPS1 = 1 - EPS
EPS2 = sqrt(EPS)

PI2  = PI * 2  # PYCHOK expected
PI_2 = PI / 2

_3rd = 1.0 / 3.0  # float!


def cbrt(x):
    '''Return the cubic root.
    '''
    # simpler and more accurate than Ken Turkowski's CubeRoot, see
    # <http://people.freebsd.org/~lstewart/references/apple_tr_kt32_cuberoot.pdf>
    if x < 0:
        return -(float(-x) ** _3rd)
    else:
        return float(x) ** _3rd


def degrees90(rad):
    '''Convert radians to degrees -90..+90.

       @param {radians} rad - Angle in radians.

       @returns {degrees90} Degrees -90..+90.
    '''
    return _drap(degrees(rad), 90)


def degrees180(rad):
    '''Convert radians to degrees -180..+180.

       @param {radians} rad - Angle in radians.

       @returns {degrees180} Degrees -180..+180.
    '''
    return _drap(degrees(rad), 180)


def degrees360(rad):
    '''Convert radians to degrees 0..360.

       @param {radians} rad - Angle in radians.

       @returns {degrees360} Degrees 0..360.
    '''
    return _drap(degrees(rad), 360)


def _drap(deg, wrap):
    d = deg % 360  # -1.5 % 360 == 358.5
    if d > wrap:
        d -= 360
    return d


def fdot(a, *b):
    '''Precision dot product.

       @param {numbers} a - List or tuple of numbers.
       @param {numbers} b - List or tuple of numbers.

       @returns {number} Dot product.
    '''
    assert len(a) == len(b)
    return fsum(map(mul, a, b))


def fdot3(a, b, c, start=0):
    '''Precision dot product.

       @param {numbers} a - List or tuple of numbers.
       @param {numbers} b - List or tuple of numbers.
       @param {numbers} c - List or tuple of numbers.
       @param {number} [start=0] - Optional start number.

       @returns {number} Dot product.
    '''
    def mul3(a, b, c):  # map function
        return a * b * c

    assert len(a) == len(b) == len(c)
    m3 = map(mul3, a, b, c)
    if start:
        m3 = (start,) + tuple(m3)
    return fsum(m3)


def fStr(floats, prec=6, sep=', ', fmt='%.*f', ints=False):
    '''Convert floats to string with zero decimals stripped.

       @param {float[]} floats - List of floating point numbers.
       @param {number} [prec=6] - Number of decimals, unstripped.
       @param {string} [sep=', '] - Separator to join.
       @param {string} [fmt='%.*f'] - Float format.

       @returns {string} Floats as '[f, f, ... f]' string.
    '''
    def _fstr(p, f):
        t = fmt % (abs(p), f)
        if p > 1 and t.endswith('0'):
            z = len(t) - p + 1
            t = t[:z] + t[z:].rstrip('0')
        if ints and isint(f):
            t = t.split('.')[0]
        return t

    if isscalar(floats):
        return _fstr(prec, floats)
    else:
        return sep.join(_fstr(prec, f) for f in floats)


def halfs(str2):
    '''Split string in 2 halfs.
    '''
    h, odd = divmod(len(str2), 2)
    if odd or not h:
        raise ValueError('%s invalid: %r' % ('str2', str2))
    return str2[:h], str2[h:]


def hypot1(x):
    '''Compute the sqrt(1 + x^2).

       @param {number} x - X coordinate.

       @returns {number} Length, norm.
    '''
    h = abs(x)
    if h > 1:
        x = 1.0 / h
        if x > EPS2:
            h *= sqrt(1 + x * x)
    elif h > EPS2:
        h = sqrt(1 + x * x)
    else:
        h = 1.0
    return h


def hypot3(x, y, z):
    '''Compute the sqrt(x^2 + y^2 + z^2).

       @param {number} x - X coordinate.
       @param {number} y - Y coordinate.
       @param {number} z - Z coordinate.

       @returns {number} Length, norm.
    '''
    x, y, z = abs(x), abs(y), abs(z)
    if x < y:
        x, y = y, x
    if x < z:
        x, z = z, x
    x = float(x)
    if x > EPS2:
        y /= x
        z /= x
        h = x * sqrt(1 + y * y + z * z)
    elif x:
        h = EPS2
    else:
        h = 0
    return h


def isint(inst):
    '''Check for integer.

       @param {object} inst - any object.

       @returns {bool} True if inst is integer.
    '''
    return isinstance(inst, _Ints)


def isscalar(inst):
    '''Check for scalar.

       @param {object} inst - any object.

       @returns {bool} True if inst is scalar.
    '''
    return isinstance(inst, _Scalars)


def len2(geniter):
    '''Make len() work for generators, iterators, etc.
       plus return a list of items (since generators,
       iterators, etc. can not be restarted).
    '''
    if not isinstance(geniter, (list, tuple)):
        geniter = list(geniter)
    return len(geniter), geniter


def map2(func, *args):
    '''Apply a function to arguments, like built-in
       map and return a tuple with the results.

       However, Python 3+ map returns a map object, an
       iterator- like object which can generate the map
       result only once.  By converting the map object
       into a tuple, the results can be used any number
       of times.
    '''
    r = map(func, *args)
    if not isinstance(r, tuple):
        r = tuple(r)
    return r


def radiansPI(deg):
    '''Convert degrees to radians -PI..PI.

       @param {degrees} rad - Angle in degrees.

       @returns {radiansPI} Radians -PI..PI.
    '''
    return _wrap(radians(deg), PI)


def radiansPI_2(deg):
    '''Convert degrees to radians -PI_2..PI_2.

       @param {degrees} rad - Angle in degrees.

       @returns {radiansPI_2} Radians -PI_2..PI_2.
    '''
    return _wrap(radians(deg), PI_2)


def sin_2(rad):
    '''Return sin(rad / 2).
    '''
    return sin(rad * 0.5)


def tanPI_2_2(rad):
    '''Return tan((rad + PI/2) / 2).
    '''
    return tan((rad + PI_2) * 0.5)


def wrap90(deg):
    '''Wrap degrees to -90..+90.

       @param {degrees} deg - Angle in degrees.

       @returns {degrees90} Degrees -90..+90.
    '''
    return _drap(deg, 90)


def wrap180(deg):
    '''Wrap degrees to -180..+180.

       @param {degrees} deg - Angle in degrees.

       @returns {degrees180} Degrees -180..+180.
    '''
    return _drap(deg, 180)


def _wrap(rad, wrap):
    r = rad % PI2  # -1.5 % 3.14 == 1.64
    if r > wrap:
        r -= PI2
    return r


def wrapPI(rad):
    '''Wrap radians to -PI..PI.

       @param {number} rad - Angle in radians.

       @returns {radiansPI} Radians -PI..PI.
    '''
    return _wrap(rad, PI)


def wrapPI2(rad):
    '''Wrap radians to -PI2..PI2.

       @param {number} rad - Angle in radians.

       @returns {radiansPI2} Radians -PI2..PI2.
    '''
    return _wrap(rad, PI2)


def wrapPI_2(rad):
    '''Wrap radians to -PI_2..PI_2.

       @param {number} rad - Angle in radians.

       @returns {radiansPI_2} Radians -PI_2..PI_2.
    '''
    return _wrap(rad, PI_2)

# **) MIT License
#
# Copyright (c) 2016-2017 -- mrJean1@Gmail.com
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
