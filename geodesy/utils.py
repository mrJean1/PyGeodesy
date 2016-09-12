
# -*- coding: utf-8 -*-

# Common base classes and functions.

# After (C) Chris Veness 2011-2015 published under the same MIT Licence,
# see <http://www.movable-type.co.uk/scripts/latlong.html>
# and <http://www.movable-type.co.uk/scripts/latlong-vectors.html>

from math import degrees, pi as PI, radians, sin, sqrt, tan
from operator import mul
import sys

# all public contants, classes and functions
__all__ = ('EPS', 'EPS1', 'EPS2', 'PI', 'PI2', 'PI_2',  # math constants
           'degrees', 'degrees90', 'degrees180', 'degrees360',
           'cbrt', 'fdot', 'fsum', 'fStr', 'hypot3', 'isscalar',
           'len2', 'sin_2', 'tanPI_2_2', 'radians', 'radiansPI',
           'radiansPI_2', 'wrapPI', 'wrapPI2', 'wrapPI_2')
__version__ = '16.09.03'

try:
    from math import fsum  # precision sum
except ImportError:
    fsum = sum  # use standard, built-in sum
try:
    _Scalars = int, long, float
except NameError:  # Python 3+
    _Scalars = int, float
try:
    EPS = sys.float_info.epsilon
except AttributeError:
    EPS = 2.2204460492503131e-16
EPS1 = 1 - EPS
EPS2 = sqrt(EPS)
PI2  = PI * 2  # PYCHOK expected
PI_2 = PI * 0.5
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


def _degrees(rad, wrap):
    d = degrees(rad) % 360  # -1.5 % 360 == 358.5
    if d > wrap:
        d -= 360
    return d


def degrees90(rad):
    '''Convert radians to degrees -90..+90.

       @param {radians} rad - Angle in radians.

       @returns {degrees90} Degrees -90..+90.
    '''
    return _degrees(rad, 90)


def degrees180(rad):
    '''Convert radians to degrees -180..+180.

       @param {radians} rad - Angle in radians.

       @returns {degrees180} Degrees -180..+180.
    '''
    return _degrees(rad, 180)


def degrees360(rad):
    '''Convert radians to degrees 0..360.

       @param {radians} rad - Angle in radians.

       @returns {degrees360} Degrees 0..360.
    '''
    return _degrees(rad, 360)


def fdot(a, *b):
    '''Precision dot product.

       @param {numbers} a - List or tuple of numbers.
       @param {numbers} b - List or tuple of numbers.

       @returns {number} Dot product.
    '''
    assert len(a) == len(b)
    return fsum(map(mul, a, b))


def fStr(floats, prec=6, sep=', '):
    '''Convert floats to string with zero decimals stripped.

       @param {float[]} floats - List of floating point numbers.
       @param {number} [prec=6] - Number of decimals, unstripped.
       @param {string} [sep=', '] - Separator to join.

       @returns {string} Floats as '[f, f, ... f]' string.
    '''
    def _fstr(p, f):
        t = '%.*f' % (abs(p), f)
        if p > 1 and t.endswith('0'):
            z = len(t) - p + 1
            t = t[:z] + t[z:].rstrip('0')
        return t

    if isscalar(floats):
        floats = (floats,)

    return sep.join(_fstr(prec, f) for f in floats)


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


def isscalar(inst):
    '''Check for scalar.

       @param {object} inst - any object.

       @returns {bool} True if inst is scalar.
    '''
    return isinstance(inst, _Scalars)


def len2(geniter):
    '''Make len() work for generators, iterators, etc.
       and return a list of the items since it is not
       possible to restart a generator, iterator, etc.
    '''
    if not isinstance(geniter, (list, tuple)):
        geniter = list(geniter)
    return len(geniter), geniter


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
