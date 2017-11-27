
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

from math import atan2, cos, degrees, hypot, pi as PI, \
                 radians, sin, sqrt, tan  # pow
from operator import mul
import sys

# all public contants, classes and functions
__all__ = ('EPS', 'EPS1', 'EPS2', 'PI', 'PI2', 'PI_2', 'R_M',  # constants
           'CrossError',  # classes
           'cbrt', 'cbrt2', 'classname', 'crosserrors',
           'degrees', 'degrees90', 'degrees180', 'degrees360',
           'equirectangular3',
           'false2f', 'favg', 'fdot', 'fdot3', 'fmean', 'fpolynomial',
           'fStr', 'fStrzs', 'fsum', 'ft2m',
           'halfs',
           'haversine', 'haversine_',  # XXX removed 'hsin', 'hsin3',
           'hypot', 'hypot1', 'hypot3',
           'inStr',
           'isfinite', 'isint', 'isscalar', 'issequence',
           'isNumpy2', 'isTuple2',
           'iterNumpy2', 'iterNumpy2over',
           'len2',
           'm2ft', 'm2km', 'm2NM', 'm2SM', 'map1', 'map2',
           'polygon',
           'radians', 'radiansPI_2', 'radiansPI', 'radiansPI2',
           'scalar',
           'tan_2', 'tanPI_2_2',
           'wrap90', 'wrap180', 'wrap360',
           'wrapPI_2', 'wrapPI', 'wrapPI2')
__version__ = '17.11.26'

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

PI2  = PI * 2  #: Two PI, M{PI * 2} (float)  # PYCHOK expected
PI_2 = PI / 2  #: Half PI, M{PI / 2} (float)

# R_M moved here to avoid circular import for bases and datum
R_M = 6371008.771415  #: Mean, spherical earth radius (meter).

_1_3rd = 1.0 / 3.0  #: (INTERNAL) One third (float)
_2_3rd = 2.0 / 3.0  #: (INTERNAL) Two third (float)

_crosserrors = True


class CrossError(ValueError):
    '''Error for zero cross product or coincident or colinear
       points or paths.
    '''
    pass


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


def classname(obj):
    '''Build module.class name of this object.

       @param obj: The object (any type).

       @return: Name of module and class (string).
    '''
    n = obj.__class__.__name__
    try:
        m = obj.__module__
        n = '.'.join(m.split('.')[-1:] + [n])
    except AttributeError:
        pass
    return n


def crosserrors(raiser=None):
    '''Get/set cross product exceptions.

       @param raiser: New on or off setting (bool).

       @return: Previous setting (bool).
    '''
    global _crosserrors
    t = _crosserrors
    if raiser in (True, False):
        _crosserrors = raiser
    return t


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


def equirectangular3(lat1, lon1, lat2, lon2, adjust=True, wrap=False):
    '''Compute the distance between two points using
       the U{Equirectangular Approximation/Projection
       <http://www.movable-type.co.uk/scripts/latlong.html>}.

       @param lat1: Latitude1 (degrees).
       @param lon1: Longitude1 (degrees).
       @param lat2: Latitude2 (degrees).
       @param lon2: Longitude2 (degrees).
       @keyword adjust: Optionally, adjust longitudinal delta by the
                        cosine of the mean of the latitudes (bool).
       @keyword wrap: Optionally, keep the longitudinal delta within
                      the -180..+180 range (bool).

       @return: 3-Tuple (distance2, delta_lat, delta_lon) with
                the distance in degrees squared, the latitudinal
                delta lat2-lat1 and the I{adjusted}, I{wrapped}
                longitudinal delta lon2-lon1.  To convert distance2
                to meter, use M{radians(sqrt(distance2)) * radius}
                where radius is the mean earth radius in the desired
                units, for example L{R_M} in meter.

       @see: Function L{haversine} for an accurate distance.
    '''
    d_lat = lat2 - lat1
    d_lon = lon2 - lon1

    if wrap:
        if d_lon > 180:
            d_lon -= 360
        elif d_lon < -180:
            d_lon += 360

    if adjust:  # scale lon
        d_lon *= cos(radians(lat1 + lat2) * 0.5)

    d2 = d_lat**2 + d_lon**2  # degrees squared!
    return d2, d_lat, d_lon


def false2f(value, name='value', false=True):
    '''Convert a false east-/northing to non-negative float.

       @param value: Value to convert (scalar).
       @keyword name: Optional name of the value (string).
       @keyword false: Optionally, value includes false origin (bool).

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
       @keyword f: Optional fraction (scalar).

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


def fmean(floats):
    '''Compute the mean of float values.

       @param floats: Values (float).

       @return: Mean value (float).

       @raise ValueError: No floats.
    '''
    n, floats = len2(floats)
    if n > 0:
        return fsum(floats) / n
    raise ValueError('%s missing: %r' % ('floats', floats))


def fpolynomial(x, *cs):
    '''Evaluate the polynomial M{sum(cs[i] * x**i), i=0..len(cs))}.

       @param x: Polynomial argument (scalar).
       @param cs: Polynomial coeffients (scalars).

       @return: Polynomial value (float).

       @raise TypeError: Argument not scalar.

       @raise ValueError: No coefficients.
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

       @param fstr: Float (string).

       @return: Float (string).
    '''
    if fstr.endswith('0'):
        z = fstr.find('.') + 2  # keep 1st zero decimal
        if z > 1:
            fstr = fstr[:z] + fstr[z:].rstrip('0')
    return fstr


try:  # MCCABE 21
    from math import fsum  # precision IEEE-754 sum, Python 2.6+

    # make sure fsum works as expected (XXX check
    # float.__getformat__('float')[:4] == 'IEEE'?)
    if fsum((1, 1e101, 1, -1e101)) != 2:
        raise ImportError  # no, use fsum below

except ImportError:

    def fsum(iterable):
        '''Precision summation similar to I{math.fsum}.

           @param iterable: Sequence, list, tuple, etc. (scalars).

           @return: Precision sum (float).

           @raise OverflowError: Intermediate sum overflow.

           @raise ValueError: Iterable not finite or otherwise invalid.

           @note: Exception and non-finite handling differs from I{math.fsum}.

           @see: U{Hettinger<http://code.activestate.com/recipes/393090/>},
                 U{Klein<http://link.springer.com/article/10.1007/s00607-005-0139-x>},
                 U{Kahan<http://wikipedia.org/wiki/Kahan_summation_algorithm>},
                 Python 2.6+ file I{Modules/mathmodule.c} and the issue log
                 U{Full precision summation<https://bugs.python.org/issue2819>}.
        '''
        def _signof(x):
            return +1 if x > 0 else (-1 if x < 0 else 0)

        def _2sum(a, b):
            if abs(b) > abs(a):
                a, b = b, a
            s = a + b
            if not isfinite(s):
                raise OverflowError('intermediate %s: %r' % ('fsum', s))
            t = s - a
            return s, b - t

        ps = []
        for a in iterable:
            try:
                a = float(a)
                if not isfinite(a):
                    raise ValueError
            except (TypeError, ValueError):
                raise ValueError('%s invalid: %r' % ('fsum', a))

            i = 0
            for b in ps:
                a, p = _2sum(a, b)
                if p:
                    ps[i] = p
                    i += 1
            ps[i:] = [a]

        if ps:  # sum_exact(ps)
            s = ps.pop()
            while ps:
                s, p = _2sum(s, ps.pop())
                if p:
                    break

            if ps:  # half-even round
                if _signof(p) == _signof(ps.pop()):
                    a, p = _2sum(s, p * 2)
                    if not p:
                        s = a
        else:
            s = 0.0
        return s


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


def haversine(lat1, lon1, lat2, lon2, radius=R_M):
    '''Compute the distance between two points using the U{Haversine
       <http://www.movable-type.co.uk/scripts/latlong.html>} formula.

       @param lat1: Latitude1 (degrees).
       @param lon1: Longitude1 (degrees).
       @param lat2: Latitude2 (degrees).
       @param lon2: Longitude2 (degrees).
       @keyword radius: Optional, mean earth radius (meter).

       @return: Distance (meter, same units as I{radius}).

       @see: U{Distance between two points
             <http://www.edwilliams.org/avform.htm#Dist>}
             and function L{equirectangular3} for an approximation.
    '''
    r = haversine_(radians(lat2), radians(lat1), radians(lon2 - lon1))
    return r * float(radius)


def haversine_(a2, a1, b21):
    '''Compute the I{angular} distance using the U{Haversine
       <http://www.movable-type.co.uk/scripts/latlong.html>} formula.

       @param a2: Latitude2 (radians).
       @param a1: Latitude1 (radians).
       @param b21: Longitudinal delta (radians).

       @return: Angular distance (radians).

       @see: This U{Distance between two points
             <http://www.edwilliams.org/avform.htm#Dist>},
             function L{haversine}.
    '''
    def _hsin(rad):
        return sin(rad / 2)**2

    h = _hsin(a2 - a1) + cos(a1) * cos(a2) * _hsin(b21)  # haversine
    try:
        r = atan2(sqrt(h), sqrt(1 - h)) * 2  # == asin(sqrt(h)) * 2
    except ValueError:
        r = 0 if h < 0.5 else PI
    return r


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
            t = (y / h)**2 + (z / h)**2
            if t > EPS:
                h *= sqrt(1.0 + t)
    else:
        h = hypot(x, y)
    return h


def inStr(inst, *args, **kwds):
    '''Return the string representation of an instance.

       @param inst: The instance (any type).
       @param args: Optional positional arguments (tuple).
       @keyword kwds: Optional keyword arguments (dict).

       @return: Representation (string).
    '''
    t = tuple('%s=%s' % t for t in sorted(kwds.items()))
    if args:
        t = map2(str, args) + t
    return '%s(%s)' % (classname(inst), ', '.join(t))


try:
    from math import isfinite  # new in Python 3+
except ImportError:
    from math import isinf, isnan

    def isfinite(obj):
        '''Check for I{Inf} and I{NaN} values.

           @param obj: Value (scalar).

           @return: False if I{Inf} or I{NaN}, True otherwise (bool).

           @raise TypeError: Value not scalar.
        '''
        if isscalar(obj):
            return not (isinf(obj) or isnan(obj))
        raise TypeError('%s invalid: %r' % ('isfinite', obj))


def isint(obj, both=False):
    '''Check for integer type or integer value.

       @param obj: The object (any).
       @keyword both: Optionally, check both type and value (bool).

       @return: True if obj is integer (bool).
    '''
    if both and isinstance(obj, float):  # NOT _Scalars!
        try:
            return obj.is_integer()
        except AttributeError:
            return False  # XXX float(int(obj)) == obj?
    return isinstance(obj, _Ints)


def isNumpy2(obj):
    '''Check for I{Numpy2LatLon} points wrapper.

       @param obj: The object (any).

       @return: True if obj is Numpy2 (bool).
    '''
    # isinstance(self, (Numpy2LatLon, ...))
    return getattr(obj, 'isNumpy2', False)


def isscalar(obj):
    '''Check for scalar types.

       @param obj: The object (any).

       @return: True if obj is scalar (bool).
    '''
    return isinstance(obj, _Scalars)


def issequence(obj, *excluded):
    '''Check for sequence types.

       @param obj: The object (any).
       @param excluded: Optional, exclusions (types).

       @note: Excluding tuple implies namedtuple.

       @return: True if obj is a sequence (bool).
    '''
    if excluded:
        return isinstance(obj, _Seqs) and not \
               isinstance(obj, excluded)
    else:
        return isinstance(obj, _Seqs)


def isTuple2(obj):
    '''Check for I{Tuple2LatLon} points wrapper.

       @param obj: The object (any).

       @return: True if obj is Tuple2 (bool).
    '''
    # isinstance(self, (Tuple2LatLon, ...))
    return getattr(obj, 'isTuple2', False)


def iterNumpy2(obj):
    '''Iterate over Numpy2 wrappers or other sequences exceeding
       the threshold.

       @param obj: Points array, list, sequence, set, etc. (any).

       @return: True, do iterate (bool).
    '''
    try:
        return isNumpy2(obj) or len(obj) > _iterNumpy2len
    except TypeError:
        return False


_iterNumpy2len = 1  # adjustable for testing purposes


def iterNumpy2over(n=None):
    '''Get or set the L{iterNumpy2} threshold.

       @keyword n: Optional, new threshold (integer).

       @return: Previous threshold (integer).
    '''
    global _iterNumpy2len
    p = _iterNumpy2len
    if n is not None:
        try:
            i = int(n)
            if i > 0:
                _iterNumpy2len = i
            else:
                raise ValueError
        except (TypeError, ValueError):
            raise ValueError('%s invalid: %r' % ('n', n))
    return p


def len2(seq):
    '''Make built-in function L{len} work for generators,
       iterators, etc. since those can only be started once.

       @param seq: Generator, iterator, list, range, tuple, etc.

       @return: 2-Tuple (number, list) of items (int, list).
    '''
    if not isinstance(seq, _Seqs):
        seq = list(seq)
    return len(seq), seq


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
    return meter * 6.21369949e-4  # XXX 6.213712e-4 == 1.0 / 1609.344


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

       Unlike Python 2 built-in L{map}, Python 3+ L{map} returns a L{map}
       object, an iterator-like object which generates the results only
       once.  Converting the L{map} object to a tuple maintains Python 2
       behavior.

       @param func: Function to apply (callable).
       @param args: Arguments to apply (list, tuple, ...).

       @return: N-Tuple of function results (tuple).
    '''
    return tuple(map(func, *args))


def polygon(points, closed=True, base=None):
    '''Check a polygon given as an array, list, sequence, set or
       tuple of points.

       @param points: The points of the polygon (I{LatLon}[])
       @keyword closed: Optionally, treat polygon as closed and remove
                        any duplicate or closing final points (bool).
       @keyword base: Optional points base class (None).

       @return: 2-Tuple (number, sequence) of points (int, sequence).

       @raise TypeError: Some points are not I{LatLon}.

       @raise ValueError: Too few points.
    '''
    n, points = len2(points)

    if closed:
        # remove duplicate or closing final points
        while n > 1 and (points[n-1] == points[0] or
                         points[n-1] == points[n-2]):
            n -= 1
        # XXX following line is unneeded if points
        # are always indexed as ... i in range(n)
        points = points[:n]  # XXX numpy.array slice is a view!

    if n < (3 if closed else 1):
        raise ValueError('too few points: %s' % (n,))

    if base and not (isNumpy2(points) or isTuple2(points)):
        for i in range(n):
            base.others(points[i], name='points[%s]' % (i,))

    return n, points


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


def scalar(value, low=EPS, high=1.0, name='scalar'):
    '''Validate a scalar.

       @param value: The value (scalar).
       @keyword low: Optional lower bound (scalar).
       @keyword high: Optional upper bound (scalar).
       @keyword name: Optional name of value (string).

       @return: New value (type(low)).

       @raise TypeError: Value not scalar.

       @raise ValueError: Value out of bounds.
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


def tan_2(rad):
    '''Compute the tangent of half angle.

       @param rad: Angle (radians).

       @return: M{tan(rad / 2)} (float).
    '''
    return tan(rad * 0.5)


def tanPI_2_2(rad):
    '''Compute the tangent of half angle, 90 degrees rotated.

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
