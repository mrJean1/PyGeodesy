
# -*- coding: utf-8 -*-

u'''Utility, geodetic/geometric functions and constants.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**, see
U{Latitude/Longitude<http://www.movable-type.co.uk/scripts/latlong.html>} and
U{Vector-based geodesy<http://www.movable-type.co.uk/scripts/latlong-vectors.html>}.

@newfield example: Example, Examples
'''

# make sure int division yields float quotient
from __future__ import division

from fmath import _Seqs, EPS, fStr, fsum_, hypot, len2, map2

from math import atan2, cos, degrees, pi as PI, \
                 radians, sin, sqrt, tan  # pow

# all public contants, classes and functions
__all__ = ('PI', 'PI2', 'PI_2', 'R_M',  # constants
           'CrossError',  'LimitError',  # classes
           'antipode',
           'classname', 'crosserrors',
           'degrees', 'degrees90', 'degrees180', 'degrees360',
           'enStr2',
           'equirectangular', 'equirectangular_',
           'false2f', 'ft2m',
           'halfs',
           'haversine', 'haversine_',  # XXX removed 'hsin', 'hsin3',
           'heightOf', 'horizon',
           'inStr', 'isantipode', 'issequence',
           'isNumpy2', 'isPoints2', 'isTuple2',
           'iterNumpy2', 'iterNumpy2over',
           'limiterrors',
           'm2ft', 'm2km', 'm2NM', 'm2SM',
           'polygon',
           'radians', 'radiansPI_2', 'radiansPI', 'radiansPI2',
           'tan_2', 'tanPI_2_2',
           'unroll180', 'unrollPI', 'unStr',
           'wrap90', 'wrap180', 'wrap360',
           'wrapPI_2', 'wrapPI', 'wrapPI2')
__version__ = '18.03.04'

division = 1 / 2  # double check int division, see datum.py
if not division:
    raise ImportError('%s 1/2 == %d' % ('division', division))
del division

PI2  = PI * 2  #: Two PI, M{PI * 2} (float)  # PYCHOK expected
PI_2 = PI / 2  #: Half PI, M{PI / 2} (float)

# R_M moved here to avoid circular import for bases and datum
R_M = 6371008.771415  #: Mean, spherical earth radius (meter).

_crosserrors   = True
_iterNumpy2len = 1  # adjustable for testing purposes
_limiterrors   = True


class CrossError(ValueError):
    '''Error raised for zero or near-zero cross products or for
       coincident or colinear points or paths.
    '''
    pass


class LimitError(ValueError):
    '''Error raised for lat- and/or longitudinal deltas exceeding
       the I{limit} in functions L{equirectangular} and
       L{equirectangular_}.
    '''
    pass


def antipode(lat, lon):
    '''Return the antipode, the point diametrically opposite to the
       given lat-/longitude.

       @param lat: Latitude (degrees).
       @param lon: Longitude (degrees).

       @return: 2-Tuple (lat, lon) of the antipodal points in degrees.

       @see: U{Geosphere<http://cran.r-project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return -lat, _drap(lon + 180, 180)


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
    '''Get/set raising of cross product errors.

       @keyword raiser: Choose True to raise or False to not raise
                        L{CrossError} exceptions.  Use None to
                        leave the setting unchanged.

       @return: Previous setting (bool).
    '''
    global _crosserrors
    t = _crosserrors
    if raiser in (True, False):
        _crosserrors = raiser
    return t


def degrees90(rad):
    '''Convert and wrap radians to degrees M{(-270..+90]}.

       @param rad: Angle (radians).

       @return: Degrees, wrapped (degrees90).
    '''
    return _drap(degrees(rad), 90)


def degrees180(rad):
    '''Convert and wrap radians to degrees M{(-180..+180]}.

       @param rad: Angle (radians).

       @return: Degrees, wrapped (degrees180).
    '''
    return _drap(degrees(rad), 180)


def degrees360(rad):
    '''Convert and wrap radians to degrees M{(0..+360]}.

       @param rad: Angle (radians).

       @return: Degrees, wrapped (degrees360).
    '''
    return _drap(degrees(rad), 360)


def _drap(deg, wrap):
    '''(INTERNAL) Degree wrapper M{((wrap-360)..+wrap]}.

       @param deg: Angle (degrees).
       @param wrap: Limit (degrees).

       @return: Degrees, wrapped (degrees).
    '''
    if not wrap >= deg > (wrap - 360):
        # math.fmod(-1.5, 360) == -1.5, but ...
        deg %= 360  # -1.5 % 360 == 358.5
        if deg > wrap:
            deg -= 360
    return deg


def enStr2(easting, northing, prec, *extras):
    '''Return easting, northing string representations.

       @param easting: Easting from false easting (meter).
       @param northing: Northing from from false northing (meter).
       @param prec: Precision in number of digits (int).
       @param extras: Optional leading items (strings).

       @return: I{extras} + 2-Tuple (eastingStr, northingStr).

       @raise ValueError: Invalid I{prec}.
    '''
    w = prec // 2
    try:
        p10 = (1e-4, 1e-3, 1e-2, 1e-1, 1)[w - 1]  # 10**(5 - w)
    except IndexError:
        raise ValueError('%s invalid: %r' % ('prec', prec))
    return extras + ('%0*d' % (w, int(easting * p10)),
                     '%0*d' % (w, int(northing * p10)))


def equirectangular(lat1, lon1, lat2, lon2, radius=R_M, **options):
    '''Compute the distance between two points using
       the U{Equirectangular Approximation / Projection
       <http://www.movable-type.co.uk/scripts/latlong.html>}.

       See function L{equirectangular_} for more details, the
       available I{options} and errors raised.

       @param lat1: Start latitude (degrees).
       @param lon1: Start longitude (degrees).
       @param lat2: End latitude (degrees).
       @param lon2: End longitude (degrees).
       @keyword radius: Optional, mean earth radius (meter).
       @keyword options: Optional keyword arguments for function
                         L{equirectangular_}.

       @return: Distance (meter, same units as I{radius}).

       @see: U{Local, Flat Earth<http://www.edwilliams.org/avform.htm#flat>},
             method L{Ellipsoid.distance2} or function L{haversine} for
             more accurate and/or larger distances.
    '''
    _, dy, dx, _ = equirectangular_(lat1, lon1, lat2, lon2, **options)
    return radians(hypot(dx, dy)) * radius


def equirectangular_(lat1, lon1, lat2, lon2,
                     adjust=True, limit=45, wrap=False):
    '''Compute the distance between two points using
       the U{Equirectangular Approximation / Projection
       <http://www.movable-type.co.uk/scripts/latlong.html>}.

       This approximation is valid for smaller distance of several
       hundred Km or Miles, see the I{limit} keyword argument and
       the L{LimitError}.

       @param lat1: Start latitude (degrees).
       @param lon1: Start longitude (degrees).
       @param lat2: End latitude (degrees).
       @param lon2: End longitude (degrees).
       @keyword adjust: Adjust the wrapped, unrolled longitudinal
                        delta by the cosine of the mean latitude (bool).
       @keyword limit: Optional limit for the lat- and longitudinal
                       deltas (degrees) or None or 0 for unlimited.
       @keyword wrap: Wrap and L{unroll180} longitudes and longitudinal
                      delta (bool).

       @return: 4-Tuple (distance2, delta_lat, delta_lon, lon2_unroll)
                with the distance in degrees squared, the latitudinal
                delta I{lat2}-I{lat1}, the wrapped, unrolled, and
                adjusted longitudinal delta I{lon2}-I{lon1} and the
                unrollment for I{lon2}.  To convert I{distance2} to
                meter, use M{radians(sqrt(distance2)) * radius} where
                I{radius} is the mean earth radius in the desired units,
                for example L{R_M} meter.

       @raise LimitError: If the lat- and/or longitudinal delta exceeds
                          the I{-limit..+limit} range and I{limiterrors}
                          set to True.

       @see: U{Local, Flat Earth<http://www.edwilliams.org/avform.htm#flat>},
             method L{Ellipsoid.distance2} or function L{equirectangular}
             for distance only and function L{haversine} for accurate
             and/or larger distances.
    '''
    d_lat = lat2 - lat1
    d_lon, ulon2 = unroll180(lon1, lon2, wrap=wrap)

    if limit and _limiterrors \
             and max(abs(d_lat), abs(d_lon)) > limit > 0:
        t = fStr((lat1, lon1, lat2, lon2), prec=4)
        raise LimitError('%s(%s, limit=%s) delta exceeds limit' %
                        ('equirectangular_', t, fStr(limit, prec=2)))

    if adjust:  # scale delta lon
        d_lon *= cos(radians(lat1 + lat2) * 0.5)

    d2 = d_lat**2 + d_lon**2  # degrees squared!
    return d2, d_lat, d_lon, ulon2 - lon2


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


def ft2m(feet):
    '''Convert feet to meter (m).

       @param feet: Value in feet (scalar).

       @return: Value in m (float).
    '''
    return feet * 0.3048


def halfs(str2):
    '''Split a string in 2 halfs.

       @param str2: String to split (string).

       @return: 2-Tuple (1st, 2nd) half (strings).

       @raise ValueError: Zero or odd I{len(str2)}.
    '''
    h, r = divmod(len(str2), 2)
    if r or not h:
        raise ValueError('%s invalid: %r' % ('str2', str2))
    return str2[:h], str2[h:]


def haversine(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two points using the U{Haversine
       <http://www.movable-type.co.uk/scripts/latlong.html>} formula.

       @param lat1: Start latitude (degrees).
       @param lon1: Start longitude (degrees).
       @param lat2: End latitude (degrees).
       @param lon2: End longitude (degrees).
       @keyword radius: Optional, mean earth radius (meter).
       @keyword wrap: Wrap and L{unrollPI} longitudes (bool).

       @return: Distance (meter, same units as I{radius}).

       @see: U{Distance between two points
             <http://www.edwilliams.org/avform.htm#Dist>} and
             function L{equirectangular} for an approximation.
    '''
    d, lon2 = unroll180(lon1, lon2, wrap=wrap)
    r = haversine_(radians(lat2), radians(lat1), radians(d))
    return r * float(radius)


def haversine_(a2, a1, b21):
    '''Compute the I{angular} distance between two points using the
       U{Haversine<http://www.movable-type.co.uk/scripts/latlong.html>}
       formula.

       @param a2: End latitude (radians).
       @param a1: Start latitude (radians).
       @param b21: Longitudinal delta, M{end-start} (radians).

       @return: Angular distance (radians).

       @see: Function L{haversine}.
    '''
    def _hsin(rad):
        return sin(rad * 0.5)**2

    h = _hsin(a2 - a1) + cos(a1) * cos(a2) * _hsin(b21)  # haversine
    try:
        r = atan2(sqrt(h), sqrt(1 - h)) * 2  # == asin(sqrt(h)) * 2
    except ValueError:
        r = 0 if h < 0.5 else PI
    return r


def heightOf(angle, distance, radius=R_M):
    '''Determine the height above the (spherical) earth after
       traveling along a straight line at a given tilt.

       @param angle: Tilt angle above horizontal (degrees).
       @param distance: Distance along the line (meter or same units as I{radius}).
       @keyword radius: Optional mean earth radius (meter).

       @return: Height (meter, same units as I{distance} and I{radius}).

       @raise ValueError: Invalid I{angle}, I{distance} or I{radius}.

       @see: U{MultiDop GeogBeamHt<http://github.com/nasa/MultiDop/>}
             (U{Shapiro et al. 2009, JTECH
             <http://journals.ametsoc.org/doi/abs/10.1175/2009JTECHA1256.1>}
             and U{Potvin et al. 2012, JTECH
             <http://journals.ametsoc.org/doi/abs/10.1175/JTECH-D-11-00019.1>}).
    '''
    d, r = distance, radius
    if d > r:
        d, r = r, d

    if d > EPS:
        d = d / float(r)
        s = sin(radians(angle))
        s = fsum_(1, 2 * s * d, d**2)
        if s > 0:
            return r * sqrt(s) - float(radius)

    raise ValueError('%s%r' % ('heightOf', (angle, distance, radius)))


def horizon(height, radius=R_M, refraction=False):
    '''Determine the distance to the horizon from a given altitude
       above the (spherical) earth.

       @param height: Altitude (meter or same units as I{radius}).
       @keyword radius: Optional mean earth radius (meter).
       @keyword refraction: Consider atmospheric refraction (bool).

       @return: Distance (meter, same units as I{height} and I{radius}).

       @raise ValueError: Invalid I{height} or I{radius}.

       @see: U{Distance to horizon<http://www.edwilliams.org/avform.htm#Horizon>}.
    '''
    if min(height, radius) < 0:
        raise ValueError('%s%r' % ('horizon', (height, radius)))

    if refraction:
        d2 = 2.415750694528 * height * radius  # 2.0 / 0.8279
    else:
        d2 = height * fsum_(radius, radius, height)
    return sqrt(d2)


def inStr(inst, *args, **kwds):
    '''Return the string representation of an instance.

       @param inst: The instance (any type).
       @param args: Optional positional arguments (tuple).
       @keyword kwds: Optional keyword arguments (dict).

       @return: Representation (string).
    '''
    return unStr(classname(inst), *args, **kwds)


def isantipode(lat1, lon1, lat2, lon2, eps=EPS):
    '''Check whether two points are anitpodal, on diametrically
       opposite sides of the earth.

       @param lat1: Latitude of one point (degrees).
       @param lon1: Longitude of one point (degrees).
       @param lat2: Latitude of the other point (degrees).
       @param lon2: Longitude of the other point (degrees).
       @keyword eps: Tolerance for near-equality (degrees).

       @return: True if points are antipodal within the given
                tolerance, False otherwise.

       @see: U{Geosphere<http://cran.r-project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return abs( lat1 + lat2) < eps and \
           abs((lon1 - lon2) % 360 - 180) < eps


def isNumpy2(obj):
    '''Check for an I{Numpy2LatLon} points wrapper.

       @param obj: The object (any).

       @return: True if obj is an I{Numpy2LatLon}
                instance, False otherwise (bool).
    '''
    # isinstance(self, (Numpy2LatLon, ...))
    return getattr(obj, 'isNumpy2', False)


def isPoints2(obj):
    '''Check for an I{LatLon2psxy} points wrapper.

       @param obj: The object (any).

       @return: True if obj is an I{LatLon2psxy}
                instance, False otherwise (bool).
    '''
    # isinstance(self, (LatLon2psxy, ...))
    return getattr(obj, 'isPoints2', False)


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
    '''Check for an I{Tuple2LatLon} points wrapper.

       @param obj: The object (any).

       @return: True if obj is an I{Tuple2LatLon}
                instance, False otherwise (bool).
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


def limiterrors(raiser=None):
    '''Get/set the raising of limit errors.

       @keyword raiser: Choose True to raise or False to not raise
                        L{LimitError} exceptions.  Use None to leave
                        the setting unchanged.

       @return: Previous setting (bool).
    '''
    global _limiterrors
    t = _limiterrors
    if raiser in (True, False):
        _limiterrors = raiser
    return t


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


def polygon(points, closed=True, base=None):
    '''Check a polygon given as an array, list, sequence, set or
       tuple of points.

       @param points: The points of the polygon (I{LatLon}[])
       @keyword closed: Optionally, treat polygon as closed and remove
                        any duplicate or closing final I{points} (bool).
       @keyword base: Optional I{points} base class (None).

       @return: 2-Tuple (number, sequence) of points (int, sequence).

       @raise TypeError: Some I{points} are not I{LatLon}.

       @raise ValueError: Insufficient number of I{points}.
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
    '''Convert and wrap degrees to radians M{(-PI..+PI]}.

       @param deg: Angle (degrees).

       @return: Radians, wrapped (radiansPI)
    '''
    return _wrap(radians(deg), PI)


def radiansPI2(deg):
    '''Convert and wrap degrees to radians M{(0..+2PI]}.

       @param deg: Angle (degrees).

       @return: Radians, wrapped (radiansPI2)
    '''
    return _wrap(radians(deg), PI2)


def radiansPI_2(deg):
    '''Convert and wrap degrees to radians M{(-3PI/2..+PI/2]}.

       @param deg: Angle (degrees).

       @return: Radians, wrapped (radiansPI_2)
    '''
    return _wrap(radians(deg), PI_2)


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


def unroll180(lon1, lon2, wrap=True):
    '''Unroll longitudinal delta and wrap longitude in degrees.

       @param lon1: Start longitude (degrees).
       @param lon2: End longitude (degrees).
       @keyword wrap: Wrap and unroll to the M{(-180..+180]} range (bool).

       @return: 2-Tuple (delta I{lon2}-I{lon1}, I{lon2}) unrolled
                (degrees, degrees).

       @see: Capability I{LONG_UNROLL} in U{GeographicLib
       <http://geographiclib.sourceforge.io/html/python/interface.html#outmask>}.
    '''
    d = lon2 - lon1
    if wrap and abs(d) > 180:
        u = _drap(d, 180)
        if u != d:
            return u, lon1 + u
    return d, lon2


def unrollPI(rad1, rad2, wrap=True):
    '''Unroll longitudinal delta and wrap longitude in radians.

       @param rad1: Start longitude (radians).
       @param rad2: End longitude (radians).
       @keyword wrap: Wrap and unroll to the M{(-PI..+PI]} range (bool).

       @return: 2-Tuple (delta I{rad2}-I{rad1}, I{rad2}) unrolled
                (radians, radians).

       @see: Capability I{LONG_UNROLL} in U{GeographicLib
       <http://geographiclib.sourceforge.io/html/python/interface.html#outmask>}.
    '''
    r = rad2 - rad1
    if wrap and abs(r) > PI:
        u = _wrap(r, PI)
        if u != r:
            return u, rad1 + u
    return r, rad2


def unStr(name, *args, **kwds):
    '''Return the string representation of an invokation.

       @param name: Function, method or class name (string).
       @param args: Optional positional arguments (tuple).
       @keyword kwds: Optional keyword arguments (dict).

       @return: Representation (string).
    '''
    t = tuple('%s=%s' % t for t in sorted(kwds.items()))
    if args:
        t = map2(str, args) + t
    return '%s(%s)' % (name, ', '.join(t))


def wrap90(deg):
    '''Wrap degrees to M{(-270..+90]}.

       @param deg: Angle (degrees).

       @return: Degrees, wrapped (degrees90).
    '''
    return _drap(deg, 90)


def wrap180(deg):
    '''Wrap degrees to M{(-180..+180]}.

       @param deg: Angle (degrees).

       @return: Degrees, wrapped (degrees180).
    '''
    return _drap(deg, 180)


def wrap360(deg):
    '''Wrap degrees to M{(0..+360]}.

       @param deg: Angle (degrees).

       @return: Degrees, wrapped (degrees360).
    '''
    return _drap(deg, 360)


def _wrap(rad, wrap):
    '''(INTERNAL) Radians wrapper M{((wrap-2PI)..+wrap]}.

       @param rad: Angle (radians).
       @param wrap: Range (radians).

       @return: Radians, wrapped (radians).
    '''
    if not wrap >= rad > (wrap - PI2):
        # math.fmod(-1.5, 3.14) == -1.5, but ...
        rad %= PI2  # ... -1.5 % 3.14 == 1.64
        if rad > wrap:
            rad -= PI2
    return rad


def wrapPI(rad):
    '''Wrap radians to M{(-PI..+PI]}.

       @param rad: Angle (radians).

       @return: Radians, wrapped (radiansPI).
    '''
    return _wrap(rad, PI)


def wrapPI2(rad):
    '''Wrap radians to M{(0..+2PI]}.

       @param rad: Angle (radians).

       @return: Radians, wrapped (radiansPI2).
    '''
    return _wrap(rad, PI2)


def wrapPI_2(rad):
    '''Wrap radians to M{(-3PI/2..+PI/2]}.

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
