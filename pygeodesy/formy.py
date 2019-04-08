
# -*- coding: utf-8 -*-

u'''Formulary of basic geodesy functions and approximations.

@newfield example: Example, Examples
'''

from fmath import EPS, fStr, fsum_, map1
from lazily import _ALL_LAZY
from utily import PI, PI2, PI_2, R_M, degrees2m, degrees360, \
                  LimitError, _limiterrors, sincos2, \
                  unroll180, unrollPI, wrap90, wrap180

from math import atan2, cos, degrees, hypot, radians, sin, sqrt  # pow

# all public contants, classes and functions
__all__ = _ALL_LAZY.formy + ('equirectangular3',)  # DEPRECATED
__version__ = '19.04.05'


def _scaled(lat1, lat2):  # degrees
    # scale delta lon by cos(mean of lats)
    m = (lat1 + lat2) * 0.5
    return cos(radians(m)) if abs(m) < 90 else 0


def _scaler(rad1, rad2):  # radians, imported by heights.HeighIDW2
    # scale delta lon by cos(mean of lats)
    m = (rad1 + rad2) * 0.5
    return cos(m) if abs(m) < PI_2 else 0


def antipode(lat, lon):
    '''Return the antipode, the point diametrically opposite
       to a given point.

       @param lat: Latitude (C{degrees}).
       @param lon: Longitude (C{degrees}).

       @return: 2-Tuple (lat, lon) of the antipodal point
                (C{degrees}, C{degrees180}).

       @see: U{Geosphere<http://CRAN.R-Project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return -lat, wrap180(lon + 180)


def bearing(lat1, lon1, lat2, lon2, **options):
    '''Compute the initial or final bearing (forward or reverse
       azimuth) between a (spherical) start and end point.

       @param lat1: Start latitude (C{degrees}).
       @param lon1: Start longitude (C{degrees}).
       @param lat2: End latitude (C{degrees}).
       @param lon2: End longitude (C{degrees}).
       @keyword options: Optional keyword arguments for function
                         L{bearing_}.

       @return: Initial or final bearing (compass C{degrees360}) or
                zero if start and end point coincide.
    '''
    ab4 = map1(radians, lat1, lon1, lat2, lon2)
    return degrees(bearing_(*ab4, **options))


def bearing_(a1, b1, a2, b2, final=False, wrap=False):
    '''Compute the initial or final bearing (forward or reverse
       azimuth) between a (spherical) start and end point.

       @param a1: Start latitude (C{radians}).
       @param b1: Start longitude (C{radians}).
       @param a2: End latitude (C{radians}).
       @param b2: End longitude (C{radians}).
       @keyword final: Return final bearing if C{True}, initial
                       otherwise (C{bool}).
       @keyword wrap: Wrap and L{unrollPI} longitudes (C{bool}).

       @return: Initial or final bearing (compass C{radiansPI2}) or
                zero if start and end point coincide.
    '''
    if final:
        a1, b1, a2, b2 = a2, b2, a1, b1
        r = PI + PI2
    else:
        r = PI2

    db, _ = unrollPI(b1, b2, wrap=wrap)
    sa1, ca1, sa2, ca2, sdb, cdb = sincos2(a1, a2, db)

    # see <http://MathForum.org/library/drmath/view/55417.html>
    x = ca1 * sa2 - sa1 * ca2 * cdb
    y = sdb * ca2

    return (atan2(y, x) + r) % PI2  # wrapPI2


def compassAngle(lat1, lon1, lat2, lon2, adjust=True, wrap=False):
    '''Return the angle from North for the direction vector
       M{(lon2 - lon1, lat2 - lat1)} between two points.

       Suitable only for short, non-near-polar vectors up to a few
       hundred Km or Miles.  Use function L{bearing} for longer
       vectors.

       @param lat1: From latitude (C{degrees}).
       @param lon1: From longitude (C{degrees}).
       @param lat2: To latitude (C{degrees}).
       @param lon2: To longitude (C{degrees}).
       @keyword adjust: Adjust the longitudinal delta by the
                        cosine of the mean latitude (C{bool}).
       @keyword wrap: Wrap and L{unroll180} longitudes (C{bool}).

       @return: Compass angle from North (C{degrees360}).

       @note: Courtesy Martin Schultz.

       @see: U{Local, flat earth approximation
             <http://www.EdWilliams.org/avform.htm#flat>}.
    '''
    d_lon, _ = unroll180(lon1, lon2, wrap=wrap)
    if adjust:  # scale delta lon
        d_lon *= _scaled(lat1, lat2)
    return degrees360(atan2(d_lon, lat2 - lat1))


def euclidean(lat1, lon1, lat2, lon2, radius=R_M, adjust=True, wrap=False):
    '''Approximate the C{Euclidian} distance between two (spherical) points.

       @param lat1: Start latitude (C{degrees}).
       @param lon1: Start longitude (C{degrees}).
       @param lat2: End latitude (C{degrees}).
       @param lon2: End longitude (C{degrees}).
       @keyword radius: Optional, mean earth radius (C{meter}).
       @keyword adjust: Adjust the longitudinal delta by the
                        cosine of the mean latitude (C{bool}).
       @keyword wrap: Wrap and L{unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as I{radius}).

       @see: U{Distance between two (spherical) points
             <http://www.EdWilliams.org/avform.htm#Dist>}, functions
             L{equirectangular}, L{haversine} and L{vincentys} and
             methods L{Ellipsoid.distance2}, C{LatLon} I{distanceTo*}
             and I{equirectangularTo}.
    '''
    d, _ = unroll180(lon1, lon2, wrap=wrap)
    r = euclidean_(radians(lat2), radians(lat1), radians(d), adjust=adjust)
    return r * float(radius)


def euclidean_(a2, a1, b21, adjust=True):
    '''Approximate the I{angular} C{Euclidean} distance between two
       (spherical) points.

       @param a2: End latitude (C{radians}).
       @param a1: Start latitude (C{radians}).
       @param b21: Longitudinal delta, M{end-start} (C{radians}).
       @keyword adjust: Adjust the longitudinal delta by the
                        cosine of the mean latitude (C{bool}).

       @return: Angular distance (C{radians}).

       @see: Functions L{euclidean}, L{equirectangular_}, L{haversine_}
             and L{vincentys_}.
    '''
    a, b = abs(a2 - a1), abs(b21)
    if adjust:
        b *= _scaler(a2, a1)
    if a < b:
        a, b = b, a
    return a + b * 0.5


def equirectangular(lat1, lon1, lat2, lon2, radius=R_M, **options):
    '''Compute the distance between two points using
       the U{Equirectangular Approximation / Projection
       <http://www.Movable-Type.co.UK/scripts/latlong.html>}.

       @param lat1: Start latitude (C{degrees}).
       @param lon1: Start longitude (C{degrees}).
       @param lat2: End latitude (C{degrees}).
       @param lon2: End longitude (C{degrees}).
       @keyword radius: Optional, mean earth radius (C{meter}).
       @keyword options: Optional keyword arguments for function
                         L{equirectangular_}.

       @return: Distance (C{meter}, same units as I{radius}).

       @see: Function L{equirectangular_} for more details, the
             available I{options}, errors, restrictions and other,
             approximate or accurate distance functions.
    '''
    _, dy, dx, _ = equirectangular_(lat1, lon1, lat2, lon2, **options)
    return degrees2m(hypot(dx, dy), radius=radius)


def equirectangular3(lat1, lon1, lat2, lon2, **options):
    '''DEPRECATED, use function C{equirectangular_}.

       @return: 3-Tuple (distance2, delta_lat, delta_lon).
    '''
    return equirectangular_(lat1, lon1, lat2, lon2, **options)[:3]


def equirectangular_(lat1, lon1, lat2, lon2,
                     adjust=True, limit=45, wrap=False):
    '''Compute the distance between two points using
       the U{Equirectangular Approximation / Projection
       <http://www.Movable-Type.co.UK/scripts/latlong.html>}.

       This approximation is valid for short distance of several
       hundred Km or Miles, see the I{limit} keyword argument and
       the L{LimitError}.

       @param lat1: Start latitude (C{degrees}).
       @param lon1: Start longitude (C{degrees}).
       @param lat2: End latitude (C{degrees}).
       @param lon2: End longitude (C{degrees}).
       @keyword adjust: Adjust the wrapped, unrolled longitudinal
                        delta by the cosine of the mean latitude (C{bool}).
       @keyword limit: Optional limit for lat- and longitudinal deltas
                       (C{degrees}) or C{None} or C{0} for unlimited.
       @keyword wrap: Wrap and L{unroll180} longitudes (C{bool}).

       @return: 4-Tuple (distance2, delta_lat, delta_lon, unroll_lon2)
                with the distance in C{degrees squared}, the latitudinal
                delta I{lat2}-I{lat1}, the wrapped, unrolled, and
                adjusted longitudinal delta I{lon2}-I{lon1} and the
                unrollment for I{lon2}.  Use Function L{degrees2m} to
                convert C{degrees squared} to distance in C{meter} as
                M{degrees2m(sqrt(distance2), ...)} or
                M{degrees2m(hypot(delta_lat, delta_lon), ...)}.

       @raise LimitError: If the lat- and/or longitudinal delta exceeds
                          the I{-limit..+limit} range and L{limiterrors}
                          set to C{True}.

       @see: U{Local, flat earth approximation
             <http://www.EdWilliams.org/avform.htm#flat>}, functions
             L{equirectangular}, L{euclidean}, L{haversine} and
             L{vincentys} and methods L{Ellipsoid.distance2}, C{LatLon}
             I{distanceTo*} and I{equirectangularTo}.
    '''
    d_lat = lat2 - lat1
    d_lon, ulon2 = unroll180(lon1, lon2, wrap=wrap)

    if limit and _limiterrors \
             and max(abs(d_lat), abs(d_lon)) > limit > 0:
        t = fStr((lat1, lon1, lat2, lon2), prec=4)
        raise LimitError('%s(%s, limit=%s) delta exceeds limit' %
                        ('equirectangular_', t, fStr(limit, prec=2)))

    if adjust:  # scale delta lon
        d_lon *= _scaled(lat1, lat2)

    d2 = d_lat**2 + d_lon**2  # degrees squared!
    return d2, d_lat, d_lon, ulon2 - lon2


def haversine(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two (spherical) points using the
       U{Haversine<http://www.Movable-Type.co.UK/scripts/latlong.html>}
       formula.

       @param lat1: Start latitude (C{degrees}).
       @param lon1: Start longitude (C{degrees}).
       @param lat2: End latitude (C{degrees}).
       @param lon2: End longitude (C{degrees}).
       @keyword radius: Optional, mean earth radius (C{meter}).
       @keyword wrap: Wrap and L{unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as I{radius}).

       @see: U{Distance between two (spherical) points
             <http://www.EdWilliams.org/avform.htm#Dist>}, functions
             L{equirectangular}, L{euclidean} and L{vincentys} and
             methods L{Ellipsoid.distance2}, C{LatLon} I{distanceTo*}
             and I{equirectangularTo}.
    '''
    d, _ = unroll180(lon1, lon2, wrap=wrap)
    r = haversine_(radians(lat2), radians(lat1), radians(d))
    return r * float(radius)


def haversine_(a2, a1, b21):
    '''Compute the I{angular} distance between two (spherical) points
       using the U{Haversine<http://www.Movable-Type.co.UK/scripts/latlong.html>}
       formula.

       @param a2: End latitude (C{radians}).
       @param a1: Start latitude (C{radians}).
       @param b21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Angular distance (C{radians}).

       @see: Functions L{haversine}, L{equirectangular_}, L{euclidean_}
             and L{vincentys_}.
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

       @param angle: Tilt angle above horizontal (C{degrees}).
       @param distance: Distance along the line (C{meter} or same units
                        as I{radius}).
       @keyword radius: Optional mean earth radius (C{meter}).

       @return: Height (C{meter}, same units as I{distance} and I{radius}).

       @raise ValueError: Invalid I{angle}, I{distance} or I{radius}.

       @see: U{MultiDop GeogBeamHt<http://GitHub.com/NASA/MultiDop>}
             (U{Shapiro et al. 2009, JTECH
             <http://Journals.AMetSoc.org/doi/abs/10.1175/2009JTECHA1256.1>}
             and U{Potvin et al. 2012, JTECH
             <http://Journals.AMetSoc.org/doi/abs/10.1175/JTECH-D-11-00019.1>}).
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

       @param height: Altitude (C{meter} or same units as I{radius}).
       @keyword radius: Optional mean earth radius (C{meter}).
       @keyword refraction: Consider atmospheric refraction (C{bool}).

       @return: Distance (C{meter}, same units as I{height} and I{radius}).

       @raise ValueError: Invalid I{height} or I{radius}.

       @see: U{Distance to horizon<http://www.EdWilliams.org/avform.htm#Horizon>}.
    '''
    if min(height, radius) < 0:
        raise ValueError('%s%r' % ('horizon', (height, radius)))

    if refraction:
        d2 = 2.415750694528 * height * radius  # 2.0 / 0.8279
    else:
        d2 = height * fsum_(radius, radius, height)
    return sqrt(d2)


def isantipode(lat1, lon1, lat2, lon2, eps=EPS):
    '''Check whether two points are antipodal, on diametrically
       opposite sides of the earth.

       @param lat1: Latitude of one point (C{degrees}).
       @param lon1: Longitude of one point (C{degrees}).
       @param lat2: Latitude of the other point (C{degrees}).
       @param lon2: Longitude of the other point (C{degrees}).
       @keyword eps: Tolerance for near-equality (C{degrees}).

       @return: C{True} if points are antipodal within the
                I{eps} tolerance, C{False} otherwise.

       @see: U{Geosphere<http://CRAN.R-Project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return abs(wrap90(lat1) + wrap90(lat2)) < eps and \
           abs(abs(wrap180(lon1) - wrap180(lon2)) % 360 - 180) < eps


def vincentys(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two (spherical) points using
       U{Vincenty's<http://WikiPedia.org/wiki/Great-circle_distance>}
       spherical formula.

       @param lat1: Start latitude (C{degrees}).
       @param lon1: Start longitude (C{degrees}).
       @param lat2: End latitude (C{degrees}).
       @param lon2: End longitude (C{degrees}).
       @keyword radius: Optional, mean earth radius (C{meter}).
       @keyword wrap: Wrap and L{unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as I{radius}).

       @see: Functions L{equirectangular}, L{euclidean} and L{haversine}
             and methods L{Ellipsoid.distance2}, C{LatLon} I{distanceTo*}
             and I{equirectangularTo}.

       @note: See note under L{vincentys_}.
    '''
    d, _ = unroll180(lon1, lon2, wrap=wrap)
    r = vincentys_(radians(lat2), radians(lat1), radians(d))
    return r * float(radius)


def vincentys_(a2, a1, b21):
    '''Compute the I{angular} distance between two (spherical) points using
       U{Vincenty's<http://WikiPedia.org/wiki/Great-circle_distance>}
       spherical formula.

       @param a2: End latitude (C{radians}).
       @param a1: Start latitude (C{radians}).
       @param b21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Angular distance (C{radians}).

       @see: Functions L{vincentys}, L{equirectangular_}, L{euclidean_}
             and L{haversine_}.

       @note: Functions L{vincentys_} and L{haversine_} produce equivalent
              results, but L{vincentys_} is suitable for antipodal points
              and slightly more expensive than L{haversine_} (M{3 cos,
              3 sin, 1 hypot, 1 atan2, 7 mul, 2 add} versus M{2 cos, 2
              sin, 2 sqrt, 1 atan2, 5 mul, 1 add}).
    '''
    sa1, ca1, sa2, ca2, sb21, cb21 = sincos2(a1, a2, b21)

    x = sa1 * sa2 + ca1 * ca2 * cb21
    y = ca1 * sa2 - sa1 * ca2 * cb21
    y = hypot(ca2 * sb21, y)
    return atan2(y, x)

# **) MIT License
#
# Copyright (C) 2016-2019 -- mrJean1 at Gmail dot com
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
