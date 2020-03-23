
# -*- coding: utf-8 -*-

u'''Geometric and other utility functions and constants.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**, see
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>} and
U{Vector-based geodesy<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.

@newfield example: Example, Examples
'''

from pygeodesy.basics import EPS, R_M  # PYCHOK PI_4
from pygeodesy.lazily import _ALL_LAZY

from math import cos, degrees, pi as PI, radians, sin, tan  # pow

# all public contants, classes and functions
__all__ = _ALL_LAZY.utily
__version__ = '20.03.19'

PI2  = PI * 2.0  #: Two PI, M{PI * 2} aka Tau (C{float})  # PYCHOK expected
PI_2 = PI / 2.0  #: Half PI, M{PI / 2} (C{float})
PI_4 = PI / 4.0  #: Quarter PI, M{PI / 4} (C{float})

# <https://Numbers.Computation.Free.FR/Constants/Miscellaneous/digits.html>
_1_90 = 1.0 / 90  # 0.011111111111111111111111111111111111111111111111
_2_PI = 2.0 / PI  # 0.63661977236758134307553505349005744813783858296182

_iterNumpy2len = 1  # adjustable for testing purposes


def degrees90(rad):
    '''Convert radians to degrees and wrap M{[-270..+90]}.

       @arg rad: Angle (C{radians}).

       @return: Angle in degrees, wrapped (C{degrees90}).
    '''
    return _wrap(degrees(rad), 90, 360)


def degrees180(rad):
    '''Convert radians to degrees and wrap M{[-180..+180]}.

       @arg rad: Angle (C{radians}).

       @return: Angle in degrees, wrapped (C{degrees180}).
    '''
    return _wrap(degrees(rad), 180, 360)


def degrees360(rad):
    '''Convert radians to degrees and wrap M{[0..+360)}.

       @arg rad: Angle (C{radians}).

       @return: Angle in degrees, wrapped (C{degrees360}).
    '''
    return _wrap(degrees(rad), 360, 360)


def degrees2m(deg, radius=R_M, lat=0):
    '''Convert angle to distance along the equator or along a
       parallel at an other latitude.

       @arg deg: Angle (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg lat: Parallel latitude (C{degrees90}).

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @raise RangeError: Latitude B{C{lat}} outside valid range
                          and L{rangerrors} set to C{True}.
    '''
    m = radians(deg) * radius
    if lat:
        from pygeodesy.dms import clipDMS
        m *= cos(radians(clipDMS(lat, 90)))
    return m


def false2f(value, name='value', false=True, Error=ValueError):
    '''Convert a false east-/northing to non-negative float.

       @arg value: Value to convert (C{scalar}).
       @kwarg name: Optional name of the value (C{str}).
       @kwarg false: Optionally, value includes false origin (C{bool}).
       @kwarg Error: Exception to raise (C{ValueError}).

       @return: The value (C{float}).

       @raise Error: Invalid or negative B{C{value}}.
    '''
    try:
        f = float(value)
        if f < 0 and false:
            raise ValueError
    except (TypeError, ValueError):
        raise Error('%s invalid: %r' % (name, value))
    return f


def ft2m(feet, usurvey=False):
    '''Convert I{International} or I{US Survey} feet to meter.

       @arg feet: Value in feet (C{scalar}).
       @kwarg usurvery: Convert I{US Survey} feet (C{bool}),
                        I{International} feet otherwise.

       @return: Value in C{meter} (C{float}).
    '''
    # US Survey 1200./3937. == 0.3048006096012192
    return feet * (0.3048006096012192 if usurvey else 0.3048)


def isNumpy2(obj):
    '''Check for an B{C{Numpy2LatLon}} points wrapper.

       @arg obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is an B{C{Numpy2LatLon}}
                instance, C{False} otherwise.
    '''
    # isinstance(self, (Numpy2LatLon, ...))
    return getattr(obj, 'isNumpy2', False)


def isPoints2(obj):
    '''Check for an B{C{LatLon2psxy}} points wrapper.

       @arg obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is an B{C{LatLon2psxy}}
                instance, C{False} otherwise.
    '''
    # isinstance(self, (LatLon2psxy, ...))
    return getattr(obj, 'isPoints2', False)


def isTuple2(obj):
    '''Check for an B{C{Tuple2LatLon}} points wrapper.

       @arg obj: The object (any).

       @return: C{True} if B{C{obj}} is an B{C{Tuple2LatLon}}
                instance, C{False} otherwise.
    '''
    # isinstance(self, (Tuple2LatLon, ...))
    return getattr(obj, 'isTuple2', False)


def iterNumpy2(obj):
    '''Iterate over Numpy2 wrappers or other sequences exceeding
       the threshold.

       @arg obj: Points array, list, sequence, set, etc. (any).

       @return: C{True} do, C{False} don't iterate.
    '''
    try:
        return isNumpy2(obj) or len(obj) > _iterNumpy2len
    except TypeError:
        return False


def iterNumpy2over(n=None):
    '''Get or set the L{iterNumpy2} threshold.

       @kwarg n: Optional, new threshold (C{int}).

       @return: Previous threshold (C{int}).

       @raise ValueError: Invalid B{C{n}}.
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


def m2degrees(meter, radius=R_M):
    '''Convert distance to angle along equator.

       @arg meter: Distance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter}).

       @return: Angle (C{degrees}).

       @raise ValueError: Invalid B{C{radius}}.
    '''
    if radius < EPS:
        raise ValueError('%s invalid: %r' % ('radius', radius))
    return degrees(float(meter) / radius)


def m2ft(meter, usurvey=False):
    '''Convert meter to I{International} or I{US Survey} feet (C{ft}).

       @arg meter: Value in meter (C{scalar}).
       @kwarg usurvery: Convert to I{US Survey} feet (C{bool}),
                        I{International} feet otherwise.

       @return: Value in C{feet} (C{float}).
    '''
    # US Survey == 3937./1200. = 3.2808333333333333
    return meter * (3.2808333333333333 if usurvey else 3.2808399)


def m2km(meter):
    '''Convert meter to kilo meter (km).

       @arg meter: Value in meter (C{scalar}).

       @return: Value in km (C{float}).
    '''
    return meter * 1.0e-3


def m2NM(meter):
    '''Convert meter to nautical miles (NM).

       @arg meter: Value in meter (C{scalar}).

       @return: Value in NM (C{float}).
    '''
    return meter * 5.39956804e-4  # == * 1.0 / 1852.0


def m2SM(meter):
    '''Convert meter to statute miles (SM).

       @arg meter: Value in meter (C{scalar}).

       @return: Value in SM (C{float}).
    '''
    return meter * 6.21369949e-4  # XXX 6.213712e-4 == 1.0 / 1609.344


def radiansPI(deg):
    '''Convert and wrap degrees to radians M{[-PI..+PI]}.

       @arg deg: Angle (C{degrees}).

       @return: Radians, wrapped (C{radiansPI})
    '''
    return _wrap(radians(deg), PI, PI2)


def radiansPI2(deg):
    '''Convert and wrap degrees to radians M{[0..+2PI)}.

       @arg deg: Angle (C{degrees}).

       @return: Radians, wrapped (C{radiansPI2})
    '''
    return _wrap(radians(deg), PI2, PI2)


def radiansPI_2(deg):
    '''Convert and wrap degrees to radians M{[-3PI/2..+PI/2]}.

       @arg deg: Angle (C{degrees}).

       @return: Radians, wrapped (C{radiansPI_2})
    '''
    return _wrap(radians(deg), PI_2, PI2)


def _sincos2(q, r):
    '''(INTERNAL) 2-tuple (C{sin(r), cos(r)}) in quadrant C{q}.
    '''
    if r:
        s, c = sin(r), cos(r)
        t = s, c, -s, -c, s
    else:  # XXX sin(-0.0)?
        t = 0.0, 1.0, -0.0, -1.0, 0.0
    q &= 3
    return t[q], t[q + 1]


def sincos2(*rad):
    '''Return the C{sine} and C{cosine} of angle(s).

       @arg rad: One or more angles (C{radians}).

       @return: The C{sin(rad)} and C{cos(rad)} for each angle.

       @see: U{GeographicLib<https://GeographicLib.SourceForge.io/html/
             classGeographicLib_1_1Math.html#sincosd>} function U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             python/geographiclib/geomath.py#l155>} and C++ U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             include/GeographicLib/Math.hpp#l558>}.
    '''
    for r in rad:
        q = int(r * _2_PI)  # int(math.floor)
        if r < 0:
            q -= 1
        s, c = _sincos2(q, r - q * PI_2)  # 0 <= r < PI_2
        yield s
        yield c


def sincos2d(*deg):
    '''Return the C{sine} and C{cosine} of angle(s) in C{degrees}.

       @arg deg: One or more angles (C{degrees}).

       @return: The C{sin(rad)} and C{cos(rad)} for each angle.

       @see: U{GeographicLib<https://GeographicLib.SourceForge.io/html/
             classGeographicLib_1_1Math.html#sincosd>} function U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             python/geographiclib/geomath.py#l155>} and C++ U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             include/GeographicLib/Math.hpp#l558>}.
    '''
    for d in deg:
        q = int(d * _1_90)  # int(math.floor)
        if d < 0:
            q -= 1
        s, c = _sincos2(q, radians(d - q * 90))  # 0 <= r < PI_2
        yield s
        yield c


def tan_2(rad):
    '''Compute the tangent of half angle.

       @arg rad: Angle (C{radians}).

       @return: M{tan(rad / 2)} (C{float}).
    '''
    return tan(rad * 0.5)


def tanPI_2_2(rad):
    '''Compute the tangent of half angle, 90 degrees rotated.

       @arg rad: Angle (C{radians}).

       @return: M{tan((rad + PI/2) / 2)} (C{float}).
    '''
    return tan((rad + PI_2) * 0.5)


def unroll180(lon1, lon2, wrap=True):
    '''Unroll longitudinal delta and wrap longitude in degrees.

       @arg lon1: Start longitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg wrap: Wrap and unroll to the M{(-180..+180]} range (C{bool}).

       @return: 2-Tuple (B{C{lon2-lon1}}, B{C{lon2}}) unrolled (C{degrees},
                C{degrees}).

       @see: Capability C{LONG_UNROLL} in U{GeographicLib
             <https://GeographicLib.SourceForge.io/html/python/interface.html#outmask>}.
    '''
    d = lon2 - lon1
    if wrap and abs(d) > 180:
        u = _wrap(d, 180, 360)
        if u != d:
            return u, lon1 + u
    return d, lon2


def unrollPI(rad1, rad2, wrap=True):
    '''Unroll longitudinal delta and wrap longitude in radians.

       @arg rad1: Start longitude (C{radians}).
       @arg rad2: End longitude (C{radians}).
       @kwarg wrap: Wrap and unroll to the M{(-PI..+PI]} range (C{bool}).

       @return: 2-Tuple (B{C{rad2 - rad1}}, B{C{rad2}}) unrolled (C{radians},
                C{radians}).

       @see: Capability C{LONG_UNROLL} in U{GeographicLib
             <https://GeographicLib.SourceForge.io/html/python/interface.html#outmask>}.
    '''
    r = rad2 - rad1
    if wrap and abs(r) > PI:
        u = _wrap(r, PI, PI2)
        if u != r:
            return u, rad1 + u
    return r, rad2


def _wrap(angle, wrap, modulo):
    '''(INTERNAL) Angle wrapper M{((wrap-modulo)..+wrap]}.

       @arg angle: Angle (C{degrees} or C{radians}).
       @arg wrap: Range (C{degrees} or C{radians}).
       @arg modulo: Upper limit (360 C{degrees} or PI2 C{radians}).

       @return: The B{C{angle}}, wrapped (C{degrees} or C{radians}).
    '''
    if not wrap > angle >= (wrap - modulo):
        # math.fmod(-1.5, 3.14) == -1.5, but -1.5 % 3.14 == 1.64
        # math.fmod(-1.5, 360) == -1.5, but -1.5 % 360 == 358.5
        angle %= modulo
        if angle > wrap:
            angle -= modulo
    return angle


def wrap90(deg):
    '''Wrap degrees to M{[-270..+90]}.

       @arg deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees90}).
    '''
    return _wrap(deg, 90, 360)


def wrap180(deg):
    '''Wrap degrees to M{[-180..+180]}.

       @arg deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees180}).
    '''
    return _wrap(deg, 180, 360)


def wrap360(deg):
    '''Wrap degrees to M{[0..+360)}.

       @arg deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees360}).
    '''
    return _wrap(deg, 360, 360)


def wrapPI(rad):
    '''Wrap radians to M{[-PI..+PI]}.

       @arg rad: Angle (C{radians}).

       @return: Radians, wrapped (C{radiansPI}).
    '''
    return _wrap(rad, PI, PI2)


def wrapPI2(rad):
    '''Wrap radians to M{[0..+2PI)}.

       @arg rad: Angle (C{radians}).

       @return: Radians, wrapped (C{radiansPI2}).
    '''
    return _wrap(rad, PI2, PI2)


def wrapPI_2(rad):
    '''Wrap radians to M{[-3PI/2..+PI/2]}.

       @arg rad: Angle (C{radians}).

       @return: Radians, wrapped (C{radiansPI_2}).
    '''
    return _wrap(rad, PI_2, PI2)

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
