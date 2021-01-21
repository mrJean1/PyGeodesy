
# -*- coding: utf-8 -*-

u'''Geometric and other utility functions and constants.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**, see
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>} and
U{Vector-based geodesy<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.

@newfield example: Example, Examples
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division

from pygeodesy.basics import copysign, isint
from pygeodesy.errors import _xkwds_get, _TypeError, _ValueError
from pygeodesy.interns import EPS, EPS0, INF, MISSING, PI, PI2, PI_2, R_M, \
                             _0_0, _0_5, _1_0, _90_0, _180_0, _360_0
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.units import Feet, Float, Lam, Lam_, Meter

from math import acos, asin, atan2, cos, degrees, radians, sin, tan  # pow

__all__ = _ALL_LAZY.utily
__version__ = '21.01.14'

# <https://Numbers.Computation.Free.FR/Constants/Miscellaneous/digits.html>
_1_90 = _1_0 / _90_0  # 0.01111111111111111111111111111111111111111111111111
_2_PI = _1_0 /  PI_2  # 0.63661977236758134307553505349005744813783858296182
# sqrt(2) + 1 <https://WikiPedia.org/wiki/Square_root_of_2>
# _R2_1 = 2.41421356237309504880  # _16887_24209_69807_85696_71875_37694_80731_76679_73799

_iterNumpy2len = 1  # adjustable for testing purposes


def acos1(x):
    '''Return C{math.acos(max(-1, min(1, B{x})))}.
    '''
    return acos(x) if abs(x) < _1_0 else (PI if x < 0 else _0_0)


def acre2ha(acres):
    '''Convert acres to hectare.

       @arg acres: Value in acres (C{scalar}).

       @return: Value in C{hectare} (C{float}).

       @raise ValueError: Invalid B{C{acres}}.
    '''
    # 0.40468564224 == acre2m2(1) / 10_000
    return Float(acres) * 0.40468564224


def acre2m2(acres):
    '''Convert acres to I{square} meter.

       @arg acres: Value in acres (C{scalar}).

       @return: Value in C{meter^2} (C{float}).

       @raise ValueError: Invalid B{C{acres}}.
    '''
    # 4046.8564224 == chain2m(1) * furlong2m(1)
    return Float(acres) * 4046.8564224


def asin1(x):
    '''Return C{math.asin(max(-1, min(1, B{x})))}.
    '''
    return asin(x) if abs(x) < _1_0 else (-PI_2 if x < 0 else PI_2)  # -PI_2, not PI3_2!


def atand(y_x):
    '''Return C{atan(B{y_x})} angle in C{degrees}.

       @see: Function L{atan2d}.
    '''
    return atan2d(y_x, 1)


def atan2b(y, x):
    '''Return C{atan2(B{y}, B{x})} in degrees M{[0..+360]}.

       @see: Function L{atan2d}.
    '''
    d = atan2d(y, x)
    if d < 0:
        d += _360_0
    return d


def atan2d(y, x):
    '''Return C{atan2(B{y}, B{x})} in degrees M{[-180..+180]}.

       @see: I{Karney}'s C++ function U{Math.atan2d
             <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Math.html>}.
    '''
    if abs(y) > abs(x) > 0:
        if y < 0:  # q = 3
            d = degrees(atan2(x, -y)) - _90_0
        else:  # q = 2
            d = _90_0 - degrees(atan2(x, y))
    elif x < 0:  # q = 1
        d = copysign(_180_0, y) - degrees(atan2(y, -x))
    elif x > 0:  # q = 0
        d = degrees(atan2(y, x)) if y else _0_0
    else:
        d = -_90_0 if y < 0 else (_90_0 if y > 0 else _0_0)
    return d


def chain2m(chains):
    '''Convert I{UK} chains to meter.

       @arg chains: Value in chains (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{chains}}.
    '''
    # 20.1168 = 22 * yard2m(1)
    return Float(chains) * 20.1168


def circle4(earth, lat):
    '''Get the equatorial or a parallel I{circle of latitude}.

       @arg earth: The earth radius, ellipsoid or datum
                   (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                   L{Datum} or L{a_f2Tuple}).
       @arg lat: Geodetic latitude (C{degrees90}, C{str}).

       @return: A L{Circle4Tuple}C{(radius, height, lat, beta)}
                instance.

       @raise RangeError: Latitude B{C{lat}} outside valid range
                          and L{rangerrors} set to C{True}.

       @raise TypeError: Invalid B{C{earth}}.

       @raise ValueError: B{C{earth}} or B{C{lat}}.
    '''
    from pygeodesy.datums import _spherical_datum
    E = _spherical_datum(earth).ellipsoid
    return E.circle4(lat)


def degrees90(rad):
    '''Convert radians to degrees and wrap M{[-270..+90]}.

       @arg rad: Angle (C{radians}).

       @return: Angle in degrees, wrapped (C{degrees90}).
    '''
    return _wrap(degrees(rad), _90_0, _360_0)


def degrees180(rad):
    '''Convert radians to degrees and wrap M{[-180..+180]}.

       @arg rad: Angle (C{radians}).

       @return: Angle in degrees, wrapped (C{degrees180}).
    '''
    return _wrap(degrees(rad), _180_0, _360_0)


def degrees360(rad):
    '''Convert radians to degrees and wrap M{[0..+360)}.

       @arg rad: Angle (C{radians}).

       @return: Angle in degrees, wrapped (C{degrees360}).
    '''
    return _wrap(degrees(rad), _360_0, _360_0)


def degrees2m(deg, radius=R_M, lat=0):
    '''Convert an angle to a distance along the equator or
       along the parallel at an other (geodetic) latitude.

       @arg deg: The angle (C{degrees}).
       @kwarg radius: Mean earth radius, ellipsoid or datum
                      (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                      L{Datum} or L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

       @return: Distance (C{meter}, same units as B{C{radius}}
                or ellipsoidal and polar radii) or C{0} for
                near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range
                          and L{rangerrors} set to C{True}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise ValueError: Invalid B{C{deg}}, B{C{radius}} or
                          B{C{lat}}.

       @see: Function L{radians2m} and L{m2degrees}.
    '''
    return radians2m(Lam_(deg=deg, clip=0), radius=radius, lat=lat)


def fathom2m(fathoms):
    '''Convert I{UK} fathom to meter.

       @arg fathoms: Value in fathoms (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{fathoms}}.
    '''
    # 1.8288 == 2 * yard2m(1)
    return Float(fathoms) * 1.8288


def ft2m(feet, usurvey=False):
    '''Convert I{International} or I{US Survey} feet to meter.

       @arg feet: Value in feet (C{scalar}).
       @kwarg usurvery: Convert I{US Survey} feet (C{bool}),
                        I{International} feet otherwise.

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{feet}}.
    '''
    # US Survey 1200 / 3937 == 0.3048006096012192
    # Int'l 0.3048 == 254 * 12 / 10_000
    return Feet(feet) * (0.3048006096 if usurvey else 0.3048)


def furlong2m(furlongs):
    '''Convert a I{UK} furlong to meter.

       @arg furlongs: Value in furlongs (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{furlongs}}.
    '''
    # 201.168 = 220 * yard2m(1)
    return Float(furlongs) * 201.168


def isNumpy2(obj):
    '''Check for an B{C{Numpy2LatLon}} points wrapper.

       @arg obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is an B{C{Numpy2LatLon}}
                instance, C{False} otherwise.
    '''
    # isinstance(self, (Numpy2LatLon, ...))
    return getattr(obj, isNumpy2.__name__, False)


def isPoints2(obj):
    '''Check for an B{C{LatLon2psxy}} points wrapper.

       @arg obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is an B{C{LatLon2psxy}}
                instance, C{False} otherwise.
    '''
    # isinstance(self, (LatLon2psxy, ...))
    return getattr(obj, isPoints2.__name__, False)


def isTuple2(obj):
    '''Check for an B{C{Tuple2LatLon}} points wrapper.

       @arg obj: The object (any).

       @return: C{True} if B{C{obj}} is an B{C{Tuple2LatLon}}
                instance, C{False} otherwise.
    '''
    # isinstance(self, (Tuple2LatLon, ...))
    return getattr(obj, isTuple2.__name__, False)


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
            raise _ValueError(n=n)
    return p


def m2degrees(distance, radius=R_M, lat=0):
    '''Convert a distance to an angle along the equator or
       along the parallel at an other (geodetic) latitude.

       @arg distance: Distance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      an L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or
                      L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

       @return: Angle (C{degrees}) or C{INF} for near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range
                          and L{rangerrors} set to C{True}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise ValueError: Invalid B{C{distance}}, B{C{radius}}
                          or B{C{lat}}.

       @see: Function L{m2radians} and L{degrees2m}.
    '''
    return degrees(m2radians(distance, radius=radius, lat=lat))


def m2ft(meter, usurvey=False):
    '''Convert meter to I{International} or I{US Survey} feet (C{ft}).

       @arg meter: Value in meter (C{scalar}).
       @kwarg usurvery: Convert to I{US Survey} feet (C{bool}),
                        I{International} feet otherwise.

       @return: Value in C{feet} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    # US Survey == 3937 / 1200  == 3.2808333333333333
    # Int'l 10_000 / (254 * 12) == 3.2808398950131235
    return Meter(meter) * (3.280833333 if usurvey else 3.280839895)


def m2km(meter):
    '''Convert meter to kilo meter (km).

       @arg meter: Value in meter (C{scalar}).

       @return: Value in km (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Meter(meter) * 1.0e-3


def m2NM(meter):
    '''Convert meter to nautical miles (NM).

       @arg meter: Value in meter (C{scalar}).

       @return: Value in NM (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Meter(meter) * 5.39956804e-4  # == * _1_0 / 1852


def m2radians(distance, radius=R_M, lat=0):
    '''Convert a distance to an angle along the equator or
       along the parallel at an other (geodetic) latitude.

       @arg distance: Distance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      an L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or
                      L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

       @return: Angle (C{radians}) or C{INF} for near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range
                          and L{rangerrors} set to C{True}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise ValueError: Invalid B{C{distance}}, B{C{radius}}
                          or B{C{lat}}.

       @see: Function L{m2degrees} and L{radians2m}.
    '''
    m = circle4(radius, lat).radius
    return INF if m < EPS0 else (Float(distance) / m)


def m2SM(meter):
    '''Convert meter to statute miles (SM).

       @arg meter: Value in meter (C{scalar}).

       @return: Value in SM (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Meter(meter) * 6.21369949e-4  # == _1_0 / 1609.344


def m2yard(meter):
    '''Convert meter to I{UK} yards.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in yards (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    # 1.0936132983377078 == 10_000 / (254 * 12 * 3)
    return Meter(meter) * 1.09361329833771


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


def radians2m(rad, radius=R_M, lat=0):
    '''Convert an angle to a distance along the equator or
       along the parallel at an other (geodetic) latitude.

       @arg rad: The angle (C{radians}).
       @kwarg radius: Mean earth radius, ellipsoid or datum
                      (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                      L{Datum} or L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

       @return: Distance (C{meter}, same units as B{C{radius}}
                or ellipsoidal and polar radii) or C{0} for
                near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range
                          and L{rangerrors} set to C{True}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise ValueError: Invalid B{C{rad}}, B{C{radius}} or
                          B{C{lat}}.

       @see: Function L{degrees2m} and L{m2radians}.
    '''
    m = circle4(radius, lat).radius
    return _0_0 if m < EPS0 else (Lam(rad=rad, clip=0) * m)


def _sincos2(q, r):
    '''(INTERNAL) 2-tuple (C{sin(r), cos(r)}) in quadrant M{0 <= q <= 3}.
    '''
    if r < EPS:
        s, c = _0_0, _1_0
    elif r < PI_2:
        s, c = sin(r), cos(r)
    else:  # r == PI_2
        s, c = _1_0, _0_0
    t = s, c, -s, -c, s
#   q &= 3
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
        s, c = _sincos2(q & 3, r - q * PI_2)
        yield s
        yield c


def sincos2d(*deg):
    '''Return the C{sine} and C{cosine} of angle(s) in C{degrees}.

       @arg deg: One or more angles (C{degrees}).

       @return: The C{sin(deg)} and C{cos(deg)} for each angle.

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
        s, c = _sincos2(q & 3, radians(d - q * _90_0))
        yield s
        yield c


def splice(iterable, n=2, **fill):
    '''Split an iterable into C{n} slices.

       @arg iterable: Items to be spliced (C{list}, C{tuple}, ...).
       @kwarg n: Number of slices to generate (C{int}).
       @kwarg fill: Optional fill value for missing items.

       @return: A generator for each of B{C{n}} slices,
                M{iterable[i::n] for i=0..n}.

       @raise TypeError: Invalid B{C{n}}.

       @note: Each generated slice is a C{tuple} or a C{list},
              the latter only if the B{C{iterable}} is a C{list}.

       @example:

       >>> from pygeodesy import splice

       >>> a, b = splice(range(10))
       >>> a, b
       ((0, 2, 4, 6, 8), (1, 3, 5, 7, 9))

       >>> a, b, c = splice(range(10), n=3)
       >>> a, b, c
       ((0, 3, 6, 9), (1, 4, 7), (2, 5, 8))

       >>> a, b, c = splice(range(10), n=3, fill=-1)
       >>> a, b, c
       ((0, 3, 6, 9), (1, 4, 7, -1), (2, 5, 8, -1))

       >>> tuple(splice(list(range(9)), n=5))
       ([0, 5], [1, 6], [2, 7], [3, 8], [4])

       >>> splice(range(9), n=1)
       <generator object splice at 0x0...>
    '''
    if not isint(n):
        raise _TypeError(n=n)

    t = iterable
    if not isinstance(t, (list, tuple)):
        t = tuple(t)  # force tuple, also for PyPy3

    if n > 1:
        if fill:
            fill = _xkwds_get(fill, fill=MISSING)
            if fill is not MISSING:
                m = len(t) % n
                if m > 0:  # fill with same type
                    t += type(t)((fill,)) * (n - m)
        for i in range(n):
            # XXX t[i::n] chokes PyChecker
            yield t[slice(i, None, n)]
    else:
        yield t


def tan_2(rad):
    '''Compute the tangent of half angle.

       @arg rad: Angle (C{radians}).

       @return: M{tan(rad / 2)} (C{float}).
    '''
    return tan(rad * _0_5)


def tanPI_2_2(rad):
    '''Compute the tangent of half angle, 90 degrees rotated.

       @arg rad: Angle (C{radians}).

       @return: M{tan((rad + PI/2) / 2)} (C{float}).
    '''
    return tan((rad + PI_2) * _0_5)


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
    if wrap and abs(d) > _180_0:
        u = _wrap(d, _180_0, _360_0)
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
    if wrap > angle >= (wrap - modulo):
        angle = float(angle)
    else:
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
    return _wrap(deg, _90_0, _360_0)


def wrap180(deg):
    '''Wrap degrees to M{[-180..+180]}.

       @arg deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees180}).
    '''
    return _wrap(deg, _180_0, _360_0)


def wrap360(deg):
    '''Wrap degrees to M{[0..+360)}.

       @arg deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees360}).
    '''
    return _wrap(deg, _360_0, _360_0)


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


def yard2m(yards):
    '''Convert I{UK} yards to meter.

       @arg yards: Value in yards (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{yards}}.
    '''
    # 0.9144 == 254 * 12 * 3 / 10_000 == 3 * ft2m(1) Int'l
    return Float(yards) * 0.9144

# **) MIT License
#
# Copyright (C) 2016-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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
