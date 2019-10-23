
# -*- coding: utf-8 -*-

u'''Geometric and other utility functions and constants.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**, see
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>} and
U{Vector-based geodesy<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.

@newfield example: Example, Examples
'''
# make sure int division yields float quotient
from __future__ import division
division = 1 / 2  # double check int division, see .datum.py, .fmath.py
if not division:
    raise ImportError('%s 1/2 == %d' % ('division', division))
del division

from pygeodesy.fmath import _Ints, _Seqs, EPS, map2
from pygeodesy.lazily import _ALL_LAZY

from inspect import isclass
from math import cos, degrees, pi as PI, radians, sin, tan  # pow

_MISSING  = object()  # singleton, imported by .utily

# all public contants, classes and functions
__all__ = _ALL_LAZY.utily
__version__ = '19.10.19'

try:
    _Strs = basestring, str  # PYCHOK .datum.py, .geohash.py
except NameError:
    _Strs = str,

OK = 'OK'  # OK for test like I{if ... is OK: ...}

PI2  = PI * 2  #: Two PI, M{PI * 2} aka Tau (C{float})  # PYCHOK expected
PI_2 = PI / 2  #: Half PI, M{PI / 2} (C{float})
PI_4 = PI / 4  #: Quarter PI, M{PI / 4} (C{float})

# R_M moved here to avoid circular imports
R_M = 6371008.771415  #: Mean, spherical earth radius (C{meter}).

_1_90 = 1 / 90.0  # 0.011111111111111111111111111111111111111111111111
# <https://Numbers.Computation.Free.FR/Constants/Miscellaneous/digits.html>
_2_PI = 2 / PI  # 0.63661977236758134307553505349005744813783858296182

_iterNumpy2len = 1  # adjustable for testing purposes
_limiterrors   = True


def _TypeError(*Types, **pairs):
    '''(INTERNAL) Check C{Types} of all C{name=value} pairs.
    '''
    for n, v in pairs.items():
        if not isinstance(v, Types):
            t = ' or '.join(t.__name__ for t in Types)
            # first letter of Type name I{pronounced} as vowel
            a = 'an' if t[:1].lower() in 'aeinoux' else 'a'
            raise TypeError('%s not %s %s: %r' % (n, a, t, v))


class LimitError(ValueError):
    '''Error raised for lat- or longitudinal deltas exceeding
       the B{C{limit}} in functions L{equirectangular} and
       L{equirectangular_}.
    '''
    pass


def anStr(name, OKd='._-', sub='_'):
    '''Make a valid name of alphanumeric and OKd characters.

       @param name: The original name (C{str}).
       @keyword OKd: Other acceptable characters (C{str}).
       @keyword sub: Substitute for invalid charactes (C{str}).

       @return: The modified name (C{str}).

       @note: Leading and trailing whitespace characters are removed
              and intermediate whitespace characters are coalesced
              and substituted.
    '''
    s = n = str(name).strip()
    for c in n:
        if not (c.isalnum() or c in OKd or c in sub):
            s = s.replace(c, ' ')
    return sub.join(s.strip().split())


def clipStr(bstr, limit=50, white=''):
    '''Clip a string to the given length limit.

       @param bstr: String (C{bytes} or C{str}).
       @keyword limit: Length limit (C{int}).
       @keyword white: Whitespace replacement (C{str}).

       @return: Clipped C{bytes} or C{str}.
    '''
    t = type(bstr)
    if bstr and limit > 8:
        n = len(bstr)
        if n > limit:
            h = limit // 2
            bstr = bstr[:h] + t('....') + bstr[-h:]
    if white:  # replace whitespace
        bstr = t(white).join(bstr.split())
    return bstr


def degrees90(rad):
    '''Convert radians to degrees and wrap M{[-270..+90]}.

       @param rad: Angle (C{radians}).

       @return: Angle in degrees, wrapped (C{degrees90}).
    '''
    return _wrap(degrees(rad), 90, 360)


def degrees180(rad):
    '''Convert radians to degrees and wrap M{[-180..+180]}.

       @param rad: Angle (C{radians}).

       @return: Angle in degrees, wrapped (C{degrees180}).
    '''
    return _wrap(degrees(rad), 180, 360)


def degrees360(rad):
    '''Convert radians to degrees and wrap M{[0..+360)}.

       @param rad: Angle (C{radians}).

       @return: Angle in degrees, wrapped (C{degrees360}).
    '''
    return _wrap(degrees(rad), 360, 360)


def degrees2m(deg, radius=R_M, lat=0):
    '''Convert angle to distance along the equator or along a
       parallel at an other latitude.

       @param deg: Angle (C{degrees}).
       @keyword radius: Mean earth radius (C{meter}).
       @keyword lat: Parallel latitude (C{degrees90}).

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @raise RangeError: Latitude B{C{lat}} outside valid range
                          and L{rangerrors} set to C{True}.
    '''
    m = radians(deg) * radius
    if lat:
        from pygeodesy.dms import clipDMS
        m *= cos(radians(clipDMS(lat, 90)))
    return m


def enStr2(easting, northing, prec, *extras):
    '''Return easting, northing string representations.

       @param easting: Easting from false easting (C{meter}).
       @param northing: Northing from from false northing (C{meter}).
       @param prec: Precision in number of digits (C{int}).
       @param extras: Optional leading items (strings).

       @return: B{C{extras}} + 2-Tuple C{(eastingStr, northingStr)}.

       @raise ValueError: Invalid B{C{prec}}.
    '''
    w = prec // 2
    try:
        p10 = (1e-4, 1e-3, 1e-2, 1e-1, 1)[w - 1]  # 10**(5 - w)
    except IndexError:
        raise ValueError('%s invalid: %r' % ('prec', prec))
    return extras + ('%0*d' % (w, int(easting * p10)),
                     '%0*d' % (w, int(northing * p10)))


def false2f(value, name='value', false=True, Error=ValueError):
    '''Convert a false east-/northing to non-negative float.

       @param value: Value to convert (C{scalar}).
       @keyword name: Optional name of the value (C{str}).
       @keyword false: Optionally, value includes false origin (C{bool}).
       @keyword Error: Exception to raise (C{ValueError}).

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

       @param feet: Value in feet (C{scalar}).
       @keyword usurvery: Convert I{US Survey} feet (C{bool}),
                          I{International} feet otherwise.

       @return: Value in C{meter} (C{float}).
    '''
    # US Survey 1200./3937. == 0.3048006096012192
    return feet * (0.3048006096012192 if usurvey else 0.3048)


def halfs2(str2):
    '''Split a string in 2 halfs.

       @param str2: String to split (C{str}).

       @return: 2-Tuple (1st, 2nd) half (C{str}).

       @raise ValueError: Zero or odd C{len}(B{str2}).
    '''
    h, r = divmod(len(str2), 2)
    if r or not h:
        raise ValueError('%s invalid: %r' % ('str2', str2))
    return str2[:h], str2[h:]


def isNumpy2(obj):
    '''Check for an B{C{Numpy2LatLon}} points wrapper.

       @param obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is an B{C{Numpy2LatLon}}
                instance, C{False} otherwise.
    '''
    # isinstance(self, (Numpy2LatLon, ...))
    return getattr(obj, 'isNumpy2', False)


def isPoints2(obj):
    '''Check for an B{C{LatLon2psxy}} points wrapper.

       @param obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is an B{C{LatLon2psxy}}
                instance, C{False} otherwise.
    '''
    # isinstance(self, (LatLon2psxy, ...))
    return getattr(obj, 'isPoints2', False)


def issequence(obj, *excluded):
    '''Check for sequence types.

       @param obj: The object (any C{type}).
       @param excluded: Optional, exclusions (C{type}).

       @note: Excluding C{tuple} implies excluding C{namedtuple}.

       @return: C{True} if B{C{obj}} is a sequence, C{False} otherwise.
    '''
    if excluded:
        return isinstance(obj, _Seqs) and not \
               isinstance(obj, excluded)
    else:
        return isinstance(obj, _Seqs)


def issubclassof(sub, sup):
    '''Check whether a class is a subclass of a super class.

       @param sub: The sub class (C{class}).
       @param sup: The super class (C{class}).

       @return: C{True} if B{C{sub}} is a subclass of B{C{sup}}.
    '''
    return isclass(sub) and isclass(sup) and issubclass(sub, sup)


def isTuple2(obj):
    '''Check for an B{C{Tuple2LatLon}} points wrapper.

       @param obj: The object (any).

       @return: C{True} if B{C{obj}} is an B{C{Tuple2LatLon}}
                instance, C{False} otherwise.
    '''
    # isinstance(self, (Tuple2LatLon, ...))
    return getattr(obj, 'isTuple2', False)


def iterNumpy2(obj):
    '''Iterate over Numpy2 wrappers or other sequences exceeding
       the threshold.

       @param obj: Points array, list, sequence, set, etc. (any).

       @return: C{True} do, C{False} don't iterate.
    '''
    try:
        return isNumpy2(obj) or len(obj) > _iterNumpy2len
    except TypeError:
        return False


def iterNumpy2over(n=None):
    '''Get or set the L{iterNumpy2} threshold.

       @keyword n: Optional, new threshold (C{int}).

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


def limiterrors(raiser=None):
    '''Get/set the raising of limit errors.

       @keyword raiser: Choose C{True} to raise or C{False} to
                        ignore L{LimitError} exceptions.  Use
                        C{None} to leave the setting unchanged.

       @return: Previous setting (C{bool}).
    '''
    global _limiterrors
    t = _limiterrors
    if raiser in (True, False):
        _limiterrors = raiser
    return t


def m2degrees(meter, radius=R_M):
    '''Convert distance to angle along equator.

       @param meter: Distance (C{meter}, same units as B{C{radius}}).
       @keyword radius: Mean earth radius (C{meter}).

       @return: Angle (C{degrees}).

       @raise ValueError: Invalid B{C{radius}}.
    '''
    if radius < EPS:
        raise ValueError('%s invalid: %r' % ('radius', radius))
    return degrees(meter / radius)


def m2ft(meter, usurvey=False):
    '''Convert meter to I{International} or I{US Survey} feet (C{ft}).

       @param meter: Value in meter (C{scalar}).
       @keyword usurvery: Convert to I{US Survey} feet (C{bool}),
                          I{International} feet otherwise.

       @return: Value in C{feet} (C{float}).
    '''
    # US Survey == 3937./1200. = 3.2808333333333333
    return meter * (3.2808333333333333 if usurvey else 3.2808399)


def m2km(meter):
    '''Convert meter to kilo meter (km).

       @param meter: Value in meter (C{scalar}).

       @return: Value in km (C{float}).
    '''
    return meter * 1.0e-3


def m2NM(meter):
    '''Convert meter to nautical miles (NM).

       @param meter: Value in meter (C{scalar}).

       @return: Value in NM (C{float}).
    '''
    return meter * 5.39956804e-4  # == * 1.0 / 1852.0


def m2SM(meter):
    '''Convert meter to statute miles (SM).

       @param meter: Value in meter (C{scalar}).

       @return: Value in SM (C{float}).
    '''
    return meter * 6.21369949e-4  # XXX 6.213712e-4 == 1.0 / 1609.344


def property_RO(method):
    '''Decorator for C{Read_Only} property.

       @param method: The callable to be decorated as C{property.getter}.

       @note: Like standard Python C{property} without a C{property.setter}
              with a more descriptive error message when set.
    '''
    def Read_Only(self, ignored):
        '''Throws an C{AttributeError}, always.
        '''
        raise AttributeError('Read_Only property: %r.%s = %r' %
                             (self, method.__name__, ignored))

    return property(method, Read_Only, None, method.__doc__ or 'N/A')


def radiansPI(deg):
    '''Convert and wrap degrees to radians M{[-PI..+PI]}.

       @param deg: Angle (C{degrees}).

       @return: Radians, wrapped (C{radiansPI})
    '''
    return _wrap(radians(deg), PI, PI2)


def radiansPI2(deg):
    '''Convert and wrap degrees to radians M{[0..+2PI)}.

       @param deg: Angle (C{degrees}).

       @return: Radians, wrapped (C{radiansPI2})
    '''
    return _wrap(radians(deg), PI2, PI2)


def radiansPI_2(deg):
    '''Convert and wrap degrees to radians M{[-3PI/2..+PI/2]}.

       @param deg: Angle (C{degrees}).

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

       @param rad: One or more angles (C{radians}).

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
    '''Return the C{sine} and C{cosine} of an angle.

       @param deg: One or more angles (C{degrees}).

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


def splice(iterable, n=2, fill=_MISSING):
    '''Split an iterable into C{n} slices.

       @param iterable: Items to be spliced (C{list}, C{tuple}, ...).
       @keyword n: Number of slices to generate (C{int}).
       @keyword fill: Fill value for missing items.

       @return: Generator of B{C{n}} slices M{iterable[i::n] for i=0..n}.

       @note: Each generated slice is a C{tuple} or a C{list},
              the latter only if the B{C{iterable}} is a C{list}.

       @raise ValueError: Non-C{int} or non-positive B{C{n}}.

       @example:

       >>> from pygeodesy import splice

       >>> a, b = splice(range(10))
       >>> a, b
       ((0, 2, 4, 6, 8), (1, 3, 5, 7, 9))

       >>> a, b, c = splice(range(10), n=3)
       >>> a, b, c
       ((0, 3, 6, 9), (1, 4, 7], [2, 5, 8))

       >>> a, b, c = splice(range(10), n=3, fill=-1)
       >>> a, b, c
       ((0, 3, 6, 9), (1, 4, 7, -1), (2, 5, 8, -1))

       >>> list(splice(range(12), n=5))
       [(0, 5, 10), (1, 6, 11), (2, 7), (3, 8), (4, 9)]

       >>> splice(range(9), n=1)
       <generator object splice at 0x0...>
    '''
    if not (isinstance(n, _Ints) and n > 0):
        raise ValueError('%s %s=%s' % ('splice', 'n', n))

    t = iterable
    if not isinstance(t, (list, tuple)):
        t = tuple(t)  # force tuple, also for PyPy3
    if n > 1:
        if fill is not _MISSING:
            m = len(t) % n
            if m > 0:  # fill with same type
                t += type(t)((fill,)) * (n - m)
        for i in range(n):
            yield t[i::n]  # [i:None:n] pychok -Tb ...
    else:
        yield t


def tan_2(rad):
    '''Compute the tangent of half angle.

       @param rad: Angle (C{radians}).

       @return: M{tan(rad / 2)} (C{float}).
    '''
    return tan(rad * 0.5)


def tanPI_2_2(rad):
    '''Compute the tangent of half angle, 90 degrees rotated.

       @param rad: Angle (C{radians}).

       @return: M{tan((rad + PI/2) / 2)} (C{float}).
    '''
    return tan((rad + PI_2) * 0.5)


def unroll180(lon1, lon2, wrap=True):
    '''Unroll longitudinal delta and wrap longitude in degrees.

       @param lon1: Start longitude (C{degrees}).
       @param lon2: End longitude (C{degrees}).
       @keyword wrap: Wrap and unroll to the M{(-180..+180]} range (C{bool}).

       @return: 2-Tuple (delta B{C{lon2}}-B{lon1}, B{C{lon2}}) unrolled
                (C{degrees}, C{degrees}).

       @see: Capability B{C{LONG_UNROLL}} in U{GeographicLib
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

       @param rad1: Start longitude (C{radians}).
       @param rad2: End longitude (C{radians}).
       @keyword wrap: Wrap and unroll to the M{(-PI..+PI]} range (C{bool}).

       @return: 2-Tuple (delta B{C{rad2}}-B{rad1}, B{C{rad2}}) unrolled
                (C{radians}, C{radians}).

       @see: Capability B{C{LONG_UNROLL}} in U{GeographicLib
       <https://GeographicLib.SourceForge.io/html/python/interface.html#outmask>}.
    '''
    r = rad2 - rad1
    if wrap and abs(r) > PI:
        u = _wrap(r, PI, PI2)
        if u != r:
            return u, rad1 + u
    return r, rad2


def unStr(name, *args, **kwds):
    '''Return the string representation of an invokation.

       @param name: Function, method or class name (C{str}).
       @param args: Optional positional arguments.
       @keyword kwds: Optional keyword arguments.

       @return: Representation (C{str}).
    '''
    t = tuple('%s=%s' % t for t in sorted(kwds.items()))
    if args:
        t = map2(str, args) + t
    return '%s(%s)' % (name, ', '.join(t))


def _wrap(angle, wrap, modulo):
    '''(INTERNAL) Angle wrapper M{((wrap-modulo)..+wrap]}.

       @param angle: Angle (C{degrees} or C{radians}).
       @param wrap: Range (C{degrees} or C{radians}).
       @param modulo: Upper limit (360 C{degrees} or PI2 C{radians}).

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

       @param deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees90}).
    '''
    return _wrap(deg, 90, 360)


def wrap180(deg):
    '''Wrap degrees to M{[-180..+180]}.

       @param deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees180}).
    '''
    return _wrap(deg, 180, 360)


def wrap360(deg):
    '''Wrap degrees to M{[0..+360)}.

       @param deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees360}).
    '''
    return _wrap(deg, 360, 360)


def wrapPI(rad):
    '''Wrap radians to M{[-PI..+PI]}.

       @param rad: Angle (C{radians}).

       @return: Radians, wrapped (C{radiansPI}).
    '''
    return _wrap(rad, PI, PI2)


def wrapPI2(rad):
    '''Wrap radians to M{[0..+2PI)}.

       @param rad: Angle (C{radians}).

       @return: Radians, wrapped (C{radiansPI2}).
    '''
    return _wrap(rad, PI2, PI2)


def wrapPI_2(rad):
    '''Wrap radians to M{[-3PI/2..+PI/2]}.

       @param rad: Angle (C{radians}).

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
