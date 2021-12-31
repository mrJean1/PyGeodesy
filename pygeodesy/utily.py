
# -*- coding: utf-8 -*-

u'''Various utility functions.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**, see
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>} and
U{Vector-based geodesy<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import copysign0, isint, isnear0
from pygeodesy.interns import EPS, EPS0, INF, PI, PI2, PI_2, R_M, \
                             _edge_, _radians_, _semi_circular_, _SPACE_, \
                             _0_0, _0_5, _1_0, _90_0, _N_90_0, _180_0, \
                             _N_180_0, _360_0, _400_0
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.units import Degrees, Feet, Float, Lam, Lam_, Meter

from math import acos, asin, atan2, cos, degrees, radians, sin, tan  # pow

__all__ = _ALL_LAZY.utily
__version__ = '21.12.28'

# <https://Numbers.Computation.Free.Fr/Constants/Miscellaneous/digits.html>
_1__90 = _1_0 / _90_0  # 0.01111111111111111111111111111111111111111111111111
_2__PI = _1_0 /  PI_2  # 0.63661977236758134307553505349005744813783858296182
# sqrt(2) + 1 <https://WikiPedia.org/wiki/Square_root_of_2>
# _1sqrt2 = 2.41421356237309504880  # _16887_24209_69807_85696_71875_37694_80731_76679_73799


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

       @see: Function L{pygeodesy.atan2d}.
    '''
    return atan2d(y_x, _1_0)


def atan2b(y, x):
    '''Return C{atan2(B{y}, B{x})} in degrees M{[0..+360]}.

       @see: Function L{pygeodesy.atan2d}.
    '''
    d = atan2d(y, x)
    if d < 0:
        d += _360_0
    return d


def atan2d(y, x, reverse=False):
    '''Return C{atan2(B{y}, B{x})} in degrees M{[-180..+180]},
       optionally reversed (by 180 degrees for C{azi2}).

       @see: I{Karney}'s C++ function U{Math.atan2d
             <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Math.html>}.
    '''
    if abs(y) > abs(x) > 0:
        if y < 0:  # q = 3
            d = degrees(atan2(x, -y)) - _90_0
        else:  # q = 2
            d = _90_0 - degrees(atan2(x, y))
    elif x < 0:  # q = 1
        d = copysign0(_180_0, y) - degrees(atan2(y, -x))
    elif x > 0:  # q = 0
        d = degrees(atan2(y, x)) if y else _0_0
    else:  # x == 0
        d = _N_90_0 if y < 0 else (_90_0 if y > 0 else _0_0)
    if reverse:
        d += _180_0 if d < 0 else _N_180_0
    return d


def chain2m(chains):
    '''Convert I{UK} chains to meter.

       @arg chains: Value in chains (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{chains}}.
    '''
    return Float(chains=chains) * 20.1168  # == 22 * yard2m(1)


def circle4(earth, lat):
    '''Get the equatorial or a parallel I{circle of latitude}.

       @arg earth: The earth radius, ellipsoid or datum
                   (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                   L{Datum} or L{a_f2Tuple}).
       @arg lat: Geodetic latitude (C{degrees90}, C{str}).

       @return: A L{Circle4Tuple}C{(radius, height, lat, beta)}
                instance.

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{pygeodesy.rangerrors} set to C{True}.

       @raise TypeError: Invalid B{C{earth}}.

       @raise ValueError: B{C{earth}} or B{C{lat}}.
    '''
    E = _MODS.datums._spherical_datum(earth).ellipsoid
    return E.circle4(lat)


def cot(rad, **error_kwds):
    '''Return the C{cotangent} of an angle in C{radians}.

       @arg rad: Angle (C{radians}).
       @kwarg error_kwds: Error to raise (C{ValueError}).

       @return: C{cot(B{rad})}.

       @raise ValueError: L{pygeodesy.isnear0}C{(sin(B{rad})}.
    '''
    s, c = sincos2(rad)
    if isnear0(s):
        kwds = _MODS.errors._xkwds(error_kwds, cot=rad)
        raise  _MODS.errors._ValueError(**kwds)
    return c / s


def cot_(*rads, **error_kwds):
    '''Return the C{cotangent} of angle(s) in C{radiansresection}.

       @arg rads: One or more angles (C{radians}).
       @kwarg error_kwds: Error to raise (C{ValueError}).

       @return: Yield the C{cot(B{rad})} for each angle.

       @raise ValueError: See L{pygeodesy.cot}.
    '''
    for r in rads:
        yield cot(r, **error_kwds)


def cotd(deg, **error_kwds):
    '''Return the C{cotangent} of an angle in C{degrees}.

       @arg deg: Angle (C{degrees}).
       @kwarg error_kwds: Error to raise (C{ValueError}).

       @return: C{cot(B{deg})}.

       @raise ValueError: L{pygeodesy.isnear0}C{(sin(B{deg})}.
    '''
    s, c = sincos2d(deg)
    if isnear0(s):
        kwds = _MODS.errors._xkwds(error_kwds, cotd=deg)
        raise  _MODS.errors._ValueError(**kwds)
    return c / s


def cotd_(*degs, **error_kwds):
    '''Return the C{cotangent} of angle(s) in C{degrees}.

       @arg degs: One or more angles (C{degrees}).
       @kwarg error_kwds: Error to raise (C{ValueError}).

       @return: Yield the C{cot(B{deg})} for each angle.

       @raise ValueError: See L{pygeodesy.cotd}.
    '''
    for d in degs:
        yield cotd(d, **error_kwds)


def degrees90(rad):
    '''Convert radians to degrees and wrap M{[-270..+90]}.

       @arg rad: Angle (C{radians}).

       @return: Angle, wrapped (C{degrees90}).
    '''
    return _wrap(degrees(rad), _90_0, _360_0)


def degrees180(rad):
    '''Convert radians to degrees and wrap M{[-180..+180]}.

       @arg rad: Angle (C{radians}).

       @return: Angle, wrapped (C{degrees180}).
    '''
    return _wrap(degrees(rad), _180_0, _360_0)


def degrees360(rad):
    '''Convert radians to degrees and wrap M{[0..+360)}.

       @arg rad: Angle (C{radians}).

       @return: Angle, wrapped (C{degrees360}).
    '''
    return _wrap(degrees(rad), _360_0, _360_0)


def degrees2grades(deg):
    '''Convert degrees to I{grades} (aka I{gons} or I{gradians}).

       @arg deg: Angle (C{degrees}).

       @return: Angle (C{grades}).
    '''
    return Degrees(deg) * _400_0 / _360_0


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

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{pygeodesy.rangerrors} set to C{True}.

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

       @see: Function L{toise2m}.
    '''
    # 1.8288 == 2 * yard2m(1)
    return Float(fathoms=fathoms) * 1.8288


def ft2m(feet, usurvey=False):
    '''Convert I{International} or I{US Survey} feet to meter.

       @arg feet: Value in feet (C{scalar}).
       @kwarg usurvey: Convert I{US Survey} feet (C{bool}),
                       I{International} feet otherwise.

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{feet}}.
    '''
    # US Survey 1200 / 3937 == 0.3048006096012192
    # Int'l 0.3048 == 254 * 12 / 10_000
    return Feet(feet) * (0.3048006096 if usurvey else 0.3048)


def furlong2m(furlongs):
    '''Convert a furlong to meter.

       @arg furlongs: Value in furlongs (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{furlongs}}.
    '''
    # 201.168 = 220 * yard2m(1)
    return Float(furlongs=furlongs) * 201.168


def grades(rad):
    '''Convert radians to I{grades} (aka I{gons} or I{gradians}).

       @arg rad: Angle (C{radians}).

       @return: Angle (C{grades}).
    '''
    return Float(rad=rad) * _400_0 / PI2


def grades400(rad):
    '''Convert radians to I{grades} (aka I{gons} or I{gradians}) and wrap M{[0..+400)}.

       @arg rad: Angle (C{radians}).

       @return: Angle, wrapped (C{grades}).
    '''
    return _wrap(grades(rad), _400_0, _400_0)


def grades2degrees(gon):
    '''Convert I{grades} (aka I{gons} or I{gradians}) to C{degrees}.

       @arg gon: Angle (C{grades}).

       @return: Angle (C{degrees}).
    '''
    return Float(gon=gon) * _360_0 / _400_0


def grades2radians(gon):
    '''Convert I{grades} (aka I{gons} or I{gradians}) to C{radians}.

       @arg gon: Angle (C{grades}).

       @return: Angle (C{radians}).
    '''
    return Float(gon=gon) * PI2 / _400_0


def m2chain(meter):
    '''Convert meter to I{UK} chains.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{chains} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Meter(meter) * 0.049709695378986715  # == 1 / 20.1168


def m2degrees(distance, radius=R_M, lat=0):
    '''Convert a distance to an angle along the equator or
       along the parallel at an other (geodetic) latitude.

       @arg distance: Distance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      an L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or
                      L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

       @return: Angle (C{degrees}) or C{INF} for near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{pygeodesy.rangerrors} set to C{True}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise ValueError: Invalid B{C{distance}}, B{C{radius}}
                          or B{C{lat}}.

       @see: Function L{m2radians} and L{degrees2m}.
    '''
    return degrees(m2radians(distance, radius=radius, lat=lat))


def m2fathom(meter):
    '''Convert meter to I{UK} fathoms.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{fathoms} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.

       @see: Function L{m2toise}.
    '''
    return Meter(meter) * 0.546806649  # == 1 / 1.8288


def m2ft(meter, usurvey=False):
    '''Convert meter to I{International} or I{US Survey} feet (C{ft}).

       @arg meter: Value in meter (C{scalar}).
       @kwarg usurvey: Convert to I{US Survey} feet (C{bool}),
                       I{International} feet otherwise.

       @return: Value in C{feet} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    # US Survey == 3937 / 1200  == 3.2808333333333333
    # Int'l 10_000 / (254 * 12) == 3.2808398950131235
    return Meter(meter) * (3.280833333 if usurvey else 3.280839895)


def m2furlong(meter):
    '''Convert meter to furlongs.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{furlongs} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Meter(meter) * 0.00497096954  # == 1 / 201.168


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

       @return: Value in C{NM} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Meter(meter) * 5.39956804e-4  # == 1 / 1852


def m2radians(distance, radius=R_M, lat=0):
    '''Convert a distance to an angle along the equator or
       along the parallel at an other (geodetic) latitude.

       @arg distance: Distance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      an L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or
                      L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

       @return: Angle (C{radians}) or C{INF} for near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{pygeodesy.rangerrors} set to C{True}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise ValueError: Invalid B{C{distance}}, B{C{radius}}
                          or B{C{lat}}.

       @see: Function L{m2degrees} and L{radians2m}.
    '''
    m = circle4(radius, lat).radius
    return INF if m < EPS0 else (Float(distance=distance) / m)


def m2SM(meter):
    '''Convert meter to statute miles (SM).

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{SM} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Meter(meter) * 6.21369949e-4  # == _1_0 / 1609.344


def m2toise(meter):
    '''Convert meter to French U{toises<https://WikiPedia.org/wiki/Toise>}.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{toises} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.

       @see: Function L{m2fathom}.
    '''
    return Meter(meter) * 0.513083632632119  # == 1 / 1.949


def m2yard(meter):
    '''Convert meter to I{UK} yards.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{yards} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Meter(meter) * 1.0936132983377078  # == 10_000 / (254 * 12 * 3)


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

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{pygeodesy.rangerrors} set to C{True}.

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
    if r < EPS:  # XXX EPS0
        s, c = _0_0, _1_0
    elif r < PI_2:
        s, c = sin(r), cos(r)
    else:  # r == PI_2
        s, c = _1_0, _0_0
    t = s, c, -s, -c, s
#   q &= 3
    return t[q], t[q + 1]


def sincos2(rad):
    '''Return the C{sine} and C{cosine} of an angle in C{radians}.

       @arg rad: Angle (C{radians}).

       @return: 2-Tuple (C{sin(B{rad})}, C{cos(B{rad})}).

       @see: U{GeographicLib<https://GeographicLib.SourceForge.io/html/
             classGeographicLib_1_1Math.html#sincosd>} function U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             python/geographiclib/geomath.py#l155>} and C++ U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             include/GeographicLib/Math.hpp#l558>}.
    '''
    q = int(rad * _2__PI)  # int(math.floor)
    if rad < 0:
        q -= 1
    return _sincos2(q & 3, rad - q * PI_2)


def sincos2_(*rads):
    '''Return the C{sine} and C{cosine} of angle(s) in {Cradians}.

       @arg rads: One or more angles (C{radians}).

       @return: Yield the C{sin(B{rad})} and C{cos(B{rad})} for each angle.

       @see: function L{sincos2}.
    '''
    for r in rads:
        s, c = sincos2(r)
        yield s
        yield c


def sincos2d(deg):
    '''Return the C{sine} and C{cosine} of an angle in C{degrees}.

       @arg deg: Angle (C{degrees}).

       @return: 2-Tuple (C{sin(B{deg})}, C{cos(B{deg})}).

       @see: U{GeographicLib<https://GeographicLib.SourceForge.io/html/
             classGeographicLib_1_1Math.html#sincosd>} function U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             python/geographiclib/geomath.py#l155>} and C++ U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             include/GeographicLib/Math.hpp#l558>}.
    '''
    q = int(deg * _1__90)  # int(math.floor)
    if deg < 0:
        q -= 1
    return _sincos2(q & 3, radians(deg - q * _90_0))


def sincos2d_(*degs):
    '''Return the C{sine} and C{cosine} of angle(s) in C{degrees}.

       @arg degs: One or more angles (C{degrees}).

       @return: Yield the C{sin(B{deg})} and C{cos(B{deg})} for each angle.

       @see: Function L{sincos2d}.
    '''
    for d in degs:
        s, c = sincos2d(d)
        yield s
        yield c


def tan_2(rad, **semi):  # edge=1
    '''Compute the tangent of half angle.

       @arg rad: Angle (C{radians}).
       @kwarg semi: Angle or edge name and index
                    for semi-circular error.

       @return: M{tan(rad / 2)} (C{float}).

       @raise ValueError: If B{C{rad}} is semi-circular
                          and B{C{semi}} is given.
    '''
    # .formy.excessKarney_, .sphericalTrigonometry.areaOf
    if semi and isnear0(abs(rad) - PI):
        for n, v in semi.items():
            break
        n = _SPACE_(n, _radians_) if not isint(v) else \
            _SPACE_(_MODS.streprs.Fmt.SQUARE(**semi), _edge_)
        raise _MODS.errors._ValueError(n, rad, txt=_semi_circular_)

    return tan(rad * _0_5)


def tand(deg, **error_kwds):
    '''Return the C{tangent} of an angle in C{degrees}.

       @arg deg: Angle (C{degrees}).
       @kwarg error_kwds: Error to raise (C{ValueError}).

       @return: C{tan(B{deg})}.

       @raise ValueError: If L{pygeodesy.isnear0}C{(cos(B{deg})}.
    '''
    s, c = sincos2d(deg)
    if isnear0(c):
        raise _MODS.errors._ValueError(tand.__name__, deg, **error_kwds)
    return s / c


def tand_(*degs, **error_kwds):
    '''Return the C{tangent} of angle(s) in C{degrees}.

       @arg degs: One or more angles (C{degrees}).
       @kwarg error_kwds: Error to raise (C{ValueError}).

       @return: Yield the C{tan(B{deg})} for each angle.

       @raise ValueError: See L{pygeodesy.tand}.
    '''
    for d in degs:
        yield tand(d, **error_kwds)


def tanPI_2_2(rad):
    '''Compute the tangent of half angle, 90 degrees rotated.

       @arg rad: Angle (C{radians}).

       @return: M{tan((rad + PI/2) / 2)} (C{float}).
    '''
    return tan((rad + PI_2) * _0_5)


def toise2m(toises):
    '''Convert French U{toises<https://WikiPedia.org/wiki/Toise>} to meter.

       @arg toises: Value in toises (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{toises}}.

       @see: Function L{fathom2m}.
    '''
    return Float(toises=toises) * 1.949


def unroll180(lon1, lon2, wrap=True):
    '''Unroll longitudinal delta and wrap longitude in degrees.

       @arg lon1: Start longitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg wrap: Wrap and unroll to the M{(-180..+180]} range (C{bool}).

       @return: 2-Tuple C{(B{lon2}-B{lon1}, B{lon2})} unrolled (C{degrees},
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


def _unrollon(p1, p2):  # unroll180 == .karney._unroll2
    '''(INTERNAL) Wrap, unroll and replace longitude if different.
    '''
    _, lon = unroll180(p1.lon, p2.lon, wrap=True)
    if abs(lon - p2.lon) > EPS:
        p2 = p2.dup(lon=wrap180(lon))
    return p2


def unrollPI(rad1, rad2, wrap=True):
    '''Unroll longitudinal delta and wrap longitude in radians.

       @arg rad1: Start longitude (C{radians}).
       @arg rad2: End longitude (C{radians}).
       @kwarg wrap: Wrap and unroll to the M{(-PI..+PI]} range (C{bool}).

       @return: 2-Tuple C{(B{rad2}-B{rad1}, B{rad2})} unrolled (C{radians},
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

       @arg angle: Angle (C{degrees}, C{radians} or C{grades}).
       @arg wrap: Range (C{degrees}, C{radians} or C{grades}).
       @arg modulo: Upper limit (360 C{degrees}, PI2 C{radians} or 400 C{grades}).

       @return: The B{C{angle}}, wrapped (C{degrees}, C{radians} or C{grades}).
    '''
    angle = float(angle)
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
    return Float(yards=yards) * 0.9144

# **) MIT License
#
# Copyright (C) 2016-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
