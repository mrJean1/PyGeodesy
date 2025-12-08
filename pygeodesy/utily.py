
# -*- coding: utf-8 -*-

u'''Various utility functions.

After I{Karney}'s C++ U{Math<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Math.html>}
class and I{Veness}' U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>}
and U{Vector-based geodesy<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}
and published under the same MIT Licence**.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

from pygeodesy.basics import _copysign, isinstanceof, isint, isstr
from pygeodesy.constants import EPS, EPS0, NAN, PI, PI2, PI_2, PI_4, PI_6, R_M, \
                               _M_KM, _M_NM, _M_SM, _0_0, _0_5, _1_0, _N_1_0, \
                               _10_0, _90_0, _180_0, _360_0, _copysign_0_0, \
                               _copysignINF, _float, _isfinite, isnan, isnear0, \
                               _over_1, _umod_360, _umod_PI2, OVERFLOW
from pygeodesy.errors import _ValueError, _xkwds, _xkwds_get
from pygeodesy.internals import _Enum, _passargs, typename
from pygeodesy.interns import _edge_, _radians_, _semi_circular_, _SPACE_
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.units import Degrees, Degrees_, Feet, Float, Lam, Lamd, \
                            Meter, Meter2, Radians  # Radians_

from math import acos, asin, asinh, atan2 as _atan2, cos, degrees, fabs, \
                 radians, sin, sinh, tan as _tan  # pow

__all__ = _ALL_LAZY.utily
__version__ = '25.11.09'

# sqrt(3) <https://WikiPedia.org/wiki/Square_root_of_3>
_COS_30,  _SIN_30 = 0.86602540378443864676, _0_5  # sqrt(3) / 2
_COS_45 = _SIN_45 = 0.70710678118654752440  # sqrt(2) / 2

_G = _Enum(  # grades per ...
    DEG     = _float(  400.0 / _360_0),  # degree
    RAD     = _float(  400.0 /  PI2))    # radian

_M = _Enum(  # meter per ...
    ACRE    = _float( 4046.8564224),     # acre, chain2m(1) * furlong2m(1), squared
    CHAIN   = _float(   20.1168),        # yard2m(1) * 22
    FATHOM  = _float(    1.8288),        # yard2m(1) * 2 or _M.NM * 1e-3
    FOOT    = _float(    0.3048),        # Int'l foot, 1 / 3.280_839_895_0131 = 10_000 / (254 * 12)
    FOOT_GE = _float(    0.31608),       # German Fuss, 1 / 3.163_756_011_1364
    FOOT_FR = _float(    0.3248406),     # French Pied-du-Roi or pied, 1 / 3.078_432_929_8739
    FOOT_US = _float(    0.3048006096012192),  # US Survey foot, 1_200 / 3_937
    FURLONG = _float(  201.168),         # furlong, 220 * yard2m(1) = 10 * m2chain(1)
    HA      = _float(10000.0),           # hectare, 100 * 100, squared
    KM      = _M_KM,                     # kilo meter
    NM      = _M_NM,                     # nautical mile
    SM      = _M_SM,                     # statute mile
    TOISE   = _float(    1.9490436),     # French toise, 6 pieds = 6 / 3.078_432_929_8739
    YARD_UK = _float(    0.9144))        # yard, 254 * 12 * 3 / 10_000 = 3 * _M.FOOT


def _abs1nan(x):
    '''(INTERNAL) Bracket C{x}.
    '''
    return _N_1_0 < x < _1_0 or isnan(x)


def acos1(x):
    '''Return C{math.acos(max(-1, min(1, B{x})))}.
    '''
    return acos(x) if _abs1nan(x) else (PI if x < 0 else _0_0)


def acre2ha(acres):
    '''Convert acres to hectare.

       @arg acres: Value in acres (C{scalar}).

       @return: Value in C{hectare} (C{float}).

       @raise ValueError: Invalid B{C{acres}}.
    '''
    return m2ha(acre2m2(acres))


def acre2m2(acres):
    '''Convert acres to I{square} meter.

       @arg acres: Value in acres (C{scalar}).

       @return: Value in C{meter^2} (C{float}).

       @raise ValueError: Invalid B{C{acres}}.
    '''
    return Meter2(Float(acres=acres) * _M.ACRE)


def agdf(phi):
    '''Inverse U{Gudermannian function
       <https://WikiPedia.org/wiki/Gudermannian_function>}.

       @arg phi: Angle (C{radians}).

       @return: Gudermannian (psi, C{float}).

       @see: Function L{gdf}.
    '''
    return asinh(tan(phi))


def asin1(x):
    '''Return C{math.asin(max(-1, min(1, B{x})))}.
    '''
    return asin(x) if _abs1nan(x) else _copysign(PI_2, x)


def atan1(y, x=_1_0):
    '''Return C{atan(B{y} / B{x})} angle in C{radians} M{[-PI/2..+PI/2]}
       using C{atan2} for consistency and to avoid C{ZeroDivisionError}.
    '''
    return _atan1u(y, x, atan2)


def atan1d(y, x=_1_0):
    '''Return C{atan(B{y} / B{x})} angle in C{degrees} M{[-90..+90]}
       using C{atan2d} for consistency and to avoid C{ZeroDivisionError}.

       @see: Function L{pygeodesy.atan2d}.
    '''
    return _atan1u(y, x, atan2d)


def _atan1u(y, x, _2u):
    '''(INTERNAL) Helper for functions C{atan1} and C{atan1d}.
    '''
    if x < 0:
        x = -x
        y = -y
    return _2u(y, x or _0_0)


atan2 = _atan2
'''Return C{atan2(B{y}, B{x})} in radians M{[-PI..+PI]}.

   @see: I{Karney}'s C++ function U{Math.atan2d
         <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Math.html>}.
'''


def atan2b(y, x):
    '''Return C{atan2(B{y}, B{x})} in degrees M{[0..+360], counter-clockwise}.

       @see: Function L{pygeodesy.atan2d}.
    '''
    b = atan2d(y, x)
    if b < 0:
        b += _360_0
    return b or _0_0  # unsigned-0


def atan2d(y, x, reverse=False):
    '''Return C{atan2(B{y}, B{x})} in degrees M{[-180..+180]},
       optionally I{reversed} (by 180 degrees for C{azimuth}s).

       @see: I{Karney}'s C++ function U{Math.atan2d
             <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Math.html>}.
    '''
    d = degrees(_atan2(y, x))  # preserves signed-0
    return _azireversed(d) if reverse else d


def _azireversed(azi):  # in .rhumbBase
    '''(INTERNAL) Return the I{reverse} B{C{azi}} in degrees M{[-180..+180]}.
    '''
    return azi - _copysign(_180_0, azi)


def chain2m(chains):
    '''Convert I{UK} chains to meter.

       @arg chains: Value in chains (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{chains}}.
    '''
    return Meter(Float(chains=chains) * _M.CHAIN)


def circle4(earth, lat):
    '''Get the equatorial or a parallel I{circle of latitude}.

       @arg earth: The earth radius (C{meter}) or an ellipsoid or datum
                   (L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}).
       @arg lat: Geodetic latitude (C{degrees90}, C{str}).

       @return: A L{Circle4Tuple}C{(radius, height, lat, beta)}.

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @raise TypeError: Invalid B{C{earth}}.

       @raise ValueError: invalid B{C{earth}} or B{C{lat}}.
    '''
    E = _MODS.datums._earth_ellipsoid(earth)
    return E.circle4(lat)


def _circle4radius(earth, lat, **radius):
    '''(INTERNAL) Get C{circle4(earth, lat).radius}.
    '''
    e = _xkwds_get(radius, radius=earth) if radius else earth
    return e if e is R_M and not lat else circle4(e, lat).radius


def cot(rad, **raiser_kwds):
    '''Return the C{cotangent} of an angle in C{radians}.

       @arg rad: Angle (C{radians}).
       @kwarg raiser_kwds: Use C{B{raiser}=False} to avoid
                     ValueErrors and optional, additional
                     ValueError keyword argments.

       @return: C{cot(B{rad})}.

       @raise Error: If L{pygeodesy.isnear0}C{(sin(B{rad})}.
    '''
    try:
        return _cotu(*sincos2(rad))
    except ZeroDivisionError:
        return _nonfinite(cot, rad, **raiser_kwds)


def cot_(*rads, **raiser_kwds):
    '''Yield the C{cotangent} of angle(s) in C{radians}.

       @arg rads: One or more angles (each in C{radians}).

       @return: Yield C{cot(B{rad})} for each angle.

       @see: Function L{pygeodesy.cot} for further details.
    '''
    try:
        for r in rads:
            yield _cotu(*sincos2(r))
    except ZeroDivisionError:
        yield _nonfinite(cot_, r, **raiser_kwds)


def cotd(deg, **raiser_kwds):
    '''Return the C{cotangent} of an angle in C{degrees}.

       @arg deg: Angle (C{degrees}).
       @kwarg raiser_kwds: Use C{B{raiser}=False} to avoid
                     ValueErrors and optional, additional
                     ValueError keyword argments.

       @return: C{cot(B{deg})}.

       @raise Error: If L{pygeodesy.isnear0}C{(sin(B{deg})}.
    '''
    try:
        return _cotu(*sincos2d(deg))
    except ZeroDivisionError:
        return _nonfinite(cotd, deg, **raiser_kwds)


def cotd_(*degs, **raiser_kwds):
    '''Yield the C{cotangent} of angle(s) in C{degrees}.

       @arg degs: One or more angles (each in C{degrees}).

       @return: Yield C{cotd(B{deg})} for each angle.

       @see: Function L{pygeodesy.cotd} for further details.
    '''
    try:
        for d in degs:
            yield _cotu(*sincos2d(d))
    except ZeroDivisionError:
        yield _nonfinite(cotd_, d, **raiser_kwds)


def _cotu(s, c):
    '''(INTERNAL) Helper for functions C{cot}, C{cotd}, C{cot_} and C{cotd_}.
    '''
    return _tanu(c, s)


def degrees90(rad):
    '''Convert radians to degrees and wrap M{[-90..+90)}.

       @arg rad: Angle (C{radians}).

       @return: Angle, wrapped (C{degrees90}).
    '''
    return wrap90(degrees(rad))


def degrees180(rad):
    '''Convert radians to degrees and wrap M{[-180..+180)}.

       @arg rad: Angle (C{radians}).

       @return: Angle, wrapped (C{degrees180}).
    '''
    return wrap180(degrees(rad))


def degrees360(rad):
    '''Convert radians to degrees and wrap M{[0..+360)}.

       @arg rad: Angle (C{radians}).

       @return: Angle, wrapped (C{degrees360}).
    '''
    return _umod_360(degrees(rad))


def degrees2grades(deg):
    '''Convert degrees to I{grades} (aka I{gons} or I{gradians}).

       @arg deg: Angle (C{degrees}).

       @return: Angle (C{grades}).
    '''
    return Float(grades=Degrees(deg) * _G.DEG)


def degrees2m(deg, earth=R_M, lat=0, **radius):
    '''Convert an angle to a distance along the equator or along a parallel
       at (geodetic) latitude.

       @arg deg: The angle (C{degrees}).
       @arg earth: The earth radius (C{meter}) or an ellipsoid or datum
                   (L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).
       @kwarg radius: For backward compatibility C{B{radius}=B{earth}}.

       @return: Distance (C{meter}, same units as B{C{earth}} or polar and
                equatorial radii) or C{0.0} for near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @raise TypeError: Invalid B{C{earth}} or B{C{radius}}.

       @raise ValueError: Invalid B{C{deg}}, B{C{earth}} or B{C{lat}}.

       @see: Function L{radians2m} and L{m2degrees}.
    '''
    m = _circle4radius(earth, lat, **radius)
    return _Radians2m(Lamd(deg=deg, clip=0), m)


def fathom2m(fathoms):
    '''Convert I{Imperial} fathom to meter.

       @arg fathoms: Value in fathoms (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{fathoms}}.

       @see: Function L{toise2m}, U{Fathom<https://WikiPedia.org/wiki/Fathom>}
             and U{Klafter<https://WikiPedia.org/wiki/Klafter>}.
    '''
    return Meter(Float(fathoms=fathoms) * _M.FATHOM)


def ft2m(feet, usurvey=False, pied=False, fuss=False):
    '''Convert I{International}, I{US Survey}, I{French} or I{German}
       B{C{feet}} to C{meter}.

       @arg feet: Value in feet (C{scalar}).
       @kwarg usurvey: If C{True}, convert I{US Survey} foot else ...
       @kwarg pied: If C{True}, convert French I{pied-du-Roi} else ...
       @kwarg fuss: If C{True}, convert German I{Fuss}, otherwise
                    I{International} foot to C{meter}.

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{feet}}.
    '''
    return Meter(Feet(feet) * (_M.FOOT_US if usurvey else
                              (_M.FOOT_FR if pied    else
                              (_M.FOOT_GE if fuss    else _M.FOOT))))


def furlong2m(furlongs):
    '''Convert a furlong to meter.

       @arg furlongs: Value in furlongs (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{furlongs}}.
    '''
    return Meter(Float(furlongs=furlongs) * _M.FURLONG)


def gdf(psi):
    '''U{Gudermannian function
       <https://WikiPedia.org/wiki/Gudermannian_function>}.

       @arg psi: Gudermannian (C{float}).

       @return: Angle (C{radians}).

       @see: Function L{agdf}.
    '''
    return atan1(sinh(psi))


def grades(rad):
    '''Convert radians to I{grades} (aka I{gons} or I{gradians}).

       @arg rad: Angle (C{radians}).

       @return: Angle (C{grades}).
    '''
    return Float(grades=Radians(rad) * _G.RAD)


def grades400(rad):
    '''Convert radians to I{grades} (aka I{gons} or I{gradians}) and wrap M{[0..+400)}.

       @arg rad: Angle (C{radians}).

       @return: Angle, wrapped (C{grades}).
    '''
    return Float(grades400=wrapPI2(rad) * _G.RAD)


def grades2degrees(gon):
    '''Convert I{grades} (aka I{gons} or I{gradians}) to C{degrees}.

       @arg gon: Angle (C{grades}).

       @return: Angle (C{degrees}).
    '''
    return Degrees(Float(gon=gon) / _G.DEG)


def grades2radians(gon):
    '''Convert I{grades} (aka I{gons} or I{gradians}) to C{radians}.

       @arg gon: Angle (C{grades}).

       @return: Angle (C{radians}).
    '''
    return Radians(Float(gon=gon) / _G.RAD)


def ha2acre(ha):
    '''Convert hectare to acre.

       @arg ha: Value in hectare (C{scalar}).

       @return: Value in acres (C{float}).

       @raise ValueError: Invalid B{C{ha}}.
    '''
    return m2acre(ha2m2(ha))


def ha2m2(ha):
    '''Convert hectare to I{square} meter.

       @arg ha: Value in hectare (C{scalar}).

       @return: Value in C{meter^2} (C{float}).

       @raise ValueError: Invalid B{C{ha}}.
    '''
    return Meter2(Float(ha=ha) * _M.HA)


def hav(rad):
    '''Return the U{haversine<https://WikiPedia.org/wiki/Haversine_formula>} of an angle.

       @arg rad: Angle (C{radians}).

       @return: C{sin(B{rad} / 2)**2}.
    '''
    return sin(rad * _0_5)**2


def km2m(km):
    '''Convert kilo meter to meter (m).

       @arg km: Value in kilo meter (C{scalar}).

       @return: Value in meter (C{float}).

       @raise ValueError: Invalid B{C{km}}.
    '''
    return Meter(Float(km=km) * _M.KM)


def _loneg(lon):
    '''(INTERNAL) "Complement" of C{lon}.
    '''
    return _180_0 - lon


def m2acre(meter2):
    '''Convert I{square} meter to acres.

       @arg meter2: Value in C{meter^2} (C{scalar}).

       @return: Value in acres (C{float}).

       @raise ValueError: Invalid B{C{meter2}}.
    '''
    return Float(acre=Meter2(meter2) / _M.ACRE)


def m2chain(meter):
    '''Convert meter to I{UK} chains.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{chains} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Float(chain=Meter(meter) / _M.CHAIN)  # * 0.049_709_695_378_986_715


def m2degrees(distance, earth=R_M, lat=0, **radius):
    '''Convert a distance to an angle along the equator or along a parallel
       at (geodetic) latitude.

       @arg distance: Distance (C{meter}, same units as B{C{radius}}).
       @kwarg earth: Mean earth radius (C{meter}) or an ellipsoid or datum
                     (L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).
       @kwarg radius: For backward compatibility C{B{radius}=B{earth}}.

       @return: Angle (C{degrees}) or C{INF} or C{NINF} for near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @raise TypeError: Invalid B{C{earth}} or B{C{radius}}.

       @raise ValueError: Invalid B{C{distance}}, B{C{earth}} or B{C{lat}}.

       @see: Function L{m2radians} and L{degrees2m}.
    '''
    return degrees(m2radians(distance, earth=earth, lat=lat, **radius))


def m2fathom(meter):
    '''Convert meter to I{Imperial} fathoms.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{fathoms} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.

       @see: Function L{m2toise}, U{Fathom<https://WikiPedia.org/wiki/Fathom>}
             and U{Klafter<https://WikiPedia.org/wiki/Klafter>}.
    '''
    return Float(fathom=Meter(meter) / _M.FATHOM)  # * 0.546_806_649


def m2ft(meter, usurvey=False, pied=False, fuss=False):
    '''Convert meter to I{International}, I{US Survey}, I{French} or
       or I{German} feet (C{ft}).

       @arg meter: Value in meter (C{scalar}).
       @kwarg usurvey: If C{True}, convert to I{US Survey} foot else ...
       @kwarg pied: If C{True}, convert to French I{pied-du-Roi} else ...
       @kwarg fuss: If C{True}, convert to German I{Fuss}, otherwise to
                    I{International} foot.

       @return: Value in C{feet} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    # * 3.280_833_333_333_3333, US Survey 3_937 / 1_200
    # * 3.280_839_895_013_1235, Int'l 10_000 / (254 * 12)
    return Float(feet=Meter(meter) / (_M.FOOT_US if usurvey else
                                     (_M.FOOT_FR if pied    else
                                     (_M.FOOT_GE if fuss    else _M.FOOT))))


def m2furlong(meter):
    '''Convert meter to furlongs.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{furlongs} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Float(furlong=Meter(meter) / _M.FURLONG)  # * 0.004_970_969_54


def m2ha(meter2):
    '''Convert I{square} meter to hectare.

       @arg meter2: Value in C{meter^2} (C{scalar}).

       @return: Value in hectare (C{float}).

       @raise ValueError: Invalid B{C{meter2}}.
    '''
    return Float(ha=Meter2(meter2) / _M.HA)


def m2km(meter):
    '''Convert meter to kilo meter (Km).

       @arg meter: Value in meter (C{scalar}).

       @return: Value in Km (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Float(km=Meter(meter) / _M.KM)


def m2NM(meter):
    '''Convert meter to nautical miles (NM).

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{NM} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Float(NM=Meter(meter) / _M.NM)  # * 5.399_568_04e-4


def m2radians(distance, earth=R_M, lat=0, **radius):
    '''Convert a distance to an angle along the equator or along a parallel
       at (geodetic) latitude.

       @arg distance: Distance (C{meter}, same units as B{C{radius}}).
       @kwarg earth: Mean earth radius (C{meter}) or an ellipsoid or datum
                     (L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).
       @kwarg radius: For backward compatibility C{B{radius}=B{earth}}.

       @return: Angle (C{radians}) or C{INF} or C{NINF} for near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @raise TypeError: Invalid B{C{earth}} or B{C{radius}}.

       @raise ValueError: Invalid B{C{distance}}, B{C{earth}} or B{C{lat}}.

       @see: Function L{m2degrees} and L{radians2m}.
    '''
    m = _circle4radius(earth, lat, **radius)
    return _copysign_0_0(distance) if m < EPS0 else \
            Radians(Float(distance=distance) / m)


def m2SM(meter):
    '''Convert meter to statute miles (SM).

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{SM} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Float(SM=Meter(meter) / _M.SM)  # * 6.213_699_49e-4 == 1 / 1_609.344


def m2toise(meter):
    '''Convert meter to French U{toises<https://WikiPedia.org/wiki/Toise>}.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{toises} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.

       @see: Function L{m2fathom}.
    '''
    return Float(toise=Meter(meter) / _M.TOISE)  # * 0.513_083_632_632_119


def m2yard(meter):
    '''Convert meter to I{UK} yards.

       @arg meter: Value in meter (C{scalar}).

       @return: Value in C{yards} (C{float}).

       @raise ValueError: Invalid B{C{meter}}.
    '''
    return Float(yard=Meter(meter) / _M.YARD_UK)  # * 1.093_613_298_337_707_8


def NM2m(nm):
    '''Convert nautical miles to meter (m).

       @arg nm: Value in nautical miles (C{scalar}).

       @return: Value in meter (C{float}).

       @raise ValueError: Invalid B{C{nm}}.
    '''
    return Meter(Float(nm=nm) * _M.NM)


def _nonfinite(where, x, raiser=True, **kwds):  # PYCHOK no cover
    '''(INTERNAL) Raise a C{_ValueError} or return C{INF} or C{NINF}.
    '''
    if raiser:
        t =  typename(where)
        t = _MODS.streprs.Fmt.PAREN(t, x)
        raise _ValueError(t, **kwds)
    return _copysignINF(x)


def radians2m(rad, earth=R_M, lat=0, **radius):
    '''Convert an angle to a distance along the equator or along a parallel
       at (geodetic) latitude.

       @arg rad: The angle (C{radians}).
       @kwarg earth: Mean earth radius (C{meter}) or an ellipsoid or datum
                     (L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}).
       @kwarg lat: Parallel latitude (C{degrees90}, C{str}).
       @kwarg radius: For backward compatibility C{B{radius}=B{earth}}.

       @return: Distance (C{meter}, same units as B{C{earth}} or polar and
                equatorial radii) or C{0.0} for near-polar B{C{lat}}.

       @raise RangeError: Latitude B{C{lat}} outside valid range and
                          L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @raise TypeError: Invalid B{C{earth}} or B{C{radius}}.

       @raise ValueError: Invalid B{C{rad}}, B{C{earth}} or B{C{lat}}.

       @see: Function L{degrees2m} and L{m2radians}.
    '''
    m = _circle4radius(earth, lat, **radius)
    return _Radians2m(Lam(rad=rad, clip=0), m)


def _Radians2m(rad, m):
    '''(INTERNAL) Helper for C{degrees2m} and C{radians2m}.
    '''
    return _copysign_0_0(rad) if m < EPS0 else (rad * m)


def radiansPI(deg):
    '''Convert and wrap degrees to radians M{[-PI..+PI]}.

       @arg deg: Angle (C{degrees}).

       @return: Radians, wrapped (C{radiansPI})
    '''
    return wrapPI(radians(deg))


def radiansPI2(deg):
    '''Convert and wrap degrees to radians M{[0..+2PI)}.

       @arg deg: Angle (C{degrees}).

       @return: Radians, wrapped (C{radiansPI2})
    '''
    return _umod_PI2(radians(deg))


def radiansPI_2(deg):
    '''Convert and wrap degrees to radians M{[-3PI/2..+PI/2]}.

       @arg deg: Angle (C{degrees}).

       @return: Radians, wrapped (C{radiansPI_2})
    '''
    return wrapPI_2(radians(deg))


def _sin0cos2(q, r, sign):
    '''(INTERNAL) 2-tuple (C{sin(r), cos(r)}) in quadrant C{0 <= B{q} <= 3} and
       C{sin} zero I{signed} with B{C{sign}}, like Karney's U{Math.sind and .cosd
       <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Math.html>}
    '''
    q &= 3  # -1: 3, -2: 2, -3: 1, -4: 0 ...
    if r < PI_2:  # Quarter turn
        a = fabs(r)
        if a:
            s = sin(r)
            if a == PI_4:
                s, c = _copysign(_SIN_45, s), _COS_45
            elif a == PI_6:
                s, c = _copysign(_SIN_30, s), _COS_30
            else:
                c = cos(r)
        else:
            s, c = _0_0, _1_0
    else:
        s, c = _1_0, _0_0
    t = s, c, -s, -c, s
    s = t[q]     or _copysign_0_0(sign)
    c = t[q + 1] or _0_0
    return s, c


def SinCos2(x, unit=Radians):
    '''Get C{sin} and C{cos} of I{typed} angle.

       @arg x: Angle (L{Degrees}, L{Radians} or scalar C{radians}).
       @kwarg unit: The C{B{x}} unit (L{Degrees}, L{Radians}).

       @return: 2-Tuple (C{sin(B{x})}, C{cos(B{x})}).
    '''
    return sincos2d(x) if unit is Degrees or unit is degrees or isinstanceof(x, Degrees, Degrees_) else (
#          sincos2(x)  if unit is Radians or unit is radians or isinstanceof(x, Radians, Radians_) else
           sincos2(Radians(x)))  # assume C{radians}


def sincos2(rad):
    '''Return the C{sine} and C{cosine} of an angle in C{radians}.

       @arg rad: Angle (C{radians}).

       @return: 2-Tuple (C{sin(B{rad})}, C{cos(B{rad})}).

       @see: U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/
             classGeographicLib_1_1Math.html#sincosd>} function U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             python/geographiclib/geomath.py#l155>} and C++ U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             include/GeographicLib/Math.hpp#l558>}.
    '''
    if _isfinite(rad):
        q, r = divmod(rad, PI_2)
        t = _sin0cos2(int(q), r, rad)
    else:
        t =  NAN, NAN
    return t


def sincos2_(*rads):
    '''Yield the C{sine} and C{cosine} of angle(s) in C{radians}.

       @arg rads: One or more angles (C{radians}).

       @return: Yield C{sin(B{rad})} and C{cos(B{rad})} for each angle.

       @see: Function L{sincos2}.
    '''
    for r in rads:
        s, c = sincos2(r)
        yield s
        yield c


def sincos2d(deg, adeg=_0_0):
    '''Return the C{sine} and C{cosine} of an angle in C{degrees}.

       @arg deg: Angle (C{degrees}).
       @kwarg adeg: Optional correction (C{degrees}).

       @return: 2-Tuple (C{sin(B{deg_})}, C{cos(B{deg_})} with C{B{deg_} =
                B{deg} + B{adeg}}).

       @see: U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/
             classGeographicLib_1_1Math.html#sincosd>} function U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             python/geographiclib/geomath.py#l155>} and C++ U{sincosd
             <https://SourceForge.net/p/geographiclib/code/ci/release/tree/
             include/GeographicLib/Math.hpp#l558>}.
    '''
    if _isfinite(deg):
        q, d = divmod(deg, _90_0)
        if adeg:
            d = _MODS.karney._around(d + adeg)
        t = _sin0cos2(int(q), radians(d), deg)
    else:
        t =  NAN, NAN
    return t


def sincos2d_(*degs):
    '''Yield the C{sine} and C{cosine} of angle(s) in C{degrees}.

       @arg degs: One or more angles (C{degrees}).

       @return: Yield C{sind(B{deg})} and C{cosd(B{deg})} for each angle.

       @see: Function L{sincos2d}.
    '''
    for d in degs:
        s, c = sincos2d(d)
        yield s
        yield c


def sincostan3(rad):
    '''Return the C{sine}, C{cosine} and C{tangent} of an angle in C{radians}.

       @arg rad: Angle (C{radians}).

       @return: 3-Tuple (C{sin(B{rad})}, C{cos(B{rad})}, C{tan(B{rad})}).

       @see: Function L{sincos2}.
    '''
    return _sincostan3(*sincos2(float(rad)))


def _sincostan3(s, c):
    '''(INTERNAL) Helper for C{sincostan3} and C{sincostan3d}.
    '''
    return s, c, _over_1(s, c)  # _tanu(s, c, raiser=False)


def sincostan3d(deg):
    '''Return the C{sine}, C{cosine} and C{tangent} of an angle in C{degrees}.

       @arg deg: Angle (C{degrees}).

       @return: 3-Tuple (C{sind(B{deg})}, C{cosd(B{deg})}, C{tand(B{deg})}).

       @see: Function L{sincos2d}.
    '''
    return _sincostan3(*sincos2d(float(deg)))


def SM2m(sm):
    '''Convert statute miles to meter (m).

       @arg sm: Value in statute miles (C{scalar}).

       @return: Value in meter (C{float}).

       @raise ValueError: Invalid B{C{sm}}.
    '''
    return Meter(Float(sm=sm) * _M.SM)


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
    if semi and isnear0(fabs(rad) - PI):
        for n, v in semi.items():
            break
        n = _SPACE_(n, _radians_) if not isint(v) else \
            _SPACE_(_MODS.streprs.Fmt.SQUARE(**semi), _edge_)
        raise _ValueError(n, rad, txt=_semi_circular_)

    return _tan(rad * _0_5) if _isfinite(rad) else NAN


def tan(rad, **raiser_kwds):
    '''Return the C{tangent} of an angle in C{radians}.

       @arg rad: Angle (C{radians}).
       @kwarg raiser_kwds: Use C{B{raiser}=False} to avoid
                     ValueErrors and optional, additional
                     ValueError keyword argments.

       @return: C{tan(B{rad})}.

       @raise Error: If L{pygeodesy.isnear0}C{(cos(B{rad})}.
    '''
    try:
        return _tanu(*sincos2(rad))
    except ZeroDivisionError:
        return _nonfinite(tan, rad, **raiser_kwds)


def tan_(*rads, **raiser_kwds):
    '''Yield the C{tangent} of angle(s) in C{radians}.

       @arg rads: One or more angles (each in C{radians}).

       @return: Yield C{tan(B{rad})} for each angle.

       @see: Function L{pygeodesy.tan} for futher details.
    '''
    try:
        for r in rads:
            yield _tanu(*sincos2(r))
    except ZeroDivisionError:
        yield _nonfinite(tan_, r, **raiser_kwds)


def tand(deg, **raiser_clamp_kwds):
    '''Return the C{tangent} of an angle in C{degrees}.

       @arg deg: Angle (C{degrees}).
       @kwarg raiser_clamp_kwds: Use C{B{raiser}=False} to avoid
                     ValueErrors, C{B{clamp}=}L{OVERFLOW} and
                     optional, additional ValueError keyword
                     argments.

       @return: C{tan(B{deg})}.

       @raise Error: If L{pygeodesy.isnear0}C{(cos(B{deg})}.
    '''
    try:
        return _tanu(*sincos2d(deg))
    except ZeroDivisionError:
        return _nonfinite(tand, deg, **raiser_clamp_kwds)


def tand_(*degs, **raiser_clamp_kwds):
    '''Yield the C{tangent} of angle(s) in C{degrees}.

       @arg degs: One or more angles (each in C{degrees}).

       @return: Yield C{tand(B{deg})} for each angle.

       @see: Function L{pygeodesy.tand} for futher details.
    '''
    try:
        for d in degs:
            yield _tanu(*sincos2d(d), **raiser_clamp_kwds)
    except ZeroDivisionError:
        yield _nonfinite(tand_, d, **raiser_clamp_kwds)


def tanPI_2_2(rad):
    '''Compute the tangent of half angle, 90 degrees rotated.

       @arg rad: Angle (C{radians}).

       @return: M{tan((rad + PI/2) / 2)} (C{float}).
    '''
    return _tan((rad + PI_2) * _0_5) if _isfinite(rad) else (
            NAN if isnan(rad) else _copysign(_90_0, rad))


def _tanu(s, c, clamp=OVERFLOW, **unused):
    '''(INTERNAL) Helper for functions C{_cotu}, C{sincostan3},
       C{sincostan3d}, C{tan}, C{tan_}, C{tand} and C{tand_}.
    '''
    if s is NAN or isnan(s):
        t = NAN
    elif isnear0(c):
        raise ZeroDivisionError()
    else:
        t = _over_1(s, c)
        if clamp:
            t = min(clamp, max(t, -clamp))
    return t


def toise2m(toises):
    '''Convert French U{toises<https://WikiPedia.org/wiki/Toise>} to meter.

       @arg toises: Value in toises (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{toises}}.

       @see: Function L{fathom2m}.
    '''
    return Meter(Float(toises=toises) * _M.TOISE)


def truncate(x, ndigits=None):
    '''Truncate to the given number of digits.

       @arg x: Value to truncate (C{scalar}).
       @kwarg ndigits: Number of digits (C{int}),
                       aka I{precision}.

       @return: Truncated B{C{x}} (C{float}).

       @see: Python function C{round}.
    '''
    if isint(ndigits):
        p = _10_0**ndigits
        x =  int(x * p) / p
    return x


def unroll180(lon1, lon2, wrap=True):
    '''Unroll longitudinal delta and wrap longitude in degrees.

       @arg lon1: Start longitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg wrap: If C{True}, wrap and unroll to the M{(-180..+180]}
                    range (C{bool}).

       @return: 2-Tuple C{(B{lon2}-B{lon1}, B{lon2})} unrolled (C{degrees},
                C{degrees}).

       @see: Capability C{LONG_UNROLL} in U{GeographicLib
             <https://GeographicLib.SourceForge.io/html/python/interface.html#outmask>}.
    '''
    d = lon2 - lon1
    if wrap:
        u = wrap180(d)
        if u != d:
            return u, (lon1 + u)
    return d, lon2


def _unrollon(p1, p2, wrap=False):  # unroll180 == .karney._unroll2
    '''(INTERNAL) Wrap/normalize, unroll and replace longitude.
    '''
    lat, lon = p2.lat, p2.lon
    if wrap and _Wrap.normal:
        lat, lon = _Wrap.latlon(lat, lon)
    _, lon = unroll180(p1.lon, lon, wrap=True)
    if lat != p2.lat or fabs(lon - p2.lon) > EPS:
        p2 = p2.dup(lat=lat, lon=wrap180(lon))
        # p2 = p2.copy(); p2.latlon = lat, wrap180(lon)
    return p2


def _unrollon3(p1, p2, p3, wrap=False):
    '''(INTERNAL) Wrap/normalize, unroll 2 points.
    '''
    w = wrap
    if w:  # PYCHOK no cover
        w  = _Wrap.normal
        p2 = _unrollon(p1, p2, wrap=w)
        p3 = _unrollon(p1, p3, wrap=w)
        p2 = _unrollon(p2, p3)
    return p2, p3, w  # was wrapped?


def unrollPI(rad1, rad2, wrap=True):
    '''Unroll longitudinal delta and wrap longitude in radians.

       @arg rad1: Start longitude (C{radians}).
       @arg rad2: End longitude (C{radians}).
       @kwarg wrap: If C{True}, wrap and unroll to the M{(-PI..+PI]}
                    range (C{bool}).

       @return: 2-Tuple C{(B{rad2}-B{rad1}, B{rad2})} unrolled
                (C{radians}, C{radians}).

       @see: Capability C{LONG_UNROLL} in U{GeographicLib
             <https://GeographicLib.SourceForge.io/html/python/interface.html#outmask>}.
    '''
    r = rad2 - rad1
    if wrap:
        u = wrapPI(r)
        if u != r:
            return u, (rad1 + u)
    return r, rad2


class _Wrap(object):

    _normal = False  # default

    @property
    def normal(self):
        '''Get the current L{normal} setting (C{True},
           C{False} or C{None}).
        '''
        return self._normal

    @normal.setter  # PYCHOK setter!
    def normal(self, setting):  # PYCHOK no cover
        '''Set L{normal} to C{True}, C{False} or C{None}.
        '''
        m = _MODS.formy
        t = {True:  (m.normal,        m.normal_),
             False: (self.wraplatlon, self.wraphilam),
             None:  (_passargs,      _passargs)}.get(setting, ())
        if t:
            self.latlon, self.philam = t
            self._normal = setting

    def latlonDMS2(self, lat, lon, **DMS2_kwds):  # PYCHOK no cover
        if isstr(lat) or isstr(lon):
            kwds = _xkwds(DMS2_kwds, clipLon=0, clipLat=0)
            lat, lon = _MODS.dms.parseDMS2(lat, lon, **kwds)
        return self.latlon(lat, lon)

#   def normalatlon(self, *latlon):
#       return _MODS.formy.normal(*latlon)

#   def normalamphi(self, *philam):
#       return _MODS.formy.normal_(*philam)

    def wraplatlon(self, lat, lon):
        return wrap90(lat), wrap180(lon)

    latlon = wraplatlon  # default

    def latlon3(self, lon1, lat2, lon2, wrap):
        if wrap:
            lat2,  lon2 = self.latlon(lat2, lon2)
            lon21, lon2 = unroll180(lon1, lon2)
        else:
            lon21 = lon2 - lon1
        return lon21, lat2, lon2

    def _latlonop(self, wrap):
        if wrap and self._normal is not None:
            return self.latlon
        else:
            return _passargs

    def wraphilam(self, phi, lam):
        return wrapPI_2(phi), wrapPI(lam)

    philam = wraphilam  # default

    def philam3(self, lam1, phi2, lam2, wrap):
        if wrap:
            phi2,  lam2 = self.philam(phi2, lam2)
            lam21, lam2 = unrollPI(lam1, lam2)
        else:
            lam21 = lam2 - lam1
        return lam21, phi2, lam2

    def _philamop(self, wrap):
        if wrap and self._normal is not None:
            return self.philam
        else:
            return _passargs

    def point(self, ll, wrap=True):  # in .points._fractional, -.PointsIter.iterate, ...
        '''Return C{ll} or a copy, I{normalized} or I{wrap}'d.
        '''
        if wrap and self._normal is not None:
            lat, lon = ll.latlon
            if fabs(lon) > _180_0 or fabs(lat) > _90_0:
                _n = self.latlon
                ll = ll.copy(name=typename(_n))
                ll.latlon = _n(lat, lon)
        return ll

_Wrap = _Wrap()  # PYCHOK singleton


# def _wrap(angle, wrap, modulo):
#     '''(INTERNAL) Angle wrapper M{((wrap-modulo)..+wrap]}.
#
#        @arg angle: Angle (C{degrees}, C{radians} or C{grades}).
#        @arg wrap: Range (C{degrees}, C{radians} or C{grades}).
#        @arg modulo: Upper limit (360 C{degrees}, PI2 C{radians} or 400 C{grades}).
#
#        @return: The B{C{angle}}, wrapped (C{degrees}, C{radians} or C{grades}).
#     '''
#     a = float(angle)
#     if not (wrap - modulo) <= a < wrap:
#         # math.fmod(-1.5, 3.14) == -1.5, but -1.5 % 3.14 == 1.64
#         # math.fmod(-1.5, 360)  == -1.5, but -1.5 % 360 == 358.5
#         a %= modulo
#         if a > wrap:
#             a -= modulo
#     return a


def wrap90(deg):
    '''Wrap degrees to M{[-90..+90]}.

       @arg deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees90}).
    '''
    return _wrapu(wrap180(deg), _180_0, _90_0)


def wrap180(deg):
    '''Wrap degrees to M{[-180..+180]}.

       @arg deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees180}).
    '''
    d =  float(deg)
    w = _umod_360(d)
    if w > _180_0:
        w -= _360_0
    elif d < 0 and w == _180_0:
        w = -w
    return w


def wrap360(deg):  # see .streprs._umod_360
    '''Wrap degrees to M{[0..+360)}.

       @arg deg: Angle (C{degrees}).

       @return: Degrees, wrapped (C{degrees360}).
    '''
    return _umod_360(float(deg))


def wrapPI(rad):
    '''Wrap radians to M{[-PI..+PI]}.

       @arg rad: Angle (C{radians}).

       @return: Radians, wrapped (C{radiansPI}).
    '''
    r =  float(rad)
    w = _umod_PI2(r)
    if w > PI:
        w -= PI2
    elif r < 0 and w == PI:
        w = -PI
    return w


def wrapPI2(rad):
    '''Wrap radians to M{[0..+2PI)}.

       @arg rad: Angle (C{radians}).

       @return: Radians, wrapped (C{radiansPI2}).
    '''
    return _umod_PI2(float(rad))


def wrapPI_2(rad):
    '''Wrap radians to M{[-PI/2..+PI/2]}.

       @arg rad: Angle (C{radians}).

       @return: Radians, wrapped (C{radiansPI_2}).
    '''
    return _wrapu(wrapPI(rad), PI, PI_2)


# def wraplatlon(lat, lon):
#     '''Both C{wrap90(B{lat})} and C{wrap180(B{lon})}.
#     '''
#     return wrap90(lat), wrap180(lon)


def wrap_normal(*normal):
    '''Define the operation for the keyword argument C{B{wrap}=True},
       across L{pygeodesy}: I{wrap}, I{normalize} or I{no-op}.  For
       backward compatibility, the default is I{wrap}.

       @arg normal: If C{True}, I{normalize} lat- and longitude using
                    L{normal} or L{normal_}, if C{False}, I{wrap} the
                    lat- and longitude individually by L{wrap90} or
                    L{wrapPI_2} respectively L{wrap180}, L{wrapPI} or
                    if C{None}, leave lat- and longitude I{unchanged}.
                    To get the current setting, do not specify.

       @return: The previous L{wrap_normal} setting (C{bool} or C{None}).
    '''
    t = _Wrap.normal
    if normal:
        _Wrap.normal = normal[0]
    return t


# def wraphilam(phi, lam,):
#     '''Both C{wrapPI_2(B{phi})} and C{wrapPI(B{lam})}.
#     '''
#     return wrapPI_2(phi), wrapPI(lam)


def _wrapu(w, H, Q):
    '''(INTERNAL) Helper for functions C{wrap180} and C{wrapPI}.
    '''
    return (w - H) if w > Q else ((w + H) if w < (-Q) else w)


def yard2m(yards):
    '''Convert I{UK} yards to meter.

       @arg yards: Value in yards (C{scalar}).

       @return: Value in C{meter} (C{float}).

       @raise ValueError: Invalid B{C{yards}}.
    '''
    return Float(yards=yards) * _M.YARD_UK

# **) MIT License
#
# Copyright (C) 2016-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
