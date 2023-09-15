
# -*- coding: utf-8 -*-

u'''Formulary of basic geodesy functions and approximations.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

# from pygeodesy.basics import isscalar  # from .fsums
# from pygeodesy.cartesianBase import CartesianBase  # _MODS
from pygeodesy.constants import EPS, EPS0, EPS1, PI, PI2, PI3, PI_2, R_M, \
                               _umod_PI2, float0_, isnon0, remainder, \
                               _0_0, _0_125, _0_25, _0_5, _1_0, _2_0, \
                               _4_0, _32_0, _90_0, _180_0, _360_0
from pygeodesy.datums import Datum, Ellipsoid, _ellipsoidal_datum, \
                            _mean_radius, _spherical_datum, _WGS84,  _EWGS84
# from pygeodesy.ellipsoids import Ellipsoid, _EWGS84  # from .datums
from pygeodesy.errors import IntersectionError, LimitError, limiterrors, \
                            _TypeError, _ValueError, _xattr, _xError, \
                            _xkwds, _xkwds_pop
from pygeodesy.fmath import euclid, hypot, hypot2, sqrt0
from pygeodesy.fsums import fsumf_,  isscalar
from pygeodesy.interns import NN, _delta_, _distant_, _SPACE_, _too_
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _NamedTuple, _xnamed,  Fmt, unstr
from pygeodesy.namedTuples import Bearing2Tuple, Distance4Tuple, \
                                  Intersection3Tuple, LatLon2Tuple, \
                                  PhiLam2Tuple, Vector3Tuple
# from pygeodesy.streprs import Fmt, unstr  # from .named
# from pygeodesy.triaxials import _hartzell3d2  # _MODS
from pygeodesy.units import Bearing, Degrees_, Distance, Distance_, Height, \
                            Lam_, Lat, Lon, Meter_, Phi_, Radians, Radians_, \
                            Radius, Radius_, Scalar, _100km
from pygeodesy.utily import acos1, atan2b, atan2d, degrees2m, _loneg, m2degrees, \
                            tan_2, sincos2, sincos2_, sincos2d_, _Wrap
# from pygeodesy.vector3d import _otherV3d  # _MODS
# from pygeodesy import ellipsoidalExact, ellipsoidalKarney, vector3d, \
#                       sphericalNvector, sphericalTrigonometry  # _MODS

from contextlib import contextmanager
from math import asin, atan, atan2, cos, degrees, fabs, radians, sin, sqrt  # pow

__all__ = _ALL_LAZY.formy
__version__ = '23.09.06'

_D2_R2  = (PI / _180_0)**2  # degrees- to radians-squared
_ratio_ = 'ratio'
_xline_ = 'xline'


def _anti2(a, b, n_2, n, n2):
    '''(INTERNAL) Helper for C{antipode} and C{antipode_}.
    '''
    r = remainder(a, n) if fabs(a) > n_2 else a
    if r == a:
        r = -r
        b += n
    if fabs(b) > n:
        b = remainder(b, n2)
    return float0_(r, b)


def antipode(lat, lon, name=NN):
    '''Return the antipode, the point diametrically opposite
       to a given point in C{degrees}.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg name: Optional name (C{str}).

       @return: A L{LatLon2Tuple}C{(lat, lon)}.

       @see: Functions L{antipode_} and L{normal} and U{Geosphere
             <https://CRAN.R-Project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return LatLon2Tuple(*_anti2(lat, lon, _90_0, _180_0, _360_0), name=name)


def antipode_(phi, lam, name=NN):
    '''Return the antipode, the point diametrically opposite
       to a given point in C{radians}.

       @arg phi: Latitude (C{radians}).
       @arg lam: Longitude (C{radians}).
       @kwarg name: Optional name (C{str}).

       @return: A L{PhiLam2Tuple}C{(phi, lam)}.

       @see: Functions L{antipode} and L{normal_} and U{Geosphere
             <https://CRAN.R-Project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return PhiLam2Tuple(*_anti2(phi, lam, PI_2, PI, PI2), name=name)


def bearing(lat1, lon1, lat2, lon2, **final_wrap):
    '''Compute the initial or final bearing (forward or reverse
       azimuth) between a (spherical) start and end point.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg final_wrap: Optional keyword arguments for function
                          L{pygeodesy.bearing_}.

       @return: Initial or final bearing (compass C{degrees360}) or
                zero if start and end point coincide.
    '''
    r = bearing_(Phi_(lat1=lat1), Lam_(lon1=lon1),
                 Phi_(lat2=lat2), Lam_(lon2=lon2), **final_wrap)
    return degrees(r)


def bearing_(phi1, lam1, phi2, lam2, final=False, wrap=False):
    '''Compute the initial or final bearing (forward or reverse azimuth)
       between a (spherical) start and end point.

       @arg phi1: Start latitude (C{radians}).
       @arg lam1: Start longitude (C{radians}).
       @arg phi2: End latitude (C{radians}).
       @arg lam2: End longitude (C{radians}).
       @kwarg final: Return final bearing if C{True}, initial otherwise (C{bool}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{phi2}} and
                    B{C{lam2}} (C{bool}).

       @return: Initial or final bearing (compass C{radiansPI2}) or zero if start
                and end point coincide.

       @see: U{Bearing<https://www.Movable-Type.co.UK/scripts/latlong.html>}, U{Course
             between two points<https://www.EdWilliams.org/avform147.htm#Crs>} and
             U{Bearing Between Two Points<https://web.Archive.org/web/20020630205931/
             https://MathForum.org/library/drmath/view/55417.html>}.
    '''
    db, phi2, lam2 = _Wrap.philam3(lam1, phi2, lam2, wrap)
    if final:  # swap plus PI
        phi1, lam1, phi2, lam2, db = phi2, lam2, phi1, lam1, -db
        r = PI3
    else:
        r = PI2
    sa1, ca1, sa2, ca2, sdb, cdb = sincos2_(phi1, phi2, db)

    x = ca1 * sa2 - sa1 * ca2 * cdb
    y = sdb * ca2
    return _umod_PI2(atan2(y, x) + r)  # .utily.wrapPI2


def _bearingTo2(p1, p2, wrap=False):  # for points.ispolar, sphericalTrigonometry.areaOf
    '''(INTERNAL) Compute initial and final bearing.
    '''
    try:  # for LatLon_ and ellipsoidal LatLon
        return p1.bearingTo2(p2, wrap=wrap)
    except AttributeError:
        pass
    # XXX spherical version, OK for ellipsoidal ispolar?
    a1, b1 = p1.philam
    a2, b2 = p2.philam
    return Bearing2Tuple(degrees(bearing_(a1, b1, a2, b2, final=False, wrap=wrap)),
                         degrees(bearing_(a1, b1, a2, b2, final=True,  wrap=wrap)),
                         name=_bearingTo2.__name__)


def compassAngle(lat1, lon1, lat2, lon2, adjust=True, wrap=False):
    '''Return the angle from North for the direction vector M{(lon2 - lon1,
       lat2 - lat1)} between two points.

       Suitable only for short, not near-polar vectors up to a few hundred
       Km or Miles.  Use function L{pygeodesy.bearing} for longer vectors.

       @arg lat1: From latitude (C{degrees}).
       @arg lon1: From longitude (C{degrees}).
       @arg lat2: To latitude (C{degrees}).
       @arg lon2: To longitude (C{degrees}).
       @kwarg adjust: Adjust the longitudinal delta by the cosine of the
                      mean latitude (C{bool}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}}
                    and B{C{lon2}} (C{bool}).

       @return: Compass angle from North (C{degrees360}).

       @note: Courtesy of Martin Schultz.

       @see: U{Local, flat earth approximation
             <https://www.EdWilliams.org/avform.htm#flat>}.
    '''
    d_lon, lat2, lon2 = _Wrap.latlon3(lon1, lat2, lon2, wrap)
    if adjust:  # scale delta lon
        d_lon *= _scale_deg(lat1, lat2)
    return atan2b(d_lon, lat2 - lat1)


def cosineAndoyerLambert(lat1, lon1, lat2, lon2, datum=_WGS84, wrap=False):
    '''Compute the distance between two (ellipsoidal) points using the
       U{Andoyer-Lambert correction<https://NavLib.net/wp-content/uploads/2013/10/
       admiralty-manual-of-navigation-vol-1-1964-english501c.pdf>} of the U{Law of
       Cosines<https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>} formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg datum: Datum (L{Datum}) or ellipsoid (L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}) to use.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                    B{C{lat2}} and B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as the B{C{datum}}'s
                ellipsoid axes or C{radians} if B{C{datum}} is C{None}).

       @raise TypeError: Invalid B{C{datum}}.

       @see: Functions L{cosineAndoyerLambert_}, L{cosineForsytheAndoyerLambert},
             L{cosineLaw}, L{equirectangular}, L{euclidean}, L{flatLocal}/L{hubeny},
             L{flatPolar}, L{haversine}, L{thomas} and L{vincentys} and method
             L{Ellipsoid.distance2}.
    '''
    return _dE(cosineAndoyerLambert_, datum, wrap, lat1, lon1, lat2, lon2)


def cosineAndoyerLambert_(phi2, phi1, lam21, datum=_WGS84):
    '''Compute the I{angular} distance between two (ellipsoidal) points using the
       U{Andoyer-Lambert correction<https://NavLib.net/wp-content/uploads/2013/10/
       admiralty-manual-of-navigation-vol-1-1964-english501c.pdf>} of the U{Law of
       Cosines<https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>} formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg datum: Datum (L{Datum}) or ellipsoid (L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}) to use.

       @return: Angular distance (C{radians}).

       @raise TypeError: Invalid B{C{datum}}.

       @see: Functions L{cosineAndoyerLambert}, L{cosineForsytheAndoyerLambert_},
             L{cosineLaw_}, L{equirectangular_}, L{euclidean_}, L{flatLocal_}/L{hubeny_},
             L{flatPolar_}, L{haversine_}, L{thomas_} and L{vincentys_} and U{Geodesy-PHP
             <https://GitHub.com/jtejido/geodesy-php/blob/master/src/Geodesy/Distance/
             AndoyerLambert.php>}.
    '''
    s2, c2, s1, c1, r, c21 = _sincosa6(phi2, phi1, lam21)
    if isnon0(c1) and isnon0(c2):
        E = _ellipsoidal(datum, cosineAndoyerLambert_)
        if E.f:  # ellipsoidal
            r2 = atan2(E.b_a * s2, c2)
            r1 = atan2(E.b_a * s1, c1)
            s2, c2, s1, c1 = sincos2_(r2, r1)
            r = acos1(s1 * s2 + c1 * c2 * c21)
            if r:
                sr, _, sr_2, cr_2 = sincos2_(r, r * _0_5)
                if isnon0(sr_2) and isnon0(cr_2):
                    s  = (sr + r) * ((s1 - s2) / sr_2)**2
                    c  = (sr - r) * ((s1 + s2) / cr_2)**2
                    r += (c - s) * E.f * _0_125
    return r


def cosineForsytheAndoyerLambert(lat1, lon1, lat2, lon2, datum=_WGS84, wrap=False):
    '''Compute the distance between two (ellipsoidal) points using the
       U{Forsythe-Andoyer-Lambert correction<https://www2.UNB.Ca/gge/Pubs/TR77.pdf>} of
       the U{Law of Cosines<https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
       formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg datum: Datum (L{Datum}) or ellipsoid (L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}) to use.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                    B{C{lat2}} and B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as the B{C{datum}}'s
                ellipsoid axes or C{radians} if B{C{datum}} is C{None}).

       @raise TypeError: Invalid B{C{datum}}.

       @see: Functions L{cosineForsytheAndoyerLambert_}, L{cosineAndoyerLambert},
             L{cosineLaw}, L{equirectangular}, L{euclidean}, L{flatLocal}/L{hubeny},
             L{flatPolar}, L{haversine}, L{thomas} and L{vincentys} and method
             L{Ellipsoid.distance2}.
    '''
    return _dE(cosineForsytheAndoyerLambert_, datum, wrap, lat1, lon1, lat2, lon2)


def cosineForsytheAndoyerLambert_(phi2, phi1, lam21, datum=_WGS84):
    '''Compute the I{angular} distance between two (ellipsoidal) points using the
       U{Forsythe-Andoyer-Lambert correction<https://www2.UNB.Ca/gge/Pubs/TR77.pdf>} of
       the U{Law of Cosines<https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
       formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg datum: Datum (L{Datum}) or ellipsoid to use (L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}).

       @return: Angular distance (C{radians}).

       @raise TypeError: Invalid B{C{datum}}.

       @see: Functions L{cosineForsytheAndoyerLambert}, L{cosineAndoyerLambert_},
             L{cosineLaw_}, L{equirectangular_}, L{euclidean_}, L{flatLocal_}/L{hubeny_},
             L{flatPolar_}, L{haversine_}, L{thomas_} and L{vincentys_} and U{Geodesy-PHP
             <https://GitHub.com/jtejido/geodesy-php/blob/master/src/Geodesy/
             Distance/ForsytheCorrection.php>}.
    '''
    s2, c2, s1, c1, r, _ = _sincosa6(phi2, phi1, lam21)
    if r and isnon0(c1) and isnon0(c2):
        E = _ellipsoidal(datum, cosineForsytheAndoyerLambert_)
        if E.f:  # ellipsoidal
            sr, cr, s2r, _ = sincos2_(r, r * 2)
            if isnon0(sr) and fabs(cr) < EPS1:
                s = (s1 + s2)**2 / (1 + cr)
                t = (s1 - s2)**2 / (1 - cr)
                x = s + t
                y = s - t

                s =  8 * r**2 / sr
                a = 64 * r  + s * cr * 2  # 16 * r**2 / tan(r)
                d = 48 * sr + s  # 8 * r**2 / tan(r)
                b = -2 * d
                e = 30 * s2r
                c = fsumf_(30 * r, e * _0_5, s * cr)  # 8 * r**2 / tan(r)
                t = fsumf_( a * x, e * y**2, b * y, -c * x**2, d * x * y)

                r += fsumf_(-r * x, 3 * y * sr, t * E.f / _32_0) * E.f * _0_25
    return r


def cosineLaw(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two points using the U{spherical Law of
       Cosines<https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
       formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum})
                      or ellipsoid (L{Ellipsoid}, L{Ellipsoid2} or
                      L{a_f2Tuple}) to use.
       @kwarg wrap: If C{True}, wrap or I{normalize} and B{C{lat2}}
                    and B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}} or the
                ellipsoid or datum axes).

       @raise TypeError: Invalid B{C{radius}}.

       @see: Functions L{cosineLaw_}, L{cosineAndoyerLambert},
             L{cosineForsytheAndoyerLambert}, L{equirectangular}, L{euclidean},
             L{flatLocal}/L{hubeny}, L{flatPolar}, L{haversine}, L{thomas} and
             L{vincentys} and method L{Ellipsoid.distance2}.

       @note: See note at function L{vincentys_}.
    '''
    return _dS(cosineLaw_, radius, wrap, lat1, lon1, lat2, lon2)


def cosineLaw_(phi2, phi1, lam21):
    '''Compute the I{angular} distance between two points using the U{spherical
       Law of Cosines<https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
       formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Angular distance (C{radians}).

       @see: Functions L{cosineLaw}, L{cosineAndoyerLambert_},
             L{cosineForsytheAndoyerLambert_}, L{equirectangular_},
             L{euclidean_}, L{flatLocal_}/L{hubeny_}, L{flatPolar_},
             L{haversine_}, L{thomas_} and L{vincentys_}.

       @note: See note at function L{vincentys_}.
    '''
    return _sincosa6(phi2, phi1, lam21)[4]


def _d3(wrap, lat1, lon1, lat2, lon2):
    '''(INTERNAL) Helper for _dE, _dS and _eA.
    '''
    if wrap:
        d_lon, lat2, _ = _Wrap.latlon3(lon1, lat2, lon2, wrap)
        return radians(lat2), Phi_(lat1=lat1), radians(d_lon)
    else:  # for backward compaibility
        return Phi_(lat2=lat2), Phi_(lat1=lat1), Phi_(d_lon=lon2 - lon1)


def _dE(func_, earth, *wrap_lls):
    '''(INTERNAL) Helper for ellipsoidal distances.
    '''
    E = _ellipsoidal(earth, func_)
    r =  func_(*_d3(*wrap_lls), datum=E)
    return r * E.a


def _dS(func_, radius, *wrap_lls, **adjust):
    '''(INTERNAL) Helper for spherical distances.
    '''
    r = func_(*_d3(*wrap_lls), **adjust)
    if radius is not R_M:
        _, lat1, _, lat2, _ = wrap_lls
        radius = _mean_radius(radius, lat1, lat2)
    return r * radius


def _eA(excess_, radius, *wrap_lls):
    '''(INTERNAL) Helper for spherical excess or area.
    '''
    r = excess_(*_d3(*wrap_lls))
    if radius:
        _, lat1, _, lat2, _ = wrap_lls
        r *= _mean_radius(radius, lat1, lat2)**2
    return r


def _ellipsoidal(earth, where):
    '''(INTERNAL) Helper for distances.
    '''
    return _EWGS84 if earth in (_WGS84, _EWGS84)   else (
            earth  if isinstance(earth, Ellipsoid) else
           (earth  if isinstance(earth, Datum)     else  # PYCHOK indent
           _ellipsoidal_datum(earth, name=where.__name__)).ellipsoid)


def equirectangular(lat1, lon1, lat2, lon2, radius=R_M, **adjust_limit_wrap):
    '''Compute the distance between two points using
       the U{Equirectangular Approximation / Projection
       <https://www.Movable-Type.co.UK/scripts/latlong.html#equirectangular>}.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum})
                      or ellipsoid (L{Ellipsoid}, L{Ellipsoid2} or
                      L{a_f2Tuple}).
       @kwarg adjust_limit_wrap: Optional keyword arguments for
                     function L{equirectangular_}.

       @return: Distance (C{meter}, same units as B{C{radius}} or
                the ellipsoid or datum axes).

       @raise TypeError: Invalid B{C{radius}}.

       @see: Function L{equirectangular_} for more details, the
             available B{C{options}}, errors, restrictions and other,
             approximate or accurate distance functions.
    '''
    d = sqrt(equirectangular_(Lat(lat1=lat1), Lon(lon1=lon1),
                              Lat(lat2=lat2), Lon(lon2=lon2),
                            **adjust_limit_wrap).distance2)  # PYCHOK 4 vs 2-3
    return degrees2m(d, radius=_mean_radius(radius, lat1, lat2))


def _equirectangular(lat1, lon1, lat2, lon2, **adjust_limit_wrap):
    '''(INTERNAL) Helper for the L{frechet._FrecherMeterRadians}
       and L{hausdorff._HausdorffMeterRedians} classes.
    '''
    return equirectangular_(lat1, lon1, lat2, lon2, **adjust_limit_wrap).distance2 * _D2_R2


def equirectangular_(lat1, lon1, lat2, lon2, adjust=True, limit=45, wrap=False):
    '''Compute the distance between two points using the U{Equirectangular
       Approximation / Projection
       <https://www.Movable-Type.co.UK/scripts/latlong.html#equirectangular>}.

       This approximation is valid for short distance of several hundred Km
       or Miles, see the B{C{limit}} keyword argument and L{LimitError}.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg adjust: Adjust the wrapped, unrolled longitudinal delta
                      by the cosine of the mean latitude (C{bool}).
       @kwarg limit: Optional limit for lat- and longitudinal deltas
                     (C{degrees}) or C{None} or C{0} for unlimited.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}}
                    and B{C{lon2}} (C{bool}).

       @return: A L{Distance4Tuple}C{(distance2, delta_lat, delta_lon,
                unroll_lon2)} in C{degrees squared}.

       @raise LimitError: If the lat- and/or longitudinal delta exceeds the
                          B{C{-limit..limit}} range and L{pygeodesy.limiterrors}
                          set to C{True}.

       @see: U{Local, flat earth approximation
             <https://www.EdWilliams.org/avform.htm#flat>}, functions
             L{equirectangular}, L{cosineAndoyerLambert},
             L{cosineForsytheAndoyerLambert}, L{cosineLaw}, L{euclidean},
             L{flatLocal}/L{hubeny}, L{flatPolar}, L{haversine}, L{thomas}
             and L{vincentys} and methods L{Ellipsoid.distance2},
             C{LatLon.distanceTo*} and C{LatLon.equirectangularTo}.
    '''
    d_lon,  lat2, ulon2 = _Wrap.latlon3(lon1, lat2, lon2, wrap)
    d_lat = lat2 - lat1

    if limit and limit > 0 and limiterrors():
        d = max(fabs(d_lat), fabs(d_lon))
        if d > limit:
            t = _SPACE_(_delta_, Fmt.PAREN_g(d), Fmt.exceeds_limit(limit))
            s =  unstr(equirectangular_, lat1, lon1, lat2, lon2,
                                         limit=limit, wrap=wrap)
            raise LimitError(s, txt=t)

    if adjust:  # scale delta lon
        d_lon *= _scale_deg(lat1, lat2)

    d2 = hypot2(d_lat, d_lon)  # degrees squared!
    return Distance4Tuple(d2, d_lat, d_lon, ulon2 - lon2)


def euclidean(lat1, lon1, lat2, lon2, radius=R_M, adjust=True, wrap=False):
    '''Approximate the C{Euclidean} distance between two (spherical) points.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum})
                      or ellipsoid (L{Ellipsoid}, L{Ellipsoid2} or
                      L{a_f2Tuple}) to use.
       @kwarg adjust: Adjust the longitudinal delta by the cosine of
                      the mean latitude (C{bool}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}}
                    and B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}} or the
                ellipsoid or datum axes).

       @raise TypeError: Invalid B{C{radius}}.

       @see: U{Distance between two (spherical) points
             <https://www.EdWilliams.org/avform.htm#Dist>}, functions L{euclid},
             L{euclidean_}, L{cosineAndoyerLambert}, L{cosineForsytheAndoyerLambert},
             L{cosineLaw}, L{equirectangular}, L{flatLocal}/L{hubeny}, L{flatPolar},
             L{haversine}, L{thomas} and L{vincentys} and methods L{Ellipsoid.distance2},
             C{LatLon.distanceTo*} and C{LatLon.equirectangularTo}.
    '''
    return _dS(euclidean_, radius, wrap, lat1, lon1, lat2, lon2, adjust=adjust)


def euclidean_(phi2, phi1, lam21, adjust=True):
    '''Approximate the I{angular} C{Euclidean} distance between two
       (spherical) points.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg adjust: Adjust the longitudinal delta by the cosine
                      of the mean latitude (C{bool}).

       @return: Angular distance (C{radians}).

       @see: Functions L{euclid}, L{euclidean}, L{cosineAndoyerLambert_},
             L{cosineForsytheAndoyerLambert_}, L{cosineLaw_}, L{equirectangular_},
             L{flatLocal_}/L{hubeny_}, L{flatPolar_}, L{haversine_}, L{thomas_}
             and L{vincentys_}.
    '''
    if adjust:
        lam21 *= _scale_rad(phi2, phi1)
    return euclid(phi2 - phi1, lam21)


def excessAbc_(A, b, c):
    '''Compute the I{spherical excess} C{E} of a (spherical) triangle
       from two sides and the included (small) angle.

       @arg A: An interior triangle angle (C{radians}).
       @arg b: Frist adjacent triangle side (C{radians}).
       @arg c: Second adjacent triangle side (C{radians}).

       @return: Spherical excess (C{radians}).

       @raise UnitError: Invalid B{C{A}}, B{C{b}} or B{C{c}}.

       @see: Functions L{excessGirard_}, L{excessLHuilier_} and U{Spherical
             trigonometry<https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    A = Radians_(A=A)
    b = Radians_(b=b) * _0_5
    c = Radians_(c=c) * _0_5

    sA, cA, sb, cb, sc, cc = sincos2_(A, b, c)
    return atan2(sA * sb * sc, cb * cc + cA * sb * sc) * _2_0


def excessCagnoli_(a, b, c):
    '''Compute the I{spherical excess} C{E} of a (spherical) triangle using
       U{Cagnoli's<https://Zenodo.org/record/35392>} (D.34) formula.

       @arg a: First triangle side (C{radians}).
       @arg b: Second triangle side (C{radians}).
       @arg c: Third triangle side (C{radians}).

       @return: Spherical excess (C{radians}).

       @raise UnitError: Invalid B{C{a}}, B{C{b}} or B{C{c}}.

       @see: Function L{excessLHuilier_} and U{Spherical trigonometry
             <https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    a = Radians_(a=a)
    b = Radians_(b=b)
    c = Radians_(c=c)

    s = fsumf_(a, b, c) * _0_5
    r = sin(s) * sin(s - a) * sin(s - b) * sin(s - c)
    c = cos(a * _0_5) * cos(b * _0_5) * cos(c * _0_5)
    r = asin(sqrt(r) * _0_5 / c) if c and r > 0 else _0_0
    return Radians(Cagnoli=r * _2_0)


def excessGirard_(A, B, C):
    '''Compute the I{spherical excess} C{E} of a (spherical) triangle using
       U{Girard's<https://MathWorld.Wolfram.com/GirardsSphericalExcessFormula.html>}
       formula.

       @arg A: First interior triangle angle (C{radians}).
       @arg B: Second interior triangle angle (C{radians}).
       @arg C: Third interior triangle angle (C{radians}).

       @return: Spherical excess (C{radians}).

       @raise UnitError: Invalid B{C{A}}, B{C{B}} or B{C{C}}.

       @see: Function L{excessLHuilier_} and U{Spherical trigonometry
             <https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    return Radians(Girard=fsumf_(Radians_(A=A),
                                 Radians_(B=B),
                                 Radians_(C=C), -PI))


def excessLHuilier_(a, b, c):
    '''Compute the I{spherical excess} C{E} of a (spherical) triangle using
       U{L'Huilier's<https://MathWorld.Wolfram.com/LHuiliersTheorem.html>}
       Theorem.

       @arg a: First triangle side (C{radians}).
       @arg b: Second triangle side (C{radians}).
       @arg c: Third triangle side (C{radians}).

       @return: Spherical excess (C{radians}).

       @raise UnitError: Invalid B{C{a}}, B{C{b}} or B{C{c}}.

       @see: Function L{excessCagnoli_}, L{excessGirard_} and U{Spherical
             trigonometry<https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    a = Radians_(a=a)
    b = Radians_(b=b)
    c = Radians_(c=c)

    s = fsumf_(a, b, c) * _0_5
    r = tan_2(s) * tan_2(s - a) * tan_2(s - b) * tan_2(s - c)
    r = atan(sqrt(r)) if r > 0 else _0_0
    return Radians(LHuilier=r * _4_0)


def excessKarney(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the surface area of a (spherical) quadrilateral bounded by a
       segment of a great circle, two meridians and the equator using U{Karney's
       <https://MathOverflow.net/questions/97711/the-area-of-spherical-polygons>}
       method.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum})
                      or ellipsoid (L{Ellipsoid}, L{Ellipsoid2} or
                      L{a_f2Tuple}) or C{None}.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                    B{C{lat2}} and B{C{lon2}} (C{bool}).

       @return: Surface area, I{signed} (I{square} C{meter} or the same units as
                B{C{radius}} I{squared}) or the I{spherical excess} (C{radians})
                if C{B{radius}=0} or C{None}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise UnitError: Invalid B{C{lat2}} or B{C{lat1}}.

       @raise ValueError: Semi-circular longitudinal delta.

       @see: Functions L{excessKarney_} and L{excessQuad}.
    '''
    return _eA(excessKarney_, radius, wrap, lat1, lon1, lat2, lon2)


def excessKarney_(phi2, phi1, lam21):
    '''Compute the I{spherical excess} C{E} of a (spherical) quadrilateral bounded
       by a segment of a great circle, two meridians and the equator using U{Karney's
       <https://MathOverflow.net/questions/97711/the-area-of-spherical-polygons>}
       method.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Spherical excess, I{signed} (C{radians}).

       @raise ValueError: Semi-circular longitudinal delta B{C{lam21}}.

       @see: Function L{excessKarney} and U{Area of a spherical polygon
       <https://MathOverflow.net/questions/97711/the-area-of-spherical-polygons>}.
    '''
    # from: Veness <https://www.Movable-Type.co.UK/scripts/latlong.html>  Area
    # method due to Karney: for each edge of the polygon,
    #
    #                 tan(Δλ / 2) · (tan(φ1 / 2) + tan(φ2 / 2))
    #    tan(E / 2) = -----------------------------------------
    #                           1 +  tan(φ1 / 2) · tan(φ2 / 2)
    #
    # where E is the spherical excess of the trapezium obtained by extending
    # the edge to the equator-circle vector for each edge (see also ***).
    t2 = tan_2(phi2)
    t1 = tan_2(phi1)
    t  = tan_2(lam21, lam21=None)
    return Radians(Karney=atan2(t * (t1 + t2),
                             _1_0 + (t1 * t2)) * _2_0)


# ***) Original post no longer available, following is a copy of the main part
# <http://OSGeo-org.1560.x6.Nabble.com/Area-of-a-spherical-polygon-td3841625.html>
#
# The area of a polygon on a (unit) sphere is given by the spherical excess
#
#    A = 2 * pi - sum(exterior angles)
#
# However this is badly conditioned if the polygon is small.  In this case, use
#
#    A = sum(S12{i, i+1}) over the edges of the polygon
#
# where S12 is the area of the quadrilateral bounded by an edge of the polygon,
# two meridians and the equator, i.e. with vertices (phi1, lambda1), (phi2,
# lambda2), (0, lambda1) and (0, lambda2).  S12 is given by
#
#    tan(S12 / 2) = tan(lambda21 / 2) * (tan(phi1 / 2) + tan(phi2 / 2)) /
#                                       (tan(phi1 / 2) * tan(phi2 / 2) + 1)
#
#                 = tan(lambda21 / 2) * tanh((Lamb(phi1) + Lamb(phi2)) / 2)
#
# where lambda21 = lambda2 - lambda1 and Lamb(x) is the Lambertian (or the
# inverse Gudermannian) function
#
#    Lambertian(x) = asinh(tan(x)) = atanh(sin(x)) = 2 * atanh(tan(x / 2))
#
# Notes: The formula for S12 is exact, except that...
#      - it is indeterminate if an edge is a semi-circle
#      - the formula for A applies only if the polygon does not include a pole
#        (if it does, then add +/- 2 * pi to the result)
#      - in the limit of small phi and lambda, S12 reduces to the trapezoidal
#        formula, S12 = (lambda2 - lambda1) * (phi1 + phi2) / 2
#      - I derived this result from the equation for the area of a spherical
#        triangle in terms of two edges and the included angle given by, e.g.
#        U{Todhunter, I. - Spherical Trigonometry (1871), Sec. 103, Eq. (2)
#        <http://Books.Google.com/books?id=3uBHAAAAIAAJ&pg=PA71>}
#      - I would be interested to know if this formula for S12 is already known
#      - Charles Karney


def excessQuad(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the surface area of a (spherical) quadrilateral bounded by a segment
       of a great circle, two meridians and the equator.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum})
                      or ellipsoid (L{Ellipsoid}, L{Ellipsoid2} or
                      L{a_f2Tuple}) or C{None}.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                    B{C{lat2}} and B{C{lon2}} (C{bool}).

       @return: Surface area, I{signed} (I{square} C{meter} or the same units as
                B{C{radius}} I{squared}) or the I{spherical excess} (C{radians})
                if C{B{radius}=0} or C{None}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise UnitError: Invalid B{C{lat2}} or B{C{lat1}}.

       @see: Function L{excessQuad_} and L{excessKarney}.
    '''
    return _eA(excessQuad_, radius, wrap, lat1, lon1, lat2, lon2)


def excessQuad_(phi2, phi1, lam21):
    '''Compute the I{spherical excess} C{E} of a (spherical) quadrilateral bounded
       by a segment of a great circle, two meridians and the equator.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Spherical excess, I{signed} (C{radians}).

       @see: Function L{excessQuad} and U{Spherical trigonometry
             <https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    s = sin((phi2 + phi1) * _0_5)
    c = cos((phi2 - phi1) * _0_5)
    return Radians(Quad=atan2(tan_2(lam21) * s, c) * _2_0)


def flatLocal(lat1, lon1, lat2, lon2, datum=_WGS84, scaled=True, wrap=False):
    '''Compute the distance between two (ellipsoidal) points using
       the U{ellipsoidal Earth to plane projection<https://WikiPedia.org/
       wiki/Geographical_distance#Ellipsoidal_Earth_projected_to_a_plane>}
       aka U{Hubeny<https://www.OVG.AT/de/vgi/files/pdf/3781/>} formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg datum: Datum (L{Datum}) or ellipsoid (L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}) to use.
       @kwarg scaled: Scale prime_vertical by C{cos(B{phi})} (C{bool}),
                      see method L{pygeodesy.Ellipsoid.roc2_}.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                    B{C{lat2}} and B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as the B{C{datum}}'s
                ellipsoid axes).

       @raise TypeError: Invalid B{C{datum}}.

       @note: The meridional and prime_vertical radii of curvature
              are taken and scaled at the mean of both latitude.

       @see: Functions L{flatLocal_} or L{hubeny_}, L{cosineLaw}, L{flatPolar},
             L{cosineAndoyerLambert}, L{cosineForsytheAndoyerLambert},
             L{equirectangular}, L{euclidean}, L{haversine}, L{thomas},
             L{vincentys}, method L{Ellipsoid.distance2} and U{local, flat
             earth approximation<https://www.EdWilliams.org/avform.htm#flat>}.
    '''
    E = _ellipsoidal(datum, flatLocal)
    return E._hubeny_2(*_d3(wrap, lat1, lon1, lat2, lon2),
                       scaled=scaled, squared=False) * E.a

hubeny = flatLocal  # PYCHOK for Karl Hubeny


def flatLocal_(phi2, phi1, lam21, datum=_WGS84, scaled=True):
    '''Compute the I{angular} distance between two (ellipsoidal) points using
       the U{ellipsoidal Earth to plane projection<https://WikiPedia.org/
       wiki/Geographical_distance#Ellipsoidal_Earth_projected_to_a_plane>}
       aka U{Hubeny<https://www.OVG.AT/de/vgi/files/pdf/3781/>} formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg datum: Datum (L{Datum}) or ellipsoid (L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}) to use.
       @kwarg scaled: Scale prime_vertical by C{cos(B{phi})} (C{bool}),
                      see method L{pygeodesy.Ellipsoid.roc2_}.

       @return: Angular distance (C{radians}).

       @raise TypeError: Invalid B{C{datum}}.

       @note: The meridional and prime_vertical radii of curvature
              are taken and scaled I{at the mean of both latitude}.

       @see: Functions L{flatLocal} or L{hubeny}, L{cosineAndoyerLambert_},
             L{cosineForsytheAndoyerLambert_}, L{cosineLaw_}, L{flatPolar_},
             L{equirectangular_}, L{euclidean_}, L{haversine_}, L{thomas_}
             and L{vincentys_} and U{local, flat earth approximation
             <https://www.EdWilliams.org/avform.htm#flat>}.
    '''
    E = _ellipsoidal(datum, flatLocal_)
    return E._hubeny_2(phi2, phi1, lam21, scaled=scaled, squared=False)

hubeny_ = flatLocal_  # PYCHOK for Karl Hubeny


def flatPolar(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two (spherical) points using
       the U{polar coordinate flat-Earth <https://WikiPedia.org/wiki/
       Geographical_distance#Polar_coordinate_flat-Earth_formula>}
       formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum})
                      or ellipsoid (L{Ellipsoid}, L{Ellipsoid2} or
                      L{a_f2Tuple}) to use.
       @kwarg wrap: If C{True}, wrap or I{normalize} and B{C{lat2}}
                    and B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}} or the
                ellipsoid or datum axes).

       @raise TypeError: Invalid B{C{radius}}.

       @see: Functions L{flatPolar_}, L{cosineAndoyerLambert},
             L{cosineForsytheAndoyerLambert},L{cosineLaw},
             L{flatLocal}/L{hubeny}, L{equirectangular},
             L{euclidean}, L{haversine}, L{thomas} and
             L{vincentys}.
    '''
    return _dS(flatPolar_, radius, wrap, lat1, lon1, lat2, lon2)


def flatPolar_(phi2, phi1, lam21):
    '''Compute the I{angular} distance between two (spherical) points
       using the U{polar coordinate flat-Earth<https://WikiPedia.org/wiki/
       Geographical_distance#Polar_coordinate_flat-Earth_formula>}
       formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Angular distance (C{radians}).

       @see: Functions L{flatPolar}, L{cosineAndoyerLambert_},
             L{cosineForsytheAndoyerLambert_}, L{cosineLaw_},
             L{equirectangular_}, L{euclidean_}, L{flatLocal_}/L{hubeny_},
             L{haversine_}, L{thomas_} and L{vincentys_}.
    '''
    a = fabs(PI_2 - phi1)  # co-latitude
    b = fabs(PI_2 - phi2)  # co-latitude
    if a < b:
        a, b = b, a
    if a < EPS0:
        a = _0_0
    elif b > 0:
        b  = b / a  # /= chokes PyChecker
        c  = b * cos(lam21) * _2_0
        c  = fsumf_(_1_0, b**2, -fabs(c))
        a *= sqrt0(c)
    return a


def hartzell(pov, los=None, earth=_WGS84, name=NN, **LatLon_and_kwds):
    '''Compute the intersection of the earth's surface and a Line-Of-Sight
       from a Point-Of-View in space.

       @arg pov: Point-Of-View outside the earth (C{Cartesian}, L{Ecef9Tuple}
                 or L{Vector3d}).
       @kwarg los: Line-Of-Sight, I{direction} to earth (L{Vector3d}) or
                   C{None} to point to the earth' center.
       @kwarg earth: The earth model (L{Datum}, L{Ellipsoid}, L{Ellipsoid2},
                     L{a_f2Tuple} or C{scalar} radius in C{meter}).
       @kwarg name: Optional name (C{str}).
       @kwarg LatLon_and_kwds: Optional C{LatLon} class for the intersection
                               point plus C{LatLon} keyword arguments, include
                               B{C{datum}} if different from B{C{earth}}.

       @return: The intersection point (L{Vector3d}, C{Cartesian type} of
                B{C{pov}} or B{C{LatLon}}).

       @raise IntersectionError: Null B{C{pov}} or B{C{los}} vector, B{C{pov}}
                                 is inside the earth or B{C{los}} points outside
                                 the earth or points in an opposite direction.

       @raise TypeError: Invalid B{C{pov}}, B{C{los}} or B{C{earth}}.

       @see: Function L{pygeodesy.hartzell4}, L{pygeodesy.tyr3d} for B{C{los}},
             method L{Ellipsoid.hartzell4} and U{I{Satellite Line-of-Sight
             Intersection with Earth}<https://StephenHartzell.Medium.com/
             satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6>}.
    '''
    D = earth if isinstance(earth, Datum) else \
           _spherical_datum(earth, name=hartzell.__name__)
    try:
        r, _ = _MODS.triaxials._hartzell3d2(pov, los, D.ellipsoid._triaxial)
    except Exception as x:
        raise IntersectionError(pov=pov, los=los, earth=earth, cause=x)

#   else:
#       E  = D.ellipsoid
#       # Triaxial(a, b, c) == (E.a, E.a, E.b)
#
#       def _Error(txt):
#           return IntersectionError(pov=pov, los=los, earth=earth, txt=txt)
#
#       a2 = b2 = E.a2  # earth' x, y, ...
#       c2 = E.b2  # ... z semi-axis squared
#       q2 = E.b2_a2  # == c2 / a2
#       bc = E.a * E.b  # == b * c
#
#       V3 = _MODS.vector3d._otherV3d
#       p3 =  V3(pov=pov)
#       u3 =  V3(los=los) if los else p3.negate()
#       u3 =  u3.unit()  # unit vector, opposing signs
#
#       x2, y2, z2 = p3.x2y2z2  # p3.times_(p3).xyz
#       ux, vy, wz = u3.times_(p3).xyz
#       u2, v2, w2 = u3.x2y2z2  # u3.times_(u3).xyz
#
#       t = c2, c2, b2
#       m = fdot(t, u2, v2, w2)  # a2 factored out
#       if m < EPS0:  # zero or near-null LOS vector
#           raise _Error(_near_(_null_))
#
#       # a2 and b2 factored out, b2 == a2 and b2 / a2 == 1
#       r = fsumf_(b2 * w2,  c2 * v2,      -v2 * z2,      vy * wz * 2,
#                  c2 * u2, -u2 * z2,      -w2 * x2,      ux * wz * 2,
#                 -w2 * y2, -u2 * y2 * q2, -v2 * x2 * q2, ux * vy * 2 * q2)
#       if r > 0:
#           r = sqrt(r) * bc  # == a * a * b * c / a2
#       elif r < 0:  # LOS pointing away from or missing the earth
#           raise _Error(_opposite_ if max(ux, vy, wz) > 0 else _outside_)
#
#       d = Fdot(t, ux, vy, wz).fadd_(r).fover(m)  # -r for antipode, a2 factored out
#       if d > 0:  # POV inside or LOS missing, outside the earth
#           s = fsumf_(_1_0, x2 / a2, y2 / b2, z2 / c2, _N_2_0)  # like _sideOf
#           raise _Error(_outside_ if s > 0 else _inside_)
#       elif fsumf_(x2, y2, z2) < d**2:  # d past earth center
#           raise _Error(_too_(_distant_))
#
#       r = p3.minus(u3.times(d))
# #     h = p3.minus(r).length  # distance to ellipsoid

    r = _xnamed(r, name or hartzell.__name__)
    if LatLon_and_kwds:
        c = _MODS.cartesianBase.CartesianBase(r, datum=D, name=r.name)
        r =  c.toLatLon(**LatLon_and_kwds)
    return r


def haversine(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two (spherical) points using the
       U{Haversine<https://www.Movable-Type.co.UK/scripts/latlong.html>}
       formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum})
                      or ellipsoid (L{Ellipsoid}, L{Ellipsoid2} or
                      L{a_f2Tuple}) to use.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                    B{C{lat2}} and B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @raise TypeError: Invalid B{C{radius}}.

       @see: U{Distance between two (spherical) points
             <https://www.EdWilliams.org/avform.htm#Dist>}, functions
             L{cosineLaw}, L{cosineAndoyerLambert}, L{cosineForsytheAndoyerLambert},
             L{equirectangular}, L{euclidean}, L{flatLocal}/L{hubeny}, L{flatPolar},
             L{thomas} and L{vincentys} and methods L{Ellipsoid.distance2},
             C{LatLon.distanceTo*} and C{LatLon.equirectangularTo}.

       @note: See note at function L{vincentys_}.
    '''
    return _dS(haversine_, radius, wrap, lat1, lon1, lat2, lon2)


def haversine_(phi2, phi1, lam21):
    '''Compute the I{angular} distance between two (spherical) points
       using the U{Haversine<https://www.Movable-Type.co.UK/scripts/latlong.html>}
       formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Angular distance (C{radians}).

       @see: Functions L{haversine}, L{cosineAndoyerLambert_},
             L{cosineForsytheAndoyerLambert_}, L{cosineLaw_},
             L{equirectangular_}, L{euclidean_}, L{flatLocal_}/L{hubeny_},
             L{flatPolar_}, L{thomas_} and L{vincentys_}.

       @note: See note at function L{vincentys_}.
    '''
    def _hsin(rad):
        return sin(rad * _0_5)**2

    h = _hsin(phi2 - phi1) + cos(phi1) * cos(phi2) * _hsin(lam21)  # haversine
    return atan2(sqrt0(h), sqrt0(_1_0 - h)) * _2_0  # == asin(sqrt(h)) * 2


def heightOf(angle, distance, radius=R_M):
    '''Determine the height above the (spherical) earth' surface after
       traveling along a straight line at a given tilt.

       @arg angle: Tilt angle above horizontal (C{degrees}).
       @arg distance: Distance along the line (C{meter} or same units as
                      B{C{radius}}).
       @kwarg radius: Optional mean earth radius (C{meter}).

       @return: Height (C{meter}, same units as B{C{distance}} and B{C{radius}}).

       @raise ValueError: Invalid B{C{angle}}, B{C{distance}} or B{C{radius}}.

       @see: U{MultiDop geog_lib.GeogBeamHt<https://GitHub.com/NASA/MultiDop>}
             (U{Shapiro et al. 2009, JTECH
             <https://Journals.AMetSoc.org/doi/abs/10.1175/2009JTECHA1256.1>}
             and U{Potvin et al. 2012, JTECH
             <https://Journals.AMetSoc.org/doi/abs/10.1175/JTECH-D-11-00019.1>}).
    '''
    r = h = Radius(radius)
    d = fabs(Distance(distance))
    if d > h:
        d, h = h, d

    if d > EPS0:  # and h > EPS0
        d = d / h  # /= h chokes PyChecker
        s = sin(Phi_(angle=angle, clip=_180_0))
        s = fsumf_(_1_0, _2_0 * s * d, d**2)
        if s > 0:
            return h * sqrt(s) - r

    raise _ValueError(angle=angle, distance=distance, radius=radius)


def heightOrthometric(h_ll, N):
    '''Get the I{orthometric} height B{H}, the height above the geoid, earth surface.

       @arg h_ll: The height above the ellipsoid (C{meter}) or an I{ellipsoidal}
                  location (C{LatLon} with a C{height} or C{h} attribute).
       @arg N: The I{geoid} height (C{meter}), the height of the geoid above the
               ellipsoid at the same B{C{h_ll}} location.

       @return: I{Orthometric} height C{B{H} = B{h} - B{N}} (C{meter}, same units
                as B{C{h}} and B{C{N}}).

       @see: U{Ellipsoid, Geoid, and Othometric Heights<https://www.NGS.NOAA.gov/
             GEOID/PRESENTATIONS/2007_02_24_CCPS/Roman_A_PLSC2007notes.pdf>}, page
             6 and module L{pygeodesy.geoids}.
    '''
    h = h_ll if isscalar(h_ll) else _xattr(h_ll, height=_xattr(h_ll, h=0))
    return Height(H=Height(h=h) - Height(N=N))


def horizon(height, radius=R_M, refraction=False):
    '''Determine the distance to the horizon from a given altitude
       above the (spherical) earth.

       @arg height: Altitude (C{meter} or same units as B{C{radius}}).
       @kwarg radius: Optional mean earth radius (C{meter}).
       @kwarg refraction: Consider atmospheric refraction (C{bool}).

       @return: Distance (C{meter}, same units as B{C{height}} and B{C{radius}}).

       @raise ValueError: Invalid B{C{height}} or B{C{radius}}.

       @see: U{Distance to horizon<https://www.EdWilliams.org/avform.htm#Horizon>}.
    '''
    h, r = Height(height), Radius(radius)
    if min(h, r) < 0:
        raise _ValueError(height=height, radius=radius)

    if refraction:
        d2 = 2.415750694528 * h * r  # 2.0 / 0.8279
    else:
        d2 = h * fsumf_(r, r, h)
    return sqrt0(d2)


class _idllmn6(object):  # see also .geodesicw._wargs, .vector2d._numpy
    '''(INTERNAL) Helper for C{intersection2} and C{intersections2}.
    '''
    @contextmanager  # <https://www.python.org/dev/peps/pep-0343/> Examples
    def __call__(self, datum, lat1, lon1, lat2, lon2, small, wrap, s, **kwds):
        try:
            if wrap:
                _, lat2, lon2 = _Wrap.latlon3(lon1, lat2, lon2, wrap)
                kwds = _xkwds(kwds, wrap=wrap)  # for _xError
            m = small if small is _100km else Meter_(small=small)
            n = (intersections2 if s else intersection2).__name__
            if datum is None or euclidean(lat1, lon1, lat2, lon2) < m:
                d, m = None, _MODS.vector3d
                _i   = m._intersects2 if s else m._intersect3d3
            elif isscalar(datum) and datum < 0 and not s:
                d = _spherical_datum(-datum, name=n)
                m = _MODS.sphericalNvector
                _i = m.intersection
            else:
                d = _spherical_datum(datum, name=n)
                if d.isSpherical:
                    m = _MODS.sphericalTrigonometry
                    _i = m._intersects2 if s else m._intersect
                elif d.isEllipsoidal:
                    try:
                        if d.ellipsoid.geodesic:
                            pass
                        m = _MODS.ellipsoidalKarney
                    except ImportError:
                        m = _MODS.ellipsoidalExact
                    _i = m._intersections2 if s else m._intersection3  # ellispoidalBaseDI
                else:
                    raise _TypeError(datum=datum)
            yield _i, d, lat2, lon2, m, n

        except (TypeError, ValueError) as x:
            raise _xError(x, lat1=lat1, lon1=lon1, datum=datum,
                             lat2=lat2, lon2=lon2, small=small, **kwds)

_idllmn6 = _idllmn6()  # PYCHOK singleton


def intersection2(lat1, lon1, bearing1,
                  lat2, lon2, bearing2, datum=None, wrap=False, small=_100km):  # was=True
    '''I{Conveniently} compute the intersection of two lines each defined
       by a (geodetic) point and a bearing from North, using either ...

       1) L{vector3d.intersection3d3} for B{C{small}} distances (below 100 Km
       or about 0.88 degrees) or if I{no} B{C{datum}} is specified, or ...

       2) L{sphericalTrigonometry.intersection} for a spherical B{C{datum}}
       or a C{scalar B{datum}} representing the earth radius, conventionally
       in C{meter} or ...

       3) L{sphericalNvector.intersection} if B{C{datum}} is a I{negative}
       C{scalar}, (negative) earth radius, conventionally in C{meter} or ...

       4) L{ellipsoidalKarney.intersection3} for an ellipsoidal B{C{datum}}
       and if I{Karney}'s U{geographiclib<https://PyPI.org/project/geographiclib>}
       is installed, otherwise ...

       5) L{ellipsoidalExact.intersection3}, provided B{C{datum}} is ellipsoidal.

       @arg lat1: Latitude of the first point (C{degrees}).
       @arg lon1: Longitude of the first point (C{degrees}).
       @arg bearing1: Bearing at the first point (compass C{degrees}).
       @arg lat2: Latitude of the second point (C{degrees}).
       @arg lon2: Longitude of the second point (C{degrees}).
       @arg bearing2: Bearing at the second point (compass C{degrees}).
       @kwarg datum: Optional datum (L{Datum}) or ellipsoid (L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}) or C{scalar} earth
                     radius (C{meter}, same units as B{C{radius1}} and
                     B{C{radius2}}) or C{None}.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}}
                    and B{C{lon2}} (C{bool}).
       @kwarg small: Upper limit for small distances (C{meter}).

       @return: A L{LatLon2Tuple}C{(lat, lon)} with the lat- and
                longitude of the intersection point.

       @raise IntersectionError: Ambiguous or infinite intersection
                                 or colinear, parallel or otherwise
                                 non-intersecting lines.

       @raise TypeError: Invalid B{C{datum}}.

       @raise UnitError: Invalid B{C{lat1}}, B{C{lon1}}, B{C{bearing1}},
                         B{C{lat2}}, B{C{lon2}} or B{C{bearing2}}.

       @see: Method L{RhumbLine.intersection2}.

       @note: The returned intersections may be near-antipodal.
    '''
    b1 = Bearing(bearing1=bearing1)
    b2 = Bearing(bearing2=bearing2)
    with _idllmn6(datum, lat1, lon1, lat2, lon2,
                  small, wrap, False, bearing1=b1, bearing2=b2) as t:
        _i, d, lat2, lon2, m, n = t
        if d is None:
            t, _, _ = _i(m.Vector3d(lon1, lat1, 0), b1,
                         m.Vector3d(lon2, lat2, 0), b2, useZ=False)
            t = LatLon2Tuple(t.y, t.x, name=n)

        else:
            t = _i(m.LatLon(lat1, lon1, datum=d), b1,
                   m.LatLon(lat2, lon2, datum=d), b2, height=0, wrap=False)
            if isinstance(t, Intersection3Tuple):  # ellipsoidal
                t, _, _ = t
            t = LatLon2Tuple(t.lat, t.lon, name=n)
    return t


def intersections2(lat1, lon1, radius1,
                   lat2, lon2, radius2, datum=None, wrap=False, small=_100km):  # was=True
    '''I{Conveniently} compute the intersections of two circles each defined
       by a (geodetic) center point and a radius, using either ...

       1) L{vector3d.intersections2} for B{C{small}} distances (below 100 Km
       or about 0.88 degrees) or if I{no} B{C{datum}} is specified, or ...

       2) L{sphericalTrigonometry.intersections2} for a spherical B{C{datum}}
       or a C{scalar B{datum}} representing the earth radius, conventionally
       in C{meter} or ...

       3) L{ellipsoidalKarney.intersections2} for an ellipsoidal B{C{datum}}
       and if I{Karney}'s U{geographiclib<https://PyPI.org/project/geographiclib>}
       is installed, otherwise ...

       4) L{ellipsoidalExact.intersections2}, provided B{C{datum}} is ellipsoidal.

       @arg lat1: Latitude of the first circle center (C{degrees}).
       @arg lon1: Longitude of the first circle center (C{degrees}).
       @arg radius1: Radius of the first circle (C{meter}, conventionally).
       @arg lat2: Latitude of the second circle center (C{degrees}).
       @arg lon2: Longitude of the second circle center (C{degrees}).
       @arg radius2: Radius of the second circle (C{meter}, same units as B{C{radius1}}).
       @kwarg datum: Optional datum (L{Datum}) or ellipsoid (L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}) or C{scalar} earth
                     radius (C{meter}, same units as B{C{radius1}} and
                     B{C{radius2}}) or C{None}.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{lat2}}
                    and B{C{lon2}} (C{bool}).
       @kwarg small: Upper limit for small distances (C{meter}).

       @return: 2-Tuple of the intersection points, each a
                L{LatLon2Tuple}C{(lat, lon)}.  For abutting circles, the
                points are the same instance, aka the I{radical center}.

       @raise IntersectionError: Concentric, antipodal, invalid or
                                 non-intersecting circles or no
                                 convergence.

       @raise TypeError: Invalid B{C{datum}}.

       @raise UnitError: Invalid B{C{lat1}}, B{C{lon1}}, B{C{radius1}},
                         B{C{lat2}}, B{C{lon2}} or B{C{radius2}}.
    '''
    r1 = Radius_(radius1=radius1)
    r2 = Radius_(radius2=radius2)
    with _idllmn6(datum, lat1, lon1, lat2, lon2,
                  small, wrap, True, radius1=r1, radius2=r2) as t:
        _i, d, lat2, lon2, m, n = t
        if d is None:
            r1 = m2degrees(r1, radius=R_M, lat=lat1)
            r2 = m2degrees(r2, radius=R_M, lat=lat2)

            def _V2T(x, y, _, **unused):  # _ == z unused
                return LatLon2Tuple(y, x, name=n)

            t = _i(m.Vector3d(lon1, lat1, 0), r1,
                   m.Vector3d(lon2, lat2, 0), r2, sphere=False,
                     Vector=_V2T)
        else:
            def _LL2T(lat, lon, **unused):
                return LatLon2Tuple(lat, lon, name=n)

            t = _i(m.LatLon(lat1, lon1, datum=d), r1,
                   m.LatLon(lat2, lon2, datum=d), r2,
                     LatLon=_LL2T, height=0, wrap=False)
    return t


def isantipode(lat1, lon1, lat2, lon2, eps=EPS):
    '''Check whether two points are I{antipodal}, on diametrically
       opposite sides of the earth.

       @arg lat1: Latitude of one point (C{degrees}).
       @arg lon1: Longitude of one point (C{degrees}).
       @arg lat2: Latitude of the other point (C{degrees}).
       @arg lon2: Longitude of the other point (C{degrees}).
       @kwarg eps: Tolerance for near-equality (C{degrees}).

       @return: C{True} if points are antipodal within the
                B{C{eps}} tolerance, C{False} otherwise.

       @see: Functions L{isantipode_} and L{antipode}.
    '''
    return (fabs(lat1 + lat2) <= eps and
            fabs(lon1 + lon2) <= eps) or _isequalTo(
            normal(lat1, lon1), antipode(lat2, lon2), eps)


def isantipode_(phi1, lam1, phi2, lam2, eps=EPS):
    '''Check whether two points are I{antipodal}, on diametrically
       opposite sides of the earth.

       @arg phi1: Latitude of one point (C{radians}).
       @arg lam1: Longitude of one point (C{radians}).
       @arg phi2: Latitude of the other point (C{radians}).
       @arg lam2: Longitude of the other point (C{radians}).
       @kwarg eps: Tolerance for near-equality (C{radians}).

       @return: C{True} if points are antipodal within the
                B{C{eps}} tolerance, C{False} otherwise.

       @see: Functions L{isantipode} and L{antipode_}.
    '''
    return (fabs(phi1 + phi2) <= eps and
            fabs(lam1 + lam2) <= eps) or _isequalTo_(
            normal_(phi1, lam1), antipode_(phi2, lam2), eps)


def _isequalTo(p1, p2, eps=EPS):
    '''Compare 2 point lat-/lons ignoring C{class}.
    '''
    return (fabs(p1.lat - p2.lat) <= eps and
            fabs(p1.lon - p2.lon) <= eps) if eps else (p1.latlon == p2.latlon)


def _isequalTo_(p1, p2, eps=EPS):
    '''(INTERNAL) Compare 2 point phi-/lams ignoring C{class}.
    '''
    return (fabs(p1.phi - p2.phi) <= eps and
            fabs(p1.lam - p2.lam) <= eps) if eps else (p1.philam == p2.philam)


def isnormal(lat, lon, eps=0):
    '''Check whether B{C{lat}} I{and} B{C{lon}} are within their
       respective I{normal} range in C{degrees}.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg eps: Optional tolerance C{degrees}).

       @return: C{True} if C{(abs(B{lat}) + B{eps}) <= 90} and
                C{(abs(B{lon}) + B{eps}) <= 180}, C{False} othwerwise.

       @see: Functions L{isnormal_} and L{normal}.
    '''
    return (_90_0 - fabs(lat)) >= eps and _loneg(fabs(lon)) >= eps


def isnormal_(phi, lam, eps=0):
    '''Check whether B{C{phi}} I{and} B{C{lam}} are within their
       respective I{normal} range in C{radians}.

       @arg phi: Latitude (C{radians}).
       @arg lam: Longitude (C{radians}).
       @kwarg eps: Optional tolerance C{radians}).

       @return: C{True} if C{(abs(B{phi}) + B{eps}) <= PI/2} and
                C{(abs(B{lam}) + B{eps}) <= PI}, C{False} othwerwise.

       @see: Functions L{isnormal} and L{normal_}.
    '''
    return (PI_2 - fabs(phi)) >= eps and (PI - fabs(lam)) >= eps


def latlon2n_xyz(lat, lon, name=NN):
    '''Convert lat-, longitude to C{n-vector} (I{normal} to the
       earth's surface) X, Y and Z components.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg name: Optional name (C{str}).

       @return: A L{Vector3Tuple}C{(x, y, z)}.

       @see: Function L{philam2n_xyz}.

       @note: These are C{n-vector} x, y and z components,
              I{NOT} geocentric ECEF x, y and z coordinates!
    '''
    return _2n_xyz(name, *sincos2d_(lat, lon))


def _normal2(a, b, n_2, n, n2):
    '''(INTERNAL) Helper for C{normal} and C{normal_}.
    '''
    if fabs(b) > n:
        b = remainder(b, n2)
    if fabs(a) > n_2:
        r = remainder(a, n)
        if r != a:
            a  = -r
            b -=  n if b > 0 else -n
    return float0_(a, b)


def normal(lat, lon, name=NN):
    '''Normalize a lat- I{and} longitude pair in C{degrees}.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg name: Optional name (C{str}).

       @return: L{LatLon2Tuple}C{(lat, lon)} with C{abs(lat) <= 90}
                and C{abs(lon) <= 180}.

       @see: Functions L{normal_} and L{isnormal}.
    '''
    return LatLon2Tuple(*_normal2(lat, lon, _90_0, _180_0, _360_0),
                          name=name or normal.__name__)


def normal_(phi, lam, name=NN):
    '''Normalize a lat- I{and} longitude pair in C{radians}.

       @arg phi: Latitude (C{radians}).
       @arg lam: Longitude (C{radians}).
       @kwarg name: Optional name (C{str}).

       @return: L{PhiLam2Tuple}C{(phi, lam)} with C{abs(phi) <= PI/2}
                and C{abs(lam) <= PI}.

       @see: Functions L{normal} and L{isnormal_}.
    '''
    return PhiLam2Tuple(*_normal2(phi, lam, PI_2, PI, PI2),
                          name=name or normal_.__name__)


def _2n_xyz(name, sa, ca, sb, cb):
    '''(INTERNAL) Helper for C{latlon2n_xyz} and C{philam2n_xyz}.
    '''
    # Kenneth Gade eqn 3, but using right-handed
    # vector x -> 0°E,0°N, y -> 90°E,0°N, z -> 90°N
    return Vector3Tuple(ca * cb, ca * sb, sa, name=name)


def n_xyz2latlon(x, y, z, name=NN):
    '''Convert C{n-vector} components to lat- and longitude in C{degrees}.

       @arg x: X component (C{scalar}).
       @arg y: Y component (C{scalar}).
       @arg z: Z component (C{scalar}).
       @kwarg name: Optional name (C{str}).

       @return: A L{LatLon2Tuple}C{(lat, lon)}.

       @see: Function L{n_xyz2philam}.
    '''
    return LatLon2Tuple(atan2d(z, hypot(x, y)), atan2d(y, x), name=name)


def n_xyz2philam(x, y, z, name=NN):
    '''Convert C{n-vector} components to lat- and longitude in C{radians}.

       @arg x: X component (C{scalar}).
       @arg y: Y component (C{scalar}).
       @arg z: Z component (C{scalar}).
       @kwarg name: Optional name (C{str}).

       @return: A L{PhiLam2Tuple}C{(phi, lam)}.

       @see: Function L{n_xyz2latlon}.
    '''
    return PhiLam2Tuple(atan2(z, hypot(x, y)), atan2(y, x), name=name)


def _opposes(d, m, n, n2):
    '''(INTERNAL) Helper for C{opposing} and C{opposing_}.
    '''
    d = d % n2  # -20 % 360 == 340, -1 % PI2 == PI2 - 1
    return False if d < m or d > (n2 - m) else (
           True if (n - m) < d < (n  + m) else None)


def opposing(bearing1, bearing2, margin=_90_0):
    '''Compare the direction of two bearings given in C{degrees}.

       @arg bearing1: First bearing (compass C{degrees}).
       @arg bearing2: Second bearing (compass C{degrees}).
       @kwarg margin: Optional, interior angle bracket (C{degrees}).

       @return: C{True} if both bearings point in opposite, C{False} if
                in similar or C{None} if in I{perpendicular} directions.

       @see: Function L{opposing_}.
    '''
    m = Degrees_(margin=margin, low=EPS0, high=_90_0)
    return _opposes(bearing2 - bearing1, m, _180_0, _360_0)


def opposing_(radians1, radians2, margin=PI_2):
    '''Compare the direction of two bearings given in C{radians}.

       @arg radians1: First bearing (C{radians}).
       @arg radians2: Second bearing (C{radians}).
       @kwarg margin: Optional, interior angle bracket (C{radians}).

       @return: C{True} if both bearings point in opposite, C{False} if
                in similar or C{None} if in perpendicular directions.

       @see: Function L{opposing}.
    '''
    m = Radians_(margin=margin, low=EPS0, high=PI_2)
    return _opposes(radians2 - radians1, m, PI, PI2)


def philam2n_xyz(phi, lam, name=NN):
    '''Convert lat-, longitude to C{n-vector} (I{normal} to the
       earth's surface) X, Y and Z components.

       @arg phi: Latitude (C{radians}).
       @arg lam: Longitude (C{radians}).
       @kwarg name: Optional name (C{str}).

       @return: A L{Vector3Tuple}C{(x, y, z)}.

       @see: Function L{latlon2n_xyz}.

       @note: These are C{n-vector} x, y and z components,
              I{NOT} geocentric ECEF x, y and z coordinates!
    '''
    return _2n_xyz(name, *sincos2_(phi, lam))


def _radical2(d, r1, r2):  # in .ellipsoidalBaseDI, .sphericalTrigonometry, .vector3d
    # (INTERNAL) See C{radical2} below
    # assert d > EPS0
    r = fsumf_(_1_0, (r1 / d)**2, -(r2 / d)**2) * _0_5
    return Radical2Tuple(max(_0_0, min(_1_0, r)), r * d)


def radical2(distance, radius1, radius2):
    '''Compute the I{radical ratio} and I{radical line} of two
       U{intersecting circles<https://MathWorld.Wolfram.com/
       Circle-CircleIntersection.html>}.

       The I{radical line} is perpendicular to the axis thru the
       centers of the circles at C{(0, 0)} and C{(B{distance}, 0)}.

       @arg distance: Distance between the circle centers (C{scalar}).
       @arg radius1: Radius of the first circle (C{scalar}).
       @arg radius2: Radius of the second circle (C{scalar}).

       @return: A L{Radical2Tuple}C{(ratio, xline)} where C{0.0 <=
                ratio <= 1.0} and C{xline} is along the B{C{distance}}.

       @raise IntersectionError: The B{C{distance}} exceeds the sum
                                 of B{C{radius1}} and B{C{radius2}}.

       @raise UnitError: Invalid B{C{distance}}, B{C{radius1}} or
                         B{C{radius2}}.

       @see: U{Circle-Circle Intersection
             <https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>}.
    '''
    d  = Distance_(distance, low=_0_0)
    r1 = Radius_(radius1=radius1)
    r2 = Radius_(radius2=radius2)
    if d > (r1 + r2):
        raise IntersectionError(distance=d, radius1=r1, radius2=r2,
                                            txt=_too_(_distant_))
    return _radical2(d, r1, r2) if d > EPS0 else \
            Radical2Tuple(_0_5, _0_0)


class Radical2Tuple(_NamedTuple):
    '''2-Tuple C{(ratio, xline)} of the I{radical} C{ratio} and
       I{radical} C{xline}, both C{scalar} and C{0.0 <= ratio <= 1.0}
    '''
    _Names_ = (_ratio_, _xline_)
    _Units_ = ( Scalar,  Scalar)


def _radistance(inst):
    '''(INTERNAL) Helper for the L{frechet._FrecherMeterRadians}
       and L{hausdorff._HausdorffMeterRedians} classes.
    '''
    kwds_ = _xkwds(inst._kwds, wrap=False)
    wrap_ = _xkwds_pop(kwds_, wrap=False)
    func_ =  inst._func_
    try:  # calling lower-overhead C{func_}
        func_(0, _0_25, _0_5, **kwds_)
        wrap_ = _Wrap._philamop(wrap_)
    except TypeError:
        return inst.distance

    def _philam(p):
        try:
            return p.phi, p.lam
        except AttributeError:  # no .phi or .lam
            return radians(p.lat), radians(p.lon)

    def _func_wrap(point1, point2):
        phi1, lam1 = wrap_(*_philam(point1))
        phi2, lam2 = wrap_(*_philam(point2))
        return func_(phi2, phi1, lam2 - lam1, **kwds_)

    inst._units = inst._units_
    return _func_wrap


def _scale_deg(lat1, lat2):  # degrees
    # scale factor cos(mean of lats) for delta lon
    m = fabs(lat1 + lat2) * _0_5
    return cos(radians(m)) if m < 90 else _0_0


def _scale_rad(phi1,  phi2):  # radians, by .frechet, .hausdorff, .heights
    # scale factor cos(mean of phis) for delta lam
    m = fabs(phi1 + phi2) * _0_5
    return cos(m) if m < PI_2 else _0_0


def _sincosa6(phi2, phi1, lam21):  # [4] in cosineLaw
    '''(INTERNAL) C{sin}es, C{cos}ines and C{acos}ine.
    '''
    s2, c2, s1, c1, _, c21 = sincos2_(phi2, phi1, lam21)
    return s2, c2, s1, c1, acos1(s1 * s2 + c1 * c2 * c21), c21


def thomas(lat1, lon1, lat2, lon2, datum=_WGS84, wrap=False):
    '''Compute the distance between two (ellipsoidal) points using
       U{Thomas'<https://apps.DTIC.mil/dtic/tr/fulltext/u2/703541.pdf>}
       formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg datum: Datum (L{Datum}) or ellipsoid (L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}) to use.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                    B{C{lat2}} and B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as the B{C{datum}}'s
                ellipsoid axes).

       @raise TypeError: Invalid B{C{datum}}.

       @see: Functions L{thomas_}, L{cosineAndoyerLambert}, L{cosineForsytheAndoyerLambert},
             L{cosineLaw}, L{equirectangular}, L{euclidean}, L{flatLocal}/L{hubeny},
             L{flatPolar}, L{haversine}, L{vincentys} and method L{Ellipsoid.distance2}.
    '''
    return _dE(thomas_, datum, wrap, lat1, lon1, lat2, lon2)


def thomas_(phi2, phi1, lam21, datum=_WGS84):
    '''Compute the I{angular} distance between two (ellipsoidal) points using
       U{Thomas'<https://apps.DTIC.mil/dtic/tr/fulltext/u2/703541.pdf>}
       formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg datum: Datum or ellipsoid to use (L{Datum}, L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}).

       @return: Angular distance (C{radians}).

       @raise TypeError: Invalid B{C{datum}}.

       @see: Functions L{thomas}, L{cosineAndoyerLambert_},
             L{cosineForsytheAndoyerLambert_}, L{cosineLaw_},
             L{equirectangular_}, L{euclidean_}, L{flatLocal_}/L{hubeny_},
             L{flatPolar_}, L{haversine_} and L{vincentys_} and U{Geodesy-PHP
             <https://GitHub.com/jtejido/geodesy-php/blob/master/src/Geodesy/
             Distance/ThomasFormula.php>}.
    '''
    s2, c2, s1, c1, r, _ = _sincosa6(phi2, phi1, lam21)
    if r and isnon0(c1) and isnon0(c2):
        E = _ellipsoidal(datum, thomas_)
        if E.f:
            r1 = atan2(E.b_a * s1, c1)
            r2 = atan2(E.b_a * s2, c2)

            j = (r2 + r1) * _0_5
            k = (r2 - r1) * _0_5
            sj, cj, sk, ck, h, _ = sincos2_(j, k, lam21 * _0_5)

            h =  fsumf_(sk**2, (ck * h)**2, -(sj * h)**2)
            u = _1_0 - h
            if isnon0(u) and isnon0(h):
                r = atan(sqrt0(h / u)) * 2  # == acos(1 - 2 * h)
                sr, cr = sincos2(r)
                if isnon0(sr):
                    u = 2 * (sj * ck)**2 / u
                    h = 2 * (sk * cj)**2 / h
                    x = u + h
                    y = u - h

                    s = r / sr
                    e = 4 * s**2
                    d = 2 * cr
                    a = e * d
                    b = 2 * r
                    c = s - (a - d) * _0_5
                    f = E.f * _0_25

                    t  = fsumf_(a * x, -b * y, c * x**2, -d * y**2, e * x * y)
                    r -= fsumf_(s * x, -y, -t * f * _0_25) * f * sr
    return r


def vincentys(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two (spherical) points using
       U{Vincenty's<https://WikiPedia.org/wiki/Great-circle_distance>}
       spherical formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}), datum (L{Datum})
                      or ellipsoid (L{Ellipsoid}, L{Ellipsoid2} or
                      L{a_f2Tuple}) to use.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                    B{C{lat2}} and B{C{lon2}} (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @raise UnitError: Invalid B{C{radius}}.

       @see: Functions L{vincentys_}, L{cosineAndoyerLambert},
             L{cosineForsytheAndoyerLambert},L{cosineLaw}, L{equirectangular},
             L{euclidean}, L{flatLocal}/L{hubeny}, L{flatPolar},
             L{haversine} and L{thomas} and methods L{Ellipsoid.distance2},
             C{LatLon.distanceTo*} and C{LatLon.equirectangularTo}.

       @note: See note at function L{vincentys_}.
    '''
    return _dS(vincentys_, radius, wrap, lat1, lon1, lat2, lon2)


def vincentys_(phi2, phi1, lam21):
    '''Compute the I{angular} distance between two (spherical) points using
       U{Vincenty's<https://WikiPedia.org/wiki/Great-circle_distance>}
       spherical formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Angular distance (C{radians}).

       @see: Functions L{vincentys}, L{cosineAndoyerLambert_},
             L{cosineForsytheAndoyerLambert_}, L{cosineLaw_},
             L{equirectangular_}, L{euclidean_}, L{flatLocal_}/L{hubeny_},
             L{flatPolar_}, L{haversine_} and L{thomas_}.

       @note: Functions L{vincentys_}, L{haversine_} and L{cosineLaw_}
              produce equivalent results, but L{vincentys_} is suitable
              for antipodal points and slightly more expensive (M{3 cos,
              3 sin, 1 hypot, 1 atan2, 6 mul, 2 add}) than L{haversine_}
              (M{2 cos, 2 sin, 2 sqrt, 1 atan2, 5 mul, 1 add}) and
              L{cosineLaw_} (M{3 cos, 3 sin, 1 acos, 3 mul, 1 add}).
    '''
    s1, c1, s2, c2, s21, c21 = sincos2_(phi1, phi2, lam21)

    c = c2 * c21
    x = s1 * s2 + c1 * c
    y = c1 * s2 - s1 * c
    return atan2(hypot(c2 * s21, y), x)

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
