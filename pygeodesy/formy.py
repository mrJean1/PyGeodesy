
# -*- coding: utf-8 -*-

u'''Formulary of basic geodesy functions and approximations.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isnon0 as _isnon0
from pygeodesy.datums import Datum, _ellipsoidal_datum, _mean_radius, \
                            _spherical_datum, _WGS84
from pygeodesy.ellipsoids import Ellipsoid
from pygeodesy.errors import _AssertionError, IntersectionError, \
                              LimitError, _limiterrors, _ValueError
from pygeodesy.fmath import euclid, fdot, fsum_, hypot, hypot2, sqrt0, unstr
from pygeodesy.interns import EPS, EPS0, EPS1, NN, PI, PI2, PI3, PI_2, R_M, \
                             _distant_, _inside_, _near_, _null_, _opposite_, \
                             _outside_, _too_, _0_0, _0_125, _0_25, _0_5, \
                             _1_0, _2_0, _4_0, _32_0, _90_0, _180_0, _360_0
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _NamedTuple, _xnamed
from pygeodesy.namedTuples import Bearing2Tuple, Distance4Tuple, \
                                  LatLon2Tuple, PhiLam2Tuple, Vector3Tuple
# from pygeodesy.streprs import unstr  # from .fmath
from pygeodesy.units import Degrees_, Distance, Distance_, Height, Lam_, Lat, \
                            Lon, Phi_, Radians, Radians_, Radius, Radius_, \
                            Scalar, _100km
from pygeodesy.utily import acos1, atan2b, degrees2m, degrees90, degrees180, \
                            m2degrees, sincos2, sincos2_, tan_2, unroll180, \
                            unrollPI, wrap90, wrap180, wrapPI, wrapPI_2

from math import atan, atan2, cos, degrees, radians, sin, sqrt  # pow

__all__ = _ALL_LAZY.formy
__version__ = '22.01.03'

_ratio_ = 'ratio'
_xline_ = 'xline'


def antipode(lat, lon):
    '''Return the antipode, the point diametrically opposite
       to a given point in C{degrees}.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).

       @return: A L{LatLon2Tuple}C{(lat, lon)}.

       @see: U{Geosphere<https://CRAN.R-Project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return LatLon2Tuple(-wrap90(lat), wrap180(lon + _180_0))


def antipode_(phi, lam):
    '''Return the antipode, the point diametrically opposite
       to a given point in C{radians}.

       @arg phi: Latitude (C{radians}).
       @arg lam: Longitude (C{radians}).

       @return: A L{PhiLam2Tuple}C{(phi, lam)}.

       @see: U{Geosphere<https://CRAN.R-Project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return PhiLam2Tuple(-wrapPI_2(phi), wrapPI(lam + PI))


def _area_or_(func_, lat1, lat2, radius, d_lon, unused):
    '''(INTERNAL) Helper for area and spherical excess.
    '''
    r = func_(Phi_(lat2=lat2),
              Phi_(lat1=lat1), radians(d_lon))
    if radius:
        r *= _mean_radius(radius, lat1, lat2)**2
    return r


def bearing(lat1, lon1, lat2, lon2, **options):
    '''Compute the initial or final bearing (forward or reverse
       azimuth) between a (spherical) start and end point.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg options: Optional keyword arguments for function
                       L{pygeodesy.bearing_}.

       @return: Initial or final bearing (compass C{degrees360}) or
                zero if start and end point coincide.
    '''
    r = bearing_(Phi_(lat1=lat1), Lam_(lon1=lon1),
                 Phi_(lat2=lat2), Lam_(lon2=lon2), **options)
    return degrees(r)


def bearing_(phi1, lam1, phi2, lam2, final=False, wrap=False):
    '''Compute the initial or final bearing (forward or reverse azimuth)
       between a (spherical) start and end point.

       @arg phi1: Start latitude (C{radians}).
       @arg lam1: Start longitude (C{radians}).
       @arg phi2: End latitude (C{radians}).
       @arg lam2: End longitude (C{radians}).
       @kwarg final: Return final bearing if C{True}, initial otherwise (C{bool}).
       @kwarg wrap: Wrap and L{pygeodesy.unrollPI} longitudes (C{bool}).

       @return: Initial or final bearing (compass C{radiansPI2}) or zero if start
                and end point coincide.
    '''
    if final:  # swap plus PI
        phi1, lam1, phi2, lam2 = phi2, lam2, phi1, lam1
        r = PI3
    else:
        r = PI2

    db, _ = unrollPI(lam1, lam2, wrap=wrap)
    sa1, ca1, sa2, ca2, sdb, cdb = sincos2_(phi1, phi2, db)

    # see <https://MathForum.org/library/drmath/view/55417.html>
    x = ca1 * sa2 - sa1 * ca2 * cdb
    y = sdb * ca2
    return (atan2(y, x) + r) % PI2  # wrapPI2


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
    '''Return the angle from North for the direction vector
       M{(lon2 - lon1, lat2 - lat1)} between two points.

       Suitable only for short, not near-polar vectors up to a few hundred
       Km or Miles.  Use function L{pygeodesy.bearing} for longer vectors.

       @arg lat1: From latitude (C{degrees}).
       @arg lon1: From longitude (C{degrees}).
       @arg lat2: To latitude (C{degrees}).
       @arg lon2: To longitude (C{degrees}).
       @kwarg adjust: Adjust the longitudinal delta by the cosine of the
                      mean latitude (C{bool}).
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

       @return: Compass angle from North (C{degrees360}).

       @note: Courtesy of Martin Schultz.

       @see: U{Local, flat earth approximation
             <https://www.EdWilliams.org/avform.htm#flat>}.
    '''
    d_lon, _ = unroll180(lon1, lon2, wrap=wrap)
    if adjust:  # scale delta lon
        d_lon *= _scale_deg(lat1, lat2)
    return atan2b(d_lon, lat2 - lat1)


def cosineAndoyerLambert(lat1, lon1, lat2, lon2, datum=_WGS84, wrap=False):
    '''Compute the distance between two (ellipsoidal) points using the
       U{Andoyer-Lambert correction<https://navlib.net/wp-content/uploads/
       2013/10/admiralty-manual-of-navigation-vol-1-1964-english501c.pdf>} of the
       U{Law of Cosines<https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
       fromula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg datum: Datum or ellipsoid to use (L{Datum}, L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}).
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as the B{C{datum}}'s
                ellipsoid axes or C{radians} if B{C{datum}} is C{None}).

       @raise TypeError: Invalid B{C{datum}}.

       @see: Functions L{cosineAndoyerLambert_}, L{cosineForsytheAndoyerLambert},
             L{cosineLaw}, L{equirectangular}, L{euclidean}, L{flatLocal}/L{hubeny},
             L{flatPolar}, L{haversine}, L{thomas} and L{vincentys} and method
             L{Ellipsoid.distance2}.
    '''
    return _distanceToE(cosineAndoyerLambert_, lat1, lat2, datum,
                                              *unroll180(lon1, lon2, wrap=wrap))


def cosineAndoyerLambert_(phi2, phi1, lam21, datum=_WGS84):
    '''Compute the I{angular} distance between two (ellipsoidal) points using the
       U{Andoyer-Lambert correction<https://navlib.net/wp-content/uploads/2013/10/
       admiralty-manual-of-navigation-vol-1-1964-english501c.pdf>} of the U{Law
       of Cosines<https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
       fromula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg datum: Datum or ellipsoid to use (L{Datum}, L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}).

       @return: Angular distance (C{radians}).

       @raise TypeError: Invalid B{C{datum}}.

       @see: Functions L{cosineAndoyerLambert}, L{cosineForsytheAndoyerLambert_},
             L{cosineLaw_}, L{equirectangular_}, L{euclidean_}, L{flatLocal_}/L{hubeny_},
             L{flatPolar_}, L{haversine_}, L{thomas_} and L{vincentys_} and U{Geodesy-PHP
             <https://GitHub.com/jtejido/geodesy-php/blob/master/src/Geodesy/
             Distance/AndoyerLambert.php>}.
    '''
    s2, c2, s1, c1, r, c21 = _sincosa6(phi2, phi1, lam21)
    if _isnon0(c1) and _isnon0(c2):
        E = _ellipsoidal(datum, cosineAndoyerLambert_)
        if E.f:  # ellipsoidal
            r2 = atan2(E.b_a * s2, c2)
            r1 = atan2(E.b_a * s1, c1)
            s2, c2, s1, c1 = sincos2_(r2, r1)
            r = acos1(s1 * s2 + c1 * c2 * c21)
            if r:
                sr, _, sr_2, cr_2 = sincos2_(r, r * _0_5)
                if _isnon0(sr_2) and _isnon0(cr_2):
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
       @kwarg datum: Datum or ellipsoid to use (L{Datum}, L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}).
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as the B{C{datum}}'s
                ellipsoid axes or C{radians} if B{C{datum}} is C{None}).

       @raise TypeError: Invalid B{C{datum}}.

       @see: Functions L{cosineForsytheAndoyerLambert_}, L{cosineAndoyerLambert},
             L{cosineLaw}, L{equirectangular}, L{euclidean}, L{flatLocal}/L{hubeny},
             L{flatPolar}, L{haversine}, L{thomas} and L{vincentys} and method
             L{Ellipsoid.distance2}.
    '''
    return _distanceToE(cosineForsytheAndoyerLambert_, lat1, lat2, datum,
                                                      *unroll180(lon1, lon2, wrap=wrap))


def cosineForsytheAndoyerLambert_(phi2, phi1, lam21, datum=_WGS84):
    '''Compute the I{angular} distance between two (ellipsoidal) points using the
       U{Forsythe-Andoyer-Lambert correction<https://www2.UNB.Ca/gge/Pubs/TR77.pdf>} of
       the U{Law of Cosines<https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
       formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg datum: Datum or ellipsoid to use (L{Datum}, L{Ellipsoid},
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
    if r and _isnon0(c1) and _isnon0(c2):
        E = _ellipsoidal(datum, cosineForsytheAndoyerLambert_)
        if E.f:  # ellipsoidal
            sr, cr, s2r, _ = sincos2_(r, r * _2_0)
            if _isnon0(sr) and abs(cr) < EPS1:
                s = (s1 + s2)**2 / (1 + cr)
                t = (s1 - s2)**2 / (1 - cr)
                x = s + t
                y = s - t

                s =  8 * r**2 / sr
                a = 64 * r + _2_0 * s * cr  # 16 * r**2 / tan(r)
                d = 48 * sr + s  # 8 * r**2 / tan(r)
                b = -2 * d
                e = 30 * s2r
                c = fsum_(30 * r, e * _0_5, s * cr)  # 8 * r**2 / tan(r)

                t  = fsum_( a * x, b * y, -c * x**2, d * x * y, e * y**2)
                r += fsum_(-r * x, 3 * y * sr, t * E.f / _32_0) * E.f * _0_25
    return r


def cosineLaw(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two points using the
       U{spherical Law of Cosines
       <https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
       formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}).
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}} or the
                ellipsoid or datum axes).

       @raise TypeError: Invalid B{C{radius}}.

       @see: Functions L{cosineLaw_}, L{cosineAndoyerLambert},
             L{cosineForsytheAndoyerLambert}, L{equirectangular}, L{euclidean},
             L{flatLocal}/L{hubeny}, L{flatPolar}, L{haversine}, L{thomas} and
             L{vincentys} and method L{Ellipsoid.distance2}.

       @note: See note at function L{vincentys_}.
    '''
    return _distanceToS(cosineLaw_, lat1, lat2, radius,
                                   *unroll180(lon1, lon2, wrap=wrap))


def cosineLaw_(phi2, phi1, lam21):
    '''Compute the I{angular} distance between two points using the
       U{spherical Law of Cosines
       <https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
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


def _distanceToE(func_, lat1, lat2, earth, d_lon, unused):
    '''(INTERNAL) Helper for ellipsoidal distances.
    '''
    E = _ellipsoidal(earth, func_)
    r =  func_(Phi_(lat2=lat2),
               Phi_(lat1=lat1), radians(d_lon), datum=E)
    return r * E.a


def _distanceToS(func_, lat1, lat2, earth, d_lon, unused, **adjust):
    '''(INTERNAL) Helper for spherical distances.
    '''
    r = func_(Phi_(lat2=lat2),
              Phi_(lat1=lat1), radians(d_lon), **adjust)
    return r * _mean_radius(earth, lat1, lat2)


def _ellipsoidal(earth, where):
    '''(INTERNAL) Helper for distances.
    '''
    return earth if isinstance(earth, Ellipsoid) else (
           earth if isinstance(earth, Datum) else
          _ellipsoidal_datum(earth, name=where.__name__)).ellipsoid  # PYCHOK indent


def equirectangular(lat1, lon1, lat2, lon2, radius=R_M, **options):
    '''Compute the distance between two points using
       the U{Equirectangular Approximation / Projection
       <https://www.Movable-Type.co.UK/scripts/latlong.html#equirectangular>}.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius, ellipsoid or datum
                      (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                      L{Datum} or L{a_f2Tuple}).
       @kwarg options: Optional keyword arguments for function
                       L{equirectangular_}.

       @return: Distance (C{meter}, same units as B{C{radius}} or
                the ellipsoid or datum axes).

       @raise TypeError: Invalid B{C{radius}}.

       @see: Function L{equirectangular_} for more details, the
             available B{C{options}}, errors, restrictions and other,
             approximate or accurate distance functions.
    '''
    d = sqrt(equirectangular_(Lat(lat1=lat1), Lon(lon1=lon1),
                              Lat(lat2=lat2), Lon(lon2=lon2),
                            **options).distance2)  # PYCHOK 4 vs 2-3
    return degrees2m(d, radius=_mean_radius(radius, lat1, lat2))


def equirectangular_(lat1, lon1, lat2, lon2,
                     adjust=True, limit=45, wrap=False):
    '''Compute the distance between two points using
       the U{Equirectangular Approximation / Projection
       <https://www.Movable-Type.co.UK/scripts/latlong.html#equirectangular>}.

       This approximation is valid for short distance of several
       hundred Km or Miles, see the B{C{limit}} keyword argument and
       the L{LimitError}.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg adjust: Adjust the wrapped, unrolled longitudinal delta
                      by the cosine of the mean latitude (C{bool}).
       @kwarg limit: Optional limit for lat- and longitudinal deltas
                     (C{degrees}) or C{None} or C{0} for unlimited.
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

       @return: A L{Distance4Tuple}C{(distance2, delta_lat, delta_lon,
                unroll_lon2)}.

       @raise LimitError: If the lat- and/or longitudinal delta exceeds the
                          B{C{-limit..+limit}} range and L{pygeodesy.limiterrors}
                          set to C{True}.

       @see: U{Local, flat earth approximation
             <https://www.EdWilliams.org/avform.htm#flat>}, functions
             L{equirectangular}, L{cosineAndoyerLambert},
             L{cosineForsytheAndoyerLambert}, L{cosineLaw}, L{euclidean},
             L{flatLocal}/L{hubeny}, L{flatPolar}, L{haversine}, L{thomas}
             and L{vincentys} and methods L{Ellipsoid.distance2},
             C{LatLon.distanceTo*} and C{LatLon.equirectangularTo}.
    '''
    d_lat = lat2 - lat1
    d_lon, ulon2 = unroll180(lon1, lon2, wrap=wrap)

    if limit and _limiterrors \
             and max(abs(d_lat), abs(d_lon)) > limit > 0:
        t = unstr(equirectangular_.__name__,
                  lat1, lon1, lat2, lon2, limit=limit)
        raise LimitError('delta exceeds limit', txt=t)

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
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}).
       @kwarg adjust: Adjust the longitudinal delta by the cosine of the
                      mean latitude (C{bool}).
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

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
    return _distanceToS(euclidean_, lat1, lat2, radius,
                                   *unroll180(lon1, lon2, wrap=wrap),
                                    adjust=adjust)


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


def excessAbc(A, b, c):
    '''Compute the I{spherical excess} C{E} of a (spherical) triangle
       from two sides and the included angle.

       @arg A: An interior triangle angle (C{radians}).
       @arg b: Frist adjacent triangle side (C{radians}).
       @arg c: Second adjacent triangle side (C{radians}).

       @return: Spherical excess ({radians}).

       @raise UnitError: Invalid B{C{A}}, B{C{b}} or B{C{c}}.

       @see: Function L{excessGirard}, L{excessLHuilier}, U{Spherical
             trigonometry<https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    sA, cA, sb, cb, sc, cc = sincos2_(Radians_(A=A), Radians_(b=b) * _0_5,
                                                     Radians_(c=c) * _0_5)
    return atan2(sA * sb * sc, cb * cc + cA * sb * sc) * _2_0


def excessGirard(A, B, C):
    '''Compute the I{spherical excess} C{E} of a (spherical) triangle using
       U{Girard's<https://MathWorld.Wolfram.com/GirardsSphericalExcessFormula.html>}
       formula.

       @arg A: First interior triangle angle (C{radians}).
       @arg B: Second interior triangle angle (C{radians}).
       @arg C: Third interior triangle angle (C{radians}).

       @return: Spherical excess ({radians}).

       @raise UnitError: Invalid B{C{A}}, B{C{B}} or B{C{C}}.

       @see: Function L{excessLHuilier}, U{Spherical trigonometry
             <https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    return Radians(Girard=fsum_(Radians_(A=A),
                                Radians_(B=B),
                                Radians_(C=C), -PI))


def excessLHuilier(a, b, c):
    '''Compute the I{spherical excess} C{E} of a (spherical) triangle using
       U{L'Huilier's<https://MathWorld.Wolfram.com/LHuiliersTheorem.html>}
       Theorem.

       @arg a: First triangle side (C{radians}).
       @arg b: Second triangle side (C{radians}).
       @arg c: Third triangle side (C{radians}).

       @return: Spherical excess ({radians}).

       @raise UnitError: Invalid B{C{a}}, B{C{b}} or B{C{c}}.

       @see: Function L{excessGirard}, U{Spherical trigonometry
             <https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    a = Radians_(a=a)
    b = Radians_(b=b)
    c = Radians_(c=c)

    s = fsum_(a, b, c) * _0_5
    r = tan_2(s) * tan_2(s - a) * tan_2(s - b) * tan_2(s - c)
    return Radians(LHuilier=atan(sqrt(r)) * _4_0)


def excessKarney(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the surface area of a (spherical) quadrilateral bounded by
       a segment of a great circle, two meridians and the equator using U{Karney's
       <http://OSGeo-org.1560.x6.Nabble.com/Area-of-a-spherical-polygon-td3841625.html>}
       method.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter}, L{Ellipsoid},
                      L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}) or C{None}.
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

       @return: Surface area, I{signed} (I{square} C{meter}, or units of B{C{radius}}
                I{squared}) or I{spherical excess} (C{radians}) if B{C{radius}} is
                C{None} or C{0}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise UnitError: Invalid B{C{lat2}} or B{C{lat1}}.

       @raise ValueError: Semi-circular longitudinal delta.

       @see: Function L{excessKarney_} and L{excessQuad}.
    '''
    return _area_or_(excessKarney_, lat1, lat2, radius,
                                   *unroll180(lon1, lon2, wrap=wrap))


def excessKarney_(phi2, phi1, lam21):
    '''Compute the I{spherical excess} C{E} of a (spherical) quadrilateral bounded
       by a segment of a great circle, two meridians and the equator using U{Karney's
       <http://OSGeo-org.1560.x6.Nabble.com/Area-of-a-spherical-polygon-td3841625.html>}
       method.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Spherical excess, I{signed} (C{radians}).

       @raise ValueError: Semi-circular longitudinal delta B{C{lam21}}.

       @see: Function L{excessKarney}, U{Area of a spherical polygon
       <http://OSGeo-org.1560.x6.Nabble.com/Area-of-a-spherical-polygon-td3841625.html>}.
    '''
    # from: Veness <https://www.Movable-Type.co.UK/scripts/latlong.html>
    # Area method due to Karney: for each edge of the polygon,
    #
    #               tan(Δλ/2) · (tan(φ1/2) + tan(φ2/2))
    #    tan(E/2) = ------------------------------------
    #                    1 + tan(φ1/2) · tan(φ2/2)
    #
    # where E is the spherical excess of the trapezium obtained by
    # extending the edge to the equator-circle vector for each edge.
    t2 = tan_2(phi2)
    t1 = tan_2(phi1)
    t  = tan_2(lam21, lam21=None)
    return Radians(Karney=atan2(t * (t1 + t2),
                             _1_0 + (t1 * t2)) * _2_0)


def excessQuad(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the surface area of a (spherical) quadrilateral bounded by a segment
       of a great circle, two meridians and the equator.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}) or C{None}.
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

       @return: Surface area, I{signed} (I{square} C{meter}, or units of B{C{radius}}
                I{squared}) or I{spherical excess} (C{radians}) if B{C{radius}} is
                C{None} or C{0}.

       @raise TypeError: Invalid B{C{radius}}.

       @raise UnitError: Invalid B{C{lat2}} or B{C{lat1}}.

       @see: Function L{excessQuad_} and L{excessKarney}.
    '''
    return _area_or_(excessQuad_, lat1, lat2, radius,
                                 *unroll180(lon1, lon2, wrap=wrap))


def excessQuad_(phi2, phi1, lam21):
    '''Compute the I{spherical excess} C{E} of a (spherical) quadrilateral bounded
       by a segment of a great circle, two meridians and the equator.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Spherical excess, I{signed} (C{radians}).

       @see: Function L{excessQuad}, U{Spherical trigonometry
             <https://WikiPedia.org/wiki/Spherical_trigonometry>}.
    '''
    s = sin((phi2 + phi1) * _0_5)
    c = cos((phi2 - phi1) * _0_5)
    return Radians(Quad=atan2(tan_2(lam21) * s, c) * _2_0)


def flatLocal(lat1, lon1, lat2, lon2, datum=_WGS84, wrap=False):
    '''Compute the distance between two (ellipsoidal) points using
       the U{ellipsoidal Earth to plane projection<https://WikiPedia.org/
       wiki/Geographical_distance#Ellipsoidal_Earth_projected_to_a_plane>}
       aka U{Hubeny<https://www.OVG.AT/de/vgi/files/pdf/3781/>} formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg datum: Datum or ellipsoid to use (L{Datum}, L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}).
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as the B{C{datum}}'s
                ellipsoid axes).

       @raise TypeError: Invalid B{C{datum}}.

       @note: The meridional and prime_vertical radii of curvature
              are taken and scaled at the mean of both latitude.

       @see: Functions L{flatLocal_}/L{hubeny_}, L{cosineLaw}, L{flatPolar},
             L{cosineAndoyerLambert}, L{cosineForsytheAndoyerLambert},
             L{equirectangular}, L{euclidean}, L{haversine}, L{thomas},
             L{vincentys}, method L{Ellipsoid.distance2} and U{local, flat
             earth approximation<https://www.EdWilliams.org/avform.htm#flat>}.
    '''
    d, _ = unroll180(lon1, lon2, wrap=wrap)
    return flatLocal_(Phi_(lat2=lat2),
                      Phi_(lat1=lat1), radians(d), datum=datum)

hubeny = flatLocal  # PYCHOK for Karl Hubeny


def flatLocal_(phi2, phi1, lam21, datum=_WGS84):
    '''Compute the I{angular} distance between two (ellipsoidal) points using
       the U{ellipsoidal Earth to plane projection<https://WikiPedia.org/
       wiki/Geographical_distance#Ellipsoidal_Earth_projected_to_a_plane>}
       aka U{Hubeny<https://www.OVG.AT/de/vgi/files/pdf/3781/>} formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg datum: Datum or ellipsoid to use (L{Datum}, L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}).

       @return: Angular distance (C{radians}).

       @raise TypeError: Invalid B{C{datum}}.

       @note: The meridional and prime_vertical radii of curvature
              are taken and scaled I{at the mean of both latitude}.

       @see: Functions L{flatLocal}/L{hubeny}, L{cosineAndoyerLambert_},
             L{cosineForsytheAndoyerLambert_}, L{cosineLaw_},
             L{flatPolar_}, L{equirectangular_}, L{euclidean_},
             L{haversine_}, L{thomas_} and L{vincentys_} and U{local, flat
             earth approximation <https://www.EdWilliams.org/avform.htm#flat>}.
    '''
    E = _ellipsoidal(datum, flatLocal_)
    m, n = E.roc2_((phi2 + phi1) * _0_5, scaled=True)
    return hypot(m * (phi2 - phi1), n * lam21)

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
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}).
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}} or the
                ellipsoid or datum axes).

       @raise TypeError: Invalid B{C{radius}}.

       @see: Functions L{flatPolar_}, L{cosineAndoyerLambert},
             L{cosineForsytheAndoyerLambert},L{cosineLaw},
             L{flatLocal}/L{hubeny}, L{equirectangular},
             L{euclidean}, L{haversine}, L{thomas} and
             L{vincentys}.
    '''
    return _distanceToS(flatPolar_, lat1, lat2, radius,
                                   *unroll180(lon1, lon2, wrap=wrap))


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
    a1 = PI_2 - phi1  # co-latitude
    a2 = PI_2 - phi2  # co-latitude
    ab = _2_0 * a1 * a2 * cos(lam21)
    r2 = fsum_(a1**2, a2**2, -abs(ab))
    return sqrt0(r2)


def hartzell(pov, los=None, earth=_WGS84, **LatLon_and_kwds):
    '''Compute the intersection of a Line-Of-Sight from a Point-Of-View in
       space with the surface of the earth.

       @arg pov: Point-Of-View outside the earth (C{Cartesian}, L{Ecef9Tuple}
                 or L{Vector3d}).
       @kwarg los: Line-Of-Sight, I{direction} to earth (L{Vector3d}) or
                   C{None} to point to the earth' center.
       @kwarg earth: The earth model (L{Datum}, L{Ellipsoid}, L{Ellipsoid2},
                     L{a_f2Tuple} or C{scalar} radius in C{meter}).
       @kwarg LatLon_and_kwds: Optional C{LatLon} class to convert insection
                               point plus C{LatLon} keyword arguments, include
                               B{C{datum}} if different from B{C{earth}}.

       @return: The earth intersection (L{Vector3d}, C{Cartesian type} of
                B{C{pov}} or B{C{LatLon}}).

       @raise IntersectionError: Null B{C{pov}} or B{C{los}} vector, B{C{pov}}
                                 is inside the earth or B{C{los}} points outside
                                 the earth or points in an opposite direction.

       @raise TypeError: Invalid B{C{pov}}, B{C{los}} or B{C{earth}}.

       @see: U{Stephen Hartzell<https://StephenHartzell.medium.com/
             satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6>}
             and function L{pygeodesy.tyr3d} for B{C{los}}.
    '''
    D = earth if isinstance(earth, Datum) else \
           _spherical_datum(earth, name=hartzell.__name__)
    E = D.ellipsoid

    a2 = b2 = E.a2  # earth x, y, ...
    c2 = E.b2  # ... z half-axis squared
    q2 = E.b2_a2  # == c2 / a2
    bc = E.a * E.b  # == b * c

    p3 = _MODS.vector3d._otherV3d(pov=pov)
    u3 = _MODS.vector3d._otherV3d(los=los) if los else p3.negate()
    u3 =  u3.unit()  # unit vector, opposing signs

    x2, y2, z2 = p3.times_(p3).xyz  # == p3.x2y2z2
    ux, vy, wz = u3.times_(p3).xyz
    u2, v2, w2 = u3.times_(u3).xyz  # == u3.x2y2z2

    t = c2, c2, b2  # a2 factored out
    m = fdot(t, u2, v2, w2)
    if m < EPS0:  # zero or near-null LOS vector
        raise IntersectionError(pov=pov, los=los, earth=earth, txt=_near_(_null_))

    # a2 and b2 factored out, b2 == a2 and b2 / a2 == 1
    r = fsum_(b2 * w2,  c2 * v2,      -v2 * z2,      vy * wz * 2,
              c2 * u2, -u2 * z2,      -w2 * x2,      ux * wz * 2,
             -w2 * y2, -u2 * y2 * q2, -v2 * x2 * q2, ux * vy * 2 * q2)
    if r > 0:
        r = bc * sqrt(r)
    elif r < 0:  # LOS pointing away from or missing the earth
        t = _opposite_ if max(ux, vy, wz) > 0 else _outside_
        raise IntersectionError(pov=pov, los=los, earth=earth, txt=t)

    n = fdot(t, ux, vy, wz)
    d = (n + r) / m  # (n - r) / m for antipode
    if d > 0:  # POV inside or LOS missing the earth
        t = _inside_ if min(x2 - a2, y2 - b2, z2 - c2) < EPS else _outside_
        raise IntersectionError(pov=pov, los=los, earth=earth, txt=t)

    if fsum_(x2, y2, z2) < d**2:  # d beyond earth center
        raise IntersectionError(pov=pov, los=los, earth=earth, txt=_too_(_distant_))

    r = _xnamed(p3.minus(u3.times(d)), hartzell.__name__)
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
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}).
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

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
    return _distanceToS(haversine_, lat1, lat2, radius,
                                   *unroll180(lon1, lon2, wrap=wrap))


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
    d = abs(Distance(distance))
    if d > h:
        d, h = h, d

    if d > EPS0:
        d = d / h  # /= h chokes PyChecker
        s = sin(Phi_(angle=angle, clip=_180_0))
        s = fsum_(_1_0, _2_0 * s * d, d**2)
        if s > 0:
            return h * sqrt(s) - r

    raise _ValueError(angle=angle, distance=distance, radius=radius)


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
        d2 = h * fsum_(r, r, h)
    return sqrt0(d2)


def intersections2(lat1, lon1, radius1,
                   lat2, lon2, radius2, datum=None, wrap=True):
    '''Conveniently compute the intersections of two circles each defined
       by a (geodetic) center point and a radius, using either ...

       1) L{vector3d.intersections2} for small distances (below 100 KM or
       about 0.9 degrees) or if no B{C{datum}} is specified, or ...

       2) L{sphericalTrigonometry.intersections2} for a spherical B{C{datum}}
       or if B{C{datum}} is a C{scalar} representing the earth radius,
       conventionally in C{meter} or ...

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
       @kwarg datum: Optional ellipsoidal or spherical datum (L{Datum},
                     L{Ellipsoid}, L{Ellipsoid2}, L{a_f2Tuple} or
                     C{scalar} earth radius in C{meter}, same units as
                     B{C{radius1}} and B{C{radius2}}) or C{None}.
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: 2-Tuple of the intersection points, each a
                L{LatLon2Tuple}C{(lat, lon)}.  For abutting circles,
                the points are the same instance, aka I{radical center}.

       @raise IntersectionError: Concentric, antipodal, invalid or
                                 non-intersecting circles or no
                                 convergence.

       @raise TypeError: Invalid B{C{datum}}.

       @raise UnitError: Invalid B{C{lat1}}, B{C{lon1}}, B{C{radius1}}
                         B{C{lat2}}, B{C{lon2}} or B{C{radius2}}.
    '''
    if datum is None or euclidean(lat1, lon1, lat2, lon2, radius=R_M,
                                  adjust=True, wrap=wrap) < _100km:
        def _V2T(x, y, _, **unused):  # _ == z unused
            return LatLon2Tuple(y, x, name=intersections2.__name__)

        r1 = m2degrees(Radius_(radius1=radius1), radius=R_M, lat=lat1)
        r2 = m2degrees(Radius_(radius2=radius2), radius=R_M, lat=lat2)

        _, lon2 = unroll180(lon1, lon2, wrap=wrap)
        m = _MODS.vector3d
        t = m.intersections2(m.Vector3d(lon1, lat1, 0), r1,
                             m.Vector3d(lon2, lat2, 0), r2, sphere=False,
                               Vector=_V2T)
    else:
        def _LL2T(lat, lon, **unused):
            return LatLon2Tuple(lat, lon, name=intersections2.__name__)

        d = _spherical_datum(datum, name=intersections2.__name__)
        if d.isSpherical:
            m = _MODS.sphericalTrigonometry
        elif d.isEllipsoidal:
            try:
                if d.ellipsoid.geodesic:
                    pass
                m = _MODS.ellipsoidalKarney
            except ImportError:
                m = _MODS.ellipsoidalExact
        else:
            raise _AssertionError(datum=d)

        t = m.intersections2(m.LatLon(lat1, lon1, datum=d), radius1,
                             m.LatLon(lat2, lon2, datum=d), radius2,
                               LatLon=_LL2T, height=0, wrap=wrap)
    return t


def isantipode(lat1, lon1, lat2, lon2, eps=EPS):
    '''Check whether two points are antipodal, on diametrically
       opposite sides of the earth.

       @arg lat1: Latitude of one point (C{degrees}).
       @arg lon1: Longitude of one point (C{degrees}).
       @arg lat2: Latitude of the other point (C{degrees}).
       @arg lon2: Longitude of the other point (C{degrees}).
       @kwarg eps: Tolerance for near-equality (C{degrees}).

       @return: C{True} if points are antipodal within the
                B{C{eps}} tolerance, C{False} otherwise.

       @see: U{Geosphere<https://CRAN.R-Project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return abs(wrap90(lat1) + wrap90(lat2)) < eps and \
           abs(abs(wrap180(lon1) - wrap180(lon2)) % _360_0 - _180_0) < eps


def isantipode_(phi1, lam1, phi2, lam2, eps=EPS):
    '''Check whether two points are antipodal, on diametrically
       opposite sides of the earth.

       @arg phi1: Latitude of one point (C{radians}).
       @arg lam1: Longitude of one point (C{radians}).
       @arg phi2: Latitude of the other point (C{radians}).
       @arg lam2: Longitude of the other point (C{radians}).
       @kwarg eps: Tolerance for near-equality (C{radians}).

       @return: C{True} if points are antipodal within the
                B{C{eps}} tolerance, C{False} otherwise.

       @see: U{Geosphere<https://CRAN.R-Project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return abs(wrapPI_2(phi1) + wrapPI_2(phi2)) < eps and abs(
           abs(wrapPI(lam1)   - wrapPI(lam2)) % PI2 - PI) < eps


def latlon2n_xyz(lat, lon, name=NN):
    '''Convert lat-, longitude to C{n-vector} (normal to the
       earth's surface) X, Y and Z components.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg name: Optional name (C{str}).

       @return: A L{Vector3Tuple}C{(x, y, z)}.

       @see: Function L{philam2n_xyz}.

       @note: These are C{n-vector} x, y and z components,
              I{NOT} geocentric ECEF x, y and z coordinates!
    '''
    return philam2n_xyz(radians(lat), radians(lon), name=name)


def n_xyz2latlon(x, y, z, name=NN):
    '''Convert C{n-vector} components to lat- and longitude in C{degrees}.

       @arg x: X component (C{scalar}).
       @arg y: Y component (C{scalar}).
       @arg z: Z component (C{scalar}).
       @kwarg name: Optional name (C{str}).

       @return: A L{LatLon2Tuple}C{(lat, lon)}.

       @see: Function L{n_xyz2philam}.
    '''
    a, b = n_xyz2philam(x, y, z)  # PYCHOK PhiLam2Tuple
    return LatLon2Tuple(degrees90(a), degrees180(b), name=name)


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


def opposing(bearing1, bearing2, margin=None):
    '''Compare the direction of two bearings given in C{degrees}.

       @arg bearing1: First bearing (compass C{degrees}).
       @arg bearing2: Second bearing (compass C{degrees}).
       @kwarg margin: Optional, interior angle bracket (C{degrees}),
                      default C{90}.

       @return: C{True} if both bearings point in opposite, C{False} if
                in similar or C{None} if in perpendicular directions.

       @see: Function L{opposing_}.
    '''
    m =  Degrees_(margin=margin, low=EPS0, high=_90_0) if margin else _90_0
    d = (bearing2 - bearing1) % _360_0  # note -20 % 360 == 340
    return False if      d < m or d > (_360_0 - m) else (
           True if (_180_0 - m) < d < (_180_0 + m) else None)


def opposing_(radians1, radians2, margin=None):
    '''Compare the direction of two bearings given in C{radians}.

       @arg radians1: First bearing (C{radians}).
       @arg radians2: Second bearing (C{radians}).
       @kwarg margin: Optional, interior angle bracket (C{radians}),
                      default C{PI_2}.

       @return: C{True} if both bearings point in opposite, C{False} if
                in similar or C{None} if in perpendicular directions.

       @see: Function L{opposing}.
    '''
    m =  Radians_(margin=margin, low=EPS0, high=PI_2) if margin else PI_2
    r = (radians2 - radians1) % PI2  # note -1 % PI2 == PI2 - 1
    return False if  r < m or r > (PI2 - m) else (
           True if (PI - m) < r < (PI  + m) else None)


def philam2n_xyz(phi, lam, name=NN):
    '''Convert lat-, longitude to C{n-vector} (normal to the
       earth's surface) X, Y and Z components.

       @arg phi: Latitude (C{radians}).
       @arg lam: Longitude (C{radians}).
       @kwarg name: Optional name (C{str}).

       @return: A L{Vector3Tuple}C{(x, y, z)}.

       @see: Function L{latlon2n_xyz}.

       @note: These are C{n-vector} x, y and z components,
              I{NOT} geocentric ECEF x, y and z coordinates!
    '''
    # Kenneth Gade eqn 3, but using right-handed
    # vector x -> 0°E,0°N, y -> 90°E,0°N, z -> 90°N
    sa, ca, sb, cb = sincos2_(phi, lam)
    return Vector3Tuple(ca * cb, ca * sb, sa, name=name)


def _radical2(d, r1, r2):  # in .ellipsoidalBaseDI, .sphericalTrigonometry, .vector3d
    # (INTERNAL) See C{radical2} below
    # assert d > EPS0
    r = fsum_(_1_0, (r1 / d)**2, -(r2 / d)**2) * _0_5
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


def _scale_deg(lat1, lat2):  # degrees
    # scale factor cos(mean of lats) for delta lon
    m = abs(lat1 + lat2) * _0_5
    return cos(radians(m)) if m < _90_0 else _0_0


def _scale_rad(phi1,  phi2):  # radians, by .frechet, .hausdorff, .heights
    # scale factor cos(mean of phis) for delta lam
    m = abs(phi1 + phi2) * _0_5
    return cos(m) if m < PI_2 else _0_0


def _sincosa6(phi2, phi1, lam21):
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
       @kwarg datum: Datum or ellipsoid to use (L{Datum}, L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}).
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as the B{C{datum}}'s
                ellipsoid axes).

       @raise TypeError: Invalid B{C{datum}}.

       @see: Functions L{thomas_}, L{cosineAndoyerLambert}, L{cosineForsytheAndoyerLambert},
             L{cosineLaw}, L{equirectangular}, L{euclidean}, L{flatLocal}/L{hubeny},
             L{flatPolar}, L{haversine}, L{vincentys} and method L{Ellipsoid.distance2}.
    '''
    return _distanceToE(thomas_, lat1, lat2, datum,
                                *unroll180(lon1, lon2, wrap=wrap))


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
    if r and _isnon0(c1) and _isnon0(c2):
        E = _ellipsoidal(datum, thomas_)
        if E.f:
            r1 = atan2(E.b_a * s1, c1)
            r2 = atan2(E.b_a * s2, c2)

            j = (r2 + r1) * _0_5
            k = (r2 - r1) * _0_5
            sj, cj, sk, ck, h, _ = sincos2_(j, k, lam21 * _0_5)

            h =  fsum_(sk**2, (ck * h)**2, -(sj * h)**2)
            u = _1_0 - h
            if _isnon0(u) and _isnon0(h):
                r = atan(sqrt0(h / u)) * _2_0  # == acos(1 - 2 * h)
                sr, cr = sincos2(r)
                if _isnon0(sr):
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

                    t  = fsum_(a * x, -b * y, c * x**2, -d * y**2, e * x * y)
                    r -= fsum_(s * x, -y, -t * f * _0_25) * f * sr
    return r


def vincentys(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two (spherical) points using
       U{Vincenty's<https://WikiPedia.org/wiki/Great-circle_distance>}
       spherical formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}).
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @raise UnitError: Invalid B{C{radius}}.

       @see: Functions L{vincentys_}, L{cosineAndoyerLambert},
             L{cosineForsytheAndoyerLambert},L{cosineLaw}, L{equirectangular},
             L{euclidean}, L{flatLocal}/L{hubeny}, L{flatPolar},
             L{haversine} and L{thomas} and methods L{Ellipsoid.distance2},
             C{LatLon.distanceTo*} and C{LatLon.equirectangularTo}.

       @note: See note at function L{vincentys_}.
    '''
    return _distanceToS(vincentys_, lat1, lat2, radius,
                                   *unroll180(lon1, lon2, wrap=wrap))


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
