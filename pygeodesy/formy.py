
# -*- coding: utf-8 -*-

u'''Formulary of basic geodesy functions and approximations.

@newfield example: Example, Examples
'''
from pygeodesy.basics import EPS, R_M, len2, LimitError, \
                            _limiterrors, map1, _TypeError
from pygeodesy.datum import Datum, Datums
from pygeodesy.fmath import fsum_, hypot, hypot2
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import Distance4Tuple, LatLon2Tuple, PhiLam2Tuple, \
                            Points2Tuple, Vector3Tuple
from pygeodesy.streprs import fstr
from pygeodesy.utily import PI, PI2, PI_2, degrees2m, \
                            degrees90, degrees180, degrees360, \
                            isNumpy2, isTuple2, sincos2, unroll180, unrollPI, \
                            wrap90, wrap180, wrapPI, wrapPI_2

from math import acos, atan2, cos, degrees, radians, sin, sqrt  # pow

# all public contants, classes and functions
__all__ = _ALL_LAZY.formy
__version__ = '20.03.30'


def _scaled(lat1, lat2):  # degrees
    # scale factor cos(mean of lats) for delta lon
    m = abs(lat1 + lat2) * 0.5
    return cos(radians(m)) if m < 90 else 0


def _scaler(phi1,  phi2):  # radians, imported by .frechet, .hausdorff, .heights
    # scale factor cos(mean of phis) for delta lam
    m = abs(phi1 + phi2) * 0.5
    return cos(m) if m < PI_2 else 0


def antipode(lat, lon):
    '''Return the antipode, the point diametrically opposite
       to a given point in C{degrees}.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).

       @return: A L{LatLon2Tuple}C{(lat, lon)}.

       @see: U{Geosphere<https://CRAN.R-Project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return LatLon2Tuple(-wrap90(lat), wrap180(lon + 180))


def antipode_(phi, lam):
    '''Return the antipode, the point diametrically opposite
       to a given point in C{radians}.

       @arg phi: Latitude (C{radians}).
       @arg lam: Longitude (C{radians}).

       @return: A L{PhiLam2Tuple}C{(phi, lam)}.

       @see: U{Geosphere<https://CRAN.R-Project.org/web/packages/geosphere/geosphere.pdf>}.
    '''
    return PhiLam2Tuple(-wrapPI_2(phi), wrapPI(lam + PI))


def bearing(lat1, lon1, lat2, lon2, **options):
    '''Compute the initial or final bearing (forward or reverse
       azimuth) between a (spherical) start and end point.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg options: Optional keyword arguments for function
                       L{bearing_}.

       @return: Initial or final bearing (compass C{degrees360}) or
                zero if start and end point coincide.
    '''
    ab4 = map1(radians, lat1, lon1, lat2, lon2)
    return degrees(bearing_(*ab4, **options))


def bearing_(phi1, lam1, phi2, lam2, final=False, wrap=False):
    '''Compute the initial or final bearing (forward or reverse
       azimuth) between a (spherical) start and end point.

       @arg phi1: Start latitude (C{radians}).
       @arg lam1: Start longitude (C{radians}).
       @arg phi2: End latitude (C{radians}).
       @arg lam2: End longitude (C{radians}).
       @kwarg final: Return final bearing if C{True}, initial
                     otherwise (C{bool}).
       @kwarg wrap: Wrap and L{unrollPI} longitudes (C{bool}).

       @return: Initial or final bearing (compass C{radiansPI2}) or
                zero if start and end point coincide.
    '''
    if final:
        phi1, lam1, phi2, lam2 = phi2, lam2, phi1, lam1
        r = PI2 + PI
    else:
        r = PI2

    db, _ = unrollPI(lam1, lam2, wrap=wrap)
    sa1, ca1, sa2, ca2, sdb, cdb = sincos2(phi1, phi2, db)

    # see <https://MathForum.org/library/drmath/view/55417.html>
    x = ca1 * sa2 - sa1 * ca2 * cdb
    y = sdb * ca2
    return (atan2(y, x) + r) % PI2  # wrapPI2


def compassAngle(lat1, lon1, lat2, lon2, adjust=True, wrap=False):
    '''Return the angle from North for the direction vector
       M{(lon2 - lon1, lat2 - lat1)} between two points.

       Suitable only for short, non-near-polar vectors up to a few
       hundred Km or Miles.  Use function L{bearing} for longer
       vectors.

       @arg lat1: From latitude (C{degrees}).
       @arg lon1: From longitude (C{degrees}).
       @arg lat2: To latitude (C{degrees}).
       @arg lon2: To longitude (C{degrees}).
       @kwarg adjust: Adjust the longitudinal delta by the
                      cosine of the mean latitude (C{bool}).
       @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

       @return: Compass angle from North (C{degrees360}).

       @note: Courtesy Martin Schultz.

       @see: U{Local, flat earth approximation
             <https://www.EdWilliams.org/avform.htm#flat>}.
    '''
    d_lon, _ = unroll180(lon1, lon2, wrap=wrap)
    if adjust:  # scale delta lon
        d_lon *= _scaled(lat1, lat2)
    return degrees360(atan2(d_lon, lat2 - lat1))


def cosineLaw(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two points using the
       U{spherical Law of Cosines
       <https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
       fromula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @raise TypeError: Invalid B{C{datum}}.

       @see: Functions L{cosineLaw_}, L{equirectangular}, L{euclidean},
             L{flatLocal}, L{flatPolar}, L{haversine}, L{vincentys} and
             method L{Ellipsoid.distance2}.

       @note: See note at function L{vincentys_}.
    '''
    r = float(radius)
    if r:
        d, _ = unroll180(lon1, lon2, wrap=wrap)
        r *= cosineLaw_(radians(lat2), radians(lat1), radians(d))
    return r


def cosineLaw_(phi2, phi1, lam21):
    '''Compute the I{angular} distance between two points using
       the U{spherical Law of Cosines
       <https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
       fromula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Angular distance (C{radians}).

       @see: Functions L{cosineLaw}, L{equirectangular_}, L{euclidean_},
             L{flatLocal_}, L{flatPolar_}, L{haversine_} and L{vincentys_}.

       @note: See note at function L{vincentys_}.
    '''
    s2, c2, s1, c1, _, c21 = sincos2(phi2, phi1, lam21)
    return acos(s1 * s2 + c1 * c2 * c21)


def equirectangular(lat1, lon1, lat2, lon2, radius=R_M, **options):
    '''Compute the distance between two points using
       the U{Equirectangular Approximation / Projection
       <https://www.Movable-Type.co.UK/scripts/latlong.html#equirectangular>}.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg options: Optional keyword arguments for function
                       L{equirectangular_}.

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @see: Function L{equirectangular_} for more details, the
             available B{C{options}}, errors, restrictions and other,
             approximate or accurate distance functions.
    '''
    _, dy, dx, _ = equirectangular_(lat1, lon1, lat2, lon2, **options)  # PYCHOK Distance4Tuple
    return degrees2m(hypot(dx, dy), radius=radius)


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
       @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

       @return: A L{Distance4Tuple}C{(distance2, delta_lat, delta_lon,
                unroll_lon2)}.

       @raise LimitError: If the lat- and/or longitudinal delta exceeds
                          the B{C{-limit..+limit}} range and L{limiterrors}
                          set to C{True}.

       @see: U{Local, flat earth approximation
             <https://www.EdWilliams.org/avform.htm#flat>}, functions
             L{equirectangular}, L{cosineLaw}, L{euclidean}, L{flatLocal},
             L{flatPolar}, L{haversine}, L{vincentys} and methods
             L{Ellipsoid.distance2}, C{LatLon.distanceTo*} and
             C{LatLon.equirectangularTo}.
    '''
    d_lat = lat2 - lat1
    d_lon, ulon2 = unroll180(lon1, lon2, wrap=wrap)

    if limit and _limiterrors \
             and max(abs(d_lat), abs(d_lon)) > limit > 0:
        t = fstr((lat1, lon1, lat2, lon2), prec=4)
        raise LimitError('%s(%s, limit=%s) delta exceeds limit' %
                        ('equirectangular_', t, fstr(limit, prec=2)))

    if adjust:  # scale delta lon
        d_lon *= _scaled(lat1, lat2)

    d2 = hypot2(d_lat, d_lon)  # degrees squared!
    return Distance4Tuple(d2, d_lat, d_lon, ulon2 - lon2)


def euclidean(lat1, lon1, lat2, lon2, radius=R_M, adjust=True, wrap=False):
    '''Approximate the C{Euclidian} distance between two (spherical) points.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg adjust: Adjust the longitudinal delta by the cosine
                      of the mean latitude (C{bool}).
       @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @see: U{Distance between two (spherical) points
             <https://www.EdWilliams.org/avform.htm#Dist>}, functions
             L{euclidean_}, L{cosineLaw}, L{equirectangular}, L{flatLocal},
             L{flatPolar}, L{haversine}, L{vincentys} and methods
             L{Ellipsoid.distance2}, C{LatLon.distanceTo*} and
             C{LatLon.equirectangularTo}.
    '''
    r = float(radius)
    if r:
        d, _ = unroll180(lon1, lon2, wrap=wrap)
        r *= euclidean_(radians(lat2), radians(lat1), radians(d), adjust=adjust)
    return r


def euclidean_(phi2, phi1, lam21, adjust=True):
    '''Approximate the I{angular} C{Euclidean} distance between two
       (spherical) points.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg adjust: Adjust the longitudinal delta by the cosine
                      of the mean latitude (C{bool}).

       @return: Angular distance (C{radians}).

       @see: Functions L{euclidean}, L{cosineLaw_}, L{equirectangular_},
             L{flatLocal_}, L{flatPolar_}, L{haversine_} and L{vincentys_}.
    '''
    a, b = abs(phi2 - phi1), abs(lam21)
    if adjust:
        b *= _scaler(phi2, phi1)
    if a < b:
        a, b = b, a
    return a + b * 0.5  # 0.4142135623731


def flatLocal(lat1, lon1, lat2, lon2, datum=Datums.WGS84, wrap=False):
    '''Compute the distance between two (ellipsoidal) points using
       the U{ellipsoidal Earth to plane projection
       <https://WikiPedia.org/wiki/Geographical_distance#Ellipsoidal_Earth_projected_to_a_plane>}
       fromula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg datum: Optional, (ellipsoidal) datum to use (L{Datum}).
       @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as the B{C{datum}}'s
                ellipsoid axes).

       @raise TypeError: Invalid B{C{datum}}.

       @note: The meridional and prime_vertical radii of curvature
              are taken and scaled at the mean latitude.

       @see: Functions L{flatLocal_}, L{cosineLaw}, L{flatPolar},
             L{equirectangular}, L{euclidean}, L{haversine},
             L{vincentys}, method L{Ellipsoid.distance2} and
             U{local, flat earth approximation
             <https://www.edwilliams.org/avform.htm#flat>}.
    '''
    d, _ = unroll180(lon1, lon2, wrap=wrap)
    return flatLocal_(radians(lat2), radians(lat1), radians(d), datum=datum)


def flatLocal_(phi2, phi1, lam21, datum=Datums.WGS84):
    '''Compute the distance between two (ellipsoidal) points using
       the U{ellipsoidal Earth to plane projection
       <https://WikiPedia.org/wiki/Geographical_distance#Ellipsoidal_Earth_projected_to_a_plane>}
       fromula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).
       @kwarg datum: Optional, (ellipsoidal) datum to use (L{Datum}).

       @return: Distance (C{meter}, same units as the B{C{datum}}'s
                ellipsoid axes).

       @raise TypeError: Invalid B{C{datum}}.

       @note: The meridional and prime_vertical radii of curvature
              are taken and scaled at the mean latitude.

       @see: Functions L{flatLocal}, L{cosineLaw_}, L{flatPolar_},
             L{equirectangular_}, L{euclidean_}, L{haversine_} and
             L{vincentys_} and U{local, flat earth approximation
             <https://www.edwilliams.org/avform.htm#flat>}.
    '''
    _TypeError(Datum, datum=datum)
    m, n = datum.ellipsoid.roc2_((phi2 + phi1) * 0.5, scaled=True)
    return hypot(m * (phi2 - phi1), n * lam21)


def flatPolar(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two (spherical) points using
       the U{polar coordinate flat-Earth
       <https://WikiPedia.org/wiki/Geographical_distance#Polar_coordinate_flat-Earth_formula>}
       formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @see: Functions L{flatPolar_}, L{cosineLaw}, L{flatLocal},
             L{equirectangular}, L{euclidean}, L{haversine} and
             L{vincentys}.
    '''
    r = float(radius)
    if r:
        d, _ = unroll180(lon1, lon2, wrap=wrap)
        r *= flatPolar_(radians(lat2), radians(lat1), radians(d))
    return r


def flatPolar_(phi2, phi1, lam21):
    '''Compute the I{angular} distance between two (spherical) points
       using the U{polar coordinate flat-Earth
       <https://WikiPedia.org/wiki/Geographical_distance#Polar_coordinate_flat-Earth_formula>}
       formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Angular distance (C{radians}).

       @see: Functions L{flatPolar}, L{cosineLaw_}, L{equirectangular_},
             L{euclidean_}, L{flatLocal_}, L{haversine_} and L{vincentys_}.
    '''
    a1 = abs(PI_2 - phi1)  # co-latitude
    a2 = abs(PI_2 - phi2)  # co-latitude
    ab = abs(2 * a1 * a2 * cos(lam21))
    a = max(a1, a2, ab)
    if a > EPS:
        s = fsum_((a1 / a)**2, (a2 / a)**2, -ab / a**2)
        a *= sqrt(s) if s > 0 else 0
    return a


def haversine(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two (spherical) points using the
       U{Haversine<https://www.Movable-Type.co.UK/scripts/latlong.html>}
       formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @see: U{Distance between two (spherical) points
             <https://www.EdWilliams.org/avform.htm#Dist>}, functions
             L{cosineLaw}, L{equirectangular}, L{euclidean},
             L{flatLocal}, L{flatPolar}, L{vincentys} and methods
             L{Ellipsoid.distance2}, C{LatLon.distanceTo*} and
             C{LatLon.equirectangularTo}.

       @note: See note at function L{vincentys_}.
    '''
    r = float(radius)
    if r:
        d, _ = unroll180(lon1, lon2, wrap=wrap)
        r *= haversine_(radians(lat2), radians(lat1), radians(d))
    return r


def haversine_(phi2, phi1, lam21):
    '''Compute the I{angular} distance between two (spherical) points
       using the U{Haversine<https://www.Movable-Type.co.UK/scripts/latlong.html>}
       formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Angular distance (C{radians}).

       @see: Functions L{haversine}, L{cosineLaw_}, L{equirectangular_},
             L{euclidean_}, L{flatLocal_}, L{flatPolar_} and L{vincentys_}.

       @note: See note at function L{vincentys_}.
    '''
    def _hsin(rad):
        return sin(rad * 0.5)**2

    h = _hsin(phi2 - phi1) + cos(phi1) * cos(phi2) * _hsin(lam21)  # haversine
    try:
        r = atan2(sqrt(h), sqrt(1 - h)) * 2  # == asin(sqrt(h)) * 2
    except ValueError:
        r = 0 if h < 0.5 else PI
    return r


def heightOf(angle, distance, radius=R_M):
    '''Determine the height above the (spherical) earth after
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
    d, r = distance, radius
    if d > r:
        d, r = r, d

    if d > EPS:
        d = d / float(r)
        s = sin(radians(angle))
        s = fsum_(1, 2 * s * d, d**2)
        if s > 0:
            return r * sqrt(s) - float(radius)

    raise ValueError('%s%r' % (heightOf.__name__, (angle, distance, radius)))


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
    if min(height, radius) < 0:
        raise ValueError('%s%r' % (horizon.__name__, (height, radius)))

    if refraction:
        d2 = 2.415750694528 * height * radius  # 2.0 / 0.8279
    else:
        d2 = height * fsum_(radius, radius, height)
    return sqrt(d2)


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
           abs(abs(wrap180(lon1) - wrap180(lon2)) % 360 - 180) < eps


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
    return abs(wrapPI_2(phi1) + wrapPI_2(phi2)) < eps and \
           abs(abs(wrapPI(lam1) - wrapPI(lam2)) % PI2 - PI) < eps


def latlon2n_xyz(lat, lon):
    '''Convert lat-, longitude to C{n-vector} (normal to the
       earth's surface) X, Y and Z components.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).

       @return: A L{Vector3Tuple}C{(x, y, z)}.

       @see: Function L{philam2n_xyz}.

       @note: These are C{n-vector} x, y and z components,
              I{NOT} geocentric ECEF x, y and z coordinates!
    '''
    return philam2n_xyz(radians(lat), radians(lon))


def n_xyz2latlon(x, y, z):
    '''Convert C{n-vector} components to lat- and longitude in C{degrees}.

       @arg x: X component (C{scalar}).
       @arg y: Y component (C{scalar}).
       @arg z: Z component (C{scalar}).

       @return: A L{LatLon2Tuple}C{(lat, lon)}.

       @see: Function L{n_xyz2philam}.
    '''
    a, b = n_xyz2philam(x, y, z)  # PYCHOK PhiLam2Tuple
    return LatLon2Tuple(degrees90(a), degrees180(b))


def n_xyz2philam(x, y, z):
    '''Convert C{n-vector} components to lat- and longitude in C{radians}.

       @arg x: X component (C{scalar}).
       @arg y: Y component (C{scalar}).
       @arg z: Z component (C{scalar}).

       @return: A L{PhiLam2Tuple}C{(phi, lam)}.

       @see: Function L{n_xyz2latlon}.
    '''
    return PhiLam2Tuple(atan2(z, hypot(x, y)), atan2(y, x))


def philam2n_xyz(phi, lam):
    '''Convert lat-, longitude to C{n-vector} (normal to the
       earth's surface) X, Y and Z components.

       @arg phi: Latitude (C{radians}).
       @arg lam: Longitude (C{radians}).

       @return: A L{Vector3Tuple}C{(x, y, z)}.

       @see: Function L{latlon2n_xyz}.

       @note: These are C{n-vector} x, y and z components,
              I{NOT} geocentric ECEF x, y and z coordinates!
    '''
    # Kenneth Gade eqn 3, but using right-handed
    # vector x -> 0°E,0°N, y -> 90°E,0°N, z -> 90°N
    sa, ca, sb, cb = sincos2(phi, lam)
    return Vector3Tuple(ca * cb, ca * sb, sa)


def points2(points, closed=True, base=None, Error=ValueError):
    '''Check a path or polygon represented by points.

       @arg points: The path or polygon points (C{LatLon}[])
       @kwarg closed: Optionally, consider the polygon closed,
                      ignoring any duplicate or closing final
                      B{C{points}} (C{bool}).
       @kwarg base: Optionally, check all B{C{points}} against
                    this base class, if C{None} don't check.
       @kwarg Error: Exception to raise (C{ValueError}).

       @return: A L{Points2Tuple}C{(number, points)} with the number
                of points and the points C{list} or C{tuple}.

       @raise TypeError: Some B{C{points}} are not B{C{base}}.

       @raise Error: Insufficient number of B{C{points}}.
    '''
    n, points = len2(points)

    if closed:
        # remove duplicate or closing final points
        while n > 1 and points[n-1] in (points[0], points[n-2]):
            n -= 1
        # XXX following line is unneeded if points
        # are always indexed as ... i in range(n)
        points = points[:n]  # XXX numpy.array slice is a view!

    if n < (3 if closed else 1):
        raise Error('too few %s: %s' % ('points', n))

    if base and not (isNumpy2(points) or isTuple2(points)):
        for i in range(n):
            base.others(points[i], name='%s[%s]' % ('points', i))

    return Points2Tuple(n, points)


def vincentys(lat1, lon1, lat2, lon2, radius=R_M, wrap=False):
    '''Compute the distance between two (spherical) points using
       U{Vincenty's<https://WikiPedia.org/wiki/Great-circle_distance>}
       spherical formula.

       @arg lat1: Start latitude (C{degrees}).
       @arg lon1: Start longitude (C{degrees}).
       @arg lat2: End latitude (C{degrees}).
       @arg lon2: End longitude (C{degrees}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @see: Functions L{vincentys_}, L{cosineLaw}, L{equirectangular},
             L{euclidean}, L{flatLocal}, L{flatPolar}, L{haversine} and
             methods L{Ellipsoid.distance2}, C{LatLon.distanceTo*} and
             C{LatLon.equirectangularTo}.

       @note: See note at function L{vincentys_}.
    '''
    r = float(radius)
    if r:
        d, _ = unroll180(lon1, lon2, wrap=wrap)
        r *= vincentys_(radians(lat2), radians(lat1), radians(d))
    return r


def vincentys_(phi2, phi1, lam21):
    '''Compute the I{angular} distance between two (spherical) points using
       U{Vincenty's<https://WikiPedia.org/wiki/Great-circle_distance>}
       spherical formula.

       @arg phi2: End latitude (C{radians}).
       @arg phi1: Start latitude (C{radians}).
       @arg lam21: Longitudinal delta, M{end-start} (C{radians}).

       @return: Angular distance (C{radians}).

       @see: Functions L{vincentys}, L{cosineLaw_}, L{equirectangular_},
             L{euclidean_}, L{flatLocal_}, L{flatPolar_} and L{haversine_}.

       @note: Functions L{vincentys_}, L{haversine_} and L{cosineLaw_}
              produce equivalent results, but L{vincentys_} is suitable
              for antipodal points but slightly more expensive (M{3 cos,
              3 sin, 1 hypot, 1 atan2, 6 mul, 2 add}) than L{haversine_}
              (M{2 cos, 2 sin, 2 sqrt, 1 atan2, 5 mul, 1 add}) and
              L{cosineLaw_} (M{3 cos, 3 sin, 1 acos, 3 mul, 1 add}).
    '''
    sa1, ca1, sa2, ca2, sb21, cb21 = sincos2(phi1, phi2, lam21)

    c = ca2 * cb21
    x = sa1 * sa2 + ca1 * c
    y = ca1 * sa2 - sa1 * c
    return atan2(hypot(ca2 * sb21, y), x)

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
