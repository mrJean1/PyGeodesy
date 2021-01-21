
# -*- coding: utf-8 -*-

u'''Geocentric conversions transcribed from I{Charles Karney}'s C++ classes U{Geocentric
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Geocentric.html>} and
U{LocalCartesian<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1LocalCartesian.html>}
into pure Python classes L{EcefKarney} respectively L{EcefCartesian}, class L{EcefSudano}
based on I{John Sudano}'s U{paper<https://www.ResearchGate.net/publication/
3709199_An_exact_conversion_from_an_Earth-centered_coordinate_system_to_latitude_longitude_and_altitude>},
class L{EcefVeness} transcribed from I{Chris Veness}' JavaScript classes U{LatLonEllipsoidal,
Cartesian<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}
and class L{EcefYou} implementing I{Rey-Jer You}'s U{transformations
<https://www.ResearchGate.net/publication/240359424>}.

Following is a copy of I{Karney}'s U{Detailed Description
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Geocentric.html>}.

Convert between geodetic coordinates C{lat}-, C{lon}gitude and height C{h}
(measured vertically from the surface of the ellipsoid) to geocentric C{x},
C{y} and C{z} coordinates, also known as I{Earth-Centered, Earth-Fixed}
(U{ECEF<https://WikiPedia.org/wiki/ECEF>}).

The origin of geocentric coordinates is at the center of the earth.  The C{z}
axis goes thru the North pole, C{lat} = 90°.  The C{x} axis goes thru C{lat}
= 0°, C{lon} = 0°.

The local cartesian origin is at (C{lat0}, C{lon0}, C{height0}).  The C{z}
axis is normal to the ellipsoid, the C{y} axis points due North.  The plane
C{z = -height0} is tangent to the ellipsoid.

Forward conversion from geodetic to geocentric (ECEF) coordinates is straightforward.

For the reverse transformation we use Hugues Vermeille's U{Direct transformation
from geocentric coordinates to geodetic coordinates
<https://DOI.org/10.1007/s00190-002-0273-6>}, J. Geodesy (2002) 76, 451-454.

Several changes have been made to ensure that the method returns accurate
results for all finite inputs (even if h is infinite).  The changes are
described in Appendix B of C. F. F. Karney U{Geodesics on an ellipsoid of
revolution<https://ArXiv.org/abs/1102.1215v1>}, Feb. 2011, 85, 105-117
(U{preprint<https://ArXiv.org/abs/1102.1215v1>}).  Vermeille similarly updated
his method in U{An analytical method to transform geocentric into geodetic
coordinates<https://DOI.org/10.1007/s00190-010-0419-x>}, J. Geodesy (2011) 85,
105-117.  See U{Geocentric coordinates
<https://GeographicLib.SourceForge.io/html/geocentric.html>} for more information.

The errors in these routines are close to round-off.  Specifically, for points
within 5,000 km of the surface of the ellipsoid (either inside or outside the
ellipsoid), the error is bounded by 7 nm (7 nanometers) for the WGS84 ellipsoid.
See U{Geocentric coordinates<https://GeographicLib.SourceForge.io/html/geocentric.html>}
for further information on the errors.
'''

from pygeodesy.basics import copysign, isscalar, neg, _xinstanceof, \
                            _xsubclassof
from pygeodesy.datums import Datums, _ellipsoidal_datum
from pygeodesy.ellipsoids import a_f2Tuple
from pygeodesy.errors import _datum_datum, LenError, _ValueError, _xkwds
from pygeodesy.fmath import cbrt, fdot, Fsum, fsum_, hypot1, hypot2_
from pygeodesy.interns import EPS, EPS0, EPS1, EPS_2, NN, PI, PI_2, _a_, \
                             _C_, _convergence_, _datum_, _ellipsoid_, \
                             _EPS0__2, _f_, _h_, _height_, _lat_, _lat0_, \
                             _lon_, _lon0_, _M_, _name_, _no_, _singular_, \
                             _SPACE_, _x_, _y_, _z_, _0_, _0_0, _0_5, \
                             _1_0, _2_0, _3_0, _4_0, _6_0, _90_0
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _FOR_DOCS
from pygeodesy.named import _NamedBase, _NamedTuple, notOverloaded, _Pass
from pygeodesy.namedTuples import LatLon2Tuple, LatLon3Tuple, \
                                  PhiLam2Tuple, Vector3Tuple
from pygeodesy.props import Property_RO, _update_all
from pygeodesy.streprs import unstr
from pygeodesy.units import Height, Int, Lat, Lon, Meter, Scalar
from pygeodesy.utily import atan2d, degrees90, sincos2, sincos2d
from pygeodesy.vector3d import _xyzn4

from math import asin, atan2, cos, degrees, hypot, radians, sqrt

__all__ = _ALL_LAZY.ecef
__version__ = '21.01.19'

_prolate_ = 'prolate'
_TRIPS    =  17  # 8..9 sufficient, EcefSudano.reverse


def _llhn4(latlonh, lon, height, suffix=NN):
    '''(INTERNAL) Get C{lat, lon, h, name} as C{4-tuple}.
    '''
    try:
        llh = latlonh.lat, latlonh.lon, getattr(latlonh, _height_,
                                        getattr(latlonh, _h_, height))
    except AttributeError:
        llh = latlonh, lon, height
    try:
        lat, lon, h = map(float, llh)
    except (TypeError, ValueError) as x:
        t = _lat_, _lon_, _height_
        if suffix:
            t = (_ + suffix for _ in t)
        raise EcefError(txt=str(x), **dict(zip(t, llh)))

    if abs(lat) > _90_0:  # XXX RangeError
        raise EcefError(_lat_ + suffix, lat)

    return lat, lon, h, getattr(latlonh, _name_, NN)


def _sch3(y, x):
    '''(INTERNAL) Compute sin, cos and hypotenuse.
    '''
    h = hypot(y, x)
    if h > 0:  # EPS_2
        s, c = y / h, x / h
    else:
        s, c = _0_0, _1_0
    return s, c, h


class EcefError(_ValueError):
    '''An ECEF or C{Ecef*} related issue.
    '''
    pass


class _EcefBase(_NamedBase):
    '''(INTERNAL) Base class for L{EcefKarney}, L{EcefVeness} and L{EcefYou}.
    '''
    _datum = None
    _E     = None

    def __init__(self, a_ellipsoid, f, name):
        '''(INTERNAL) New C{Ecef...}.
        '''
        if name:
            self.name = name
        try:
            E = a_ellipsoid
            if f is None:
                pass
            elif isscalar(E) and isscalar(f):
                E = a_f2Tuple(E, f)
            else:
                raise ValueError

            d = _ellipsoidal_datum(E, name=name)
            E = d.ellipsoid
            if E.a < EPS or E.f > EPS1:
                raise ValueError

        except (TypeError, ValueError) as x:
            t = unstr(self.classname, a=a_ellipsoid, f=f)
            raise EcefError(_SPACE_(t, _ellipsoid_), txt=str(x))

        self._datum = d
        self._E = E

    @Property_RO
    def equatoradius(self):
        '''Get the I{equatorial} radius, semi-axis (C{meter}).
        '''
        return self._E.a

    equatorialRadius = a = equatoradius  # Karney property

    @Property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @Property_RO
    def ellipsoid(self):
        '''Get the ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).
        '''
        return self._E

    @Property_RO
    def flattening(self):  # Karney property
        '''Get the I{flattening} (C{float}), M{(a - b) / a}, positive for
           I{oblate}, negative for I{prolate} or C{0} for I{near-spherical}.
        '''
        return self._E.f

    f = flattening

    def forward(self, latlonh, lon=None, height=0, M=False):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.
        '''
        notOverloaded(self, self.forward, latlonh, lon=lon,
                                      height=height, M=M)

    def _forward(self, latlonh, lon, height, M, You=False):
        '''(INTERNAL) Common C{EcefKarney}, C{EcefVeness},
           C{EcefSudano} and C{EcefYou} forward.
        '''
        lat, lon, h, name = _llhn4(latlonh, lon, height)
        sa, ca, sb, cb = sincos2d(lat, lon)

        E = self.ellipsoid
        n = E.roc1_(sa, ca) if You else E.roc1_(sa)
        z = (h + n * E.e12) * sa
        x = (h + n) * ca

        m = self._Matrix(sa, ca, sb, cb) if M else None
        r = Ecef9Tuple(x * cb, x * sb, z, lat, lon, h, 0, m, self.datum)
        return self._xnamed(r, name=name)

    def _Matrix(self, sa, ca, sb, cb):
        '''Creation a rotation matrix.

           @arg sa: C{sin(phi)} (C{float}).
           @arg ca: C{cos(phi)} (C{float}).
           @arg sb: C{sin(lambda)} (C{float}).
           @arg cb: C{cos(lambda)} (C{float}).

           @return: An L{EcefMatrix}.
        '''
        return self._xnamed(EcefMatrix(sa, ca, sb, cb))

    def reverse(self, xyz, y=None, z=None, M=False):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.
        '''
        notOverloaded(self, self.reverse, xyz, y=y, z=z, M=M)

    def toStr(self, prec=9, **unused):  # PYCHOK signature
        '''Return this C{Ecef*} as a string.

           @kwarg prec: Optional precision, number of decimal digits (0..9).

           @return: This C{Ecef*} representation (C{str}).
        '''
        return self.attrs(_a_, _f_, _datum_, _ellipsoid_, _name_, prec=prec)


class EcefKarney(_EcefBase):
    '''Conversion between geodetic and geocentric, aka I{Earth-Centered,
       Earth-Fixed} (ECEF) coordinates based on I{Karney}'s U{Geocentric
       <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Geocentric.html>}
       methods.
    '''
    _hmax = 0  # max height, 12M lightyears

    def __init__(self, a_ellipsoid, f=None, name=NN):
        '''New L{EcefKarney} converter.

           @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}, L{Ellipsoid2}, L{Datum}
                             or L{a_f2Tuple}) or C{scalar} for the equatorial
                             radius of the ellipsoid (C{meter}).
           @kwarg f: C{None} or the ellipsoid flattening (C{scalar}), required
                     for C{scalar} B{C{a_ellipsoid}}, B{C{f=0}} represents a
                     sphere, negative B{C{f}} a prolate ellipsoid.
           @kwarg name: Optional name (C{str}).

           @raise EcefError: If B{C{a_ellipsoid}} not L{Ellipsoid}, L{Ellipsoid2},
                             L{Datum} or L{a_f2Tuple} or C{scalar} or B{C{f}} not
                             C{scalar} or if C{scalar} B{C{a_ellipsoid}} not positive
                             or B{C{f}} not less than 1.0.
        '''
        _EcefBase.__init__(self, a_ellipsoid, f, name)
        self._hmax = self.equatoradius / EPS_2  # self.equatoradius * _2_EPS

    def forward(self, latlonh, lon=None, height=0, M=False):
        '''Convert from geodetic C{(lat, lon, height)} to geocentric C{(x, y, z)}.

           @arg latlonh: Either a C{LatLon}, an L{Ecef9Tuple} or C{scalar}
                         latitude in C{degrees}.
           @kwarg lon: Optional C{scalar} longitude in C{degrees} for C{scalar}
                       B{C{latlonh}}.
           @kwarg height: Optional height in C{meter}, vertically above (or
                          below) the surface of the ellipsoid.
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with geocentric C{(x, y, z)} coordinates for the given
                    geodetic ones C{(lat, lon, height)}, case C{C} 0, optional
                    L{EcefMatrix} C{M} and C{datum} if available.

           @raise EcefError: If B{C{latlonh}} not C{LatLon}, L{Ecef9Tuple} or
                             C{scalar} or B{C{lon}} not C{scalar} for C{scalar}
                             B{C{latlonh}} or C{abs(lat)} exceeds 90°.

           @note: Let C{v} be a unit vector located at C{(lat, lon, h)}.  We can
                  express C{v} as column vectors in one of two ways, C{v1} in east,
                  north, up coordinates (where the components are relative
                  to a local coordinate system at C{C(lat0, lon0, h0)}) or as
                  C{v0} in geocentric C{x, y, z} coordinates.  Then, M{v0 =
                  M ⋅ v1} where C{M} is the rotation matrix.
        '''
        return _EcefBase._forward(self, latlonh, lon, height, M)

    @Property_RO
    def hmax(self):
        '''Get the distance limit (C{float}).
        '''
        return self._hmax

    def reverse(self, xyz, y=None, z=None, M=False):
        '''Convert from geocentric C{(x, y, z)} to geodetic C{(lat, lon, height)}.

           @arg xyz: Either an L{Ecef9Tuple}, an C{(x, y, z)} 3-tuple or C{scalar}
                     ECEF C{x} coordinate in C{meter}.
           @kwarg y: ECEF C{y} coordinate in C{meter} for C{scalar} B{C{xyz}}
                     and B{C{z}}.
           @kwarg z: ECEF C{z} coordinate in C{meter} for C{scalar} B{C{xyz}}
                     and B{C{y}}.
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with geodetic coordinates C{(lat, lon, height)} for the given
                    geocentric ones C{(x, y, z)}, case C{C}, optional L{EcefMatrix}
                    C{M} and C{datum} if available.

           @raise EcefError: If B{C{xyz}} not L{Ecef9Tuple} or C{scalar} C{x}
                             or B{C{y}} and/or B{C{z}} not C{scalar} for C{scalar}
                             B{C{xyz}}.

           @note: In general, there are multiple solutions and the result
                  which minimizes C{height} is returned, i.e., C{(lat, lon)}
                  corresponds to the closest point on the ellipsoid.  If
                  there are still multiple solutions with different latitudes
                  (applies only if C{z} = 0), then the solution with C{lat} > 0
                  is returned.  If there are still multiple solutions with
                  different longitudes (applies only if C{x} = C{y} = 0) then
                  C{lon} = 0 is returned.  The returned C{height} value is not
                  below M{−E.a * (1 − E.e2) / sqrt(1 − E.e2 * sin(lat)**2)}.
                  The returned C{lon} is in the range [−180°, 180°].  Like
                  C{forward} above, M{v1 = Transpose(M) ⋅ v0}.
        '''
        E = self.ellipsoid

        x, y, z, name = _xyzn4(xyz, y, z, Error=EcefError)

        sb, cb, R = _sch3(y, x)
        h = hypot(R, z)  # distance to earth center
        if h > self._hmax:  # PYCHOK no cover
            # We are really far away (> 12M light years).  Treat the earth
            # as a point and h, above as an acceptable approximation to the
            # height.  This avoids overflow, e.g., in the computation of disc
            # below.  It's possible that h has overflowed to INF, that's OK.
            # Treat finite x, y, but R overflows to +INF by scaling by 2.
            sb, cb, R = _sch3(y * _0_5, x * _0_5)
            sa, ca, _ = _sch3(z * _0_5, R)
            C = 1

        elif E.e4:  # E.isEllipsoidal
            prolate = E.isProlate
            # Treat prolate spheroids by swapping R and Z here and by
            # switching the arguments to phi = atan2(...) at the end.
            p = (R / E.a)**2
            q = E.e12 * (z / E.a)**2
            if prolate:
                p, q = q, p
            r = p + q - E.e4
            e = E.e4 * q
            if e or r > 0:
                # Avoid possible division by zero when r = 0 by multiplying
                # equations for s and t by r^3 and r, respectively.
                s = e * p / _4_0  # s = r^3 * s
                u = r = r / _6_0
                r2 = r**2
                r3 = r * r2
                t3 = s + r3
                disc = s * (r3 + t3)
                if disc < 0:
                    # t is complex, but the way u is defined, the result is real.
                    # There are three possible cube roots.  We choose the root
                    # which avoids cancellation.  Note, disc < 0 implies r < 0.
                    u += _2_0 * r * cos(atan2(sqrt(-disc), -t3) / _3_0)
                else:
                    # Pick the sign on the sqrt to maximize abs(T3).  This
                    # minimizes loss of precision due to cancellation.  The
                    # result is unchanged because of the way the t is used
                    # in definition of u.
                    if disc > 0:
                        t3 += copysign(sqrt(disc), t3)  # t3 = (r * t)^3
                    # N.B. cbrt always returns the real root, cbrt(-8) = -2.
                    t = cbrt(t3)  # t = r * t
                    # t can be zero; but then r2 / t -> 0.
                    if t:
                        u = fsum_(u, t, r2 / t)
                v = sqrt(e + u**2)  # guaranteed positive
                # Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
                # E.e4 * q / (v - u) because u ~ e^4 when q is small and u < 0.
                uv = (e / (v - u)) if u < 0 else (u + v)  # u+v, guaranteed positive
                # Need to guard against w going negative due to roundoff in uv - q.
                w = max(_0_0, E.e2abs * (uv - q) / (_2_0 * v))
                # Rearrange expression for k to avoid loss of accuracy due to
                # subtraction.  Division by 0 not possible because uv > 0, w >= 0.
                k1 = k2 = uv / (sqrt(uv + w**2) + w)
                if prolate:
                    k1 -= E.e2
                else:
                    k2 += E.e2
                sa, ca, h = _sch3(z / k1, R / k2)
                h *= k1 - E.e12
                C  = 2

            else:  # e = E.e4 * q == 0 and r <= 0
                # This leads to k = 0 (oblate, equatorial plane) and k + E.e^2 = 0
                # (prolate, rotation axis) and the generation of 0/0 in the general
                # formulas for phi and h, using the general formula and division
                # by 0 in formula for h.  Handle this case by taking the limits:
                #   f > 0: z -> 0, k        ->  E.e2 * sqrt(q) / sqrt(E.e4 - p)
                #   f < 0: r -> 0, k + E.e2 -> -E.e2 * sqrt(q) / sqrt(E.e4 - p)
                q = E.e4 - p
                if prolate:
                    p, q = q, p
                    e = E.a
                else:
                    e = E.b2_a
                sa, ca, h = _sch3(sqrt(q / E.e12), sqrt(p))
                if z < 0:
                    sa = neg(sa)  # for tiny negative z, not for prolate
                h *= neg(e / E.e2abs)
                C  = 3

        else:  # E.e4 == 0, spherical case
            # Dealing with underflow in the general case with E.e2 = 0 is
            # difficult.  Origin maps to North pole, same as with ellipsoid.
            sa, ca, _ = _sch3(z if h else _1_0, R)
            h -= E.a
            C  = 4

        r = Ecef9Tuple(x, y, z, atan2d(sa, ca),
                                atan2d(sb, cb), h, C,
                                self._Matrix(sa, ca, sb, cb) if M else None,
                                self.datum)
        return self._xnamed(r, name=name)


class EcefCartesian(_NamedBase):
    '''Conversion between geodetic C{(lat, lon, height)} and local cartesian
       C{(x, y, z)} coordinates with a local cartesian origin at C{(lat0, lon0,
       height0)}, transcribed from I{Karney}'s C++ class U{LocalCartesian
       <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1LocalCartesian.html>}.

       The C{z} axis is normal to the ellipsoid, the C{y} axis points due
       North.  The plane C{z = -heighth0} is tangent to the ellipsoid.

       The conversions all take place via geocentric coordinates using a
       geocentric L{EcefKarney}, by default the WGS84 datum/ellipsoid.
    '''
    _ecef = EcefKarney(Datums.WGS84)
    _t0   = None  # forward(lat0, lon0, height0) L{Ecef9Tuple}

    def __init__(self, latlonh0=0, lon0=0, height0=0, ecef=None, name=NN):
        '''New L{EcefCartesian} converter.

           @kwarg latlonh0: Either a C{LatLon}, an L{Ecef9Tuple} or C{scalar}
                            latitude in C{degrees} of the cartesian origin.
           @kwarg lon0: Optional C{scalar} longitude of the cartesian origin
                        in C{degrees} for C{scalar} B{C{latlonh0}}.
           @kwarg height0: Optional height of the cartesian origin in C{meter},
                           vertically above (or below) the surface of the ellipsoid.
           @kwarg ecef: An ECEF converter (L{EcefKarney}).
           @kwarg name: Optional name (C{str}).

           @raise EcefError: If B{C{latlonh0}} not C{LatLon}, L{Ecef9Tuple} or
                             C{scalar} or B{C{lon0}} not C{scalar} for C{scalar}
                             B{C{latlonh0}} or C{abs(lat)} exceeds 90°.

           @raise TypeError: Invalid B{C{ecef}}, not L{EcefKarney}.
        '''
        if ecef:
            _xinstanceof(EcefKarney, ecef=ecef)
            self._ecef = ecef
        self.reset(latlonh0, lon0, height0, name=name)

    @Property_RO
    def ecef(self):
        '''Get the ECEF converter (L{EcefKarney}).
        '''
        return self._ecef

    def forward(self, latlonh, lon=None, height=0, M=False):
        '''Convert from geodetic C{(lat, lon, height)} to local cartesian
           C{(x, y, z)}.

           @arg latlonh: Either a C{LatLon}, an L{Ecef9Tuple} or C{scalar}
                         latitude in C{degrees}.
           @kwarg lon: Optional C{scalar} longitude in C{degrees} for C{scalar}
                       B{C{latlonh}}.
           @kwarg height: Optional height in C{meter}, vertically above (or
                          below) the surface of the ellipsoid.
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with geocentric C{(x, y, z)} coordinates for the given
                    geodetic ones C{(lat, lon, height)}, case C{C} 0, optional
                    L{EcefMatrix} C{M} and C{datum} if available.

           @raise EcefError: If B{C{latlonh}} not C{LatLon}, L{Ecef9Tuple} or
                             C{scalar} or B{C{lon}} not C{scalar} for C{scalar}
                             B{C{latlonh}} or C{abs(lat)} exceeds 90°.

           @see: Note at method L{EcefKarney.forward}.
        '''
        lat, lon, h, name = _llhn4(latlonh, lon, height)
        t = self.ecef.forward(lat, lon, h, M=M)
        x, y, z = self.M.rotate(t[:3], *self._t0[:3])  # .x, .y, .z

        m = self.M.multiply(t.M) if M else None
        r = Ecef9Tuple(x, y, z, t.lat, t.lon, t.height, 0, m, self.ecef.datum)
        return self._xnamed(r, name=name)

    @Property_RO
    def height0(self):
        '''Get origin's height (C{meter}).
        '''
        return self._t0.height

    @Property_RO
    def lat0(self):
        '''Get origin's latitude (C{degrees}).
        '''
        return self._t0.lat

    @Property_RO
    def lon0(self):
        '''Get origin's longitude (C{degrees}).
        '''
        return self._t0.lon

    @Property_RO
    def M(self):
        '''Get the rotation matrix (C{EcefMatrix}).
        '''
        return self._t0.M

    def reset(self, latlonh0=0, lon0=0, height0=0, name=NN):
        '''Reset the local cartesian origin.

           @kwarg latlonh0: Either a C{LatLon}, an L{Ecef9Tuple} or C{scalar}
                            latitude in C{degrees} of the cartesian origin.
           @kwarg lon0: Optional, C{scalar} longitude of the cartesian origin
                        in C{degrees} for C{scalar} B{C{latlonh0}}.
           @kwarg height0: Optional, height of the cartesian origin in C{meter},
                           vertically above (or below) the surface of the ellipsoid.
           @kwarg name: Optional, new name (C{str}).

           @raise EcefError: If B{C{latlonh0}} not C{LatLon}, L{Ecef9Tuple} or
                             C{scalar} or B{C{lon0}} not C{scalar} for C{scalar}
                             B{C{latlonh0}} or C{abs(lat)} exceeds 90°.
        '''
        _update_all(self)

        lat0, lon0, height0, n = _llhn4(latlonh0, lon0, height0, suffix=_0_)
        if name:
            n = name
        if n:
            self.rename(n)
            self.ecef.rename(n)
        self._t0 = self.ecef.forward(lat0, lon0, height0, M=True)

    def reverse(self, xyz, y=None, z=None, M=False):
        '''Convert from local cartesian C{(x, y, z)} to geodetic C{(lat, lon, height)}.

           @arg xyz: Either an L{Ecef9Tuple}, an C{(x, y, z)} 3-tuple or C{scalar}
                     local cartesian C{x} coordinate in C{meter}.
           @kwarg y: Local cartesian C{y} coordinate in C{meter} for C{scalar}
                     B{C{xyz}} and B{C{z}}.
           @kwarg z: Local cartesian  C{z} coordinate in C{meter} for C{scalar}
                     B{C{xyz}} and B{C{y}}.
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with geodetic coordinates C{(lat, lon, height)} for the given
                    geocentric ones C{(x, y, z)}, case C{C}, optional L{EcefMatrix}
                    C{M} and C{datum} if available.

           @raise EcefError: If B{C{xyz}} not L{Ecef9Tuple} or C{scalar} C{x}
                             or B{C{y}} and/or B{C{z}} not C{scalar} for C{scalar}
                             B{C{xyz}}.

           @see: Note at method L{EcefKarney.reverse}.
        '''
        xyz_n = _xyzn4(xyz, y, z, Error=EcefError)
        x, y, z = self.M.unrotate(xyz_n[:3], *self._t0[:3])  # .x, .y, .z
        t = self.ecef.reverse(x, y, z, M=M)

        m = self.M.multiply(t.M) if M else None
        r = Ecef9Tuple(x, y, z, t.lat, t.lon, t.height, t.C, m, self.ecef.datum)
        return self._xnamed(r, name=xyz_n[3])

    def toStr(self, prec=9):  # PYCHOK signature
        '''Return this L{EcefCartesian} as a string.

           @kwarg prec: Optional precision, number of decimal digits (0..9).

           @return: This L{EcefCartesian} representation (C{str}).
        '''
        return self.attrs(_lat0_, _lon0_, 'height0', _M_, 'ecef', _name_, prec=prec)


class EcefMatrix(_NamedTuple):
    '''A rotation matrix.
    '''
    _Names_ = ('_0_0_', '_0_1_', '_0_2_',  # row-order
               '_1_0_', '_1_1_', '_1_2_',
               '_2_0_', '_2_1_', '_2_2_')
    _Units_ = (Scalar,) * len(_Names_)

    def _validate(self, **_OK):  # PYCHOK unused
        '''(INTERNAL) Allow C{_Names_} with leading underscore.
        '''
        _NamedTuple._validate(self, _OK=True)

    def __new__(cls, sa, ca, sb, cb, *_m):
        '''New L{EcefMatrix} matrix.

           @arg sa: C{sin(phi)} (C{float}).
           @arg ca: C{cos(phi)} (C{float}).
           @arg sb: C{sin(lambda)} (C{float}).
           @arg cb: C{cos(lambda)} (C{float}).
           @arg _m: (INTERNAL) from C{.multiply}.

           @raise EcefError: If B{C{sa}}, B{C{ca}}, B{C{sb}} or
                             B{C{cb}} outside M{[-1.0, +1.0]}.
        '''
        t = sa, ca, sb, cb
        if _m:  # all 9 matrix elements ...
            t += _m  # ... from .multiply

        elif max(map(abs, t)) > 1:
            raise EcefError(unstr(EcefMatrix.__name__, *t))

        else:  # build matrix from the following quaternion operations
            #   qrot(lam, [0,0,1]) * qrot(phi, [0,-1,0]) * [1,1,1,1]/2
            # or
            #   qrot(pi/2 + lam, [0,0,1]) * qrot(-pi/2 + phi, [-1,0,0])
            # where
            #   qrot(t,v) = [cos(t/2), sin(t/2)*v[1], sin(t/2)*v[2], sin(t/2)*v[3]]

            # Local X axis (east) in geocentric coords
            #  M[0] = -slam;        M[3] =  clam;        M[6] = 0;
            # Local Y axis (north) in geocentric coords
            #  M[1] = -clam * sphi; M[4] = -slam * sphi; M[7] = cphi;
            # Local Z axis (up) in geocentric coords
            #  M[2] =  clam * cphi; M[5] =  slam * cphi; M[8] = sphi;
            t = (-sb, -cb * sa, cb * ca,
                  cb, -sb * sa, sb * ca,
                _0_0,       ca,      sa)

        return _NamedTuple.__new__(cls, *t)

    def copy(self, **unused):  # PYCHOK signature
        '''Make a shallow or deep copy of this instance.

           @return: The copy (C{This class} or subclass thereof).
        '''
        return self.classof(*self)

    __copy__ = __deepcopy__ = copy

    def multiply(self, other):
        '''Matrix multiply M{M0' ⋅ M} this matrix transposed with
           an other matrix.

           @arg other: The other matrix (L{EcefMatrix}).

           @return: The matrix product (L{EcefMatrix}).

           @raise TypeError: If B{C{other}} is not L{EcefMatrix}.
        '''
        _xinstanceof(EcefMatrix, other=other)

        # like LocalCartesian.MatrixMultiply, transposed(self) x other
        # <https://GeographicLib.SourceForge.io/html/LocalCartesian_8cpp_source.html>
        M = (fdot(self[r::3], *other[c::3]) for r in range(3) for c in range(3))
        return EcefMatrix(*M)

    def rotate(self, xyz, *xyz0):
        '''Forward rotation M{M0' ⋅ ([x, y, z] - [x0, y0, z0])'}.

           @arg xyz: Local C{(x, y, z)} coordinates (C{3-tuple}).
           @arg xyz0: Optional, local C{(x0, y0, z0)} origin (C{3-tuple}).

           @return: Rotated C{(x, y, z)} location (C{3-tuple}).

           @raise LenError: Unequal C{len(B{xyz})} and C{len(B{xyz0})}.
        '''
        if xyz0:
            if len(xyz0) != len(xyz):
                raise LenError(self.rotate, xyz0=len(xyz0), xyz=len(xyz))

            xyz = tuple(c - c0 for c, c0 in zip(xyz, xyz0))

        # x' = M[0] * x + M[3] * y + M[6] * z
        # y' = M[1] * x + M[4] * y + M[7] * z
        # z' = M[2] * x + M[5] * y + M[8] * z
        return (fdot(xyz, *self[0::3]),
                fdot(xyz, *self[1::3]),
                fdot(xyz, *self[2::3]))

    def unrotate(self, xyz, *xyz0):
        '''Inverse rotation M{[x0, y0, z0] + M0 ⋅ [x,y,z]'}.

           @arg xyz: Local C{(x, y, z)} coordinates (C{3-tuple}).
           @arg xyz0: Optional, local C{(x0, y0, z0)} origin (C{3-tuple}).

           @return: Unrotated C{(x, y, z)} location (C{3-tuple}).

           @raise LenError: Unequal C{len(B{xyz})} and C{len(B{xyz0})}.
        '''
        if xyz0:
            if len(xyz0) != len(xyz):
                raise LenError(self.unrotate, xyz0=len(xyz0), xyz=len(xyz))

            _xyz = (_1_0,) + xyz
            # x' = x0 + M[0] * x + M[1] * y + M[2] * z
            # y' = y0 + M[3] * x + M[4] * y + M[5] * z
            # z' = z0 + M[6] * x + M[7] * y + M[8] * z
            xyz_ = (fdot(_xyz, xyz0[0], *self[0:3]),
                    fdot(_xyz, xyz0[1], *self[3:6]),
                    fdot(_xyz, xyz0[2], *self[6:]))
        else:
            # x' = M[0] * x + M[1] * y + M[2] * z
            # y' = M[3] * x + M[4] * y + M[5] * z
            # z' = M[6] * x + M[7] * y + M[8] * z
            xyz_ = (fdot(xyz, *self[0:3]),
                    fdot(xyz, *self[3:6]),
                    fdot(xyz, *self[6:]))
        return xyz_


class Ecef9Tuple(_NamedTuple):  # .ecef.py
    '''9-Tuple C{(x, y, z, lat, lon, height, C, M, datum)} with geocentric
       coordinates C{x}, C{y} and C{z}, geodetic coordinates C{lat}, C{lon}
       and C{height}, case C{C} (see the C{Ecef*.reverse} methods) and
       optionally, the L{EcefMatrix} C{M} and C{datum}, with C{lat} and
       C{lon} in C{degrees} and C{x}, C{y}, C{z} and C{height} in C{meter},
       conventionally.
    '''
    _Names_ = (_x_,   _y_,   _z_,   _lat_, _lon_, _height_, _C_,  _M_,   _datum_)
    _Units_ = ( Meter, Meter, Meter, Lat,   Lon,   Height,   Int, _Pass, _Pass)

    @Property_RO
    def lam(self):
        '''Get the longitude in C{radians} (C{float}).
        '''
        return self.philam.lam

    @Property_RO
    def lamVermeille(self):
        '''Get the longitude in C{[-PI*3/2..+PI*3/2] radians} after U{Vermeille
           <https://Search.ProQuest.com/docview/639493848>} (2004), p 95.

           @see: U{Karney<https://GeographicLib.SourceForge.io/html/geocentric.html>},
                 U{Vermeille<https://Search.ProQuest.com/docview/847292978>} 2011, pp 112-113, 116
                 and U{Featherstone, et.al.<https://Search.ProQuest.com/docview/872827242>}, p 7.
        '''
        x, y = self.x, self.y
        if y > 0:
            r = -_2_0 * atan2(x, hypot(y, x) + y) + PI_2
        elif y < 0:
            r =  _2_0 * atan2(x, hypot(y, x) - y) - PI_2
        else:  # y == 0
            r = PI if x < 0 else _0_0
        return r

    @Property_RO
    def latlon(self):
        '''Get the lat-, longitude in C{degrees} (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return self._xnamed(LatLon2Tuple(self.lat, self.lon))

    @Property_RO
    def latlonheight(self):
        '''Get the lat-, longitude in C{degrees} and height (L{LatLon3Tuple}C{(lat, lon, height)}).
        '''
        return self.latlon.to3Tuple(self.height)

    @Property_RO
    def latlonheightdatum(self):
        '''Get the lat-, longitude in C{degrees} with height and datum (L{LatLon4Tuple}C{(lat, lon, height, datum)}).
        '''
        return self.latlonheight.to4Tuple(self.datum)

    @Property_RO
    def latlonVermeille(self):
        '''Get the latitude and I{Vermeille} longitude in C{[-225..+225] degrees} (L{LatLon2Tuple}C{(lat, lon)}).

           @see: Property C{lonVermeille}.
        '''
        return self._xnamed(LatLon2Tuple(self.lat, degrees(self.lamVermeille)))

    @Property_RO
    def lonVermeille(self):
        '''Get the longitude in C{[-225..+225] degrees} after U{Vermeille
           <https://Search.ProQuest.com/docview/639493848>} (2004), p 95.

           @see: Property C{lamVermeille}.
        '''
        return degrees(self.lamVermeille)

    @Property_RO
    def phi(self):
        '''Get the latitude in C{radians} (C{float}).
        '''
        return self.philam.phi

    @Property_RO
    def philam(self):
        '''Get the lat-, longitude in C{radians} (L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        return self._xnamed(PhiLam2Tuple(radians(self.lat), radians(self.lon)))

    @Property_RO
    def philamheight(self):
        '''Get the lat-, longitude in C{radians} and height (L{PhiLam3Tuple}C{(phi, lam, height)}).
        '''
        return self.philam.to3Tuple(self.height)

    @Property_RO
    def philamheightdatum(self):
        '''Get the lat-, longitude in C{radians} with height and datum (L{PhiLam4Tuple}C{(phi, lam, height, datum)}).
        '''
        return self.philamheight.to4Tuple(self.datum)

    @Property_RO
    def philamVermeille(self):
        '''Get the latitude and I{Vermeille} longitude in C{[-PI*3/2..+PI*3/2] radians} (L{PhiLam2Tuple}C{(phi, lam)}).

           @see: Property C{lamVermeille}.
        '''
        return self._xnamed(PhiLam2Tuple(radians(self.lat), self.lamVermeille))

    def toCartesian(self, Cartesian, **Cartesian_kwds):
        '''Return the geocentric C{(x, y, z)} coordinates as an ellipsoidal or spherical
           C{Cartesian}.

           @arg Cartesian: L{ellipsoidalKarney.Cartesian}, L{ellipsoidalNvector.Cartesian},
                           L{ellipsoidalVincenty.Cartesian}, L{sphericalNvector.Cartesian} or
                           L{sphericalTrigonometry.Cartesian} class to return the C{(x, y, z)}
                           coordinates.
           @kwarg Cartesian_kwds: Optional B{C{Cartesian}} keyword arguments.

           @return: A B{C{Cartesian}}C{(x, y, z)} instance.

           @raise TypeError: Invalid B{C{Cartesian}} or B{C{Cartesian_kwds}}.
        '''
        from pygeodesy.cartesianBase import CartesianBase
        _xsubclassof(CartesianBase, Cartesian=Cartesian)
        r = Cartesian(self, **Cartesian_kwds)
        return self._xnamed(r)

    def toDatum(self, datum2):
        '''Convert this C{Ecef9Tuple} to an other datum.

           @arg datum2: Datum to convert I{to} (L{Datum}).

           @return: The converted 9-Tuple (C{Ecef9Tuple}).

           @raise TypeError: The B{C{datum2}} is not a L{Datum}.
        '''
        if self.datum in (None, datum2):  # PYCHOK _Names_
            r = self.copy()
        else:
            from pygeodesy.cartesianBase import CartesianBase
            c = CartesianBase(self, datum=self.datum)  # PYCHOK _Names_
            # c.toLatLon converts datum, x, y, z, lat, lon, etc.
            # and returns another Ecef9Tuple iff LatLon is None
            r = c.toLatLon(datum=datum2, LatLon=None)
        return self._xnamed(r)

    convertDatum = toDatum  # for backward compatibility

    def toLatLon(self, LatLon=None, **LatLon_height_datum_kwds):
        '''Return the geodetic C{(lat, lon, height[, datum])} coordinates.

           @kwarg LatLon: Optional class to return C{(lat, lon, height[,
                          datum])} or C{None}.
           @kwarg LatLon_height_datum_kwds: Optional B{C{height}}, B{C{datum}}
                                            and other B{C{LatLon}} keyword
                                            arguments.

           @return: An instance of C{LatLon}C{(lat, lon, **_height_datum_kwds)}
                    or if B{C{LatLon}} is C{None}, a L{LatLon3Tuple}C{(lat, lon,
                    height)} respectively L{LatLon4Tuple}C{(lat, lon, height,
                    datum)} depending on whether C{datum} is un-/available.

           @raise TypeError: Invalid B{C{LatLon}} or B{C{LatLon_height_datum_kwds}}.
        '''
        kwds = _xkwds(LatLon_height_datum_kwds, height=self.height, datum=self.datum)  # PYCHOK Ecef9Tuple
        d = kwds[_datum_]
        if LatLon is None:
            r = LatLon3Tuple(self.lat, self.lon, kwds[_height_])  # PYCHOK Ecef9Tuple
            if d:
                r = r.to4Tuple(d)  # checks type(d)
        else:
            if d is None:  # remove datum
                _ = kwds.pop[_datum_]
            r = LatLon(self.lat, self.lon, **kwds)  # PYCHOK Ecef9Tuple
        _datum_datum(getattr(r, _datum_, self.datum), self.datum)  # PYCHOK Ecef9Tuple
        return self._xnamed(r)

    def toVector(self, Vector=None, **Vector_kwds):
        '''Return the geocentric C{(x, y, z)} coordinates as vector.

           @kwarg Vector: Optional vector class to return C{(x, y, z)} or
                          C{None}.
           @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword
                               arguments, ignored if C{B{Vector}=None}.

           @return: A C{Vector}C{(x, y, z, **Vector_kwds)} instance or a
                    L{Vector3Tuple}C{(x, y, z)} if B{C{Vector}} is C{None}.

           @see: Propertes C{xyz} and C{xyzh}
        '''
        return self.xyz if Vector is None else \
              self._xnamed(Vector(self.x, self.y, self.z, **Vector_kwds))  # PYCHOK Ecef9Tuple

    @Property_RO
    def xyz(self):
        '''Get the geocentric C{(x, y, z)} coordinates (L{Vector3Tuple}C{(x, y, z)}).
        '''
        return self._xnamed(Vector3Tuple(self.x, self.y, self.z))

    @Property_RO
    def xyzh(self):
        '''Get the geocentric C{(x, y, z)} coordinates and height (L{Vector4Tuple}C{(x, y, z, h)})
        '''
        return self._xnamed(self.xyz.to4Tuple(self.height))


class EcefVeness(_EcefBase):
    '''Conversion between geodetic and geocentric, aka I{Earth-Centered,
       Earth-Fixed} (ECEF) coordinates transcribed from I{Chris Veness}'
       JavaScript classes U{LatLonEllipsoidal, Cartesian
       <https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}.

       @see: U{A Guide to Coordinate Systems in Great Britain
             <https://www.OrdnanceSurvey.co.UK/documents/resources/guide-coordinate-systems-great-britain.pdf>},
             section I{B) Converting between 3D Cartesian and ellipsoidal
             latitude, longitude and height coordinates}.
    '''

    def __init__(self, a_ellipsoid, f=None, name=NN):
        '''New L{EcefVeness}/L{EcefSudano} converter.

           @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}, L{Ellipsoid2}, L{Datum}
                             or L{a_f2Tuple}) or C{scalar} for the equatorial
                             radius of the ellipsoid (C{meter}).
           @kwarg f: C{None} or the ellipsoid flattening (C{scalar}), required
                     for C{scalar} B{C{a_ellipsoid}}, B{C{f=0}} represents a
                     sphere, negative B{C{f}} a prolate ellipsoid.
           @kwarg name: Optional name (C{str}).

           @raise EcefError: If B{C{a_ellipsoid}} not L{Ellipsoid}, L{Ellipsoid2},
                             L{Datum} or L{a_f2Tuple} or C{scalar} or B{C{f}} not
                             C{scalar} or if C{scalar} B{C{a_ellipsoid}} not positive
                             or B{C{f}} not less than 1.0.
        '''
        _EcefBase.__init__(self, a_ellipsoid, f, name)

    def forward(self, latlonh, lon=None, height=0, M=False):
        '''Convert from geodetic C{(lat, lon, height)} to geocentric C{(x, y, z)}.

           @arg latlonh: Either a C{LatLon}, an L{Ecef9Tuple} or C{scalar}
                         latitude in C{degrees}.
           @kwarg lon: Optional C{scalar} longitude in C{degrees} for C{scalar}
                       B{C{latlonh}}.
           @kwarg height: Optional height in C{meter}, vertically above (or
                          below) the surface of the ellipsoid.
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with geocentric C{(x, y, z)} coordinates for the given
                    geodetic ones C{(lat, lon, height)}, case C{C} 0,
                    L{EcefMatrix} C{M} and C{datum} if available.

           @raise EcefError: If B{C{latlonh}} not C{LatLon}, L{Ecef9Tuple} or
                             C{scalar} or B{C{lon}} not C{scalar} for C{scalar}
                             B{C{latlonh}} or C{abs(lat)} exceeds 90°.
        '''
        return _EcefBase._forward(self, latlonh, lon, height, M)

    def reverse(self, xyz, y=None, z=None, **no_M):  # PYCHOK unused M
        '''Convert from geocentric C{(x, y, z)} to geodetic C{(lat, lon, height)}
           transcribed from I{Chris Veness}' U{JavaScript
           <https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}.

           Uses B. R. Bowring’s formulation for μm precision in concise
           form: U{'The accuracy of geodetic latitude and height equations'
           <https://www.ResearchGate.net/publication/
           233668213_The_Accuracy_of_Geodetic_Latitude_and_Height_Equations>},
           Survey Review, Vol 28, 218, Oct 1985.

           @arg xyz: Either an L{Ecef9Tuple}, an C{(x, y, z)} 3-tuple or C{scalar}
                     ECEF C{x} coordinate in C{meter}.
           @kwarg y: ECEF C{y} coordinate in C{meter} for C{scalar} B{C{xyz}} and B{C{z}}.
           @kwarg z: ECEF C{z} coordinate in C{meter} for C{scalar} B{C{xyz}} and B{C{y}}.
           @kwarg no_M: Rotation matrix C{M} not available.

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with geodetic coordinates C{(lat, lon, height)} for the given
                    geocentric ones C{(x, y, z)}, case C{C}, L{EcefMatrix} C{M}
                    always C{None} and C{datum} if available.

           @raise EcefError: If B{C{xyz}} not L{Ecef9Tuple} or C{scalar} C{x}
                             or B{C{y}} and/or B{C{z}} not C{scalar} for C{scalar}
                             B{C{xyz}}.

           @see: Ralph M. Toms U{'An Efficient Algorithm for Geocentric to Geodetic
                 Coordinate Conversion'<https://www.OSTI.gov/scitech/biblio/110235>},
                 Sept 1995 and U{'An Improved Algorithm for Geocentric to Geodetic
                 Coordinate Conversion'<https://www.OSTI.gov/scitech/servlets/purl/231228>},
                 Apr 1996, both from Lawrence Livermore National Laboratory (LLNL) and
                 John J. Sudano U{An exact conversion from an Earth-centered coordinate
                 system to latitude longitude and altitude<https://www.ResearchGate.net/
                 publication/3709199_An_exact_conversion_from_an_Earth-centered_coordinate_system_to_latitude_longitude_and_altitude>}.
        '''
        x, y, z, name = _xyzn4(xyz, y, z, Error=EcefError)

        E = self.ellipsoid

        p = hypot(x, y)  # distance from minor axis
        r = hypot(p, z)  # polar radius
        if min(p, r) > EPS0:
            # parametric latitude (Bowring eqn 17, replaced)
            t = (E.b * z) / (E.a * p) * (_1_0 + E.e22 * E.b / r)
            c = _1_0 / hypot1(t)
            s = t * c

            # geodetic latitude (Bowring eqn 18)
            a = atan2(z + E.e22 * E.b * s**3,
                      p - E.e2  * E.a * c**3)

            # height above ellipsoid (Bowring eqn 7)
            sa, ca = sincos2(a)
#           r = E.a / E.e2s(sa)  # length of normal terminated by minor axis
#           h = p * ca + z * sa - (E.a * E.a / r)
            h = fsum_(p * ca, z * sa, -E.a * E.e2s(sa))

            C, lat, lon = 1, degrees90(a), atan2d(y, x)

        # see <https://GIS.StackExchange.com/questions/28446>
        elif p > EPS:  # lat arbitrarily zero
            C, lat, lon, h = 2, _0_0, atan2d(y, x), p - E.a

        else:  # polar lat, lon arbitrarily zero
            C, lat, lon, h = 3, copysign(_90_0, z), _0_0, abs(z) - E.b

        r = Ecef9Tuple(x, y, z, lat, lon, h, C, None, self.datum)
        return self._xnamed(r, name=name)


class EcefSudano(EcefVeness):
    '''Conversion between geodetic and geocentric, aka I{Earth-Centered,
       Earth-Fixed} (ECEF) coordinates based on I{John J. Sudano}'s U{paper
       <https://www.ResearchGate.net/publication/
       3709199_An_exact_conversion_from_an_Earth-centered_coordinate_system_to_latitude_longitude_and_altitude>}.
    '''

    if _FOR_DOCS:
        __init__ = EcefVeness.__init__
        forward  = EcefVeness.forward

    def reverse(self, xyz, y=None, z=None, **no_M):  # PYCHOK unused M
        '''Convert from geocentric C{(x, y, z)} to geodetic C{(lat, lon, height)} using
           I{Sudano}'s U{iterative method<https://www.ResearchGate.net/publication/
           3709199_An_exact_conversion_from_an_Earth-centered_coordinate_system_to_latitude_longitude_and_altitude>}.

           @arg xyz: Either an L{Ecef9Tuple}, an C{(x, y, z)} 3-tuple or C{scalar}
                     ECEF C{x} coordinate in C{meter}.
           @kwarg y: ECEF C{y} coordinate in C{meter} for C{scalar} B{C{xyz}} and B{C{z}}.
           @kwarg z: ECEF C{z} coordinate in C{meter} for C{scalar} B{C{xyz}} and B{C{y}}.
           @kwarg no_M: Rotation matrix C{M} not available.

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with geodetic coordinates C{(lat, lon, height)} for the given
                    geocentric ones C{(x, y, z)}, iteration C{C}, L{EcefMatrix} C{M}
                    always C{None} and C{datum} if available.

           @raise EcefError: If B{C{xyz}} not L{Ecef9Tuple} or C{scalar} C{x}
                             or B{C{y}} and/or B{C{z}} not C{scalar} for C{scalar}
                             B{C{xyz}} or no convergence.
        '''
        x, y, z, name = _xyzn4(xyz, y, z, Error=EcefError)

        E = self.ellipsoid
        e = E.e2 * E.a
        h = hypot(x, y)  # Rh
        d = e - h

        a = atan2(z, h * E.e12)
        sa, ca = sincos2(abs(a))
        # Sudano's Eq (A-6) and (A-7) refactored/reduced,
        # replacing Rn from Eq (A-4) with n = E.a / ca:
        # N = ca**2 * ((z + E.e2 * n * sa) * ca - h * sa)
        #   = ca**2 * (z * ca + E.e2 * E.a * sa - h * sa)
        #   = ca**2 * (z * ca + (E.e2 * E.a - h) * sa)
        # D = ca**3 * (E.e2 * n / E.e2s2(sa)) - h
        #   = ca**2 * (E.e2 * E.a / E.e2s2(sa) - h / ca**2)
        # N / D = (z * ca + (E.e2 * E.a - h) * sa) /
        #         (E.e2 * E.a / E.e2s2(sa) - h / ca**2)
        S = Fsum(sa)
        for C in range(1, _TRIPS):
            ca2 = _1_0 - sa**2
            if ca2 < EPS_2:  # PYCHOK no cover
                ca = _0_0
                break
            ca = sqrt(ca2)
            t = e / E.e2s2(sa) - h / ca2
            if abs(t) < EPS_2:
                break
            a = None
            sa, t = S.fsum2_(-(z * ca + d * sa) / t)
            if abs(t) < EPS:
                break
        else:
            raise EcefError(unstr(self.reverse.__name__, x=x, y=y, z=z), txt=_no_(_convergence_))

        if a is None:
            a = copysign(asin(sa), z)
        h = fsum_(h * ca, abs(z * sa), -E.a * E.e2s(sa))  # use Veness',
        # since Sudano's Eq (7) doesn't provide the correct height
        # h = (abs(z) + h - E.a * cos(a + E.e12) * sa / ca) / (ca + sa)

        r = Ecef9Tuple(x, y, z, degrees90(a), atan2d(y, x), h, C, None, self.datum)
        r._iteration = C
        return self._xnamed(r, name=name)


class EcefYou(_EcefBase):
    '''Conversion between geodetic and geocentric, aka I{Earth-Centered,
       Earth-Fixed} (ECEF) coordinates using I{Rey-Jer You}'s U{transformations
       <https://www.ResearchGate.net/publication/240359424>}.

       @see: W.E. Featherstone, S.J. (Sten) Claessens U{Closed-form transformation
             between geodetic and ellipsoidal coordinates
             <https://espace.Curtin.edu.AU/bitstream/handle/20.500.11937/11589/115114_9021_geod2ellip_final.pdf>}
             Studia Geophysica et Geodaetica, 2008, 52, 1-18 and U{PyMap3D
             <https://PyPI.org/project/pymap3d>}.
    '''

    def __init__(self, a_ellipsoid, f=None, name=NN):
        '''New L{EcefYou} converter.

           @arg a_ellipsoid: A (non-prolate) ellipsoid (L{Ellipsoid}, L{Ellipsoid2},
                             L{Datum} or L{a_f2Tuple}) or C{scalar} for the equatorial
                             radius of the ellipsoid (C{meter}).
           @kwarg f: C{None} or the ellipsoid flattening (C{scalar}), required for
                     C{scalar} B{C{a_ellipsoid}}, use B{C{f=0}} to represent a sphere.
           @kwarg name: Optional name (C{str}).

           @raise EcefError: If B{C{a_ellipsoid}} not L{Ellipsoid}, L{Ellipsoid2},
                             L{Datum} or L{a_f2Tuple} or C{scalar} or B{C{f}} not
                             C{scalar} or if C{scalar} B{C{a_ellipsoid}} not positive
                             or B{C{f}} negative or not less than 1.0.
        '''
        _EcefBase.__init__(self, a_ellipsoid, f, name)
        E = self.ellipsoid
        if E.isProlate or (E.a2 - E.b2) < 0:
            raise EcefError(ellipsoid=E, txt=_prolate_)

    def forward(self, latlonh, lon=None, height=0, M=False):
        '''Convert from geodetic C{(lat, lon, height)} to geocentric C{(x, y, z)}.

           @arg latlonh: Either a C{LatLon}, an L{Ecef9Tuple} or C{scalar}
                         latitude in C{degrees}.
           @kwarg lon: Optional C{scalar} longitude in C{degrees} for C{scalar}
                       B{C{latlonh}}.
           @kwarg height: Optional height in C{meter}, vertically above (or
                          below) the surface of the ellipsoid.
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with geocentric C{(x, y, z)} coordinates for the given
                    geodetic ones C{(lat, lon, height)}, case C{C} 0,
                    L{EcefMatrix} C{M} and C{datum} if available.

           @raise EcefError: If B{C{latlonh}} not C{LatLon}, L{Ecef9Tuple} or
                             C{scalar} or B{C{lon}} not C{scalar} for C{scalar}
                             B{C{latlonh}} or C{abs(lat)} exceeds 90°.
        '''
        return _EcefBase._forward(self, latlonh, lon, height, M, You=True)

    def reverse(self, xyz, y=None, z=None, **no_M):  # PYCHOK unused M
        '''Convert from geocentric C{(x, y, z)} to geodetic C{(lat, lon, height)}
           using I{Rey-Jer You}'s transformation.

           @arg xyz: Either an L{Ecef9Tuple}, an C{(x, y, z)} 3-tuple or C{scalar}
                     ECEF C{x} coordinate in C{meter}.
           @kwarg y: ECEF C{y} coordinate in C{meter} for C{scalar} B{C{xyz}}
                     and B{C{z}}.
           @kwarg z: ECEF C{z} coordinate in C{meter} for C{scalar} B{C{xyz}}
                     and B{C{y}}.
           @kwarg no_M: Rotation matrix C{M} not available.

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with geodetic coordinates C{(lat, lon, height)} for the given
                    geocentric ones C{(x, y, z)}, case C{C} 1, L{EcefMatrix} C{M}
                    always C{None} and C{datum} if available.

           @raise EcefError: If B{C{xyz}} not L{Ecef9Tuple} or C{scalar} C{x}
                             or B{C{y}} and/or B{C{z}} not C{scalar} for C{scalar}
                             B{C{xyz}}.
        '''
        x, y, z, name = _xyzn4(xyz, y, z, Error=EcefError)

        r2 = hypot2_(x, y, z)

        E  = self.ellipsoid
        e2 = E.a2 - E.b2  # == E.e2 * E.a2
        if e2 < 0:
            raise EcefError(ellipsoid=E, txt=_prolate_)
        e = sqrt(e2)  # XXX sqrt0(e2)?

        q = hypot(x, y)
        u = fsum_(r2, -e2, hypot(r2 - e2, 2 * e * z)) * _0_5
        if u > _EPS0__2:
            u = sqrt(u)
            p = hypot(u, e)
            B = atan2(p * z, u * q)  # beta0 = atan(p / u * z / q)
            sB, cB = sincos2(B)
            if cB and sB:
                p *= E.a
                d  = (p / cB - e2 * cB) / sB
                if abs(d) > EPS0:
                    B += fsum_(u * E.b, -p, e2) / d
                    sB, cB = sincos2(B)
        elif u < 0:
            raise EcefError(x=x, y=y, z=z, txt=_singular_)
        else:
            sB, cB = copysign(_1_0, z), _0_0

        h = hypot(z - E.b * sB, q - E.a * cB)
        if hypot2_(x, y, z * E.a_b) < E.a2:
            h = neg(h)  # inside ellipsoid

        r = Ecef9Tuple(x, y, z, atan2d(E.a * sB, E.b * cB),  # atan(E.a_b * tan(B))
                                atan2d(y, x), h, 1,  # C=1
                                None,  # M=None
                                self.datum)
        return self._xnamed(r, name=name)


__all__ += _ALL_DOCS(_EcefBase)

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
