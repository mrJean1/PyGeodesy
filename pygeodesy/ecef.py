
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

Convert between geodetic coordinates C{lat}-, C{lon}gitude and height C{h}
(measured vertically from the surface of the ellipsoid) to geocentric C{x},
C{y} and C{z} coordinates, also known as I{Earth-Centered, Earth-Fixed}
(U{ECEF<https://WikiPedia.org/wiki/ECEF>}).

The origin of geocentric coordinates is at the center of the earth.  The C{z}
axis goes thru the North pole, C{lat} = 90°.  The C{x} axis goes thru C{lat}
= 0°, C{lon} = 0°.

The origin of local cartesian coordinate system is at C{lat0}, C{lon0} and
C{height0}.  The C{z} axis is normal to the ellipsoid, the C{y} axis points due
North.  The plane C{z = -height0} is tangent to the ellipsoid.

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

from pygeodesy.basics import EPS, EPS1, EPS_2, isscalar, map1, property_RO, \
                            _xinstanceof, _xkwds, _xsubclassof
from pygeodesy.datum import Datum, Datums, Ellipsoid
from pygeodesy.errors import _datum_datum, LenError, _ValueError
from pygeodesy.fmath import cbrt, fdot, fsum_, hypot1
from pygeodesy.interns import _C_, _datum_, _ellipsoid_, _h_, _height_, _lat_, \
                              _lat0_, _lon_, _lon0_, _name_, NN, _no_convergence_, \
                              _SPACE_, _UNDERSCORE_, _x_, _y_, _z_, _0_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _FOR_DOCS
from pygeodesy.named import LatLon2Tuple, LatLon3Tuple, _NamedBase, _NamedTuple, \
                            notOverloaded, PhiLam2Tuple, Vector3Tuple, _xnamed
from pygeodesy.streprs import unstr
from pygeodesy.utily import degrees90, degrees180, sincos2, sincos2d
from pygeodesy.vector3d import _xyzn4

from math import asin, atan2, copysign, cos, degrees, hypot, radians, sqrt

__all__ = _ALL_LAZY.ecef
__version__ = '20.08.04'

_M_ = 'M'


def _llhn4(latlonh, lon, height, suffix=NN):
    '''(INTERNAL) Get C{lat, lon, h, name} as C{4-tuple}.
    '''
    try:
        llh = latlonh.lat, latlonh.lon, getattr(latlonh, _height_,
                                        getattr(latlonh, _h_, height))
    except AttributeError:
        llh = latlonh, lon, height
    try:
        lat, lon, h = map1(float, *llh)
    except (TypeError, ValueError) as x:
        t = _lat_, _lon_, _height_
        if suffix:
            t = (_ + suffix for _ in t)
        raise EcefError(txt=str(x), **dict(zip(t, llh)))

    if abs(lat) > 90:  # XXX RangeError
        raise EcefError(_lat_ + suffix, lat)

    return lat, lon, h, getattr(latlonh, _name_, NN)


def _sch3(y, x):
    '''(INTERNAL) Compute sin, cos and hypotenuse.
    '''
    h = hypot(y, x)
    if h > 0:  # EPS_2
        s, c = y / h, x / h
    else:
        s, c = 0, 1
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
        try:
            E = a_ellipsoid
            if f is None:
                if isinstance(E, Datum):
                    self._datum = E
                    E = E.ellipsoid
                elif not isinstance(E, Ellipsoid):
                    raise TypeError
                if not name:
                    name = E.name

            elif isscalar(E) and isscalar(f):
                a  = float(E)
                f_ = (1.0 / f) if f else 0  # sphere
                b  = None if f_ else a
                E  = Ellipsoid(a, b, f_, name=_UNDERSCORE_ + name)

            else:
                raise ValueError

            if not (E.a > 0 and E.f < 1):
                raise ValueError

        except (TypeError, ValueError) as x:
            t = unstr(self.classname, a=a_ellipsoid, f=f)
            raise EcefError(t + _SPACE_ + _ellipsoid_, txt=str(x))

        self._E = E
        if name:
            self.name = name

    @property_RO
    def a(self):
        '''Get C{E.a}, the major, equatorial radius (C{meter}).
        '''
        return self._E.a

    equatorialRadius = a  # Karney property

    @property_RO
    def datum(self):
        '''Get the datum (L{Datum} or C{None} if not available).
        '''
        return self._datum

    @property_RO
    def ellipsoid(self):
        '''Get C{E}, the ellipsoid (L{Ellipsoid}).
        '''
        return self._E

    @property_RO
    def f(self):
        '''Get C{E.f}, the flattening (C{float}), zero for sphere.
        '''
        return self._E.f

    flattening = f  # Karney property

    def forward(self, latlonh, lon=None, height=0, M=False):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.
        '''
        notOverloaded(self, self.forward, latlonh, lon=lon,
                                      height=height, M=M)

    def _Matrix(self, sa, ca, sb, cb):
        '''Creation a rotation matrix.

           @arg sa: C{sin(phi)} (C{float}).
           @arg ca: C{cos(phi)} (C{float}).
           @arg sb: C{sin(lambda)} (C{float}).
           @arg cb: C{cos(lambda)} (C{float}).

           @return: A L{EcefMatrix}.
        '''
        return self._xnamed(EcefMatrix(sa, ca, sb, cb))

    def reverse(self, xyz, y=None, z=None, M=False):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.
        '''
        notOverloaded(self, self.reverse, xyz, y=y, z=z, M=M)

    def toStr(self, prec=9):  # PYCHOK signature
        '''Return this C{Ecef*} as a string.

           @kwarg prec: Optional precision, number of decimal digits (0..9).

           @return: This C{Ecef*} representation (C{str}).
        '''
        return self.attrs('a', 'f', _datum_, _ellipsoid_, _name_, prec=prec)

#   @property_RO
#   def xyz(self):
#       '''Get the geocentric C{(x, y, z)} coordinates (L{Vector3Tuple}C{(x, y, z)}).
#       '''
#       return self._xnamed(Vector3Tuple(self.x, self.y, self.z))


class EcefKarney(_EcefBase):
    '''Conversion between geodetic and geocentric, aka I{Earth-Centered,
       Earth-Fixed} (ECEF) coordinates based on I{Karney}'s U{Geocentric
       <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Geocentric.html>}
       methods.
    '''
    _hmax = 0  # 12M light years

    def __init__(self, a_ellipsoid, f=None, name=NN):
        '''New L{EcefKarney} converter.

           @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}), a datum (L{Datum}) or
                             C{scalar} for the major, equatorial radius of the
                             ellipsoid (C{meter}).
           @kwarg f: C{None} or the ellipsoid flattening (C{scalar}), required
                     for C{scalar} B{C{a_datum_ellipsoid}}, B{C{f=0}} represents
                     a sphere, negative B{C{f}} a prolate ellipsoid.
           @kwarg name: Optional name (C{str}).

           @raise EcefError: If B{C{a_ellipsoid}} not L{Ellipsoid}, L{Datum} or
                             C{scalar} or B{C{f}} not C{scalar} or if C{scalar}
                             B{C{a_ellipsoid}} not positive or B{C{f}} not less
                             than 1.0.
        '''
        _EcefBase.__init__(self, a_ellipsoid, f, name)
        self._hmax = self.a / EPS_2

    @property_RO
    def e2(self):
        '''Get C{E.f * (2 - E.f)}, 1st eccentricty squared (C{float}).
        '''
        return self._E.e2

    @property_RO
    def e2a(self):
        '''Get C{abs(E.e2)} (C{float}).
        '''
        return abs(self._E.e2)

    @property_RO
    def e2m(self):
        '''Get C{1 - E.e2a} == C{(1 - E.f)**2} (C{float}).
        '''
        return self._E.e12  # == (1 - E.f)**2

    @property_RO
    def e4a(self):
        '''Get C{E.e2a**2} (C{float}).
        '''
        return self._E.e4

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
        lat, lon, h, name = _llhn4(latlonh, lon, height)
        sa, ca, sb, cb = sincos2d(lat, lon)

        n = self.a / self.ellipsoid.e2s(sa)  # ... / sqrt(1 - self.e2 * sa**2)
        z = (self.e2m * n + h) * sa
        x = (n + h) * ca

        m = self._Matrix(sa, ca, sb, cb) if M else None
        r = Ecef9Tuple(x * cb, x * sb, z, lat, lon, h, 0, m, self.datum)
        return self._xnamed(r, name)

    @property_RO
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
        x, y, z, name = _xyzn4(xyz, y, z, Error=EcefError)

        sb, cb, d = _sch3(y, x)
        h = hypot(d, z)  # distance to earth center
        if h > self._hmax:  # PYCHOK no cover
            # We are really far away (> 12M light years).  Treat the earth
            # as a point and h, above, is an acceptable approximation to the
            # height.  This avoids overflow, e.g., in the computation of d
            # below.  It's possible that h has overflowed to INF, that's OK.
            # Treat finite x, y, but r overflows to +INF by scaling by 2.
            sb, cb, r = _sch3(y / 2, x / 2)
            sa, ca, _ = _sch3(z / 2, r)
            C = 1

        elif self.e4a:
            # Treat prolate spheroids by swapping R and Z here and by switching
            # the arguments to phi = atan2(...) at the end.
            p = (d / self.a)**2
            q = self.e2m * (z / self.a)**2
            if self.f < 0:
                p, q = q, p
            r = p + q - self.e4a
            e = self.e4a * q
            if e or r > 0:
                # Avoid possible division by zero when r = 0 by multiplying
                # equations for s and t by r^3 and r, resp.
                s = e * p / 4  # s = r^3 * s
                u = r = r / 6
                r2 = r**2
                r3 = r * r2
                t3 = s + r3
                disc = s * (r3 + t3)
                if disc < 0:
                    # t is complex, but the way u is defined the result is real.
                    # There are three possible cube roots.  We choose the root
                    # which avoids cancellation.  Note, disc < 0 implies r < 0.
                    u += 2 * r * cos(atan2(sqrt(-disc), -t3) / 3)
                else:
                    # Pick the sign on the sqrt to maximize abs(T3).  This
                    # minimizes loss of precision due to cancellation.  The
                    # result is unchanged because of the way the T is used
                    # in definition of u.
                    if disc > 0:
                        t3 += copysign(sqrt(disc), t3)  # t3 = (r * t)^3
                    # N.B. cbrt always returns the real root, cbrt(-8) = -2.
                    t = cbrt(t3)  # t = r * t
                    # t can be zero; but then r2 / t -> 0.
                    if t:
                        u += t + r2 / t
                v = sqrt(u**2 + e)  # guaranteed positive
                # Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
                # E.e4 * q / (v - u) because u ~ e^4 when q is small and u < 0.
                uv = (e / (v - u)) if u < 0 else (u + v)  # u+v, guaranteed positive
                # Need to guard against w going negative due to roundoff in uv - q.
                w = max(0.0, self.e2a * (uv - q) / (2 * v))
                # Rearrange expression for k to avoid loss of accuracy due to
                # subtraction.  Division by 0 not possible because uv > 0, w >= 0.
                k1 = k2 = uv / (sqrt(uv + w**2) + w)
                if self.f < 0:
                    k1 -= self.e2
                else:
                    k2 += self.e2
                sa, ca, h = _sch3(z / k1, d / k2)
                h *= k1 - self.e2m
                C = 2

            else:  # e = self.e4a * q == 0 and r <= 0
                # This leads to k = 0 (oblate, equatorial plane) and k + E.e^2 = 0
                # (prolate, rotation axis) and the generation of 0/0 in the general
                # formulas for phi and h, using the general formula and division
                # by 0 in formula for h.  Handle this case by taking the limits:
                #   f > 0: z -> 0, k        ->  E.e2 * sqrt(q) / sqrt(E.e4 - p)
                #   f < 0: r -> 0, k + E.e2 -> -E.e2 * sqrt(q) / sqrt(E.e4 - p)
                q = self.e4a - p
                if self.f < 0:
                    p, q, e = q, p, -1
                else:
                    e = -self.e2m
                sa, ca, h = _sch3(sqrt(q / self.e2m), sqrt(p))
                h *= self.a * e / self.e2a
                if z < 0:
                    sa = -sa  # for tiny negative z, not for prolate
                C = 3

        else:  # self.e4a == 0
            # Treat the spherical case.  Dealing with underflow in the general
            # case with E.e2 = 0 is difficult.  Origin maps to North pole, same
            # as with ellipsoid.
            sa, ca, _ = _sch3(z if h else 1, d)
            h -= self.a
            C = 4

        r = Ecef9Tuple(x, y, z, degrees(atan2(sa, ca)),
                                degrees(atan2(sb, cb)), h, C,
                                self._Matrix(sa, ca, sb, cb) if M else None,
                                self.datum)
        return self._xnamed(r, name)


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
    _t0   = None

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

    @property_RO
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
        return self._xnamed(r, name)

    @property_RO
    def height0(self):
        '''Get origin's height (C{meter}).
        '''
        return self._t0.height

    @property_RO
    def lat0(self):
        '''Get origin's latitude (C{degrees}).
        '''
        return self._t0.lat

    @property_RO
    def lon0(self):
        '''Get origin's longitude (C{degrees}).
        '''
        return self._t0.lon

    @property_RO
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
        lat0, lon0, height0, n = _llhn4(latlonh0, lon0, height0, suffix=_0_)
        if name or n:
            self.ecef.name = self.name = name or n
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
        return self._xnamed(r, xyz_n[3])

    def toStr(self, prec=9):  # PYCHOK signature
        '''Return this L{EcefCartesian} as a string.

           @kwarg prec: Optional precision, number of decimal digits (0..9).

           @return: This L{EcefCartesian} representation (C{str}).
        '''
        return self.attrs(_lat0_, _lon0_, 'height0', _M_, 'ecef', _name_, prec=prec)


class EcefMatrix(_NamedTuple):
    '''A rotation matrix.
    '''
    _Names_ = ('_0_0', '_0_1', '_0_2',  # row-order
               '_1_0', '_1_1', '_1_2',
               '_2_0', '_2_1', '_2_2')

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
            raise EcefError(EcefMatrix.__name__, t)

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
                   0,       ca,      sa)

        return _NamedTuple.__new__(cls, *t)

    def copy(self, **unused):  # PYCHOK signature
        '''Make a shallow or deep copy of this instance.

           @return: The copy (C{This class} or subclass thereof).
        '''
        return _xnamed(self.classof(*self), self.name)

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
        '''
        if xyz0:
            if len(xyz0) != len(xyz):
                raise LenError(self.unrotate, xyz0=len(xyz0), xyz=len(xyz))

            _xyz = (1.0,) + xyz
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
    _Names_ = (_x_, _y_, _z_, _lat_, _lon_, _height_, _C_, _M_, _datum_)

    def convertDatum(self, datum2):
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

    @property_RO
    def lam(self):
        '''Get the longitude in C{radians} (C{float}).
        '''
        return self.philam.lam

    @property_RO
    def latlon(self):
        '''Get the lat-, longitude in C{degrees} (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return self._xnamed(LatLon2Tuple(self.lat, self.lon))

    @property_RO
    def latlonheight(self):
        '''Get the lat-, longitude in C{degrees} and height (L{LatLon3Tuple}C{(lat, lon, height)}).
        '''
        return self._xnamed(self.latlon.to3Tuple(self.height))

    @property_RO
    def latlonheightdatum(self):
        '''Get the lat-, longitude in C{degrees} with height and datum (L{LatLon4Tuple}C{(lat, lon, height, datum)}).
        '''
        return self._xnamed(self.latlon.to3Tuple(self.height).to4Tuple(self.datum))

    @property_RO
    def phi(self):
        '''Get the latitude in C{radians} (C{float}).
        '''
        return self.philam.phi

    @property_RO
    def philam(self):
        '''Get the lat-, longitude in C{radians} (L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        return self._xnamed(PhiLam2Tuple(radians(self.lat), radians(self.lon)))

    @property_RO
    def philamheight(self):
        '''Get the lat-, longitude in C{radians} and height (L{PhiLam3Tuple}C{(phi, lam, height)}).
        '''
        return self._xnamed(self.philam.to3Tuple(self.height))

    @property_RO
    def philamheightdatum(self):
        '''Get the lat-, longitude in C{radians} with height and datum (L{PhiLam4Tuple}C{(phi, lam, height, datum)}).
        '''
        return self._xnamed(self.philam.to3Tuple(self.height).to4Tuple(self.datum))

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
                               arguments, ignored if B{C{Vector=None}}.

           @return: A C{Vector}C{(x, y, z, **Vector_kwds)} instance or a
                    L{Vector3Tuple}C{(x, y, z)} if B{C{Vector}} is C{None}.

           @see: Propertes C{xyz} and C{xyzh}
        '''
        return self.xyz if Vector is None else \
              self._xnamed(Vector(self.x, self.y, self.z, **Vector_kwds))  # PYCHOK Ecef9Tuple

    @property_RO
    def xyz(self):
        '''Get the geocentric C{(x, y, z)} coordinates (L{Vector3Tuple}C{(x, y, z)}).
        '''
        return self._xnamed(Vector3Tuple(self.x, self.y, self.z))

    @property_RO
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

           @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}), a datum (L{Datum}) or
                             C{scalar} for the major, equatorial radius of the
                             ellipsoid (C{meter}).
           @kwarg f: C{None} or the ellipsoid flattening (C{scalar}), required
                     for C{scalar} B{C{a_datum_ellipsoid}}, B{C{f=0}} represents
                     a sphere, negative B{C{f}} a prolate ellipsoid.
           @kwarg name: Optional name (C{str}).

           @raise EcefError: If B{C{a_ellipsoid}} not L{Ellipsoid}, L{Datum} or
                             C{scalar} or B{C{f}} not C{scalar} or if C{scalar}
                             B{C{a_ellipsoid}} not positive or B{C{f}} not less
                             than 1.0.
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
        lat, lon, h, name = _llhn4(latlonh, lon, height)
        sa, ca, sb, cb = sincos2d(lat, lon)

        E = self.ellipsoid

        # radius of curvature in prime vertical
        t = E.e2s2(sa)  # r, _ = E.roc2_(sa, 1)
        r = 0 if t < EPS else (E.a if t > EPS1 else (E.a / sqrt(t)))
        x = (h + r) * ca
        z = (h + r * E.e12) * sa

        m = self._Matrix(sa, ca, sb, cb) if M else None
        r = Ecef9Tuple(x * cb, x * sb, z, lat, lon, h, 0, m, self.datum)
        return self._xnamed(r, name)

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
        if min(p, r) > EPS:
            # parametric latitude (Bowring eqn 17, replaced)
            t = (E.b * z) / (E.a * p) * (1 + E.e22 * E.b / r)
            c = 1 / hypot1(t)
            s = t * c

            # geodetic latitude (Bowring eqn 18)
            a = atan2(z + E.e22 * E.b * s**3,
                      p - E.e2  * E.a * c**3)
            b = atan2(y, x)  # ... and longitude

            # height above ellipsoid (Bowring eqn 7)
            sa, ca = sincos2(a)
#           r = E.a / E.e2s(sa)  # length of normal terminated by minor axis
#           h = p * ca + z * sa - (E.a * E.a / r)
            h = fsum_(p * ca, z * sa, -E.a * E.e2s(sa))

            C, lat, lon = 1, degrees90(a), degrees180(b)

        # see <https://GIS.StackExchange.com/questions/28446>
        elif p > EPS:  # lat arbitrarily zero
            C, lat, lon, h = 2, 0.0, degrees180(atan2(y, x)), p - E.a

        else:  # polar lat, lon arbitrarily zero
            C, lat, lon, h = 3, copysign(90.0, z), 0.0, abs(z) - E.b

        r = Ecef9Tuple(x, y, z, lat, lon, h, C, None, self.datum)
        return self._xnamed(r, name)


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
           I{John Sudano}'s U{iterative method<https://www.ResearchGate.net/publication/
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
        for i in range(16):  # 8..9 is sufficient
            ca2 = 1 - sa**2
            if ca2 < EPS_2:
                ca = 0
                break
            ca = sqrt(ca2)
            t = e / E.e2s2(sa) - h / ca2
            if abs(t) < EPS_2:
                break
            t = (z * ca + d * sa) / t
            if abs(t) < EPS:
                break
            sa -= t
        else:
            raise EcefError(_no_convergence_, txt=unstr(self.reverse.__name__, x=x, y=y, z=z))

        if i:
            a = copysign(asin(sa), z)
        b = atan2(y, x)
        h = fsum_(h * ca, abs(z * sa), -E.a * E.e2s(sa))  # use Veness',
        # Sudano's Eq (7) doesn't seem to provide the correct height
        # h = (abs(z) + h - E.a * cos(a + E.e12) * sa / ca) / (ca + sa)

        r = Ecef9Tuple(x, y, z, degrees90(a), degrees180(b), h, i, None, self.datum)
        return self._xnamed(r, name)


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

           @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}), a datum (L{Datum}) or
                             C{scalar} for the major, equatorial radius of the
                             ellipsoid (C{meter}).
           @kwarg f: C{None} or the ellipsoid flattening (C{scalar}), required
                     for C{scalar} B{C{a_datum_ellipsoid}}, B{C{f=0}} represents
                     a sphere, negative B{C{f}} a prolate ellipsoid.
           @kwarg name: Optional name (C{str}).

           @raise EcefError: If B{C{a_ellipsoid}} not L{Ellipsoid}, L{Datum} or
                             C{scalar} or B{C{f}} not C{scalar} or if C{scalar}
                             B{C{a_ellipsoid}} not positive or B{C{f}} not less
                             than 1.0.
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
        lat, lon, h, name = _llhn4(latlonh, lon, height)
        sa, ca, sb, cb = sincos2d(lat, lon)

        E = self.ellipsoid
        n = 1.0 / hypot(E.a * ca, E.b * sa)
        z = (n * E.b2 + h) * sa
        x = (n * E.a2 + h) * ca

        m = self._Matrix(sa, ca, sb, cb) if M else None
        r = Ecef9Tuple(x * cb, x * sb, z, lat, lon, h, 0, m, self.datum)
        return self._xnamed(r, name)

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

        x2, y2, z2 = x**2, y**2, z**2
        r2 = fsum_(x2, y2, z2)  # = hypot3(x2, y2, z2)**2

        E = self.ellipsoid
        e = sqrt(E.a2 - E.b2)
        e2 = e**2

        u = sqrt(fsum_(r2, -e2, hypot(r2 - e2, 2 * e * z)) / 2)
        p = hypot(u, e)
        q = hypot(x, y)
        B = atan2(p * z, u * q)  # beta0 = atan(p / u * z / q)
        sB, cB = sincos2(B)
        p *= E.a
        B += fsum_(u * E.b, -p, e2) * sB / (p / cB - e2 * cB)
        sB, cB = sincos2(B)

        h = hypot(z - E.b * sB, q - E.a * cB)
        if fsum_(x2, y2, z2 * E.a_b**2) < E.a2:
            h = -h  # inside ellipsoid

        r = Ecef9Tuple(x, y, z, degrees(atan2(E.a * sB, E.b * cB)),  # atan(E.a_b * tan(B))
                                degrees(atan2(y, x)), h, 1,  # C=1
                                None,  # M=None
                                self.datum)
        return self._xnamed(r, name)


__all__ += _ALL_DOCS(_EcefBase, Ecef9Tuple)

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
