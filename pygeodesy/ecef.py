
# -*- coding: utf-8 -*-

u'''I{Geocentric} Earth-Centered, Earth-Fixed (ECEF) coordinates.

Geocentric conversions transcoded from I{Charles Karney}'s C++ class U{Geocentric
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Geocentric.html>}
into pure Python class L{EcefKarney}, class L{EcefSudano} based on I{John Sudano}'s
U{paper<https://www.ResearchGate.net/publication/
3709199_An_exact_conversion_from_an_Earth-centered_coordinate_system_to_latitude_longitude_and_altitude>},
class L{EcefVeness} transcoded from I{Chris Veness}' JavaScript classes U{LatLonEllipsoidal,
Cartesian<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}, class L{EcefYou}
implementing I{Rey-Jer You}'s U{transformations <https://www.ResearchGate.net/publication/240359424>} and
classes L{EcefFarrell22} and L{EcefFarrell22} from I{Jay A. Farrell}'s U{Table 2.1 and 2.2
<https://Books.Google.com/books?id=fW4foWASY6wC>}, page 29-30.

Following is a copy of I{Karney}'s U{Detailed Description
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Geocentric.html>}.

Convert between geodetic coordinates C{lat}-, C{lon}gitude and height C{h}
(measured vertically from the surface of the ellipsoid) to geocentric C{x},
C{y} and C{z} coordinates, also known as I{Earth-Centered, Earth-Fixed}
(U{ECEF<https://WikiPedia.org/wiki/ECEF>}).

The origin of geocentric coordinates is at the center of the earth.  The C{z}
axis goes thru the North pole, C{lat} = 90°.  The C{x} axis goes thru C{lat}
= 0°, C{lon} = 0°.

The I{local (cartesian) origin} is at (C{lat0}, C{lon0}, C{height0}).  The I{local} C{x}
axis points East, the I{local} C{y} axis points North and the I{local} C{z} axis is
normal to the ellipsoid.  The plane C{z = -height0} is tangent to the ellipsoid, hence
the alternate name I{local tangent plane}.

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

@see: Module L{ltp} and class L{LocalCartesian}, a transcription of I{Charles Karney}'s
C++ class U{LocalCartesian
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1LocalCartesian.html>},
providing conversion to and from I{local} cartesian cordinates in a I{local tangent
plane} as opposed to I{geocentric} (ECEF) ones.
'''

from pygeodesy.basics import copysign0, isnon0, isscalar, issubclassof, \
                             neg, map1, _xinstanceof, _xsubclassof
from pygeodesy.datums import _ellipsoidal_datum
from pygeodesy.ellipsoids import a_f2Tuple
from pygeodesy.errors import _datum_datum, _IndexError, LenError, \
                             _ValueError, _TypesError, _xkwds
from pygeodesy.fmath import cbrt, fdot, Fsum, fsum_, hypot, hypot1, hypot2_
from pygeodesy.interns import EPS, EPS0, EPS02, EPS1, EPS_2, NN, PI, PI_2, \
                             _a_, _C_, _convergence_, _datum_, _ellipsoid_, \
                             _f_, _h_, _height_, _lat_, _lon_, _M_, _name_, \
                             _no_, _singular_, _SPACE_, _x_, _xyz_, _y_, _z_, \
                             _0_0, _0_5, _1_0, _1_0_T, _2_0, _3_0, _4_0, \
                             _6_0, _90_0
from pygeodesy.interns import _N_2_0  # PYCHOK used!
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import _NamedBase, _NamedTuple, notOverloaded, \
                            _Pass, _xnamed
from pygeodesy.namedTuples import LatLon2Tuple, LatLon3Tuple, \
                                  PhiLam2Tuple, Vector3Tuple, Vector4Tuple
from pygeodesy.props import deprecated_method, Property_RO
from pygeodesy.streprs import unstr
from pygeodesy.units import Height, Int, Lam, Lat, Lon, Meter, Phi, Scalar
from pygeodesy.utily import atan2d, degrees90, degrees180, \
                            sincos2, sincos2_, sincos2d_

from math import asin, atan2, cos, degrees, radians, sqrt

__all__ = _ALL_LAZY.ecef
__version__ = '21.11.18'

_Ecef_    = 'Ecef'
_prolate_ = 'prolate'
_TRIPS    =  17  # 8..9 sufficient, EcefSudano.reverse


class EcefError(_ValueError):
    '''An ECEF or C{Ecef*} related issue.
    '''
    pass


def _llhn4(latlonh, lon, height, suffix=NN, Error=EcefError, name=NN):  # in .ltp.LocalCartesian.forward and -.reset
    '''(INTERNAL) Get C{lat, lon, h, name} as C{4-tuple}.
    '''
    try:
        lat, lon = latlonh.lat, latlonh.lon
        h = getattr(latlonh, _height_,
            getattr(latlonh, _h_, height))
        n = getattr(latlonh, _name_, NN)
    except AttributeError:
        lat, h, n = latlonh, height, NN

    try:
        llhn = Lat(lat), Lon(lon), Height(h), (name or n)
    except (TypeError, ValueError) as x:
        t = _lat_, _lon_, _height_
        if suffix:
            t = (_ + suffix for _ in t)
        d = dict(zip(t, (lat, lon, h)))
        raise Error(txt=str(x), **d)
    return llhn


def _sch3(y, x):
    '''(INTERNAL) Compute sin, cos and hypotenuse.
    '''
    h = hypot(y, x)
    if h > 0:  # EPS_2
        s, c = y / h, x / h
    else:
        s, c = _0_0, _1_0
    return s, c, h


def _xyzn4(xyz, y, z, Types, Error=EcefError, name=NN):  # in .ltp, @see: .vector3d
    '''(INTERNAL) Get an C{(x, y, z, name)} 4-tuple.
    '''
    try:
        try:
            t = xyz.x, xyz.y, xyz.z, getattr(xyz, _name_, name)
            if not isinstance(xyz, Types):
                raise _TypesError(_xyz_, xyz, *Types)
        except AttributeError:
            t = map1(float, xyz, y, z) + (name,)

    except (TypeError, ValueError) as x:
        d = dict(zip((_xyz_, _y_, _z_), (xyz, y, z)))
        raise Error(txt=str(x), **d)

    return t


class _EcefBase(_NamedBase):
    '''(INTERNAL) Base class for L{EcefKarney}, L{EcefVeness} and L{EcefYou}.
    '''
    _datum = None
    _E     = None

    def __init__(self, a_ellipsoid, f=None, name=NN):
        '''New C{Ecef*} converter.

           @arg a_ellipsoid: A (non-prolate) ellipsoid (L{Ellipsoid}, L{Ellipsoid2},
                             L{Datum} or L{a_f2Tuple}) or C{scalar} ellipsoid's
                             equatorial radius (C{meter}).
           @kwarg f: C{None} or the ellipsoid flattening (C{scalar}), required
                     for C{scalar} B{C{a_ellipsoid}}, C{B{f}=0} represents a
                     sphere, negative B{C{f}} a prolate ellipsoid.
           @kwarg name: Optional name (C{str}).

           @raise EcefError: If B{C{a_ellipsoid}} not L{Ellipsoid}, L{Ellipsoid2},
                             L{Datum} or L{a_f2Tuple} or C{scalar} or B{C{f}} not
                             C{scalar} or if C{scalar} B{C{a_ellipsoid}} not positive
                             or B{C{f}} not less than 1.0.
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
            E =  d.ellipsoid
            if E.a < EPS or E.f > EPS1:
                raise ValueError

        except (TypeError, ValueError) as x:
            t = unstr(self.classname, a=a_ellipsoid, f=f)
            raise EcefError(_SPACE_(t, _ellipsoid_), txt=str(x))

        self._datum = d
        self._E = E

    def __eq__(self, other):
        '''Compare this and an other Ecef.

           @arg other: The other ecef (C{Ecef*}).

           @return: C{True} if equal, C{False} otherwise.
        '''
        return other is self or (isinstance(other, _EcefBase) and
                                 other.ellipsoid == self.ellipsoid)

    @Property_RO
    def equatoradius(self):
        '''Get the I{equatorial} radius, semi-axis (C{meter}).
        '''
        return self._E.a

    a = equatorialRadius = equatoradius  # Karney property

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

    def _forward(self, lat, lon, h, name, M=False, _philam=False):  # in .ltp.LocalCartesian.forward and -.reset
        '''(INTERNAL) Common for all C{Ecef*}.
        '''
        E = self.ellipsoid

        if _philam:
            sa, ca, sb, cb = sincos2_(lat, lon)
            lat = Lat(degrees90( lat))
            lon = Lon(degrees180(lon))
        else:
            sa, ca, sb, cb = sincos2d_(lat, lon)

        n = E.roc1_(sa, ca) if self._isYou else E.roc1_(sa)
        z = (h + n * E.e12) * sa
        x = (h + n) * ca

        m = self._Matrix(sa, ca, sb, cb) if M else None
        return Ecef9Tuple(x * cb, x * sb, z, lat, lon, h,
                                             0, m, self.datum,
                                             name=name or self.name)

    def forward(self, latlonh, lon=None, height=0, M=False, name=NN):
        '''Convert from geodetic C{(lat, lon, height)} to geocentric C{(x, y, z)}.

           @arg latlonh: Either a C{LatLon}, an L{Ecef9Tuple} or C{scalar}
                         latitude (C{degrees}).
           @kwarg lon: Optional C{scalar} longitude for C{scalar} B{C{latlonh}}
                       (C{degrees}).
           @kwarg height: Optional height (C{meter}), vertically above (or below)
                          the surface of the ellipsoid.
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).
           @kwarg name: Optional name (C{str}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)} with
                    geocentric C{(x, y, z)} coordinates for the given geodetic ones
                    C{(lat, lon, height)}, case C{C} 0, optional C{M} (L{EcefMatrix})
                    and C{datum} if available.

           @raise EcefError: If B{C{latlonh}} not C{LatLon}, L{Ecef9Tuple} or
                             C{scalar} or B{C{lon}} not C{scalar} for C{scalar}
                             B{C{latlonh}} or C{abs(lat)} exceeds 90°.

           @note: Use method C{.forward_} to specify C{lat} and C{lon} in C{radians}
                  and avoid double angle conversions.
        '''
        llhn = _llhn4(latlonh, lon, height, name=name)
        return _EcefBase._forward(self, *llhn, M=M)

    def forward_(self, phi, lam, height=0, M=False, name=NN):
        '''Like method C{.forward} except with geodetic lat- and longitude given
           in I{radians}.

           @arg phi: Latitude in I{radians} (C{scalar}).
           @arg lam: Longitude in I{radians} (C{scalar}).
           @kwarg height: Optional height (C{meter}), vertically above (or below)
                          the surface of the ellipsoid.
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).
           @kwarg name: Optional name (C{str}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with C{lat} set to C{degrees90(B{phi})} and C{lon} to
                    C{degrees180(B{lam})}.

           @raise EcefError: If B{C{phi}} or B{C{lam}} invalid or not C{scalar}.
        '''
        try:  # like function C{_llhn4} above
            plhn = Phi(phi), Lam(lam), Height(height), name
        except (TypeError, ValueError) as x:
            raise EcefError(phi=phi, lam=lam, height=height, txt=str(x))
        return self._forward(*plhn, M=M, _philam=True)

    @Property_RO
    def _Geocentrics(self):
        '''(INTERNAL) Valid geocentric classes.
        '''
        from pygeodesy.cartesianBase import CartesianBase
        return Ecef9Tuple, CartesianBase

    @Property_RO
    def _isYou(self):
        '''(INTERNAL) Is this an C{EcefYou}?.
        '''
        return isinstance(self, EcefYou)

    def _Matrix(self, sa, ca, sb, cb):
        '''Creation a rotation matrix.

           @arg sa: C{sin(phi)} (C{float}).
           @arg ca: C{cos(phi)} (C{float}).
           @arg sb: C{sin(lambda)} (C{float}).
           @arg cb: C{cos(lambda)} (C{float}).

           @return: An L{EcefMatrix}.
        '''
        return self._xnamed(EcefMatrix(sa, ca, sb, cb))

    def reverse(self, xyz, y=None, z=None, M=False, name=NN):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self, xyz, y=y, z=z, M=M, name=name)

    def toStr(self, prec=9, **unused):  # PYCHOK signature
        '''Return this C{Ecef*} as a string.

           @kwarg prec: Optional precision, number of decimal digits (0..9).

           @return: This C{Ecef*} representation (C{str}).
        '''
        return self.attrs(_a_, _f_, _datum_, _ellipsoid_, _name_, prec=prec)


class EcefKarney(_EcefBase):
    '''Conversion between geodetic and geocentric, aka I{Earth-Centered,
       Earth-Fixed} (ECEF) coordinates transcoded from I{Karney}'s C++ U{Geocentric
       <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Geocentric.html>}
       methods.

       @note: On methods C{.forward} and C{.forwar_}, let C{v} be a unit vector located
              at C{(lat, lon, h)}.  We can express C{v} as column vectors in one of two
              ways, C{v1} in east, north, up coordinates (where the components are
              relative to a local coordinate system at C{C(lat0, lon0, h0)}) or as C{v0}
              in geocentric C{x, y, z} coordinates.  Then, M{v0 = M ⋅ v1} where C{M} is
              the rotation matrix.
    '''

    @Property_RO
    def hmax(self):
        '''Get the distance or height limit (C{meter}, conventionally).
        '''
        return self.equatoradius / EPS_2  # self.equatoradius * _2_EPS, 12M lighyears

    def reverse(self, xyz, y=None, z=None, M=False, name=NN):
        '''Convert from geocentric C{(x, y, z)} to geodetic C{(lat, lon, height)}.

           @arg xyz: A geocentric (C{Cartesian}, L{Ecef9Tuple}) or C{scalar} ECEF C{x}
                     coordinate (C{meter}).
           @kwarg y: ECEF C{y} coordinate for C{scalar} B{C{xyz}} and B{C{z}} (C{meter}).
           @kwarg z: ECEF C{z} coordinate for C{scalar} B{C{xyz}} and B{C{y}} (C{meter}).
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).
           @kwarg name: Optional name (C{str}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)} with
                    geodetic coordinates C{(lat, lon, height)} for the given geocentric
                    ones C{(x, y, z)}, case C{C}, optional C{M} (L{EcefMatrix}) and
                    C{datum} if available.

           @raise EcefError: Invalid B{C{xyz}} or C{scalar} C{x} or B{C{y}} and/or B{C{z}}
                             not C{scalar} for C{scalar} B{C{xyz}}.

           @note: In general, there are multiple solutions and the result which minimizes
                  C{height} is returned, i.e., C{(lat, lon)} corresponds to the closest
                  point on the ellipsoid.  If there are still multiple solutions with
                  different latitudes (applies only if C{z} = 0), then the solution with
                  C{lat} > 0 is returned.  If there are still multiple solutions with
                  different longitudes (applies only if C{x} = C{y} = 0) then C{lon} = 0
                  is returned.  The returned C{height} value is not below M{−E.a * (1 −
                  E.e2) / sqrt(1 − E.e2 * sin(lat)**2)}.  The returned C{lon} is in the
                  range [−180°, 180°].  Like C{forward} above, M{v1 = Transpose(M) ⋅ v0}.
        '''
        x, y, z, name = _xyzn4(xyz, y, z, self._Geocentrics, name=name)

        E = self.ellipsoid

        sb, cb, R = _sch3(y, x)
        h = hypot(R, z)  # distance to earth center
        if h > self.hmax:  # PYCHOK no cover
            # We are really far away (> 12M light years).  Treat the earth
            # as a point and h, above as an acceptable approximation to the
            # height.  This avoids overflow, e.g., in the computation of disc
            # below.  It's possible that h has overflowed to INF, that's OK.
            # Treat finite x, y, but R overflows to +INF by scaling by 2.
            sb, cb, R = _sch3(y * _0_5, x * _0_5)
            sa, ca, _ = _sch3(z * _0_5, R)
            C = 1

        elif E.e4:  # E.isEllipsoidal
            # Treat prolate spheroids by swapping R and Z here and by
            # switching the arguments to phi = atan2(...) at the end.
            p = (R / E.a)**2
            q = E.e12 * (z / E.a)**2
            if E.isProlate:
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
                        t3 += copysign0(sqrt(disc), t3)  # t3 = (r * t)^3
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
                if E.isProlate:
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
                if E.isProlate:
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
            sa, ca, _ = _sch3((z if h else _1_0), R)
            h -= E.a
            C  = 4

        m = self._Matrix(sa, ca, sb, cb) if M else None
        return Ecef9Tuple(x, y, z, atan2d(sa, ca),
                                   atan2d(sb, cb), h,
                                   C, m, self.datum,
                                   name=name or self.name)


class EcefFarrell21(_EcefBase):
    '''Conversion between geodetic and geocentric, aka I{Earth-Centered,
       Earth-Fixed} (ECEF) coordinates based on I{Jay A. Farrell}'s U{Table
       2.1<https://Books.Google.com/books?id=fW4foWASY6wC>}, page 29.
    '''

    def reverse(self, xyz, y=None, z=None, M=None, name=NN):  # PYCHOK unused M
        '''Convert from geocentric C{(x, y, z)} to geodetic C{(lat, lon, height)} using
           I{Farrell}'s U{Table 2.1<https://Books.Google.com/books?id=fW4foWASY6wC>},
           page 29.

           @arg xyz: A geocentric (C{Cartesian}, L{Ecef9Tuple}) or C{scalar} ECEF C{x}
                     coordinate (C{meter}).
           @kwarg y: ECEF C{y} coordinate for C{scalar} B{C{xyz}} and B{C{z}} (C{meter}).
           @kwarg z: ECEF C{z} coordinate for C{scalar} B{C{xyz}} and B{C{y}} (C{meter}).
           @kwarg M: I{Ignored}, rotation matrix C{M} not available.
           @kwarg name: Optional name (C{str}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)} with
                    geodetic coordinates C{(lat, lon, height)} for the given geocentric
                    ones C{(x, y, z)}, case C{C=1}, C{M=None} always and C{datum}
                    if available.

           @raise EcefError: Invalid B{C{xyz}} or C{scalar} C{x} or B{C{y}} and/or B{C{z}}
                             not C{scalar} for C{scalar} B{C{xyz}} or C{sqrt} domain or
                             zero division error.

           @see: L{EcefFarrell22} and L{EcefVeness}.
        '''
        x, y, z, name = _xyzn4(xyz, y, z, self._Geocentrics, name=name)

        E  = self.ellipsoid
        a  = E.a
        a2 = E.a2
        b2 = E.b2
        e_ = E.a_b * E.e  # 0.0820944... WGS84
        e2 = E.e2
        e4 = E.e4

        try:
            z2 = z**2
            ez = (_1_0 - e2) * z2  # E.e2s2(z)

            p  = hypot(x, y)
            p2 = p**2
            F  = 54 * b2 * z2
            G  = p2 + ez - e2 * (a2 - b2)  # p2 + ez - e4 * a2
            c  = e4 * F * p2 / G**3
            s  = cbrt(_1_0 + c + sqrt(c**2 + c + c))
            P  = F / (_3_0 * fsum_(_1_0, s, _1_0 / s)**2 * G**2)
            Q  = sqrt(_1_0 + _2_0 * e4 * P)
            Q1 = Q + _1_0
            r0 = P * e2 * p / Q1 - sqrt(fsum_(a2 * (Q1 / Q) * _0_5,
                                              -P * ez / (Q * Q1),
                                              -P * p2 * _0_5))
            r = p + e2 * r0
            v = b2 / (a * sqrt(r**2 + ez))

            h = hypot(r, z) * (_1_0 - v)
            t = atan2(z + e_**2 * v * z, p)

        except (ValueError, ZeroDivisionError) as e:
            raise EcefError(x=x, y=y, z=z, txt=str(e))

        return Ecef9Tuple(x, y, z, degrees90(t), atan2d(y, x), h,
                                   1, None, self.datum,
                                   name=name or self.name)


class EcefFarrell22(_EcefBase):
    '''Conversion between geodetic and geocentric, aka I{Earth-Centered,
       Earth-Fixed} (ECEF) coordinates based on I{Jay A. Farrell}'s U{Table
       2.2<https://Books.Google.com/books?id=fW4foWASY6wC>}, page 30.
    '''

    def reverse(self, xyz, y=None, z=None, M=None, name=NN):  # PYCHOK unused M
        '''Convert from geocentric C{(x, y, z)} to geodetic C{(lat, lon, height)} using
           I{Farrell}'s U{Table 2.2<https://Books.Google.com/books?id=fW4foWASY6wC>},
           page 30.

           @arg xyz: A geocentric (C{Cartesian}, L{Ecef9Tuple}) or C{scalar} ECEF C{x}
                     coordinate (C{meter}).
           @kwarg y: ECEF C{y} coordinate for C{scalar} B{C{xyz}} and B{C{z}} (C{meter}).
           @kwarg z: ECEF C{z} coordinate for C{scalar} B{C{xyz}} and B{C{y}} (C{meter}).
           @kwarg M: I{Ignored}, rotation matrix C{M} not available.
           @kwarg name: Optional name (C{str}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)} with
                    geodetic coordinates C{(lat, lon, height)} for the given geocentric
                    ones C{(x, y, z)}, case C{C=1}, C{M=None} always and C{datum}
                    if available.

           @raise EcefError: Invalid B{C{xyz}} or C{scalar} C{x} or B{C{y}} and/or B{C{z}}
                             not C{scalar} for C{scalar} B{C{xyz}} or C{sqrt} domain or
                             zero division error.

           @see: L{EcefFarrell21} and L{EcefVeness}.
        '''
        x, y, z, name = _xyzn4(xyz, y, z, self._Geocentrics, name=name)

        E = self.ellipsoid
        a = E.a
        b = E.b

        try:  # see EcefVeness.reverse
            p = hypot(x, y)
            s, c = sincos2(atan2(z * a, p * b))

            t = atan2(z + E.e22 * b * s**3,
                      p - E.e2  * a * c**3)

            s, c = sincos2(t)
            h = p / c - E.roc1_(s)  # E.a / sqrt(1 - e2 * s**2)

        except (ValueError, ZeroDivisionError) as e:
            raise EcefError(x=x, y=y, z=z, txt=str(e))

        return Ecef9Tuple(x, y, z, degrees90(t), atan2d(y, x), h,
                                   1, None, self.datum,
                                   name=name or self.name)


class EcefSudano(_EcefBase):
    '''Conversion between geodetic and geocentric, aka I{Earth-Centered,
       Earth-Fixed} (ECEF) coordinates based on I{John J. Sudano}'s U{paper
       <https://www.ResearchGate.net/publication/
       3709199_An_exact_conversion_from_an_Earth-centered_coordinate_system_to_latitude_longitude_and_altitude>}.
    '''

    def reverse(self, xyz, y=None, z=None, M=None, name=NN):  # PYCHOK unused M
        '''Convert from geocentric C{(x, y, z)} to geodetic C{(lat, lon, height)} using
           I{Sudano}'s U{iterative method<https://www.ResearchGate.net/publication/
           3709199_An_exact_conversion_from_an_Earth-centered_coordinate_system_to_latitude_longitude_and_altitude>}.

           @arg xyz: A geocentric (C{Cartesian}, L{Ecef9Tuple}) or C{scalar} ECEF C{x}
                     coordinate (C{meter}).
           @kwarg y: ECEF C{y} coordinate for C{scalar} B{C{xyz}} and B{C{z}} (C{meter}).
           @kwarg z: ECEF C{z} coordinate for C{scalar} B{C{xyz}} and B{C{y}} (C{meter}).
           @kwarg M: I{Ignored}, rotation matrix C{M} not available.
           @kwarg name: Optional name (C{str}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)} with geodetic
                    coordinates C{(lat, lon, height)} for the given geocentric ones C{(x, y, z)},
                    iteration C{C}, C{M=None} always and C{datum} if available.

           @raise EcefError: Invalid B{C{xyz}} or C{scalar} C{x} or B{C{y}} and/or B{C{z}}
                             not C{scalar} for C{scalar} B{C{xyz}} or no convergence.
        '''
        x, y, z, name = _xyzn4(xyz, y, z, self._Geocentrics, name=name)

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
            t = unstr(self.reverse.__name__, x=x, y=y, z=z)
            raise EcefError(t, txt=_no_(_convergence_))

        if a is None:
            a = copysign0(asin(sa), z)
        h = fsum_(h * ca, abs(z * sa), -E.a * E.e2s(sa))  # use Veness',
        # since Sudano's Eq (7) doesn't provide the correct height
        # h = (abs(z) + h - E.a * cos(a + E.e12) * sa / ca) / (ca + sa)

        r = Ecef9Tuple(x, y, z, degrees90(a), atan2d(y, x), h,
                                C, None, self.datum,
                                name=name or self.name)
        r._iteration = C
        return r


class EcefVeness(_EcefBase):
    '''Conversion between geodetic and geocentric, aka I{Earth-Centered,
       Earth-Fixed} (ECEF) coordinates transcoded from I{Chris Veness}'
       JavaScript classes U{LatLonEllipsoidal, Cartesian
       <https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}.

       @see: U{A Guide to Coordinate Systems in Great Britain
             <https://www.OrdnanceSurvey.co.UK/documents/resources/guide-coordinate-systems-great-britain.pdf>},
             section I{B) Converting between 3D Cartesian and ellipsoidal
             latitude, longitude and height coordinates}.
    '''

    def reverse(self, xyz, y=None, z=None, M=None, name=NN):  # PYCHOK unused M
        '''Convert from geocentric C{(x, y, z)} to geodetic C{(lat, lon, height)}
           transcoded from I{Chris Veness}' U{JavaScript
           <https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}.

           Uses B. R. Bowring’s formulation for μm precision in concise
           form: U{'The accuracy of geodetic latitude and height equations'
           <https://www.ResearchGate.net/publication/
           233668213_The_Accuracy_of_Geodetic_Latitude_and_Height_Equations>},
           Survey Review, Vol 28, 218, Oct 1985.

           @arg xyz: A geocentric (C{Cartesian}, L{Ecef9Tuple}) or C{scalar} ECEF C{x}
                     coordinate (C{meter}).
           @kwarg y: ECEF C{y} coordinate for C{scalar} B{C{xyz}} and B{C{z}} (C{meter}).
           @kwarg z: ECEF C{z} coordinate for C{scalar} B{C{xyz}} and B{C{y}} (C{meter}).
           @kwarg M: I{Ignored}, rotation matrix C{M} not available.
           @kwarg name: Optional name (C{str}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)} with
                    geodetic coordinates C{(lat, lon, height)} for the given geocentric
                    ones C{(x, y, z)}, case C{C}, C{M=None} always and C{datum} if available.

           @raise EcefError: Invalid B{C{xyz}} or C{scalar} C{x} or B{C{y}} and/or B{C{z}}
                             not C{scalar} for C{scalar} B{C{xyz}}.

           @see: Ralph M. Toms U{'An Efficient Algorithm for Geocentric to Geodetic
                 Coordinate Conversion'<https://www.OSTI.gov/scitech/biblio/110235>},
                 Sept 1995 and U{'An Improved Algorithm for Geocentric to Geodetic
                 Coordinate Conversion'<https://www.OSTI.gov/scitech/servlets/purl/231228>},
                 Apr 1996, both from Lawrence Livermore National Laboratory (LLNL) and
                 John J. Sudano U{An exact conversion from an Earth-centered coordinate
                 system to latitude longitude and altitude<https://www.ResearchGate.net/
                 publication/3709199_An_exact_conversion_from_an_Earth-centered_coordinate_system_to_latitude_longitude_and_altitude>}.
        '''
        x, y, z, name = _xyzn4(xyz, y, z, self._Geocentrics, name=name)

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
            C, lat, lon, h = 3, copysign0(_90_0, z), _0_0, abs(z) - E.b

        return Ecef9Tuple(x, y, z, lat, lon, h,
                                   C, None,  # M=None
                                   self.datum, name=name or self.name)


class EcefYou(_EcefBase):
    '''Conversion between geodetic and geocentric, aka I{Earth-Centered,
       Earth-Fixed} (ECEF) coordinates using I{Rey-Jer You}'s U{transformation
       <https://www.ResearchGate.net/publication/240359424>}.

       @see: W.E. Featherstone, S.J. (Sten) Claessens U{Closed-form transformation
             between geodetic and ellipsoidal coordinates
             <https://Espace.Curtin.edu.AU/bitstream/handle/20.500.11937/11589/115114_9021_geod2ellip_final.pdf>}
             Studia Geophysica et Geodaetica, 2008, 52, 1-18 and U{PyMap3D
             <https://PyPI.org/project/pymap3d>}.
    '''

    def __init__(self, a_ellipsoid, f=None, name=NN):
        _EcefBase.__init__(self, a_ellipsoid, f=f, name=name)  # inherited documentation
        E = self.ellipsoid
        if E.isProlate or (E.a2 - E.b2) < 0:
            raise EcefError(ellipsoid=E, txt=_prolate_)

    def reverse(self, xyz, y=None, z=None, M=None, name=NN):  # PYCHOK unused M
        '''Convert geocentric C{(x, y, z)} to geodetic C{(lat, lon, height)}
           using I{Rey-Jer You}'s transformation.

           @arg xyz: A geocentric (C{Cartesian}, L{Ecef9Tuple}) or C{scalar} ECEF C{x}
                     coordinate (C{meter}).
           @kwarg y: ECEF C{y} coordinate for C{scalar} B{C{xyz}} and B{C{z}} (C{meter}).
           @kwarg z: ECEF C{z} coordinate for C{scalar} B{C{xyz}} and B{C{y}} (C{meter}).
           @kwarg M: I{Ignored}, rotation matrix C{M} not available.
           @kwarg name: Optional name (C{str}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)} with
                    geodetic coordinates C{(lat, lon, height)} for the given geocentric
                    ones C{(x, y, z)}, case C{C=1}, C{M=None} always and C{datum} if
                    available.

           @raise EcefError: Invalid B{C{xyz}} or C{scalar} C{x} or B{C{y}} and/or B{C{z}}
                             not C{scalar} for C{scalar} B{C{xyz}}.
        '''
        x, y, z, name = _xyzn4(xyz, y, z, self._Geocentrics, name=name)

        r2 = hypot2_(x, y, z)

        E  = self.ellipsoid
        e2 = E.a2 - E.b2  # == E.e2 * E.a2
        if e2 < 0:
            raise EcefError(ellipsoid=E, txt=_prolate_)
        e = sqrt(e2)  # XXX sqrt0(e2)?

        q = hypot(x, y)
        u = fsum_(r2, -e2, hypot(r2 - e2, 2 * e * z)) * _0_5
        if u > EPS02:
            u = sqrt(u)
            p = hypot(u, e)
            B = atan2(p * z, u * q)  # beta0 = atan(p / u * z / q)
            sB, cB = sincos2(B)
            if cB and sB:
                p *= E.a
                d  = (p / cB - e2 * cB) / sB
                if isnon0(d):
                    B += fsum_(u * E.b, -p, e2) / d
                    sB, cB = sincos2(B)
        elif u < 0:
            raise EcefError(x=x, y=y, z=z, txt=_singular_)
        else:
            sB, cB = copysign0(_1_0, z), _0_0

        h = hypot(z - E.b * sB, q - E.a * cB)
        if hypot2_(x, y, z * E.a_b) < E.a2:
            h = neg(h)  # inside ellipsoid

        return Ecef9Tuple(x, y, z, atan2d(E.a * sB, E.b * cB),  # atan(E.a_b * tan(B))
                                   atan2d(y, x), h,
                                   1, None,  # C=1, M=None
                                   self.datum, name=name or self.name)


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

    def __new__(cls, sa, ca, sb, cb, *_more):
        '''New L{EcefMatrix} matrix.

           @arg sa: C{sin(phi)} (C{float}).
           @arg ca: C{cos(phi)} (C{float}).
           @arg sb: C{sin(lambda)} (C{float}).
           @arg cb: C{cos(lambda)} (C{float}).
           @arg _more: (INTERNAL) from C{.multiply}.

           @raise EcefError: If B{C{sa}}, B{C{ca}}, B{C{sb}} or
                             B{C{cb}} outside M{[-1.0, +1.0]}.
        '''
        t = sa, ca, sb, cb
        if _more:  # all 9 matrix elements ...
            t += _more  # ... from .multiply

        elif max(map(abs, t)) > _1_0:
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

    def column(self, column):
        '''Get matrix B{C{column}} as 3-tuple.
        '''
        if 0 <= column < 3:
            return self[column::3]
        raise _IndexError(column=column)

    @Property_RO
    def _column_0(self):
        return self.column(0)

    @Property_RO
    def _column_1(self):
        return self.column(1)

    @Property_RO
    def _column_2(self):
        return self.column(2)

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

        # like LocalCartesian.MatrixMultiply, transposed(self) X other
        # <https://GeographicLib.SourceForge.io/html/LocalCartesian_8cpp_source.html>
        M = (fdot(self[r::3], *other[c::3]) for r in range(3) for c in range(3))
        return _xnamed(EcefMatrix(*M), EcefMatrix.multiply.__name__)

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
        return (fdot(xyz, *self._column_0),
                fdot(xyz, *self._column_1),
                fdot(xyz, *self._column_2))

    def row(self, row):
        '''Get matrix B{C{row}} as 3-tuple.
        '''
        if 0 <= row < 3:
            r = row * 3
            return self[r:r+3]
        raise _IndexError(row=row)

    @Property_RO
    def _row_0(self):
        return self.row(0)

    @Property_RO
    def _row_1(self):
        return self.row(1)

    @Property_RO
    def _row_2(self):
        return self.row(2)

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

            _xyz = _1_0_T + xyz
            # x' = x0 + M[0] * x + M[1] * y + M[2] * z
            # y' = y0 + M[3] * x + M[4] * y + M[5] * z
            # z' = z0 + M[6] * x + M[7] * y + M[8] * z
            xyz_ = (fdot(_xyz, xyz0[0], *self._row_0),
                    fdot(_xyz, xyz0[1], *self._row_1),
                    fdot(_xyz, xyz0[2], *self._row_2))
        else:
            # x' = M[0] * x + M[1] * y + M[2] * z
            # y' = M[3] * x + M[4] * y + M[5] * z
            # z' = M[6] * x + M[7] * y + M[8] * z
            xyz_ = (fdot(xyz, *self._row_0),
                    fdot(xyz, *self._row_1),
                    fdot(xyz, *self._row_2))
        return xyz_


class Ecef9Tuple(_NamedTuple):
    '''9-Tuple C{(x, y, z, lat, lon, height, C, M, datum)} with I{geocentric}
       C{x}, C{y} and C{z} plus I{geodetic} C{lat}, C{lon} and C{height}, case
       C{C} (see the C{Ecef*.reverse} methods) and optionally, the rotation
       matrix C{M} (L{EcefMatrix}) and C{datum}, with C{lat} and C{lon} in
       C{degrees} and C{x}, C{y}, C{z} and C{height} in C{meter}, conventionally.
    '''
    _Names_ = (_x_,   _y_,   _z_,   _lat_, _lon_, _height_, _C_,  _M_,   _datum_)
    _Units_ = ( Meter, Meter, Meter, Lat,   Lon,   Height,   Int, _Pass, _Pass)

    _IndexM = _Names_.index(_M_)  # for ._M_x_M

    @deprecated_method
    def convertDatum(self, datum2):  # for backward compatibility
        '''DEPRECATED, use method L{toDatum}.'''
        return self.toDatum(datum2)

    @Property_RO
    def lam(self):
        '''Get the longitude in C{radians} (C{float}).
        '''
        return self.philam.lam

    @Property_RO
    def lamVermeille(self):
        '''Get the longitude in C{radians [-PI*3/2..+PI*3/2]} after U{Vermeille
           <https://Search.ProQuest.com/docview/639493848>} (2004), p 95.

           @see: U{Karney<https://GeographicLib.SourceForge.io/html/geocentric.html>},
                 U{Vermeille<https://Search.ProQuest.com/docview/847292978>} 2011, pp 112-113, 116
                 and U{Featherstone, et.al.<https://Search.ProQuest.com/docview/872827242>}, p 7.
        '''
        x, y = self.x, self.y
        if y > EPS0:
            r = _N_2_0 * atan2(x, hypot(y, x) + y) + PI_2
        elif y < -EPS0:
            r = _2_0   * atan2(x, hypot(y, x) - y) - PI_2
        else:  # y == 0
            r = PI if x < 0 else _0_0
        return Lam(Vermeille=r)

    @Property_RO
    def latlon(self):
        '''Get the lat-, longitude in C{degrees} (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return LatLon2Tuple(self.lat, self.lon, name=self.name)

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
        '''Get the latitude and I{Vermeille} longitude in C{degrees [-225..+225]} (L{LatLon2Tuple}C{(lat, lon)}).

           @see: Property C{lonVermeille}.
        '''
        return LatLon2Tuple(self.lat, self.lonVermeille, name=self.name)

    @Property_RO
    def lonVermeille(self):
        '''Get the longitude in C{degrees [-225..+225]} after U{Vermeille
           <https://Search.ProQuest.com/docview/639493848>} (2004), p 95.

           @see: Property C{lamVermeille}.
        '''
        return Lon(Vermeille=degrees(self.lamVermeille))

    def _T_x_M(self, T):
        '''(INTERNAL) Update M{self.M = T.multiply(self.M)}.
        '''
        t = list(self)
        M = self._IndexM
        t[M] = T.multiply(t[M])
        return self.classof(*t)

    @Property_RO
    def phi(self):
        '''Get the latitude in C{radians} (C{float}).
        '''
        return self.philam.phi

    @Property_RO
    def philam(self):
        '''Get the lat-, longitude in C{radians} (L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        return PhiLam2Tuple(radians(self.lat), radians(self.lon), name=self.name)

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
        '''Get the latitude and I{Vermeille} longitude in C{radians [-PI*3/2..+PI*3/2]} (L{PhiLam2Tuple}C{(phi, lam)}).

           @see: Property C{lamVermeille}.
        '''
        return PhiLam2Tuple(radians(self.lat), self.lamVermeille, name=self.name)

    def toCartesian(self, Cartesian=None, **Cartesian_kwds):
        '''Return the geocentric C{(x, y, z)} coordinates as an ellipsoidal or spherical
           C{Cartesian}.

           @kwarg Cartesian: Optional class to return C{(x, y, z)} (L{ellipsoidalKarney.Cartesian},
                             L{ellipsoidalNvector.Cartesian}, L{ellipsoidalVincenty.Cartesian},
                             L{sphericalNvector.Cartesian} or L{sphericalTrigonometry.Cartesian})
                             or C{None}.
           @kwarg Cartesian_kwds: Optional, additional B{C{Cartesian}} keyword arguments, ignored
                                  if C{B{Cartesian} is None}.

           @return: A C{B{Cartesian}(x, y, z, **B{Cartesian_kwds})} instance or
                    a L{Vector4Tuple}C{(x, y, z, h)} if C{B{Cartesian} is None}.

           @raise TypeError: Invalid B{C{Cartesian}} or B{C{Cartesian_kwds}}.
        '''
        if Cartesian in (None, Vector4Tuple):
            r = self.xyzh
        elif Cartesian is Vector3Tuple:
            r = self.xyz
        else:
            from pygeodesy.cartesianBase import CartesianBase
            _xsubclassof(CartesianBase, Cartesian=Cartesian)
            r = Cartesian(self, **_xkwds(Cartesian_kwds, name=self.name))
        return r

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
            c = CartesianBase(self, datum=self.datum, name=self.name)  # PYCHOK _Names_
            # c.toLatLon converts datum, x, y, z, lat, lon, etc.
            # and returns another Ecef9Tuple iff LatLon is None
            r = c.toLatLon(datum=datum2, LatLon=None)
        return r

    def toLatLon(self, LatLon=None, **LatLon_kwds):
        '''Return the geodetic C{(lat, lon, height[, datum])} coordinates.

           @kwarg LatLon: Optional class to return C{(lat, lon, height[, datum])}
                          or C{None}.
           @kwarg LatLon_kwds: Optional B{C{height}}, B{C{datum}} and other
                               B{C{LatLon}} keyword arguments.

           @return: An instance of C{B{LatLon}(lat, lon, **B{LatLon_kwds})}
                    or if B{C{LatLon}} is C{None}, a L{LatLon3Tuple}C{(lat, lon,
                    height)} respectively L{LatLon4Tuple}C{(lat, lon, height,
                    datum)} depending on whether C{datum} is un-/specified.

           @raise TypeError: Invalid B{C{LatLon}} or B{C{LatLon_kwds}}.
        '''
        kwds = _xkwds(LatLon_kwds, height=self.height, datum=self.datum, name=self.name)  # PYCHOK Ecef9Tuple
        d = kwds[_datum_]
        if LatLon is None:
            r = LatLon3Tuple(self.lat, self.lon, kwds[_height_], name=kwds[_name_])  # PYCHOK Ecef9Tuple
            if d:
                r = r.to4Tuple(d)  # checks type(d)
        else:
            if d is None:  # remove None datum
                _ = kwds.pop[_datum_]
            r = LatLon(self.lat, self.lon, **kwds)  # PYCHOK Ecef9Tuple
        _datum_datum(getattr(r, _datum_, self.datum), self.datum)  # PYCHOK Ecef9Tuple
        return r

    def toLocal(self, ltp, Xyz=None, **Xyz_kwds):
        '''Convert this geocentric to I{local} C{x}, C{y} and C{z}.

           @kwarg ltp: The I{local tangent plane} (LTP) to use (L{Ltp}).
           @kwarg Xyz: Optional class to return C{x}, C{y} and C{z}
                       (L{XyzLocal}, L{Enu}, L{Ned}) or C{None}.
           @kwarg Xyz_kwds: Optional, additional B{C{Xyz}} keyword
                            arguments, ignored if C{B{Xyz} is None}.

           @return: An B{C{Xyz}} instance or if C{B{Xyz} is None},
                    a L{Local9Tuple}C{(x, y, z, lat, lon, height,
                    ltp, ecef, M)} with C{M=None}, always.

           @raise TypeError: Invalid B{C{ltp}}.
        '''
        from pygeodesy.ltp import _xLtp
        return _xLtp(ltp)._ecef2local(self, Xyz, Xyz_kwds)

    def toVector(self, Vector=None, **Vector_kwds):
        '''Return the geocentric C{(x, y, z)} coordinates as vector.

           @kwarg Vector: Optional vector class to return C{(x, y, z)} or
                          C{None}.
           @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword
                               arguments, ignored if C{B{Vector} is None}.

           @return: A C{Vector}C{(x, y, z, **Vector_kwds)} instance or a
                    L{Vector3Tuple}C{(x, y, z)} if B{C{Vector}} is C{None}.

           @see: Propertes C{xyz} and C{xyzh}
        '''
        return self.xyz if Vector is None else self._xnamed(
               Vector(self.x, self.y, self.z, **Vector_kwds))  # PYCHOK Ecef9Tuple

    @Property_RO
    def xyz(self):
        '''Get the geocentric C{(x, y, z)} coordinates (L{Vector3Tuple}C{(x, y, z)}).
        '''
        return Vector3Tuple(self.x, self.y, self.z, name=self.name)

    @Property_RO
    def xyzh(self):
        '''Get the geocentric C{(x, y, z)} coordinates and C{height} (L{Vector4Tuple}C{(x, y, z, h)})
        '''
        return self.xyz.to4Tuple(self.height)


def _xEcef(Ecef):  # PYCHOK .latlonBase.py
    '''(INTERNAL) Validate B{C{Ecef}} I{class}.
    '''
    if issubclassof(Ecef, _EcefBase):
        return Ecef
    raise _TypesError(_Ecef_, Ecef, EcefFarrell21, EcefFarrell22, EcefKarney, EcefSudano, EcefVeness, EcefYou)


__all__ += _ALL_DOCS(_EcefBase)

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
