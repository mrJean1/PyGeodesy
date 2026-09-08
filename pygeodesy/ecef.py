
# -*- coding: utf-8 -*-

u'''I{Geocentric} Earth-Centered, Earth-Fixed (ECEF) coordinates.

Geocentric conversions transcoded from I{Charles Karney}'s C++ class U{Geocentric
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Geocentric.html>}
into pure Python class L{EcefKarney}, class L{EcefFukushima} from I{Toshio Fukushima}'s U{Fortran
<https://www.ResearchGate.net/publication/277721539>} version, class L{EcefSudano} based on
I{John Sudano}'s U{paper<https://www.ResearchGate.net/publication/3709199>}, class L{EcefUPC}
using the I{Universitat Politècnica de Catalunya}'s U{method, page 186
<https://GSSC.ESA.int/navipedia/GNSS_Book/ESA_GNSS-Book_TM-23_Vol_I.pdf>}, class L{EcefVeness}
transcoded from I{Chris Veness}' JavaScript classes U{LatLonEllipsoidal, Cartesian
<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}, class L{EcefYou}
implementing I{Rey-Jer You}'s U{transformations<https://www.ResearchGate.net/publication/240359424>}
and classes L{EcefFarrell22} and L{EcefFarrell22} from I{Jay A. Farrell}'s U{Table 2.1 and 2.2
<https://Books.Google.com/books?id=fW4foWASY6wC>}, page 29-30.

Following is a copy of I{Karney}'s U{Detailed Description
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Geocentric.html>}.

Convert between geodetic coordinates C{lat}-, C{lon}gitude and height C{h} (measured vertically
from the surface of the ellipsoid) to geocentric C{x}, C{y} and C{z} coordinates, also known as
I{Earth-Centered, Earth-Fixed} (U{ECEF<https://WikiPedia.org/wiki/ECEF>}).

The origin of geocentric coordinates is at the center of the earth.  The C{z} axis goes thru
the North pole, C{lat} = 90°.  The C{x} axis goes thru C{lat} = 0°, C{lon} = 0°.

The I{local (cartesian) origin} is at (C{lat0}, C{lon0}, C{height0}).  The I{local} C{x} axis points
East, the I{local} C{y} axis points North and the I{local} C{z} axis is normal to the ellipsoid.  The
plane C{z = -height0} is tangent to the ellipsoid, hence the alternate name I{local tangent plane}.

Forward conversion from geodetic to geocentric (ECEF) coordinates is straightforward.

For the reverse transformation we use Hugues Vermeille's U{Direct transformation from geocentric
coordinates to geodetic coordinates<https://DOI.org/10.1007/s00190-002-0273-6>}, J. Geodesy
(2002) 76, page 451-454.

Several changes have been made to ensure that the method returns accurate results for all finite
inputs (even if h is infinite).  The changes are described in Appendix B of C. F. F. Karney
U{Geodesics on an ellipsoid of revolution<https://ArXiv.org/abs/1102.1215v1>}, Feb. 2011, 85,
105-117 (U{preprint<https://ArXiv.org/abs/1102.1215v1>}).  Vermeille similarly updated his method
in U{An analytical method to transform geocentric into geodetic coordinates
<https://DOI.org/10.1007/s00190-010-0419-x>}, J. Geodesy (2011) 85, page 105-117.  See U{Geocentric
coordinates<https://GeographicLib.SourceForge.io/C++/doc/geocentric.html>} for more information.

The errors in these routines are close to round-off.  Specifically, for points within 5,000 Km of
the surface of the ellipsoid (either inside or outside the ellipsoid), the error is bounded by 7
nm (7 nanometers) for the WGS84 ellipsoid.  See U{Geocentric coordinates
<https://GeographicLib.SourceForge.io/C++/doc/geocentric.html>} for further information on the errors.

@note: The C{reverse} methods of all C{Ecef...} classes return by default C{INT0} as the (geodetic)
longitude for I{polar} ECEF location C{x == y == 0}.  Use keyword argument C{lon00} or property
C{lon00} to configure that value.

@see: Module L{ltp} and class L{LocalCartesian}, a transcription of I{Charles Karney}'s C++ class
U{LocalCartesian<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1LocalCartesian.html>},
for conversion between geodetic and I{local cartesian} coordinates in a I{local tangent
plane} as opposed to I{geocentric} (ECEF) ones.
'''

from pygeodesy.basics import copysign0, _isin, isscalar, issubclassof, neg, map1, \
                            _xinstanceof, _xsubclassof,  typename  # _args_kwds_names
from pygeodesy.constants import EPS, EPS0, EPS02, EPS1, INT0, PI, PI_2, _0_0, _0_5, \
                               _1_0, _1_0_1T, _1_5, _2_0, _3_0, _4_0, _6_0, _90_0, \
                               _copysign_1_0, _isNAN, _isNAN0, _over,  isnon0  # PYCHOK used!
from pygeodesy.datums import _ellipsoidal_datum, _WGS84,  a_f2Tuple, _EWGS84
from pygeodesy.ecefLocals import _EcefLocal
# from pygeodesy.ellipsoids import a_f2Tuple, _EWGS84  # from .datums
from pygeodesy.errors import _IndexError, LenError, _ValueError, _TypesError, \
                             _xattr, _xdatum, _xkwds, _xkwds_get
from pygeodesy.fmath import cbrt, _fdotf, hypot, hypot1, hypot2_
from pygeodesy.fsums import Fsum, fsumf_,  Fmt, unstr
# from pygeodesy.internals import typename  # from .basics
from pygeodesy.interns import NN, _a_, _C_, _datum_, _ellipsoid_, _f_, _height_, \
                             _lat_, _lon_, _M_, _name_, _singular_, _SPACE_, \
                             _x_, _xyz_, _y_, _z_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _name__, _name1__, _NamedBase, _NamedTuple, _Pass, _xnamed
from pygeodesy.namedTuples import LatLon2Tuple, LatLon3Tuple, \
                                  PhiLam2Tuple, Vector3Tuple, Vector4Tuple
from pygeodesy.props import deprecated_method, deprecated_property, Property_RO, \
                            property_RO, property_ROver
# from pygeodesy.streprs import Fmt, unstr  # from .fsums
from pygeodesy.units import _isRadius, Degrees, Degrees_, Height, Int, Lam, Lat, \
                             Lon, Meter, Phi, Scalar, Scalar_
from pygeodesy.utily import atan1, atan1d, atan2, atan2d, degrees90, degrees180, \
                            sincos2, sincos2_, sincos2d_
# from pygeodesy.vector3d import Vector3d  # _MODS

from math import cos, degrees, fabs, radians, sqrt

__all__ = _ALL_LAZY.ecef
__version__ = '26.09.06'

_Ecef_    = 'Ecef'
_prolate_ = 'prolate'
_TOL      =  1.e-12  # degrees > 1.e-14
_TRIPS    =  33  # 8..9 sufficient
_xyz_y_z  = _xyz_, _y_, _z_  # _args_kwds_names(_xyzn4)[:3]


def _Degrees2Radians(tol):  # for EcefUPC
    return Degrees_(tol=tol, low=EPS, Error=EcefError).toRadians()


class _EcefBase(_NamedBase):
    '''(INTERNAL) Base class for C{Ecef*} convertor classes.
    '''
    _datum = _WGS84
    _e_e2  =  None
    _E     = _EWGS84
    _isYou =  False
    _lon00 =  INT0  # arbitrary, "polar" lon for LocalCartesian, Ltp

    def __init__(self, a_ellipsoid=_EWGS84, f=None, lon00=INT0, **name):
        '''New C{Ecef*} converter.

           @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple})
                             or a datum (L{Datum}) or the ellipsoid's equatorial
                             radius (C{meter}).
           @kwarg f: C{None} or the ellipsoid flattening (C{scalar}), required if
                     C{B{a_ellipsoid} is scalar}.
           @kwarg lon00: An arbitrary, I{"polar"} longitude (C{degrees}), see the
                         C{reverse} method.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise EcefError: If B{C{a_ellipsoid}} is not an L{Ellipsoid}, L{Ellipsoid2},
                             L{a_f2Tuple} or L{Datum} instance or a positive C{scalar}
                             or if B{C{f}} is not C{scalar} and less than C{1.0}.
        '''
        try:
            E = a_ellipsoid
            if f is None:
                pass
            elif _isRadius(E) and isscalar(f):
                E = a_f2Tuple(E, f)
            else:
                raise ValueError()  # _invalid_

            if not _isin(E, _EWGS84, _WGS84):
                d = _ellipsoidal_datum(E, **name)
                E =  d.ellipsoid
                if E.a < EPS or E.f > EPS1:
                    raise ValueError()  # _invalid_
                self._datum = d
                self._E     = E

            if self._isYou:
                E  = self.ellipsoid
                e2 = E.a2 - E.b2
                if e2 < 0 or E.f < 0:
                    raise EcefError(ellipsoid=E, txt=_prolate_)
                self._e_e2 = sqrt(e2), e2

        except (TypeError, ValueError) as x:
            t = unstr(self.classname, a=a_ellipsoid, f=f)
            raise EcefError(_SPACE_(t, _ellipsoid_), cause=x)

        if name:
            self.name = name
        if lon00 is not INT0:
            self.lon00 = lon00

    def __eq__(self, other):
        '''Compare this and an other Ecef.

           @arg other: The other ecef (C{Ecef*}).

           @return: C{True} if equal, C{False} otherwise.
        '''
        return other is self or (isinstance(other, type(self)) and
                                 other.ellipsoid == self.ellipsoid)

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
    def equatoradius(self):
        '''Get the C{ellipsoid}'s equatorial radius, semi-axis (C{meter}).
        '''
        return self.ellipsoid.a

    a = equatorialRadius = equatoradius  # Karney property

    @Property_RO
    def flattening(self):  # Karney property
        '''Get the C{ellipsoid}'s flattening (C{scalar}), positive for
           I{oblate}, negative for I{prolate} or C{0} for I{near-spherical}.
        '''
        return self.ellipsoid.f

    f = flattening

    def _forward(self, lat, lon, h, name, M=False, _philam=False):  # in .ltp.LocalCartesian.forward and -.reset
        '''(INTERNAL) Common for all C{Ecef*}.

           @note: From C{Karney}'s, let C{v} be a unit vector located at C{(lat,
                  lon, h)}.  We can express C{v} as column vectors in one of two
                  ways, C{v1} in East, North, Up (ENU) coordinates (where the
                  components are relative to a local coordinate system at C{C(lat0,
                  lon0, h0)}) or as C{v0} in geocentric C{x, y, z} coordinates.
                  Then, M{v0 = M ⋅ v1} where C{M} is the rotation matrix.
        '''
        if _philam:  # lat, lon in radians
            sa, ca, sb, cb = sincos2_(lat, lon)
            lat = Lat(degrees90( lat), Error=EcefError)
            lon = Lon(degrees180(lon), Error=EcefError)
        else:
            sa, ca, sb, cb = sincos2d_(lat, lon)

        E =  self.ellipsoid
        n =  E.roc1_(sa, ca) if self._isYou else E.roc1_(sa)
        H = _isNAN0(h)
        c = (H + n) * ca
        x = cb * c
        y = sb * c
        z = (H + n * E.e21) * sa

        m = self._Matrix(sa, ca, sb, cb) if M else None
        n = self._name__(name)
        return Ecef9Tuple(x, y, z, lat, lon, h, 0,  # C=0, forward
                                   m, self.datum, name=n)

    def forward(self, latlonh, lon=None, height=0, M=False, **name):
        '''Convert from geodetic C{(lat, lon, height)} to geocentric C{(x, y, z)}.

           @arg latlonh: Either a C{LatLon}, an L{Ecef9Tuple} or C{scalar}
                         latitude (C{degrees}).
           @kwarg lon: Optional C{scalar} longitude for C{scalar} B{C{latlonh}}
                       (C{degrees}).
           @kwarg height: Optional height (C{meter}), vertically above (or below)
                          the surface of the ellipsoid.
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)} with
                    geocentric C{(x, y, z)} coordinates for the given geodetic ones
                    C{(lat, lon, height)}, case C{C} (0, forward), rotation matrix
                    C{M} (L{EcefMatrix} or C{None}) and C{datum}.

           @raise EcefError: If B{C{latlonh}} not C{LatLon}, L{Ecef9Tuple} or
                             C{scalar} or B{C{lon}} not C{scalar} for C{scalar}
                             B{C{latlonh}} or C{abs(lat)} exceeds 90°.

           @note: Use method C{.forward_} to specify C{lat} and C{lon} in C{radians}
                  and avoid double angle conversions.
        '''
        llhn = _llhn4(latlonh, lon, height, **name)
        return self._forward(*llhn, M=M)

    def forward_(self, phi, lam, height=0, M=False, **name):
        '''Like method C{.forward} except with geodetic lat- and longitude given
           in I{radians}.

           @arg phi: Latitude in I{radians} (C{scalar}).
           @arg lam: Longitude in I{radians} (C{scalar}).
           @kwarg height: Optional height (C{meter}), vertically above (or below)
                          the surface of the ellipsoid.
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with C{lat} set to C{degrees90(B{phi})} and C{lon} to
                    C{degrees180(B{lam})}.

           @raise EcefError: If B{C{phi}} or B{C{lam}} invalid or not C{scalar}.
        '''
        try:  # like function C{_llhn4} below
            plhn = Phi(phi), Lam(lam), Height(height), _name__(name)
        except (TypeError, ValueError) as x:
            raise EcefError(phi=phi, lam=lam, height=height, cause=x)
        return self._forward(*plhn, M=M, _philam=True)

    @property_ROver
    def _Geocentrics(self):
        '''(INTERNAL) Get the valid geocentric classes. I{once}.
        '''
        return (Ecef9Tuple,  # overwrite property_ROver
               _MODS.vector3d.Vector3d)  # _MODS.cartesianBase.CartesianBase

    @property
    def lon00(self):
        '''Get the I{"polar"} longitude (C{degrees}), see method C{reverse}.
        '''
        return self._lon00

    @lon00.setter  # PYCHOK setter!
    def lon00(self, lon00):
        '''Set the I{"polar"} longitude (C{degrees}), see method C{reverse}.
        '''
        self._lon00 = Degrees(lon00=lon00)

    def _Matrix(self, *sa_ca_sb_cb):
        '''(INTERNAL) Create a rotation L{EcefMatrix}.
        '''
        return self._xnamed(EcefMatrix(*sa_ca_sb_cb))

    def _polon(self, y, x, p, **lon00_name):
        '''(INTERNAL) Handle I{"polar"} longitude.
        '''
        return atan2d(y, x) if p else _xkwds_get(lon00_name, lon00=self.lon00)

    def reverse(self, xyz, y=None, z=None, M=False, **lon00_name):
        '''Convert from geocentric C{(x, y, z)} to geodetic C{(lat, lon, height)}.

           @arg xyz: A geocentric (C{Cartesian}, L{Ecef9Tuple}) or C{scalar} ECEF X
                     coordinate (C{meter}).
           @kwarg y: ECEF Y coordinate for C{scalar} B{C{xyz}} and B{C{z}} (C{meter}).
           @kwarg z: ECEF Z coordinate for C{scalar} B{C{xyz}} and B{C{y}} (C{meter}).
           @kwarg M: Optionally, return the rotation L{EcefMatrix} (C{bool}).
           @kwarg lon00_name: Optional C{B{name}=NN} (C{str}) and optional keyword argument
                        C{B{lon00}=INT0} (C{degrees}), an arbitrary I{"polar"} longitude
                        returned if C{B{x}=0} and C{B{y}=0}, see property C{lon00}.

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)} with
                    geodetic coordinates C{(lat, lon, height)} for the given geocentric
                    ones C{(x, y, z)}, case indicator C{C} (C{int} 1..5), rotation matrix
                    C{M} (L{EcefMatrix} or C{None}) and C{datum}.

           @raise EcefError: Invalid B{C{xyz}} or C{scalar} C{x} or B{C{y}} and/or B{C{z}}
                             not C{scalar} for C{scalar} B{C{xyz}}.
        '''
        x, y, z, name = _xyzn4(xyz, y, z, self._Geocentrics, **lon00_name)

        E = self.ellipsoid
        i = None

        sa, ca, sb, cb, h, p, C = _norm7(y, x, z, E)
        if C:  # PYCHOK no cover
            pass  # too high, too far

        elif p < EPS:  # near polar
            p  =  0  # force lon00
            sa = _copysign_1_0(z)
            ca = _0_0
            h  =  fabs(z) - E.b
            C  =  2  # polar

        elif E.e4:  # E.isEllipsoidal
            s, t = _equatorial2(E, p, z)
            r = fsumf_(s, t, -E.e4)
            if t or r > 0:
                try:
                    sa, ca, h, i, C = self._reverse5(1, h, p, z, y, x, r, s, t)
                except (TypeError, ValueError) as X:
                    t = unstr(self.reverse, x=x, y=y, z=z)
                    raise EcefError(t, cause=X)

            else:  # near equatorial plane: e = E.e4 * q == 0 and r <= 0
                # This leads to k = 0 (oblate, equatorial plane) and k + E.e^2 = 0
                # (prolate, rotation axis) and the generation of 0/0 in the general
                # formulas for phi and h, using the general formula and division
                # by 0 in formula for h.  Handle this case by taking the limits:
                #   f > 0: z -> 0, k        ->  E.e2 * sqrt(q) / sqrt(E.e4 - s)
                #   f < 0: r -> 0, k + E.e2 -> -E.e2 * sqrt(q) / sqrt(E.e4 - s)
                sa, ca, h = _equatorial3(E, s, z)
                C = 3  # equatorial

        else:  # E.isSpherical: E.e4 == 0
            # Dealing with underflow in the general case with E.e2 = 0 is
            # difficult.  Origin maps to North pole, same as with ellipsoid.
            sa, ca, _ = _norm3((z if h else _1_0), p)
            h -= E.a
            C  = 4  # spherical

        lat = atan1d(sa, ca)
        # lon00 <https://GitHub.com/mrJean1/PyGeodesy/issues/77>
        lon = self._polon(sb, cb, p, **lon00_name)
        m   = self._Matrix(sa, ca, sb, cb) if M else None
        return Ecef9Tuple(x, y, z, lat, lon, h, C, m, self.datum,
                                   iteration=i, name=self._name__(name))  # PYCHOK return

    def _reverse5(self, *C_h_p_z_y_x_r_s_t):  # PYCHOK no cover
        '''I{Must be overloaded}.'''
        self._notOverloaded(*C_h_p_z_y_x_r_s_t)

    def toStr(self, prec=9, **unused):  # PYCHOK signature
        '''Return this C{Ecef*} as a string.

           @kwarg prec: Precision, number of decimal digits (0..9).

           @return: This C{Ecef*} (C{str}).
        '''
        return self.attrs(_a_, _f_, _datum_, _name_, prec=prec)  # _ellipsoid_


class EcefError(_ValueError):
    '''An ECEF or C{Ecef*} related issue.
    '''
    pass


class EcefFarrell21(_EcefBase):
    '''Conversion between geodetic and geocentric, I{Earth-Centered, Earth-Fixed} (ECEF)
       coordinates based on I{Jay A. Farrell}'s U{Table 2.1<https://Books.Google.com/
       books?id=fW4foWASY6wC>}, page 29, aka the I{Heikkinen application} of U{Ferrari's
       solution<https://WikiPedia.org/wiki/Geographic_coordinate_conversion>}.

       @see: Classes L{EcefFarrell22} and L{EcefVeness}.
    '''
    def _reverse5(self, C, h, p, z, *unused):  # PYCHOK signature
        E   = self.ellipsoid
        a   = E.a
        a2  = E.a2
        b2  = E.b2
        e2  = E.e2
        e2_ = E.e2abs * E.a2_b2  # (E.e * E.a_b)**2 = 0.0820944... WGS84
        e4  = E.e4

        z2 = z**2  # names as page 29
        ez = z2 * (_1_0 - e2)  # E.e2s2(z)

        p2 = p**2
        G  = p2 + ez - e2 * (a2 - b2)  # p2 + ez - e4 * a2
        F  = b2 * z2 * 54
        c  = e4 * p2 * F / G**3
        s  = sqrt(c * (c + _2_0))
        c  = cbrt(s +  c + _1_0)
        G *= fsumf_(c, _1_0, _1_0 / c)  # k
        P  = F / (G**2 * _3_0)
        Q  = sqrt(_2_0 * e4 * P + _1_0)
        Q1 = Q +  _1_0
        s  = fsumf_(a2 * (Q1 / Q) * _0_5,
                    -P * ez / (Q * Q1),
                    -P * p2 * _0_5)
        r = p * P * e2 / Q1 - sqrt(s)
        r = p + r * e2
        v = b2 / (sqrt(r**2 + ez) * a)  # z0 / z

        h  = hypot(r, z) * (_1_0 - v)
        z += e2_ * v * z  # lat = atan1d(z, p)
        return z, p, h, None, C
        # note, phi and lam are swapped on page 29


class EcefFarrell22(_EcefBase):
    '''Conversion between geodetic and geocentric, I{Earth-Centered, Earth-Fixed} (ECEF)
       coordinates based on I{Jay A. Farrell}'s U{Table 2.2<https://Books.Google.com/
       books?id=fW4foWASY6wC>}, page 30.

       @see: Classes L{EcefFarrell21} and L{EcefVeness}.
    '''
    def _reverse5(self, C, h, p, z, *unused):  # PYCHOK signature
        E = self.ellipsoid
        a, b    =  E.a, E.b
        s, c, _ = _norm3(z * a, p * b)  # Bowring
        s, c, _ = _norm3(z + s**3 * b * E.e22,
                         p - c**3 * a * E.e2)
        if c:
            h = p / fabs(c)
            if s:
                h -= E.roc1_(s)
            else:
                h -= a
#               C  = 3  # XXX 1?
        else:
            h = fabs(z) - b
            C = 2
        # lat = atan1d(s, c)
        return s, c, h, None, C
        # note, phi and lam are swapped on page 30


class EcefFukushima(_EcefBase):
    '''Conversion between geodetic and geocentric, I{Earth-Centered, Earth-Fixed} (ECEF) coordinates
       transcoded from I{Toshio Fukushima}'s U{Fortran<https://www.ResearchGate.net/publication/277721539>}
       implementation.

       @see: Fukushima, T. U{Transformation from Cartesian to Geodetic Coordinates Accelerated by
             Halley’s Method<https://www.researchgate.net/publication/227215135>} and Eleiche, M.
             U{A comparison between Fukushima-Halley algorithm and Trilateration algorithm for
             geodetic conversion<https://link.Springer.com/article/10.1007/s12145-022-00779-7>}.
    '''
    def _reverse5(self, C, h, p, z, *unused):  # PYCHOK signature
        E  =  self.ellipsoid
        a  =  E.a
        e2 =  E.e2
        e4 = _1_5 * E.e4  # e4T
        ec = _1_0 - E.f  # sqrt(_1_0 - e2)
#       assert (a * ec) == E.b

        za = fabs(z)
        s0 = za / a
        zc = s0 * ec
        pn = p  / a
        # Newton Correction Factors
        c0 = pn * ec
        c2 = c0**2
        c3 = c2 * c0
        s2 = s0**2
        s3 = s2 * s0
#       a2 = s2 + c2
        a0 = hypot(s0, c0)  # sqrt(a2)
        a3 = a0**3  # a0 * a2
        d0 = a3 * zc + e2 * s3
        f0 = a3 * pn - e2 * c3
        # Halley Correction Factor
        b0 = e4 * s2 * c2 * pn * (a0 - ec)
        sa = d0 * f0 - b0 * s0
        ca = (f0**2  - b0 * c0) * ec

        # lat = atan1d(sa, ca)
        h =  hypot(ec * sa, ca)
        h =  fsumf_(p * ca, za * sa, -h * a)
        h = _over(h, hypot(sa, ca))
        return sa, ca, h, None, C


class EcefKarney(_EcefBase):
    '''Conversion between geodetic and geocentric, I{Earth-Centered, Earth-Fixed} (ECEF) coordinates
       transcoded from I{Karney}'s C++ U{Geocentric<https://GeographicLib.SourceForge.io/C++/doc/
       classGeographicLib_1_1Geocentric.html>} methods.

       @note: In general, there are multiple solutions and the result which minimizes C{height} is
              returned, i.e., the C{(lat, lon)} corresponding to the closest point on the ellipsoid.
              If there are still multiple solutions with different latitudes (applies only if C{z}
              = 0), then the solution with C{lat} > 0 is returned.  If there are still multiple
              solutions with different longitudes (applies only if C{x} = C{y} = 0), then C{lon00}
              is returned.  The returned C{lon} is in the range [−180°, 180°] and C{height} is not
              below M{−E.a * (1 − E.e2) / sqrt(1 − E.e2 * sin(lat)**2)}.  Like C{forward} above,
              M{v1 = Transpose(M) ⋅ v0}.
    '''
    def _reverse5(self, C, h, p, z, y, x, r, s, q):  # PYCHOK unused y, x
        E = self.ellipsoid
        e = E.e4 * q  # p renamed to s
        # Avoid possible division by zero when r = 0 by multiplying
        # equations for s and t by r^3 and r, respectively.
        d  = s = e * s / _4_0  # s = r^3 * s
        u  = r = r / _6_0
        r2 = r**2
        r3 = r2 * r
        t3 = r3 + s
        d *= t3 + r3
        if d < 0:
            # t is complex, but the way u is defined, the result is real.
            # There are three possible cube roots.  We choose the root
            # which avoids cancellation.  Note, d < 0 implies r < 0.
            u += cos(atan2(sqrt(-d), -t3) / _3_0) * r * _2_0
        else:
            # Pick the sign on the sqrt to maximize abs(t3).  This
            # minimizes loss of precision due to cancellation.  The
            # result is unchanged because of the way the t is used
            # in definition of u.
            if d > 0:
                t3 += copysign0(sqrt(d), t3)  # t3 = (r * t)^3
            # N.B. cbrt always returns the real root, cbrt(-8) = -2.
            t = cbrt(t3)  # t = r * t
            if t:  # t can be zero; but then r2 / t -> 0.
                u = fsumf_(u, t, r2 / t)
        v = sqrt(u**2 + e)  # guaranteed positive
        # Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
        # E.e4 * q / (v - u) because u ~ e^4 when q is small and u < 0.
        u = (e / (v - u)) if u < 0 else (u + v)  # u+v, guaranteed positive
        # Need to guard against w going negative due to roundoff in u - q.
        w = E.e2abs * (u - q) / (_2_0 * v)
        # Rearrange expression for k to avoid loss of accuracy due to
        # subtraction.  Division by 0 not possible because u > 0, w >= 0.
        k1 = k2 = (u / (sqrt(w**2 + u) + w)) if w > 0 else sqrt(u)
        if E.f < 0:
            k1 -= E.e2
        else:
            k2 += E.e2
        sa, ca, h = _norm3(z / k1, p / k2)
        h *= k1 - E.e21
        return sa, ca, h, None, C


class EcefSudano(_EcefBase):
    '''Conversion between geodetic and geocentric, I{Earth-Centered, Earth-Fixed} (ECEF) coordinates
       based on I{John J. Sudano}'s U{paper<https://www.ResearchGate.net/publication/3709199>}.
    '''
    _TOL = \
    _tol = EPS

    def reverse(self, xyz, y=None, z=None, M=False, tol=EPS, **lon00_name):  # PYCHOK tol
        '''Convert from geocentric C{(x, y, z)} to geodetic C{(lat, lon, height)} using
           I{Sudano}'s U{iterative method<https://www.ResearchGate.net/publication/3709199>}.

           @kwarg tol: Convergence tolerance for C{sin(latitude)} (C{scalar}).

           @see: L{Parent method<_EcefBase.reverse>} for all other information.

           @raise EcefError: Invalid B{C{xyz}} or C{scalar} C{x} or B{C{y}} and/or B{C{z}} not
                             C{scalar} for C{scalar} B{C{xyz}} or no convergence for C{B{tol}}.
        '''
        if tol != self._TOL:
            self._tol = Scalar_(tol=tol, low=EPS, Error=EcefError)
        return _EcefBase.reverse(self, xyz, y=y, z=z, M=M, **lon00_name)

    def _reverse5(self, C, h, p, z, *unused):  # PYCHOK signature
        E = self.ellipsoid
        e = E.e2 * E.a
        d = e - p

        sa, ca, _ = _norm3(fabs(z), p * E.e21)
        # Sudano's Eq (A-6) and (A-7) refactored/reduced,
        # replacing Rn from Eq (A-4) with n = E.a / ca:
        # N = ca**2 * ((z + E.e2 * n * sa) * ca - p * sa)
        #   = ca**2 * (z * ca + E.e2 * E.a * sa - p * sa)
        #   = ca**2 * (z * ca + (E.e2 * E.a - p) * sa)
        # D = ca**3 * (E.e2 * n / E.e2s2(sa)) - p
        #   = ca**2 * (E.e2 * E.a / E.e2s2(sa) - p / ca**2)
        # N / D = (z * ca + (E.e2 * E.a - p) * sa) /
        #         (E.e2 * E.a / E.e2s2(sa) - p / ca**2)
        tol = self._tol
        _S2 = Fsum(sa).fsum2f_
        for i in range(1, _TRIPS):  # 6+ max
            ca2 = _1_0 - sa**2
            if ca2 < EPS02:
                break
            D = p / ca2 - e / E.e2s2(sa)
            if fabs(D) < EPS0:
                break
            ca = sqrt(ca2)
            sa, D = _S2(z * ca / D, d * sa / D)
            if fabs(D) < tol:
                break
        else:  # PYCHOK no cover
            raise ValueError(Fmt.no_convergence(fabs(D), tol))

        sa = copysign0(sa, z)
        # lat = atan1d(sa, ca)
        # h = (fabs(z) + p - E.a * cos(a + E.e21) * sa / ca) / (ca + sa)
        # Sudano's Eq (7) doesn't produce the correct height, ...
        h = E._heightB(sa, ca, z, p)  # ... use Veness' (Bowring eqn 7)
        return sa, ca, h, i, C

    @deprecated_property
    def tolerance(self):
        '''DEPRECATED on 2025.08.22, use keyword argument C{tol}.'''
        return self._tol

    @tolerance.setter  # PYCHOK setter!
    def tolerance(self, tol):
        self._tol = Scalar_(tolerance=tol, low=EPS, Error=EcefError)


class EcefUPC(_EcefBase):
    '''Conversion between geodetic and geocentric, I{Earth-Centered, Earth-Fixed} (ECEF) coordinates based on
       I{UPC}'s U{method<https://GSSC.ESA.int/navipedia/index.php/Ellipsoidal_and_Cartesian_Coordinates_Conversion>}.
    '''
    _TOL = _TOL
    _tol = _Degrees2Radians(_TOL)

    def reverse(self, xyz, y=None, z=None, M=False, tol=_TOL, **lon00_name):  # PYCHOK tol
        '''Convert from geocentric C{(x, y, z)} to geodetic C{(lat, lon, height)} using I{UPC}'s
           U{iterative method<https://GSSC.ESA.int/navipedia/GNSS_Book/ESA_GNSS-Book_TM-23_Vol_I.pdf>}, page 186.

           @kwarg tol: Convergence tolerance for the C{latitude} (C{degrees}).

           @see: L{Parent method<_EcefBase.reverse>} for all other information.

           @raise EcefError: Invalid B{C{xyz}} or C{scalar} C{x} or B{C{y}} and/or B{C{z}} not
                             C{scalar} for C{scalar} B{C{xyz}} or no convergence for C{B{tol}}.
        '''
        if tol != _TOL:
            self._tol = _Degrees2Radians(tol)
        return _EcefBase.reverse(self, xyz, y=y, z=z, M=M, **lon00_name)

    def _reverse5(self, C, h, p, z, *unused):  # PYCHOK signature
        E  = self.ellipsoid
        a  = E.a
        e2 = E.e2  # signed

        za  = fabs(z)
        ph_ = atan1(za, E.e21 * p)
        tol = self._tol
        for i in range(1, _TRIPS):  # 5..6 max
            s, c = sincos2(ph_)
            N  = a / sqrt(_1_0 - s**2 * e2)  # N + h == N + p / c - N == p / c
            ca = p - N * c * e2  # == p * (1 - N * e2 / (N + h)) == p * (1 - N * e2 * c / p)
            ph = atan1(za, ca)  # atan1(z / p, 1 - N * e2 / (N + h)) == atan1(z, ca)
            r  = fabs(ph - ph_)
            if r < tol:
                # lat = copysign0(degrees(ph), z)
                #    == atan1d(z, ca)
                h = p / c - N
                break
            ph_ = ph
        else:  # PYCHOK no cover
            r, tol = map1(degrees, r, tol)
            raise ValueError(Fmt.no_convergence(r, tol))
        return z, ca, h, i, C


class EcefVeness(_EcefBase):
    '''Conversion between geodetic and geocentric, I{Earth-Centered, Earth-Fixed} (ECEF) coordinates
       transcoded from I{Chris Veness}' JavaScript classes U{LatLonEllipsoidal, Cartesian<https://
       www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}.

       @note: Uses B. R. Bowring’s formulation for μm precision in concise form U{The accuracy of
              geodetic latitude and height equations<https://www.ResearchGate.net/publication/233668213>},
              Survey Review, Vol 28, 218, Oct 1985.

       @see: U{A Guide to Coordinate Systems in Great Britain<https://www.OrdnanceSurvey.co.UK/documents/
             resources/guide-coordinate-systems-great-britain.pdf>}, section I{B) Converting between 3D
             Cartesian and ellipsoidal latitude, longitude and height coordinates}.

       @see: Toms, Ralph M. U{An Efficient Algorithm for Geocentric to Geodetic Coordinate Conversion
             <https://www.OSTI.gov/scitech/biblio/110235>}, Sept 1995 and U{An Improved Algorithm for
             Geocentric to Geodetic Coordinate Conversion<https://www.OSTI.gov/scitech/servlets/purl/231228>},
             Apr 1996, both from Lawrence Livermore National Laboratory (LLNL).
    '''
    def _reverse5(self, C, h, p, z, *unused):  # PYCHOK signature
        # assert h >= p > 0  # h = hypot(z, p)
        E =  self.ellipsoid
        a =  E.a
        B =  E.b * E.e22
        # parametric latitude (Bowring eqn 17, replaced)
        t = (E.b * z) / (a * p) * (B / h + _1_0)  # theta
        c = _1_0 / hypot1(t)  # t == atan2(z * a, p * E.b)
        s =  c * t  # s, c == sincos2(t)
        # geodetic latitude (Bowring eqn 18)
        sa, ca, _ = _norm3(z + s**3 * B,
                           p - c**3 * a * E.e2)
        h = E._heightB(sa, ca, z, p)  # height (Bowring eqn 7)
        # lat = atan1d(sa, ca)
        return sa, ca, h, None, C


class EcefYou(_EcefBase):
    '''Conversion between geodetic and geocentric, I{Earth-Centered, Earth-Fixed} (ECEF) coordinates
       using I{Rey-Jer You}'s U{transformation<https://www.ResearchGate.net/publication/240359424>}
       for I{non-prolate} ellipsoids.

       @see: Featherstone, W.E., Claessens, S.J. U{Closed-form transformation between geodetic and
             ellipsoidal coordinates<https://Espace.Curtin.edu.AU/bitstream/handle/20.500.11937/11589/
             115114_9021_geod2ellip_final.pdf>} Studia Geophysica et Geodaetica, 2008, 52, pages 1-18
             and U{PyMap3D<https://PyPI.org/project/pymap3d>}.
    '''
    _isYou = True

    def _reverse5(self, C, h, p, z, y, x, *unused):  # PYCHOK signature
        E = self.ellipsoid
        a, b  = E.a, E.b
        e, e2 = self._e_e2

        u  =  hypot2_(x, y, z) - e2
        u +=  hypot(u, e * z * _2_0)
        u *= _0_5
        if u > EPS02:
            u = sqrt(u)
            q = hypot(u, e)
            B = atan1(q * z, u * p)  # beta0 = atan(q / u * z / p)
            sB, cB = sincos2(B)
            if cB and sB:
                q *= a
                d  = (q / cB - e2 * cB) / sB
                if isnon0(d):
                    B += fsumf_(u * b, -q, e2) / d
                    sB, cB = sincos2(B)
        elif u < (-EPS02):
            raise EcefError(u=u, txt=_singular_)
        else:  # near polar  # PYCHOK no cover
            sB, cB, C = _copysign_1_0(z), _0_0, 2

        h = hypot(p - a * cB, z - b * sB)
        if hypot2_(x, y, z * E.a_b) < E.a2:  # or lat < 0 or z < 0
            h = neg(h)  # inside ellipsoid
        # lat = atand(E.a_b * tan(B)) == atan1d(a * sB, b * cB)
        return (a * sB), (b * cB), h, None, C


class EcefMatrix(_NamedTuple):
    '''A rotation matrix known as I{East-North-Up (ENU) to ECEF}.

       @see: U{From ENU to ECEF<https://WikiPedia.org/wiki/
             Geographic_coordinate_conversion#From_ECEF_to_ENU>} and
             U{Issue #74<https://Github.com/mrJean1/PyGeodesy/issues/74>}.
    '''
    _Names_ = ('_0_0_', '_0_1_', '_0_2_',  # row-major order
               '_1_0_', '_1_1_', '_1_2_',
               '_2_0_', '_2_1_', '_2_2_')
    _Units_ = (Scalar,) * len(_Names_)

    def _validate(self, **unused):  # PYCHOK unused
        '''(INTERNAL) Allow C{_Names_} with leading underscore.
        '''
        _NamedTuple._validate(self, underOK=True)

    def __new__(cls, sa, ca, sb, cb, *_more, **name):
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

        elif max(map(fabs, t)) > _1_0:
            raise EcefError(unstr(EcefMatrix, *t))

        else:  # build matrix from the following quaternion operations
            #   qrot(lam, [0,0,1]) * qrot(phi, [0,-1,0]) * [1,1,1,1]/2
            # or
            #   qrot(pi/2 + lam, [0,0,1]) * qrot(-pi/2 + phi, [-1,0,0])
            # where
            #   qrot(t,v) = [cos(t/2), sin(t/2)*v[1], sin(t/2)*v[2], sin(t/2)*v[3]]

            # Local X axis (East) in geocentric coords
            #  M[0] = -slam;        M[3] =  clam;        M[6] = 0;
            # Local Y axis (North) in geocentric coords
            #  M[1] = -clam * sphi; M[4] = -slam * sphi; M[7] = cphi;
            # Local Z axis (Up) in geocentric coords
            #  M[2] =  clam * cphi; M[5] =  slam * cphi; M[8] = sphi;
            t = (-sb, -cb * sa, cb * ca,
                  cb, -sb * sa, sb * ca,
                _0_0,       ca,      sa)

        return _NamedTuple.__new__(cls, *t, **name)

    def column(self, column):
        '''Get this matrix' B{C{column}} 0, 1 or 2 as C{3-tuple}.
        '''
        if 0 <= column < 3:
            return self[column::3]
        raise _IndexError(column=column)

    @property_RO
    def _columns(self):
        for c in range(3):
            yield self[c::3]

    def copy(self, **unused):  # PYCHOK signature
        '''Make a shallow or deep copy of this instance.

           @return: The copy (C{This class} or subclass thereof).
        '''
        return self.classof(*self)

    __copy__ = __deepcopy__ = copy

    @Property_RO
    def matrix3(self):
        '''Get this matrix' rows (C{3-tuple} of 3 C{3-tuple}s).
        '''
        return tuple(self._rows)

    @Property_RO
    def matrixTransposed3(self):
        '''Get this matrix' I{Transposed} rows (C{3-tuple} of 3 C{3-tuple}s).
        '''
        return tuple(self._columns)

    def multiply(self, other):
        '''Matrix multiply M{M0' ⋅ M} this matrix I{Transposed} with an other matrix.

           @arg other: The other matrix (L{EcefMatrix}).

           @return: The matrix product (L{EcefMatrix}).

           @raise TypeError: If B{C{other}} is not an L{EcefMatrix}.
        '''
        _xinstanceof(EcefMatrix, other=other)
        # like LocalCartesian.MatrixMultiply, C{self.matrixTransposed3 X other.matrix3}
        # <https://GeographicLib.SourceForge.io/C++/doc/LocalCartesian_8cpp_source.html>
        X = (_fdotf(t, *c) for t in self._columns for c in other._columns)
        return _xnamed(EcefMatrix(*X), typename(EcefMatrix.multiply))

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
            xyz = tuple(s - s0 for s, s0 in zip(xyz, xyz0))

        # x' = M[0] * x + M[3] * y + M[6] * z
        # y' = M[1] * x + M[4] * y + M[7] * z
        # z' = M[2] * x + M[5] * y + M[8] * z
        return tuple(_fdotf(xyz, *c) for c in self._columns)

    def row(self, row):
        '''Get this matrix' B{C{row}} 0, 1 or 2 as C{3-tuple}.
        '''
        if 0 <= row < 3:
            r = row * 3
            return self[r:r+3]
        raise _IndexError(row=row)

    @property_RO
    def _rows(self):
        for r in (0, 3, 6):
            yield self[r:r+3]

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
            _xyz = _1_0_1T + xyz
            # x' = x0 + M[0] * x + M[1] * y + M[2] * z
            # y' = y0 + M[3] * x + M[4] * y + M[5] * z
            # z' = z0 + M[6] * x + M[7] * y + M[8] * z
            xyz_ = (_fdotf(_xyz, s, *r) for s, r in zip(xyz0, self._rows))
        else:
            # x' = M[0] * x + M[1] * y + M[2] * z
            # y' = M[3] * x + M[4] * y + M[5] * z
            # z' = M[6] * x + M[7] * y + M[8] * z
            xyz_ = (_fdotf(xyz, *r) for r in self._rows)
        return tuple(xyz_)


class Ecef9Tuple(_NamedTuple, _EcefLocal):
    '''9-Tuple C{(x, y, z, lat, lon, height, C, M, datum)} with I{geocentric} C{x},
       C{y} and C{z} plus I{geodetic} C{lat}, C{lon} and C{height}, case C{C} and
       optionally, rotation matrix C{M} (L{EcefMatrix} or C{None}) and C{datum},
       with C{lat} and C{lon} in C{degrees} and C{x}, C{y}, C{z} and C{height} in
       C{meter}, conventionally.  Case C{C=0} means C{x, y,z} from foward, C{C=1}
       C{lat, lon, height} from reverse, C{C=2} near-polar C{lat, lon}, C{C=3}
       near-equatorial C{lat, lon}, C{C=4} spherical C{lat, lon} and C{C=5} means
       the C{height} exceeds C{datum}'s C{ellipsoid.heightMax}.
    '''
    _Names_ = (_x_,   _y_,   _z_,   _lat_, _lon_, _height_, _C_,  _M_,   _datum_)
    _Units_ = ( Meter, Meter, Meter, Lat,   Lon,   Height,   Int, _Pass, _Pass)

    @property_ROver
    def _CartesianBase(self):
        '''(INTERNAL) Get class C{CartesianBase}, I{once}.
        '''
        return _MODS.cartesianBase.CartesianBase  # overwrite property_ROver

    @deprecated_method
    def convertDatum(self, datum2):  # for backward compatibility
        '''DEPRECATED, use method L{toDatum}.'''
        return self.toDatum(datum2)

    @property_RO
    def _ecef9(self):  # in ._EcefLocal._Ltp_ecef2local
        return self

    @property_RO
    def ellipsoid(self):
        '''Get the ellipsoid (L{Ellipsoid}).
        '''
        return (self.datum or _WGS84).ellipsoid

    @Property_RO
    def lam(self):
        '''Get the longitude in C{radians} (C{float}).
        '''
        return self.philam.lam

    @Property_RO
    def lamVermeille(self):
        '''Get the longitude in C{radians} M{[-PI*3/2..+PI*3/2]} after U{Vermeille
           <https://Search.ProQuest.com/docview/639493848>} (2004), page 95.

           @see: U{Karney<https://GeographicLib.SourceForge.io/C++/doc/geocentric.html>},
                 U{Vermeille<https://Search.ProQuest.com/docview/847292978>} 2011, pp 112-113, 116
                 and U{Featherstone, et.al.<https://Search.ProQuest.com/docview/872827242>}, page 7.
        '''
        x, y = self.x, self.y
        a = fabs(y)
        if a > EPS0:
            r = PI_2 - atan2(x, hypot(x, a) + a) * _2_0
            if y < 0:
                r = -r
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
           <https://Search.ProQuest.com/docview/639493848>} 2004, page 95.

           @see: Property C{lamVermeille}.
        '''
        return Lon(Vermeille=degrees(self.lamVermeille))

    @Property_RO
    def Mx(self):
        '''Compute rotation matrix (L{EcefMatrix}), seperate from C{M}.
        '''
        sa, ca, sb, cb, _, _, _ = _norm7(self.y, self.x, self.z, self.ellipsoid)
        return EcefMatrix(sa, ca, sb, cb, name=self.name)

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

    phiVermeille = phi

    def toCartesian(self, Cartesian=None, **Cartesian_kwds):
        '''Return the geocentric C{(x, y, z)} coordinates as an ellipsoidal or spherical
           C{Cartesian}.

           @kwarg Cartesian: Optional class to return C{(x, y, z)} (L{ellipsoidalKarney.Cartesian},
                             L{ellipsoidalNvector.Cartesian}, L{ellipsoidalVincenty.Cartesian},
                             L{sphericalNvector.Cartesian} or L{sphericalTrigonometry.Cartesian})
                             or C{None}.
           @kwarg Cartesian_kwds: Optionally, additional B{C{Cartesian}} keyword arguments, ignored
                                  if C{B{Cartesian} is None}.

           @return: A B{C{Cartesian}} instance or a L{Vector4Tuple}C{(x, y, z, h)} if C{B{Cartesian}
                    is None}.

           @raise TypeError: Invalid B{C{Cartesian}} or B{C{Cartesian_kwds}} item.
        '''
        if _isin(Cartesian, None, Vector4Tuple):
            r = self.xyzh
        elif Cartesian is Vector3Tuple:
            r = self.xyz
        else:
            _xsubclassof(self._CartesianBase, Cartesian=Cartesian)
            r = Cartesian(self, **_name1__(Cartesian_kwds, _or_nameof=self))
        return r

    def toDatum(self, datum2, **name):
        '''Convert this C{Ecef9Tuple} to an other datum.

           @arg datum2: Datum to convert I{to} (L{Datum}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: The converted 9-Tuple (C{Ecef9Tuple}).

           @raise TypeError: The B{C{datum2}} is not a L{Datum}.
        '''
        n = _name__(name, _or_nameof=self)
        if _isin(self.datum, None, datum2):  # PYCHOK _Names_
            r = self.copy(name=n)
        else:
            c = self._CartesianBase(self, datum=self.datum, name=n)  # PYCHOK _Names_
            # c.toLatLon converts datum, x, y, z, lat, lon, etc.
            # and returns another Ecef9Tuple iff LatLon is None
            r = c.toLatLon(datum=datum2, LatLon=None)
        return r

    def toLatLon(self, LatLon=None, **LatLon_kwds):
        '''Return the geodetic C{(lat, lon, height[, datum])} coordinates.

           @kwarg LatLon: Optional class to return C{(lat, lon, height[, datum])} or C{None}.
           @kwarg LatLon_kwds: Optional B{C{height}}, B{C{datum}} and other B{C{LatLon}}
                               keyword arguments.

           @return: A B{C{LatLon}} instance or if C{B{LatLon} is None}, a L{LatLon4Tuple}C{(lat,
                    lon, height, datum)} or L{LatLon3Tuple}C{(lat, lon, height)} if C{datum} is
                    specified or not.

           @raise TypeError: Invalid B{C{LatLon}} or B{C{LatLon_kwds}} item.
        '''
        lat, lon, D = self.lat, self.lon, self.datum  # PYCHOK Ecef9Tuple
        kwds = _name1__(LatLon_kwds, _or_nameof=self)
        kwds = _xkwds(kwds, height=self.height, datum=D)  # PYCHOK Ecef9Tuple
        d    =  kwds.get(_datum_, LatLon)
        if LatLon is None:
            r = LatLon3Tuple(lat, lon, kwds[_height_], name=kwds[_name_])
            if d is not None:
                # assert d is not LatLon
                r = r.to4Tuple(d)  # checks type(d)
        else:
            if d is None:
                _ = kwds.pop(_datum_)  # remove None datum
            r = LatLon(lat, lon, **kwds)
        _xdatum(_xattr(r, datum=D), D)
        return r

    def toVector(self, Vector=None, **Vector_kwds):
        '''Return these geocentric C{(x, y, z)} coordinates as vector.

           @kwarg Vector: Optional vector class to return C{(x, y, z)} or C{None}.
           @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword arguments,
                               ignored if C{B{Vector} is None}.

           @return: A B{C{Vector}} instance or a L{Vector3Tuple}C{(x, y, z)} if
                    C{B{Vector} is None}.

           @raise TypeError: Invalid B{C{Vector}} or B{C{Vector_kwds}} item.

           @see: Propertes C{xyz} and C{xyzh}
        '''
        return self.xyz if Vector is None else Vector(
              *self.xyz, **_name1__(Vector_kwds, _or_nameof=self))  # PYCHOK Ecef9Tuple

#   def _T_x_M(self, T):
#       '''(INTERNAL) Update M{self.M = T.multiply(self.M)}.
#       '''
#       return self.dup(M=T.multiply(self.M))

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


def _4Ecef(this, Ecef):  # in .datums.Datum.ecef, .ellipsoids.Ellipsoid.ecef
    '''Return an ECEF converter for C{this} L{Datum} or L{Ellipsoid}.
    '''
    if Ecef is None:
        Ecef = EcefKarney
    else:
        _xinstanceof(*_Ecefs, Ecef=Ecef)
    return Ecef(this, name=this.name)


def _equatorial2(E, p, z):
    '''(INTERNAL) Equatorial plane from C{EcefKarney}.
    '''
    # Treat prolate spheroids by swapping p and z here and by
    # switching the arguments to phi = atan2(...) at the end
    # of method C{EcefKarney._reverse5}
    p = (p / E.a)**2
    q = (z / E.a)**2 * E.e21
    return (q, p) if E.f < 0 else (p, q)


def _equatorial3(E, s, z):
    '''(INTERNAL) Equatorial plane from C{EcefKarney}.
    '''
    t = E.e4 - s
    if E.f < 0:
        s, t = t, s
        e = E.a
    else:
        e = E.b2_a
    sa, ca, h = _norm3(*map1(sqrt, E._1_e21 * t, s))
    if z < 0:  # for tiny negative z, not for prolate
        sa = neg(sa)
    h *= neg(e / E.e2abs)
    return sa, ca, h


def _llhn4(latlonh, lon, height, suffix=NN, Error=EcefError, **name):  # in .ltp
    '''(INTERNAL) Get a C{(lat, lon, h, name)} 4-tuple.
    '''
    try:
        lat, lon = latlonh.lat, latlonh.lon
        h = _xattr(latlonh, height=_xattr(latlonh, h=height))
        n = _name__(name, _or_nameof=latlonh)  # == latlonh._name__(name)
    except AttributeError:
        lat, h, n = latlonh, height, _name__(**name)
    try:
        return Lat(lat), Lon(lon), Height(h), n
    except (TypeError, ValueError) as x:
        t = _lat_, _lon_, _height_
        if suffix:
            t = (_ + suffix for _ in t)
        d = dict(zip(t, (lat, lon, h)))
        raise Error(cause=x, **d)


def _norm3(y, x, eps=0):
    '''(INTERNAL) Return C{y, x, h} normalized.
    '''
    h = hypot(y, x)  # EPS0, EPS_2
    return (y / h, x / h, h) if h > eps else (_0_0, _1_0, h)  # copysign_1_0(x)


def _norm7(y, x, z=0, E=_EWGS84):
    '''(INTERNAL) Return C{phi, lam, h, p, C}.
    '''
    sb, cb, p = _norm3(y, x)  # lam, distance to polar axis
    sa, ca, h = _norm3(z, p)  # phi, distance to earth center
    if h > E.heightMax:
        # We are really far away (> 12M light years).  Treat the earth
        # as a point and h above as an acceptable approximation to the
        # height.  This avoids overflow, e.g., in the computation of d
        # below.  It's possible that h has overflowed to INF, that's OK.
        # Treat finite x, y, but R overflows to +INF by scaling by 2.
        sb, cb, p = _norm3(y * _0_5, x * _0_5)
        sa, ca, _ = _norm3(z * _0_5, p)
        C = 5
    else:
        C = 0
    return sa, ca, sb, cb, h, p, C


def _xEcef(Ecef):  # PYCHOK .latlonBase
    '''(INTERNAL) Validate B{C{Ecef}} I{class}.
    '''
    if issubclassof(Ecef, _EcefBase):
        return Ecef
    raise _TypesError(_Ecef_, Ecef, *_Ecefs)


# kwd lon00 unused but will throw a TypeError if misspelled, etc.
def _xyzn4(xyz, y, z, Types, Error=EcefError, lon00=0,  # PYCHOK unused, in pychlv
                     _xyz_y_z_names=_xyz_y_z, **name):  # in .ltp
    '''(INTERNAL) Get an C{(x, y, z, name)} 4-tuple.
    '''
    try:
        n = _name__(name, _or_nameof=xyz)  # == xyz._name__(name)
        try:
            t = xyz.x, xyz.y, xyz.z, n
            if not isinstance(xyz, Types):
                raise _TypesError(_xyz_y_z_names[0], xyz, *Types)
        except AttributeError:
            t = map1(float, xyz, y, z) + (n,)
    except (TypeError, ValueError) as x:
        d = dict(zip(_xyz_y_z_names, (xyz, y, z)))
        raise Error(cause=x, **d)
    return t
# assert _xyz_y_z == _args_kwds_names(_xyzn4)[:3]


_Ecefs = tuple(_ for _ in locals().values()
                       if issubclassof(_, _EcefBase) and
                     _ is not _EcefBase)
__all__ += _ALL_DOCS(_EcefBase)

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
