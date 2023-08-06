
# -*- coding: utf-8 -*-

u'''(INTERNAL) Private spherical base classes C{CartesianSphericalBase} and
C{LatLonSphericalBase} for L{sphericalNvector} and L{sphericalTrigonometry}.

A pure Python implementation of geodetic (lat-/longitude) functions,
transcoded in part from JavaScript originals by I{(C) Chris Veness 2011-2016}
and published under the same MIT Licence**, see
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isbool, isinstanceof, map1
from pygeodesy.cartesianBase import CartesianBase,  Bearing2Tuple
from pygeodesy.constants import EPS, PI, PI2, PI_2, R_M, \
                               _umod_360, isnear0, isnon0, _0_0, \
                               _0_5, _1_0, _180_0
from pygeodesy.datums import Datums, _spherical_datum
from pygeodesy.errors import IntersectionError, _ValueError, _xError
from pygeodesy.fmath import favg, fdot, hypot, sqrt_a
from pygeodesy.interns import NN, _COMMA_, _concentric_, _datum_, \
                             _distant_, _exceed_PI_radians_, _name_, \
                             _near_, _radius_, _too_
from pygeodesy.latlonBase import LatLonBase, _trilaterate5  # PYCHOK passed
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
# from pygeodesy.namedTuples import Bearing2Tuple  # from .cartesianBase
from pygeodesy.nvectorBase import NvectorBase,  Fmt, _xattrs
from pygeodesy.props import deprecated_method, property_doc_, \
                            property_RO, _update_all
# from pygeodesy.streprs import Fmt, _xattrs  # from .nvectorBase
from pygeodesy.units import Bearing, Bearing_, Radians_, Radius, \
                            Radius_, Scalar_
from pygeodesy.utily import acos1, atan2b, atan2d, degrees90, \
                            degrees180, sincos2, sincos2d, tanPI_2_2, \
                            _unrollon, wrap360, wrapPI

from math import cos, fabs, log, sin, sqrt

__all__ = _ALL_LAZY.sphericalBase
__version__ = '23.07.01'


def _angular(distance, radius, low=EPS):  # PYCHOK in .spherical*
    '''(INTERNAL) Return an angular distance in C{radians}.

       @raise UnitError: Invalid B{C{distance}} or B{C{radius}}.
    '''
    r = float(distance)
    if radius:
        r = r / Radius_(radius=radius)  # /= chokes PyChecker
    if low is not None:
        # small near0 values from .rhumbDestination not exact OK
        r = _0_0 if low < 0 and low < r < 0 else Radians_(r, low=low)
    return r


def _rads3(rad1, rad2, radius):  # in .sphericalTrigonometry
    '''(INTERNAL) Convert radii to radians.
    '''
    r1 = Radius_(rad1=rad1)
    r2 = Radius_(rad2=rad2)
    if radius is not None:  # convert radii to radians
        r1 = _angular(r1, radius)
        r2 = _angular(r2, radius)

    x = r1 < r2
    if x:
        r1, r2 = r2, r1
    if r1 > PI:
        raise IntersectionError(rad1=rad1, rad2=rad2,
                                txt=_exceed_PI_radians_)
    return r1, r2, x


class CartesianSphericalBase(CartesianBase):
    '''(INTERNAL) Base class for spherical C{Cartesian}s.
    '''
    _datum = Datums.Sphere  # L{Datum}

    def intersections2(self, rad1, other, rad2, radius=R_M):
        '''Compute the intersection points of two circles each defined
           by a center point and a radius.

           @arg rad1: Radius of the this circle (C{meter} or C{radians},
                      see B{C{radius}}).
           @arg other: Center of the other circle (C{Cartesian}).
           @arg rad2: Radius of the other circle (C{meter} or C{radians},
                      see B{C{radius}}).
           @kwarg radius: Mean earth radius (C{meter} or C{None} if both
                          B{C{rad1}} and B{C{rad2}} are given in C{radians}).

           @return: 2-Tuple of the intersection points, each C{Cartesian}.
                    For abutting circles, the intersection points are the
                    same C{Cartesian} instance, aka the I{radical center}.

           @raise IntersectionError: Concentric, antipodal, invalid or
                                     non-intersecting circles.

           @raise TypeError: If B{C{other}} is not C{Cartesian}.

           @raise ValueError: Invalid B{C{rad1}}, B{C{rad2}} or B{C{radius}}.

           @see: U{Calculating intersection of two Circles
                 <https://GIS.StackExchange.com/questions/48937/
                 calculating-intersection-of-two-circles>} and method
                 or function C{trilaterate3d2}.
        '''
        x1, x2 = self, self.others(other)
        r1, r2, x = _rads3(rad1, rad2, radius)
        if x:
            x1, x2 = x2, x1
        try:
            n, q = x1.cross(x2), x1.dot(x2)
            n2, q1 = n.length2, (_1_0 - q**2)
            if n2 < EPS or isnear0(q1):
                raise ValueError(_near_(_concentric_))
            c1, c2 = cos(r1), cos(r2)
            x0 = x1.times((c1 - q * c2) / q1).plus(
                 x2.times((c2 - q * c1) / q1))
            n1 = _1_0 - x0.length2
            if n1 < EPS:
                raise ValueError(_too_(_distant_))
        except ValueError as x:
            raise IntersectionError(center=self, rad1=rad1,
                                    other=other, rad2=rad2, cause=x)
        n = n.times(sqrt(n1 / n2))
        if n.length > EPS:
            x1 = x0.plus(n)
            x2 = x0.minus(n)
        else:  # abutting circles
            x1 = x2 = x0

        return (_xattrs(x1, self, _datum_, _name_),
                _xattrs(x2, self, _datum_, _name_))


class LatLonSphericalBase(LatLonBase):
    '''(INTERNAL) Base class for spherical C{LatLon}s.
    '''
    _datum = Datums.Sphere  # spherical L{Datum}

    def __init__(self, latlonh, lon=None, height=0, datum=None, wrap=False, name=NN):
        '''Create a spherical C{LatLon} point frome the given lat-, longitude and
           height on the given datum.

           @arg latlonh: Latitude (C{degrees} or DMS C{str} with N or S suffix) or
                         a previous C{LatLon} instance provided C{B{lon}=None}.
           @kwarg lon: Longitude (C{degrees} or DMS C{str} with E or W suffix) or
                       C(None), indicating B{C{latlonh}} is a C{LatLon}.
           @kwarg height: Optional height above (or below) the earth surface (C{meter},
                          same units as the datum's ellipsoid axes or radius).
           @kwarg datum: Optional, spherical datum to use (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2}, L{a_f2Tuple}) or earth radius in C{meter},
                         conventionally).
           @kwarg wrap: If C{True}, wrap or I{normalize} B{C{lat}} and B{C{lon}}
                        (C{bool}).
           @kwarg name: Optional name (string).

           @raise TypeError: If B{C{latlonh}} is not a C{LatLon} or B{C{datum}} not
                             spherical.
        '''
        LatLonBase.__init__(self, latlonh, lon=lon, height=height, wrap=wrap, name=name)
        if datum not in (None, self.datum):
            self.datum = datum

    def bearingTo2(self, other, wrap=False, raiser=False):
        '''Return the initial and final bearing (forward and reverse
           azimuth) from this to an other point.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: A L{Bearing2Tuple}C{(initial, final)}.

           @raise TypeError: The B{C{other}} point is not spherical.

           @see: Methods C{initialBearingTo} and C{finalBearingTo}.
        '''
        # .initialBearingTo is inside .-Nvector and .-Trigonometry
        i = self.initialBearingTo(other, wrap=wrap, raiser=raiser)  # PYCHOK .initialBearingTo
        f = self.finalBearingTo(  other, wrap=wrap, raiser=raiser)
        return Bearing2Tuple(i, f, name=self.name)

    @property_doc_(''' this point's datum (L{Datum}).''')
    def datum(self):
        '''Get this point's datum (L{Datum}).
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set this point's datum I{without conversion} (L{Datum}, L{Ellipsoid},
           L{Ellipsoid2}, L{a_f2Tuple}) or C{scalar} spherical earth radius).

           @raise TypeError: If B{C{datum}} invalid or not not spherical.
        '''
        d = _spherical_datum(datum, name=self.name, raiser=_datum_)
        if self._datum != d:
            _update_all(self)
            self._datum = d

    def finalBearingTo(self, other, wrap=False, raiser=False):
        '''Return the final bearing (reverse azimuth) from this to
           an other point.

           @arg other: The other point (spherical C{LatLon}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{other}} point (C{bool}).

           @return: Final bearing (compass C{degrees360}).

           @raise TypeError: The B{C{other}} point is not spherical.

           @example:

            >>> p = LatLon(52.205, 0.119)
            >>> q = LatLon(48.857, 2.351)
            >>> b = p.finalBearingTo(q)  # 157.9
        '''
        p = self.others(other)
        if wrap:
            p = _unrollon(self, p, wrap=wrap)
        # final bearing is the reverse of the other, initial one;
        # .initialBearingTo is inside .-Nvector and .-Trigonometry
        b = p.initialBearingTo(self, wrap=False, raiser=raiser)
        return _umod_360(b + _180_0)

    def intersecant2(self, circle, point, bearing, radius=R_M, exact=False,
                                                   height=None, wrap=False):  # was=True
        '''Compute the intersections of a circle and a line.

           @arg circle: Radius of the circle centered at this location
                       (C{meter}, same units as B{C{radius}}) or a point
                       on the circle (this C{LatLon}).
           @arg point: An other point in- or outside the circle (this C{LatLon}).
           @arg bearing: Bearing at the B{C{point}} (compass C{degrees360})
                         or a second point on the line (this C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}, conventionally).
           @kwarg exact: If C{True} use the I{exact} rhumb methods for azimuth,
                         destination and distance, if C{False} use the basic
                         rhumb methods (C{bool}) or if C{None} use the I{great
                         circle} methods.
           @kwarg height: Optional height for the intersection points (C{meter},
                          conventionally) or C{None}.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{point}}, B{C{circle}} and/or B{C{bearing}} (C{bool}).

           @return: 2-Tuple of the intersection points (representing a chord),
                    each an instance of this class.  For a tangent line, both
                    points are the same instance, the B{C{point}} or wrapped
                    or I{normalized}.

           @raise IntersectionError: The circle and line do not intersect.

           @raise TypeError: If B{C{point}} is not this C{LatLon} or B{C{circle}}
                             or B{C{bearing}} invalid.

           @raise ValueError: Invalid B{C{circle}}, B{C{bearing}}, B{C{radius}},
                              B{C{exact}} or B{C{height}}.
        '''
        p = self.others(point=point)
        try:
            return _intersecant2(self, circle, p, bearing, radius=radius, exact=exact,
                                                           height=height, wrap=wrap)
        except (TypeError, ValueError) as x:
            raise _xError(x, center=self, circle=circle, point=point, bearing=bearing,
                                                         exact=exact, wrap=wrap)

    def maxLat(self, bearing):
        '''Return the maximum latitude reached when travelling
           on a great circle on given bearing from this point
           based on Clairaut's formula.

           The maximum latitude is independent of longitude
           and the same for all points on a given latitude.

           Negate the result for the minimum latitude (on the
           Southern hemisphere).

           @arg bearing: Initial bearing (compass C{degrees360}).

           @return: Maximum latitude (C{degrees90}).

           @raise ValueError: Invalid B{C{bearing}}.
        '''
        m = acos1(fabs(sin(Bearing_(bearing)) * cos(self.phi)))
        return degrees90(m)

    def minLat(self, bearing):
        '''Return the minimum latitude reached when travelling
           on a great circle on given bearing from this point.

           @arg bearing: Initial bearing (compass C{degrees360}).

           @return: Minimum latitude (C{degrees90}).

           @see: Method L{maxLat} for more details.

           @raise ValueError: Invalid B{C{bearing}}.
        '''
        return -self.maxLat(bearing)

    def parse(self, strllh, height=0, sep=_COMMA_, name=NN):
        '''Parse a string representing a similar, spherical C{LatLon}
           point, consisting of C{"lat, lon[, height]"}.

           @arg strllh: Lat, lon and optional height (C{str}),
                        see function L{pygeodesy.parse3llh}.
           @kwarg height: Optional, default height (C{meter}).
           @kwarg sep: Optional separator (C{str}).
           @kwarg name: Optional instance name (C{str}),
                        overriding this name.

           @return: The similar point (spherical C{LatLon}).

           @raise ParseError: Invalid B{C{strllh}}.
        '''
        t = _MODS.dms.parse3llh(strllh, height=height, sep=sep)
        r =  self.classof(*t)
        if name:
            r.rename(name)
        return r

    @property_RO
    def _radius(self):
        '''(INTERNAL) Get this sphere's radius.
        '''
        return self.datum.ellipsoid.equatoradius

    def _rhumbs3(self, other, wrap, r=False):  # != .latlonBase._rhumbx3
        '''(INTERNAL) Rhumb_ helper function.

           @arg other: The other point (spherical C{LatLon}).
        '''
        p = self.others(other, up=2)
        if wrap:
            p = _unrollon(self, p, wrap=wrap)
        a2, b2 = p.philam
        a1, b1 = self.philam
        # if |db| > 180 take shorter rhumb
        # line across the anti-meridian
        db = wrapPI(b2 - b1)
        dp = log(tanPI_2_2(a2) / tanPI_2_2(a1))
        da = a2 - a1
        if r:
            # on Mercator projection, longitude distances shrink
            # by latitude; the 'stretch factor' q becomes ill-
            # conditioned along E-W line (0/0); use an empirical
            # tolerance to avoid it
            q  = (da / dp) if fabs(dp) > EPS else cos(self.phi)
            da = hypot(da, q * db)  # angular distance radians
        return da, db, dp

    def rhumbAzimuthTo(self, other, radius=R_M, exact=False, wrap=False):
        '''Return the azimuth (bearing) of a rhumb line (loxodrome)
           between this and an other (spherical) point.

           @arg other: The other point (spherical C{LatLon}).
           @kwarg radius: Earth radius (C{meter}) or earth model (L{Datum},
                          L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg exact: If C{True}, use I{Krüger} L{rhumbx} (C{bool}),
                         default C{False} for backward compatibility.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: Rhumb line azimuth (compass C{degrees180}).

           @raise TypeError: The B{C{other}} point is incompatible or
                             B{C{radius}} is invalid.

           @example:

            >>> p = LatLon(51.127, 1.338)
            >>> q = LatLon(50.964, 1.853)
            >>> b = p.rhumbBearingTo(q)  # 116.7
        '''
        if exact:  # use series, always
            z = LatLonBase.rhumbAzimuthTo(self, other, exact=False,  # Krüger
                                                radius=radius, wrap=wrap)
        else:
            _, db, dp = self._rhumbs3(other, wrap)
            z = atan2d(db, dp)  # see .rhumbx.Rhumb.Inverse
        return z

    @deprecated_method
    def rhumbBearingTo(self, other):  # unwrapped
        '''DEPRECATED, use method C{.rhumbAzimuthTo}.'''
        return wrap360(self.rhumbAzimuthTo(other))  # [0..360)

    def rhumbDestination(self, distance, bearing, radius=R_M, height=None, exact=False):
        '''Return the destination point having travelled the given distance
           from this point along a rhumb line (loxodrome) at the given bearing.

           @arg distance: Distance travelled (C{meter}, same units as B{C{radius}}),
                          may be negative if C{B{exact}=True}.
           @arg bearing: Bearing (azimuth) at this point (compass C{degrees360}).
           @kwarg radius: Earth radius (C{meter}) or earth model (L{Datum},
                          L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}) if
                          C{B{exact}=True}.
           @kwarg height: Optional height, overriding the default height
                          (C{meter}, same unit as B{C{radius}}).
           @kwarg exact: If C{True}, use I{Krüger} L{rhumbx} (C{bool}),
                         default C{False} for backward compatibility.

           @return: The destination point (spherical C{LatLon}).

           @raise ValueError: Invalid B{C{distance}}, B{C{bearing}},
                              B{C{radius}} or B{C{height}}.

           @example:

            >>> p = LatLon(51.127, 1.338)
            >>> q = p.rhumbDestination(40300, 116.7)  # 50.9642°N, 001.8530°E
        '''
        if exact:  # use series, always
            r = LatLonBase.rhumbDestination(self, distance, bearing, exact=False,  # Krüger
                                                  radius=radius, height=height)
        else:  # radius=None from .rhumbMidpointTo
            if radius in (None, self._radius):
                d, r = self.datum, radius
            else:
                d = _spherical_datum(radius, raiser=_radius_)  # spherical only
                r =  d.ellipsoid.equatoradius
            r = _angular(distance, r, low=-EPS)  # distance=0 from .rhumbMidpointTo

            a1, b1 = self.philam
            sb, cb = sincos2(Bearing_(bearing))

            da = r  * cb
            a2 = a1 + da
            # normalize latitude if past pole
            if a2 > PI_2:
                a2 =  PI - a2
            elif a2 < -PI_2:
                a2 = -PI - a2

            dp = log(tanPI_2_2(a2) / tanPI_2_2(a1))
            # q becomes ill-conditioned on E-W course 0/0
            q  = (da / dp) if fabs(dp) > EPS else cos(a1)
            b2 = (b1 + r * sb / q) if fabs(q) > EPS else b1

            h = self._heigHt(height)
            r = self.classof(degrees90(a2), degrees180(b2), datum=d, height=h)
        return r

    def rhumbDistanceTo(self, other, radius=R_M, exact=False, wrap=False):
        '''Return the distance from this to an other point along
           a rhumb line (loxodrome).

           @arg other: The other point (spherical C{LatLon}).
           @kwarg radius: Earth radius (C{meter}) or earth model (L{Datum},
                          L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}) if
                          C{B{exact}=True}.
           @kwarg exact: If C{True}, use I{Krüger} L{rhumbx} (C{bool}),
                         default C{False} for backward compatibility.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: Distance (C{meter}, the same units as B{C{radius}}
                    or C{radians} if B{C{radius}} is C{None}).

           @raise TypeError: The B{C{other}} point is incompatible.

           @raise ValueError: Invalid B{C{radius}}.

           @example:

            >>> p = LatLon(51.127, 1.338)
            >>> q = LatLon(50.964, 1.853)
            >>> d = p.rhumbDistanceTo(q)  # 403100
        '''
        if exact:  # use series, always
            r = LatLonBase.rhumbDistanceTo(self, other, exact=False,  # Krüger
                                                 radius=radius, wrap=wrap)
            if radius is None:  # angular distance in radians
                r = r / self._radius  # /= chokes PyChecker
        else:
            # see <https://www.EdWilliams.org/avform.htm#Rhumb>
            r, _, _ = self._rhumbs3(other, wrap, r=True)
            if radius is not None:
                r *= Radius(radius)
        return r

    def rhumbMidpointTo(self, other, height=None, radius=R_M, exact=False,
                                                fraction=_0_5, wrap=False):
        '''Return the (loxodromic) midpoint on the rhumb line between
           this and an other point.

           @arg other: The other point (spherical LatLon).
           @kwarg height: Optional height, overriding the mean height
                          (C{meter}).
           @kwarg radius: Earth radius (C{meter}) or earth model (L{Datum},
                          L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg exact: If C{True}, use I{Krüger} L{rhumbx} (C{bool}),
                         default C{False} for backward compatibility.
           @kwarg fraction: Midpoint location from this point (C{scalar}),
                            may be negative if C{B{exact}=True}.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: The (mid)point at the given B{C{fraction}} along
                    the rhumb line (spherical C{LatLon}).

           @raise TypeError: The B{C{other}} point is incompatible.

           @raise ValueError: Invalid B{C{height}} or B{C{fraction}}

           @example:

            >>> p = LatLon(51.127, 1.338)
            >>> q = LatLon(50.964, 1.853)
            >>> m = p.rhumb_midpointTo(q)
            >>> m.toStr()  # '51.0455°N, 001.5957°E'
        '''
        if exact:  # use series, always
            r = LatLonBase.rhumbMidpointTo(self, other, exact=False,  # Krüger
                                                 radius=radius, height=height,
                                                 fraction=fraction, wrap=wrap)
        elif fraction is not _0_5:
            f = Scalar_(fraction=fraction)  # low=_0_0
            r, db, dp = self._rhumbs3(other, wrap, r=True)  # radians
            z = atan2b(db, dp)
            h = self._havg(other, f=f, h=height)
            r = self.rhumbDestination(r * f, z, radius=None, height=h)

        else:  # for backward compatibility, unwrapped
            # see <https://MathForum.org/library/drmath/view/51822.html>
            a1, b1 = self.philam
            a2, b2 = self.others(other).philam

            if fabs(b2 - b1) > PI:
                b1 += PI2  # crossing anti-meridian

            a3 = favg(a1, a2)
            b3 = favg(b1, b2)

            f1 = tanPI_2_2(a1)
            if isnon0(f1):
                f2 = tanPI_2_2(a2)
                f  = f2 / f1
                if isnon0(f):
                    f = log(f)
                    if isnon0(f):
                        f3 = tanPI_2_2(a3)
                        b3 = fdot(map1(log, f1, f2, f3),
                                           -b2, b1, b2 - b1) / f

            d = self.datum if radius in (None, self._radius) else \
               _spherical_datum(radius, name=self.name, raiser=_radius_)
            h = self._havg(other, h=height)
            r = self.classof(degrees90(a3), degrees180(b3), datum=d, height=h)
        return r

    def toNvector(self, Nvector=NvectorBase, **Nvector_kwds):  # PYCHOK signature
        '''Convert this point to C{Nvector} components, I{including
           height}.

           @kwarg Nvector_kwds: Optional, additional B{C{Nvector}}
                                keyword arguments, ignored if
                                C{B{Nvector} is None}.

           @return: An B{C{Nvector}} or a L{Vector4Tuple}C{(x, y, z, h)}
                    if B{C{Nvector}} is C{None}.

           @raise TypeError: Invalid B{C{Nvector}} or B{C{Nvector_kwds}}.
        '''
        return LatLonBase.toNvector(self, Nvector=Nvector, **Nvector_kwds)


def _intersecant2(c, r, p, b, radius=R_M, exact=False,
                              height=None, wrap=False):
    # (INTERNAL) Intersect a circle and bearing, see L{intersecant2}
    # above, separated to allow callers to embellish any exceptions

    if wrap:
        p = _unrollon(c, p, wrap=wrap)
    nonexact = exact is None

    if not isinstanceof(r, c.__class__, p.__class__):
        r = Radius_(circle=r)
    elif nonexact:
        r = c.distanceTo(r, radius=radius, wrap=wrap)
    elif isbool(exact):
        r = c.rhumbDistanceTo(r, radius=radius, exact=exact, wrap=wrap)
    else:
        raise _ValueError(exact=exact)

    if not isinstanceof(b, c.__class__, p.__class__):
        b = Bearing(b)
    elif nonexact:
        b = p.initialBearingTo(b, wrap=wrap)
    else:
        b = p.rhumbAzimuthTo(b, radius=radius, exact=exact, wrap=wrap)

    d = p.distanceTo(c, radius=radius) if nonexact else \
        p.rhumbDistanceTo(c, radius=radius, exact=exact)
    if d > EPS:
        a = p.initialBearingTo(c) if nonexact else wrap360(
            p.rhumbAzimuthTo(c, radius=radius, exact=exact))
        s, c = sincos2d(b - a)
        s = sqrt_a(r, fabs(s * d))
        if s > r:
            raise IntersectionError(_too_(Fmt.distant(s)))
        elif (r - s) < EPS:
            return p, p  # tangent
        c *= d
    else:  # p and c coincide
        s, c = r, 0
    t = ()
    for d, b in ((s + c, b), (s - c, b + _180_0)):  # bearing direction first
        t += (p.destination(d, b, radius=radius, height=height) if nonexact else
              p.rhumbDestination(d, b, radius=radius, height=height, exact=exact)),
    return t


def _r2m(r, radius):
    '''(INTERNAL) Angular distance in C{radians} to C{meter}.
    '''
    if radius is not None:  # not in (None, _0_0)
        r *= R_M if radius is R_M else Radius(radius)
    return r


__all__ += _ALL_DOCS(CartesianSphericalBase, LatLonSphericalBase)

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
