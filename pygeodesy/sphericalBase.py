
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

from pygeodesy.basics import _copysign, isbool, isinstanceof, map1
from pygeodesy.cartesianBase import CartesianBase,  Bearing2Tuple
from pygeodesy.constants import EPS, EPS0, PI, PI2, PI_2, R_M, \
                               _0_0, _0_5, _1_0, _180_0, _360_0, \
                               _over, isnear0, isnon0
from pygeodesy.datums import Datums, _earth_ellipsoid, _spherical_datum
from pygeodesy.errors import IntersectionError, _ValueError, \
                            _xattr, _xError
from pygeodesy.fmath import favg, fdot, hypot, sqrt_a
from pygeodesy.interns import NN, _COMMA_, _concentric_, _datum_, \
                             _distant_, _exceed_PI_radians_, _name_, \
                             _near_, _radius_, _too_
from pygeodesy.latlonBase import LatLonBase,  _trilaterate5  # PYCHOK passed
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
# from pygeodesy.namedTuples import Bearing2Tuple  # from .cartesianBase
from pygeodesy.nvectorBase import NvectorBase,  Fmt, _xattrs
from pygeodesy.props import deprecated_method, property_doc_, \
                            property_RO, _update_all
# from pygeodesy.streprs import Fmt, _xattrs  # from .nvectorBase
from pygeodesy.units import _isRadius, Bearing, Bearing_, Radians_, \
                             Radius, Radius_, Scalar_, _100km
from pygeodesy.utily import acos1, asin1, atan2b, atan2d, degrees90, \
                            degrees180, sincos2, sincos2d, _unrollon, \
                            tanPI_2_2, wrapPI

from math import cos, fabs, log, sin, sqrt

__all__ = _ALL_LAZY.sphericalBase
__version__ = '23.12.18'


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

    @property_RO
    def sphericalCartesian(self):
        '''Get this C{Cartesian}'s spherical class.
        '''
        return type(self)


class LatLonSphericalBase(LatLonBase):
    '''(INTERNAL) Base class for spherical C{LatLon}s.
    '''
    _datum       =  Datums.Sphere  # spherical L{Datum}
    _napieradius = _100km

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
           @kwarg name: Optional name (C{str}).

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
        '''
        p = self.others(other)
        if wrap:
            p = _unrollon(self, p, wrap=wrap)
        # final bearing is the reverse of the other, initial one
        b = p.initialBearingTo(self, wrap=False, raiser=raiser) + _180_0
        return b if b < 360 else (b - _360_0)

    def intersecant2(self, circle, point, other, radius=R_M, exact=False,  # PYCHOK signature
                                                 height=None, wrap=False):
        '''Compute the intersections of a circle and a (great circle) line
           given as two points or as a point and bearing.

           @arg circle: Radius of the circle centered at this location (C{meter},
                        same units as B{C{radius}}) or a point on the circle
                        (this C{LatLon}).
           @arg point: A point on the (great circle) line (this C{LatLon}).
           @arg other: An other point I{on} (this {LatLon}) or the bearing at
                       B{C{point}} I{of} the (great circle) line (compass
                       C{degrees}).
           @kwarg radius: Mean earth radius (C{meter}, conventionally).
           @kwarg exact: If C{True} use the I{exact} rhumb methods for azimuth,
                         destination and distance, if C{False} use the basic
                         rhumb methods (C{bool}) or if C{None} use the I{great
                         circle} methods.
           @kwarg height: Optional height for the intersection points (C{meter},
                          conventionally) or C{None} for interpolated heights.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the points
                        B{C{circle}}, B{C{point}} and/or B{C{other}} (C{bool}).

           @return: 2-Tuple of the intersection points (representing a chord), each
                    an instance of the B{C{point}} class.  Both points are the same
                    instance if the (great circle) line is tangent to the circle.

           @raise IntersectionError: The circle and line do not intersect.

           @raise TypeError: If B{C{point}} is not this C{LatLon} or B{C{circle}}
                             or B{C{other}} invalid.

           @raise UnitError: Invalid B{C{circle}}, B{C{other}}, B{C{radius}},
                             B{C{exact}}, B{C{height}} or B{C{napieradius}}.
        '''
        p = self.others(point=point)
        try:
            return _intersecant2(self, circle, p, other, radius=radius, exact=exact,
                                                         height=height, wrap=wrap)
        except (TypeError, ValueError) as x:
            raise _xError(x, center=self, circle=circle, point=point, other=other,
                             radius=radius, exact=exact, height=height, wrap=wrap)

    def maxLat(self, bearing):
        '''Return the maximum latitude reached when travelling on a great circle
           on given bearing from this point based on Clairaut's formula.

           The maximum latitude is independent of longitude and the same for all
           points on a given latitude.

           Negate the result for the minimum latitude (on the Southern hemisphere).

           @arg bearing: Initial bearing (compass C{degrees360}).

           @return: Maximum latitude (C{degrees90}).

           @raise ValueError: Invalid B{C{bearing}}.
        '''
        r = acos1(fabs(sin(Bearing_(bearing)) * cos(self.phi)))
        return degrees90(r)

    def minLat(self, bearing):
        '''Return the minimum latitude reached when travelling on a great circle
           on given bearing from this point.

           @arg bearing: Initial bearing (compass C{degrees360}).

           @return: Minimum latitude (C{degrees90}).

           @see: Method L{maxLat} for more details.

           @raise ValueError: Invalid B{C{bearing}}.
        '''
        return -self.maxLat(bearing)

    def _mpr(self, radius=R_M, exact=None):  # meter per radian
        if exact and not _isRadius(radius):  # see .rhumb.ekx.Rhumb._mpr
            radius = _earth_ellipsoid(radius)._Lpr
        return radius

    @property_doc_(''' the I{Napier} radius to apply spherical trigonometry.''')
    def napieradius(self):
        '''Get the I{Napier} radius (C{meter}, conventionally).
        '''
        return self._napieradius

    @napieradius.setter  # PYCHOK setter!
    def napieradius(self, radius):
        '''Set this I{Napier} radius (C{meter}, conventionally) or C{0}.

           In methods L{intersecant2} and L{rhumbIntersecant2}, I{Napier}'s
           spherical trigonometry is applied if the circle radius exceeds
           the I{Napier} radius, otherwise planar trigonometry is used.

           @raise UnitError: Invalid B{C{radius}}.
        '''
        self._napieradius = Radius(napieradius=radius or 0)

#   def nearestTo(self, point, other, **radius_exact_height_wrap):  # PYCHOK signature
#       p = self.others(point=point)
#       try:
#           p, q = _intersecant2(self, p, p, other, **radius_exact_height_wrap)
#       except (TypeError, ValueError) as x:
#           raise _xError(x, this=self, point=point, other=other, **radius_exact_height_wrap)
#       return p.midpointTo(q)

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
        db =  wrapPI(b2 - b1)
        dp = _logPI_2_2(a2, a1)
        da =  a2 - a1
        if r:
            # on Mercator projection, longitude distances shrink
            # by latitude; the 'stretch factor' q becomes ill-
            # conditioned along E-W line (0/0); use an empirical
            # tolerance to avoid it
            q  = (da / dp) if fabs(dp) > EPS else cos(a1)
            da = hypot(da, q * db)  # angular distance radians
        return da, db, dp

    def rhumbAzimuthTo(self, other, radius=R_M, exact=False, wrap=False, b360=False):
        '''Return the azimuth (bearing) of a rhumb line (loxodrome) between
           this and an other (spherical) point.

           @arg other: The other point (spherical C{LatLon}).
           @kwarg radius: Earth radius (C{meter}) or earth model (L{Datum},
                          L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg exact: If C{True}, use I{Elliptic, Krüger} L{Rhumb} (C{bool}),
                         default C{False} for backward compatibility.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).
           @kwarg b360: If C{True}, return the azimuth in the bearing range.

           @return: Rhumb azimuth (compass C{degrees180} or C{degrees360}).

           @raise TypeError: The B{C{other}} point is incompatible or
                             B{C{radius}} is invalid.
        '''
        if exact:  # use series, always
            z = LatLonBase.rhumbAzimuthTo(self, other, exact=False,  # Krüger
                                                radius=radius, wrap=wrap, b360=b360)
        else:
            _, db, dp = self._rhumbs3(other, wrap)
            z = (atan2b if b360 else atan2d)(db, dp)  # see .rhumbBase.RhumbBase.Inverse
        return z

    @deprecated_method
    def rhumbBearingTo(self, other):  # unwrapped
        '''DEPRECATED, use method C{.rhumbAzimuthTo}.'''
        return self.rhumbAzimuthTo(other, b360=True)  # [0..360)

    def rhumbDestination(self, distance, azimuth, radius=R_M, height=None, exact=False):
        '''Return the destination point having travelled the given distance from
           this point along a rhumb line (loxodrome) of the given azimuth.

           @arg distance: Distance travelled (C{meter}, same units as B{C{radius}}),
                          may be negative if C{B{exact}=True}.
           @arg azimuth: Azimuth (bearing) of the rhumb line (compass C{degrees}).
           @kwarg radius: Earth radius (C{meter}) or earth model (L{Datum},
                          L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}) if
                          C{B{exact}=True}.
           @kwarg height: Optional height, overriding the default height (C{meter}.
           @kwarg exact: If C{True}, use I{Elliptic, Krüger} L{Rhumb} (C{bool}),
                         default C{False} for backward compatibility.

           @return: The destination point (spherical C{LatLon}).

           @raise ValueError: Invalid B{C{distance}}, B{C{azimuth}}, B{C{radius}}
                              or B{C{height}}.
        '''
        if exact:  # use series, always
            r = LatLonBase.rhumbDestination(self, distance, azimuth, exact=False,  # Krüger
                                                  radius=radius, height=height)
        else:  # radius=None from .rhumbMidpointTo
            if radius in (None, self._radius):
                d, r = self.datum, radius
            else:
                d = _spherical_datum(radius, raiser=_radius_)  # spherical only
                r =  d.ellipsoid.equatoradius
            r = _m2radians(distance, r, low=-EPS)  # distance=0 from .rhumbMidpointTo

            a1, b1 = self.philam
            sb, cb = sincos2(Bearing_(azimuth))  # radians

            da = r  * cb
            a2 = a1 + da
            # normalize latitude if past pole
            if fabs(a2) > PI_2:
                a2 = _copysign(PI, a2) - a2

            dp = _logPI_2_2(a2, a1)
            # q becomes ill-conditioned on E-W course 0/0
            q  = cos(a1) if isnear0(dp) else (da / dp)
            b2 = b1 if isnear0(q) else (b1 + r * sb / q)

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
           @kwarg exact: If C{True}, use I{Elliptic, Krüger} L{Rhumb} (C{bool}),
                         default C{False} for backward compatibility.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: Distance (C{meter}, the same units as B{C{radius}}
                    or C{radians} if B{C{radius}} is C{None}).

           @raise TypeError: The B{C{other}} point is incompatible.

           @raise ValueError: Invalid B{C{radius}}.
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

    def rhumbIntersecant2(self, circle, point, other, radius=R_M, exact=True,  # PYCHOK signature
                                                      height=None, wrap=False):
        '''Compute the intersections of a circle and a rhumb line given as two
           points and as a point and azimuth.

           @arg circle: Radius of the circle centered at this location (C{meter},
                        same units as B{C{radius}}) or a point on the circle
                        (this C{LatLon}).
           @arg point: The rhumb line's start point (this C{LatLon}).
           @arg other: An other point (this I{on} C{LatLon}) or the azimuth I{of}
                       (compass C{degrees}) the rhumb line.
           @kwarg radius: Mean earth radius (C{meter}, conventionally).
           @kwarg exact: If C{True} use the I{exact} rhumb methods for azimuth,
                         destination and distance, if C{False} use the basic
                         rhumb methods (C{bool}) or if C{None} use the I{great
                         circle} methods.
           @kwarg height: Optional height for the intersection points (C{meter},
                          conventionally) or C{None}.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the points
                        B{C{circle}}, B{C{point}} and/or B{C{other}} (C{bool}).

           @return: 2-Tuple of the intersection points (representing a chord),
                    each an instance of this class.  For a tangent line, both
                    points are the same instance, wrapped or I{normalized}.

           @raise IntersectionError: The circle and line do not intersect.

           @raise TypeError: If B{C{point}} is not this C{LatLon} or B{C{circle}}
                             or B{C{other}} invalid.

           @raise UnitError: Invalid B{C{circle}}, B{C{other}}, B{C{radius}},
                             B{C{exact}} or B{C{height}}.
        '''
        m = LatLonBase.rhumbIntersecant2 if exact else \
            LatLonSphericalBase.intersecant2
        return m(self, circle, point, other, radius=radius, exact=exact,
                                             height=height, wrap=wrap)

    def rhumbMidpointTo(self, other, height=None, radius=R_M, exact=False,
                                                fraction=_0_5, wrap=False):
        '''Return the (loxodromic) midpoint on the rhumb line between
           this and an other point.

           @arg other: The other point (spherical LatLon).
           @kwarg height: Optional height, overriding the mean height (C{meter}).
           @kwarg radius: Earth radius (C{meter}) or earth model (L{Datum},
                          L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg exact: If C{True}, use I{Elliptic, Krüger} L{Rhumb} (C{bool}),
                         default C{False} for backward compatibility.
           @kwarg fraction: Midpoint location from this point (C{scalar}), may
                            be negative if C{B{exact}=True}.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the B{C{other}}
                        point (C{bool}).

           @return: The (mid)point at the given B{C{fraction}} along the rhumb
                    line (spherical C{LatLon}).

           @raise TypeError: The B{C{other}} point is incompatible.

           @raise ValueError: Invalid B{C{height}} or B{C{fraction}}
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

    @property_RO
    def sphericalLatLon(self):
        '''Get this C{LatLon}'s spherical class.
        '''
        return type(self)

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


def _intersecant2(c, r, p, b, radius=R_M, exact=False, height=None, wrap=False):
    # (INTERNAL) Intersect a circle and line, see L{intersecant2}
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
        b = p.rhumbAzimuthTo(b, radius=radius, exact=exact, wrap=wrap,
                             b360=True)

    d = p.distanceTo(c, radius=radius) if nonexact else \
        p.rhumbDistanceTo(c, radius=radius, exact=exact)
    if d > EPS0:
        n = _xattr(c, napieradius=0)
        a =  p.initialBearingTo(c) if nonexact else \
             p.rhumbAzimuthTo(c, radius=radius, exact=exact, b360=True)
        s, c = sincos2d(b - a)  # Napier's sin(A), cos(A)
        if r > n:
            # Napier's right spherical triangle rules (R2) and (R1)
            # <https://WikiPedia.org/wiki/Spherical_trigonometry>
            m = p._mpr(radius=radius, exact=exact)  # meter per radian
            if fabs(c) > EPS0:
                d =  d / m  # /= chokes PyChecker
                a =  asin1(sin(d) * fabs(s))  # Napier's a
                c = _copysign(cos(a), c)
                d =  acos1(cos(d) / c) * m
                a *= m  # meter
            else:  # point and chord center coincident
                a, d = d, 0
                c = cos(a / m)
            h = (acos1(cos(r / m) / c) * m) if a < r else 0
        else:  # distance from the chord center to ...
            a  = fabs(d * s)  # ... the cicle center ...
            d *= c  # ... and to the point
            h  = sqrt_a(r, a) if a < r else 0  # half chord length
        if a > r:
            raise IntersectionError(_too_(Fmt.distant(a)))
    else:
        d, h = 0, r  # point and circle center coincident

    _intersecant1, kwds = (p.destination, {}) if nonexact else \
                          (p.rhumbDestination, dict(exact=exact))
    kwds.update(radius=radius, height=height)
    t = (_intersecant1(d + h, b, **kwds),)
    if h:
        t += (_intersecant1(d - h, b, **kwds),)
    else:  # same instance twice
        t *= 2
    return t


def _logPI_2_2(a2, a1):
    '''(INTERNAL) C{log} of C{tanPI_2_2}'s quotient.
    '''
    return log(_over(tanPI_2_2(a2), tanPI_2_2(a1)))


def _m2radians(distance, radius, low=EPS):  # PYCHOK in .spherical*
    '''(INTERNAL) Distance in C{meter} to angular distance in C{radians}.

       @raise UnitError: Invalid B{C{distance}} or B{C{radius}}.
    '''
    r = float(distance)
    if radius:
        r = r / Radius_(radius=radius)  # /= chokes PyChecker
    if low is not None:
        # small near0 values from .rhumbDestination not exact OK
        r = _0_0 if low < 0 and r < 0 else Radians_(r, low=low)
        # _0_0 if low < 0 and low < r < 0 else Radians_(r, low=low)
    return r


def _radians2m(rad, radius):
    '''(INTERNAL) Angular distance in C{radians} to distance in C{meter}.
    '''
    if radius is not None:  # not in (None, _0_0)
        rad *= R_M if radius is R_M else Radius(radius)
    return rad


def _rads3(rad1, rad2, radius):  # in .sphericalTrigonometry
    '''(INTERNAL) Convert radii to radians.
    '''
    r1 = Radius_(rad1=rad1)
    r2 = Radius_(rad2=rad2)
    if radius is not None:  # convert radii to radians
        r1 = _m2radians(r1, radius)
        r2 = _m2radians(r2, radius)

    x = r1 < r2
    if x:
        r1, r2 = r2, r1
    if r1 > PI:
        raise IntersectionError(rad1=rad1, rad2=rad2,
                                txt=_exceed_PI_radians_)
    return r1, r2, x


__all__ += _ALL_DOCS(CartesianSphericalBase, LatLonSphericalBase)

# **) MIT License
#
# Copyright (C) 2016-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
