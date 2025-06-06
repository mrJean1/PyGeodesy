
# -*- coding: utf-8 -*-

u'''Spherical, C{trigonometry}-based geodesy.

Trigonometric classes geodetic (lat-/longitude) L{LatLon} and
geocentric (ECEF) L{Cartesian} and functions L{areaOf}, L{intersection},
L{intersections2}, L{isPoleEnclosedBy}, L{meanOf}, L{nearestOn3} and
L{perimeterOf}, I{all spherical}.

Pure Python implementation of geodetic (lat-/longitude) methods using
spherical trigonometry, transcoded from JavaScript originals by
I{(C) Chris Veness 2011-2024} published under the same MIT Licence**, see
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

from pygeodesy.basics import copysign0, _isin, map1, signOf,  typename
from pygeodesy.constants import EPS, EPS1, EPS4, PI, PI2, PI_2, PI_4, R_M, \
                                isnear0, isnear1, isnon0, _0_0, _0_5, \
                                _1_0, _2_0, _90_0
from pygeodesy.datums import _ellipsoidal_datum, _mean_radius
from pygeodesy.errors import _AssertionError, CrossError, crosserrors, \
                             _TypeError, _ValueError, IntersectionError, \
                             _xError, _xkwds, _xkwds_get, _xkwds_pop2
from pygeodesy.fmath import favg, fdot, fdot_, fmean, hypot
from pygeodesy.fsums import Fsum, fsum, fsumf_
from pygeodesy.formy import antipode_, bearing_, _bearingTo2, excessAbc_, \
                            excessGirard_, excessLHuilier_, opposing_, _radical2, \
                            vincentys_
# from pygeodesy.internals import typename  # from .basics
from pygeodesy.interns import _1_, _2_, _coincident_, _composite_, _colinear_, \
                              _concentric_, _convex_, _end_, _infinite_, \
                              _invalid_, _line_, _near_, _null_, _parallel_, \
                              _point_, _SPACE_, _too_
from pygeodesy.latlonBase import _trilaterate5
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS, _ALL_OTHER
# from pygeodesy.nvectorBase import NvectorBase, sumOf  # _MODS
from pygeodesy.namedTuples import LatLon2Tuple, LatLon3Tuple, NearestOn3Tuple, \
                                  Triangle7Tuple, Triangle8Tuple
from pygeodesy.points import ispolar, nearestOn5 as _nearestOn5, \
                             Fmt as _Fmt  # XXX shadowed
from pygeodesy.props import deprecated_function, deprecated_method
from pygeodesy.sphericalBase import _m2radians, CartesianSphericalBase, \
                                    _intersecant2, LatLonSphericalBase, \
                                    _rads3, _radians2m
# from pygeodesy.streprs import Fmt as _Fmt  # from .points XXX shadowed
from pygeodesy.units import Bearing_, Height, _isDegrees, _isRadius, Lamd, \
                            Phid, Radius_, Scalar
from pygeodesy.utily import acos1, asin1, atan1d, atan2, atan2d, degrees90, \
                            degrees180, degrees2m, m2radians, radiansPI2, \
                            sincos2_, tan_2, unrollPI, _unrollon, _unrollon3, \
                            wrap180, wrapPI, _Wrap
from pygeodesy.vector3d import sumOf, Vector3d

from math import asin, cos, degrees, fabs, radians, sin

__all__ = _ALL_LAZY.sphericalTrigonometry
__version__ = '25.05.28'

_PI_EPS4 = PI - EPS4
if _PI_EPS4 >= PI:
    raise _AssertionError(EPS4=EPS4, PI=PI, PI_EPS4=_PI_EPS4)


class Cartesian(CartesianSphericalBase):
    '''Extended to convert geocentric, L{Cartesian} points to
       spherical, geodetic L{LatLon}.
    '''

    def toLatLon(self, **LatLon_and_kwds):  # PYCHOK LatLon=LatLon
        '''Convert this cartesian point to a C{spherical} geodetic point.

           @kwarg LatLon_and_kwds: Optional L{LatLon} and L{LatLon} keyword
                                   arguments.  Use C{B{LatLon}=...} to override
                                   this L{LatLon} class or specify C{B{LatLon}=None}.

           @return: The geodetic point (L{LatLon}) or if C{B{LatLon} is None},
                    an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon_and_kwds}} argument.
        '''
        kwds = _xkwds(LatLon_and_kwds, LatLon=LatLon, datum=self.datum)
        return CartesianSphericalBase.toLatLon(self, **kwds)


class LatLon(LatLonSphericalBase):
    '''New point on a spherical earth model, based on trigonometry formulae.
    '''

    def _ab1_ab2_db5(self, other, wrap):
        '''(INTERNAL) Helper for several methods.
        '''
        a1, b1 = self.philam
        a2, b2 = self.others(other, up=2).philam
        if wrap:
            a2, b2 = _Wrap.philam(a2, b2)
            db, b2 =  unrollPI(b1, b2, wrap=wrap)
        else:  # unrollPI shortcut
            db = b2 - b1
        return a1, b1, a2, b2, db

    def alongTrackDistanceTo(self, start, end, radius=R_M, wrap=False):
        '''Compute the (signed) distance from the start to the closest
           point on the great circle line defined by a start and an
           end point.

           That is, if a perpendicular is drawn from this point to the
           great circle line, the along-track distance is the distance
           from the start point to the point where the perpendicular
           crosses the line.

           @arg start: Start point of the great circle line (L{LatLon}).
           @arg end: End point of the great circle line (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None}.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{start}} and B{C{end}} point (C{bool}).

           @return: Distance along the great circle line (C{radians}
                    if C{B{radius} is None} or C{meter}, same units
                    as B{C{radius}}), positive if I{after} the
                    B{C{start}} toward the B{C{end}} point of the
                    line, I{negative} if before or C{0} if at the
                    B{C{start}} point.

           @raise TypeError: Invalid B{C{start}} or B{C{end}} point.

           @raise ValueError: Invalid B{C{radius}}.
        '''
        r, x, b = self._a_x_b3(start, end, radius, wrap)
        cx = cos(x)
        return _0_0 if isnear0(cx) else \
               _radians2m(copysign0(acos1(cos(r) / cx), cos(b)), radius)

    def _a_x_b3(self, start, end, radius, wrap):
        '''(INTERNAL) Helper for .along-/crossTrackDistanceTo.
        '''
        s = self.others(start=start)
        e = self.others(end=end)
        s, e, w = _unrollon3(self, s, e, wrap)

        r = Radius_(radius)
        r = s.distanceTo(self, r, wrap=w) / r

        b = radians(s.initialBearingTo(self, wrap=w)
                  - s.initialBearingTo(e,    wrap=w))
        x = asin(sin(r) * sin(b))
        return r, x, -b

    @deprecated_method
    def bearingTo(self, other, wrap=False, raiser=False):  # PYCHOK no cover
        '''DEPRECATED, use method L{initialBearingTo}.
        '''
        return self.initialBearingTo(other, wrap=wrap, raiser=raiser)

    def crossingParallels(self, other, lat, wrap=False):
        '''Return the pair of meridians at which a great circle defined
           by this and an other point crosses the given latitude.

           @arg other: The other point defining great circle (L{LatLon}).
           @arg lat: Latitude at the crossing (C{degrees}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: 2-Tuple C{(lon1, lon2)}, both in C{degrees180} or
                    C{None} if the great circle doesn't reach B{C{lat}}.
        '''
        a1, b1, a2, b2, db = self._ab1_ab2_db5(other, wrap)
        sa,  ca,  sa1, ca1, \
        sa2, ca2, sdb, cdb = sincos2_(radians(lat), a1, a2, db)
        sa1 *= ca2 * ca

        x = sa1 * sdb
        y = sa1 * cdb - ca1 * sa2 * ca
        z = ca1 * sdb * ca2 * sa

        h = hypot(x, y)
        if h < EPS or fabs(z) > h:  # PYCHOK no cover
            return None  # great circle doesn't reach latitude

        m = atan2(-y, x) + b1  # longitude at max latitude
        d = acos1(z / h)  # delta longitude to intersections
        return degrees180(m - d), degrees180(m + d)

    def crossTrackDistanceTo(self, start, end, radius=R_M, wrap=False):
        '''Compute the (signed) distance from this point to a great
           circle from a start to an end point.

           @arg start: Start point of the great circle line (L{LatLon}).
           @arg end: End point of the great circle line (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None}.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{start}} and B{C{end}} point (C{bool}).

           @return: Distance to the great circle (C{radians} if
                    B{C{radius}} or C{meter}, same units as
                    B{C{radius}}), I{negative} if to the left or
                    I{positive} if to the right of the line.

           @raise TypeError: If B{C{start}} or B{C{end}} is not L{LatLon}.

           @raise ValueError: Invalid B{C{radius}}.
        '''
        _, x, _ = self._a_x_b3(start, end, radius, wrap)
        return _radians2m(x, radius)

    def destination(self, distance, bearing, radius=R_M, height=None):
        '''Locate the destination from this point after having
           travelled the given distance on a bearing from North.

           @arg distance: Distance travelled (C{meter}, same units as
                          B{C{radius}}).
           @arg bearing: Bearing from this point (compass C{degrees360}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg height: Optional height at destination (C{meter}, same
                          units a B{C{radius}}).

           @return: Destination point (L{LatLon}).

           @raise ValueError: Invalid B{C{distance}}, B{C{bearing}},
                              B{C{radius}} or B{C{height}}.
        '''
        a, b =  self.philam
        r, t = _m2radians(distance, radius, low=None), Bearing_(bearing)

        a, b = _destination2(a, b, r, t)
        h = self._heigHt(height)
        return self.classof(degrees90(a), degrees180(b), height=h)

    def distanceTo(self, other, radius=R_M, wrap=False):
        '''Compute the (angular) distance from this to an other point.

           @arg other: The other point (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None}.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{other}} point (C{bool}).

           @return: Distance between this and the B{C{other}} point
                    (C{meter}, same units as B{C{radius}} or
                    C{radians} if C{B{radius} is None}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: Invalid B{C{radius}}.
        '''
        a1, _, a2, _, db = self._ab1_ab2_db5(other, wrap)
        return _radians2m(vincentys_(a2, a1, db), radius)

#   @Property_RO
#   def Ecef(self):
#       '''Get the ECEF I{class} (L{EcefVeness}), I{lazily}.
#       '''
#       return _MODS.ecef.EcefKarney

    def greatCircle(self, bearing, Vector=Vector3d, **Vector_kwds):
        '''Compute the vector normal to great circle obtained by heading
           from this point on the bearing from North.

           Direction of vector is such that initial bearing vector
           b = c × n, where n is an n-vector representing this point.

           @arg bearing: Bearing from this point (compass C{degrees360}).
           @kwarg Vector: Vector class to return the great circle,
                          overriding the default L{Vector3d}.
           @kwarg Vector_kwds: Optional, additional keyword argunents
                               for B{C{Vector}}.

           @return: Vector representing great circle (C{Vector}).

           @raise ValueError: Invalid B{C{bearing}}.
        '''
        a, b = self.philam
        sa, ca, sb, cb, st, ct = sincos2_(a, b, Bearing_(bearing))

        sa *= st
        return Vector(fdot_(sb, ct, -cb, sa),
                     -fdot_(cb, ct,  sb, sa),
                      ca * st, **Vector_kwds)  # XXX .unit()?

    def initialBearingTo(self, other, wrap=False, raiser=False):
        '''Compute the initial bearing (forward azimuth) from this
           to an other point.

           @arg other: The other point (spherical L{LatLon}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{other}} point (C{bool}).
           @kwarg raiser: Optionally, raise L{CrossError} (C{bool}),
                          use C{B{raiser}=True} for behavior like
                          C{sphericalNvector.LatLon.initialBearingTo}.

           @return: Initial bearing (compass C{degrees360}).

           @raise CrossError: If this and the B{C{other}} point coincide
                              and if B{C{raiser}} and L{crosserrors
                              <pygeodesy.crosserrors>} are both C{True}.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.
        '''
        a1, b1, a2, b2, db = self._ab1_ab2_db5(other, wrap)
        # XXX behavior like sphericalNvector.LatLon.initialBearingTo
        if raiser and crosserrors() and max(fabs(a2 - a1), fabs(db)) < EPS:
            raise CrossError(_point_, self, other=other, wrap=wrap, txt=_coincident_)

        return degrees(bearing_(a1, b1, a2, b2, final=False))

    def intermediateTo(self, other, fraction, height=None, wrap=False):
        '''Locate the point at given fraction between (or along) this
           and an other point.

           @arg other: The other point (L{LatLon}).
           @arg fraction: Fraction between both points (C{scalar},
                          0.0 at this and 1.0 at the other point).
           @kwarg height: Optional height, overriding the intermediate
                          height (C{meter}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: Intermediate point (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: Invalid B{C{fraction}} or B{C{height}}.

           @see: Methods C{midpointTo} and C{rhumbMidpointTo}.
        '''
        p = self
        f = Scalar(fraction=fraction)
        if not isnear0(f):
            p = p.others(other)
            if wrap:
                p = _Wrap.point(p)
            if not isnear1(f):  # and not near0
                a1, b1 = self.philam
                a2, b2 = p.philam
                db, b2 = unrollPI(b1, b2, wrap=wrap)
                r  = vincentys_(a2, a1, db)
                sr = sin(r)
                if isnon0(sr):
                    sa1, ca1, sa2, ca2, \
                    sb1, cb1, sb2, cb2 = sincos2_(a1, a2, b1, b2)

                    t = f * r
                    a = sin(r - t)  # / sr  superflous
                    b = sin(    t)  # / sr  superflous

                    x = fdot_(a, ca1 * cb1, b, ca2 * cb2)
                    y = fdot_(a, ca1 * sb1, b, ca2 * sb2)
                    z = fdot_(a, sa1,       b, sa2)

                    a = atan1d(z, hypot(x, y))
                    b = atan2d(y, x)

                else:  # PYCHOK no cover
                    a = degrees90( favg(a1, a2, f=f))  # coincident
                    b = degrees180(favg(b1, b2, f=f))

                h = self._havg(other, f=f, h=height)
                p = self.classof(a, b, height=h)
        return p

    def intersection(self, end1, other, end2, height=None, wrap=False):
        '''Compute the intersection point of two lines, each defined by
           two points or a start point and a bearing from North.

           @arg end1: End point of this line (L{LatLon}) or the initial
                      bearing at this point (compass C{degrees360}).
           @arg other: Start point of the other line (L{LatLon}).
           @arg end2: End point of the other line (L{LatLon}) or the
                      initial bearing at the B{C{other}} point (compass
                      C{degrees360}).
           @kwarg height: Optional height for intersection point,
                          overriding the mean height (C{meter}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        B{C{start2}} and both B{C{end*}} points (C{bool}).

           @return: The intersection point (L{LatLon}).  An alternate
                    intersection point might be the L{antipode} to
                    the returned result.

           @raise IntersectionError: Ambiguous or infinite intersection
                                     or colinear, parallel or otherwise
                                     non-intersecting lines.

           @raise TypeError: If B{C{other}} is not L{LatLon} or B{C{end1}}
                             or B{C{end2}} not C{scalar} nor L{LatLon}.

           @raise ValueError: Invalid B{C{height}} or C{null} line.
        '''
        try:
            s2 = self.others(other)
            return _intersect(self, end1, s2, end2, height=height, wrap=wrap,
                                                    LatLon=self.classof)
        except (TypeError, ValueError) as x:
            raise _xError(x, start1=self, end1=end1,
                             other=other, end2=end2, wrap=wrap)

    def intersections2(self, rad1, other, rad2, radius=R_M, eps=_0_0,
                                                height=None, wrap=True):
        '''Compute the intersection points of two circles, each defined
           by a center point and a radius.

           @arg rad1: Radius of the this circle (C{meter} or C{radians},
                      see B{C{radius}}).
           @arg other: Center point of the other circle (L{LatLon}).
           @arg rad2: Radius of the other circle (C{meter} or C{radians},
                      see B{C{radius}}).
           @kwarg radius: Mean earth radius (C{meter} or C{None} if B{C{rad1}},
                          B{C{rad2}} and B{C{eps}} are given in C{radians}).
           @kwarg eps: Required overlap (C{meter} or C{radians}, see
                       B{C{radius}}).
           @kwarg height: Optional height for the intersection points (C{meter},
                          conventionally) or C{None} for the I{"radical height"}
                          at the I{radical line} between both centers.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: 2-Tuple of the intersection points, each a L{LatLon}
                    instance.  For abutting circles, both intersection
                    points are the same instance, aka the I{radical center}.

           @raise IntersectionError: Concentric, antipodal, invalid or
                                     non-intersecting circles.

           @raise TypeError: If B{C{other}} is not L{LatLon}.

           @raise ValueError: Invalid B{C{rad1}}, B{C{rad2}}, B{C{radius}},
                              B{C{eps}} or B{C{height}}.
        '''
        try:
            c2 = self.others(other)
            return _intersects2(self, rad1, c2, rad2, radius=radius, eps=eps,
                                                      height=height, wrap=wrap,
                                                      LatLon=self.classof)
        except (TypeError, ValueError) as x:
            raise _xError(x, center=self,  rad1=rad1,
                              other=other, rad2=rad2, wrap=wrap)

    @deprecated_method
    def isEnclosedBy(self, points):  # PYCHOK no cover
        '''DEPRECATED, use method C{isenclosedBy}.'''
        return self.isenclosedBy(points)

    def isenclosedBy(self, points, wrap=False):
        '''Check whether a (convex) polygon or composite encloses this point.

           @arg points: The polygon points or composite (L{LatLon}[],
                        L{BooleanFHP} or L{BooleanGH}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{points}} (C{bool}).

           @return: C{True} if this point is inside the polygon or
                    composite, C{False} otherwise.

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} are not L{LatLon}.

           @raise ValueError: Invalid B{C{points}}, non-convex polygon.

           @see: Functions L{pygeodesy.isconvex}, L{pygeodesy.isenclosedBy}
                 and L{pygeodesy.ispolar} especially if the B{C{points}} may
                 enclose a pole or wrap around the earth I{longitudinally}.
        '''
        if _MODS.booleans.isBoolean(points):
            return points._encloses(self.lat, self.lon, wrap=wrap)

        Ps = self.PointsIter(points, loop=2, dedup=True, wrap=wrap)
        n0 = self._N_vector

        v2 = Ps[0]._N_vector
        p1 = Ps[1]
        v1 = p1._N_vector
        # check whether this point on same side of all
        # polygon edges (to the left or right depending
        # on the anti-/clockwise polygon direction)
        gc1 = v2.cross(v1)
        t0  = gc1.angleTo(n0) > PI_2
        s0  = None
        # get great-circle vector for each edge
        for i, p2 in Ps.enumerate(closed=True):
            if wrap and not Ps.looped:
                p2 = _unrollon(p1, p2)
            p1 = p2
            v2 = p2._N_vector
            gc = v1.cross(v2)
            t  = gc.angleTo(n0) > PI_2
            if t != t0:  # different sides of edge i
                return False  # outside

            # check for convex polygon: angle between
            # gc vectors, signed by direction of n0
            # (otherwise the test above is not reliable)
            s = signOf(gc1.angleTo(gc, vSign=n0))
            if s != s0:
                if s0 is None:
                    s0 = s
                else:
                    t = _Fmt.SQUARE(points=i)
                    raise _ValueError(t, p2, wrap=wrap, txt_not_=_convex_)
            gc1, v1 = gc, v2

        return True  # inside

    def midpointTo(self, other, height=None, fraction=_0_5, wrap=False):
        '''Find the midpoint between this and an other point.

           @arg other: The other point (L{LatLon}).
           @kwarg height: Optional height for midpoint, overriding
                          the mean height (C{meter}).
           @kwarg fraction: Midpoint location from this point (C{scalar}),
                            may be negative or greater than 1.0.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: Midpoint (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: Invalid B{C{height}}.

           @see: Methods C{intermediateTo} and C{rhumbMidpointTo}.
        '''
        if fraction is _0_5:
            # see <https://MathForum.org/library/drmath/view/51822.html>
            a1, b, a2, _, db = self._ab1_ab2_db5(other, wrap)
            sa1, ca1, sa2, ca2, sdb, cdb = sincos2_(a1, a2, db)

            x = ca2 * cdb + ca1
            y = ca2 * sdb

            a = atan1d(sa1 + sa2, hypot(x, y))
            b = degrees180(b + atan2(y, x))

            h = self._havg(other, h=height)
            r = self.classof(a, b, height=h)
        else:
            r = self.intermediateTo(other, fraction, height=height, wrap=wrap)
        return r

    def nearestOn(self, point1, point2, radius=R_M, **wrap_adjust_limit):
        '''Locate the point between two other points closest to this point.

           Distances are approximated by function L{pygeodesy.equirectangular4},
           subject to the supplied B{C{options}}.

           @arg point1: Start point (L{LatLon}).
           @arg point2: End point (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg wrap_adjust_limit: Optional keyword arguments for functions
                              L{sphericalTrigonometry.nearestOn3} and
                              L{pygeodesy.equirectangular4},

           @return: Closest point on the great circle line (L{LatLon}).

           @raise LimitError: Lat- and/or longitudinal delta exceeds B{C{limit}},
                              see function L{pygeodesy.equirectangular4}.

           @raise NotImplementedError: Keyword argument C{B{within}=False}
                                       is not (yet) supported.

           @raise TypeError: Invalid B{C{point1}} or B{C{point2}}.

           @raise ValueError: Invalid B{C{radius}} or B{C{options}}.

           @see: Functions L{pygeodesy.equirectangular4} and L{pygeodesy.nearestOn5}
                 and method L{sphericalTrigonometry.LatLon.nearestOn3}.
        '''
        # remove kwarg B{C{within}} if present
        w, kwds = _xkwds_pop2(wrap_adjust_limit, within=True)
        if not w:
            self._notImplemented(within=w)

#       # UNTESTED - handle C{B{within}=False} and C{B{within}=True}
#       wrap = _xkwds_get(options, wrap=False)
#       a = self.alongTrackDistanceTo(point1, point2, radius=radius, wrap=wrap)
#       if fabs(a) < EPS or (within and a < EPS):
#           return point1
#       d = point1.distanceTo(point2, radius=radius, wrap=wrap)
#       if isnear0(d):
#           return point1  # or point2
#       elif fabs(d - a) < EPS or (a + EPS) > d:
#           return point2
#       f = a / d
#       if within:
#           if f > EPS1:
#               return point2
#           elif f < EPS:
#               return point1
#       return point1.intermediateTo(point2, f, wrap=wrap)

        # without kwarg B{C{within}}, use backward compatible .nearestOn3
        return self.nearestOn3([point1, point2], closed=False, radius=radius,
                                                             **kwds)[0]

    @deprecated_method
    def nearestOn2(self, points, closed=False, radius=R_M, **options):  # PYCHOK no cover
        '''DEPRECATED, use method L{sphericalTrigonometry.LatLon.nearestOn3}.

           @return: ... 2-Tuple C{(closest, distance)} of the closest
                    point (L{LatLon}) on the polygon and the distance
                    to that point from this point in C{meter}, same
                    units of B{C{radius}}.
        '''
        r = self.nearestOn3(points, closed=closed, radius=radius, **options)
        return r.closest, r.distance

    def nearestOn3(self, points, closed=False, radius=R_M, **wrap_adjust_limit):
        '''Locate the point on a polygon closest to this point.

           Distances are approximated by function L{pygeodesy.equirectangular4},
           subject to the supplied B{C{options}}.

           @arg points: The polygon points (L{LatLon}[]).
           @kwarg closed: Optionally, close the polygon (C{bool}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg wrap_adjust_limit: Optional keyword arguments for function
                              L{sphericalTrigonometry.nearestOn3} and
                              L{pygeodesy.equirectangular4},

           @return: A L{NearestOn3Tuple}C{(closest, distance, angle)} of the
                    C{closest} point (L{LatLon}), the L{pygeodesy.equirectangular4}
                    C{distance} between this and the C{closest} point converted to
                    C{meter}, same units as B{C{radius}}.  The C{angle} from this
                    to the C{closest} point is in compass C{degrees360}, like
                    function L{pygeodesy.compassAngle}.

           @raise LimitError: Lat- and/or longitudinal delta exceeds B{C{limit}},
                              see function L{pygeodesy.equirectangular4}.

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} are not C{LatLon}.

           @raise ValueError: Invalid B{C{radius}} or B{C{options}}.

           @see: Functions L{pygeodesy.compassAngle}, L{pygeodesy.equirectangular4}
                 and L{pygeodesy.nearestOn5}.
        '''
        return nearestOn3(self, points, closed=closed, radius=radius,
                                LatLon=self.classof, **wrap_adjust_limit)

    def toCartesian(self, **Cartesian_datum_kwds):  # PYCHOK Cartesian=Cartesian, datum=None
        '''Convert this point to C{Karney}-based cartesian (ECEF) coordinates.

           @kwarg Cartesian_datum_kwds: Optional L{Cartesian}, B{C{datum}} and other
                                        keyword arguments, ignored if C{B{Cartesian} is
                                        None}.  Use C{B{Cartesian}=...} to override this
                                        L{Cartesian} class or specify C{B{Cartesian}=None}.

           @return: The cartesian point (L{Cartesian}) or if C{B{Cartesian} is None},
                    an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)} with C{C}
                    and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian_datum_kwds}} argument.
        '''
        kwds = _xkwds(Cartesian_datum_kwds, Cartesian=Cartesian, datum=self.datum)
        return LatLonSphericalBase.toCartesian(self, **kwds)

    def triangle7(self, otherB, otherC, radius=R_M, wrap=False):
        '''Compute the angles, sides and area of a spherical triangle.

           @arg otherB: Second triangle point (C{LatLon}).
           @arg otherC: Third triangle point (C{LatLon}).
           @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter}, L{Ellipsoid},
                          L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}) or C{None}.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll points B{C{otherB}}
                        and B{C{otherC}} (C{bool}).

           @return: L{Triangle7Tuple}C{(A, a, B, b, C, c, area)} or if B{C{radius} is
                    None}, a L{Triangle8Tuple}C{(A, a, B, b, C, c, D, E)}.

           @see: Function L{triangle7} and U{Spherical trigonometry
                 <https://WikiPedia.org/wiki/Spherical_trigonometry>}.
        '''
        B = self.others(otherB=otherB)
        C = self.others(otherC=otherC)
        B, C, _ = _unrollon3(self, B, C, wrap)

        r = self.philam + B.philam + C.philam
        t = triangle8_(*r, wrap=wrap)
        return self._xnamed(_t7Tuple(t, radius))

    def triangulate(self, bearing1, other, bearing2, **height_wrap):
        '''Locate a point given this, an other point and a bearing from
           North at both points.

           @arg bearing1: Bearing at this point (compass C{degrees360}).
           @arg other: The other point (C{LatLon}).
           @arg bearing2: Bearing at the other point (compass C{degrees360}).
           @kwarg height_wrap_tol: Optional keyword arguments C{B{height}=None},
                         C{B{wrap}=False}, see method L{intersection}.

           @return: Triangulated point (C{LatLon}).

           @see: Method L{intersection} for further details.
        '''
        if _isDegrees(bearing1) and _isDegrees(bearing2):
            return self.intersection(bearing1, other, bearing2, **height_wrap)
        raise _TypeError(bearing1=bearing1, bearing2=bearing2, **height_wrap)

    def trilaterate5(self, distance1, point2, distance2, point3, distance3,
                           area=True, eps=EPS1, radius=R_M, wrap=False):
        '''Trilaterate three points by I{area overlap} or I{perimeter intersection}
           of three corresponding circles.

           @arg distance1: Distance to this point (C{meter}, same units as B{C{radius}}).
           @arg point2: Second center point (C{LatLon}).
           @arg distance2: Distance to point2 (C{meter}, same units as B{C{radius}}).
           @arg point3: Third center point (C{LatLon}).
           @arg distance3: Distance to point3 (C{meter}, same units as B{C{radius}}).
           @kwarg area: If C{True}, compute the area overlap, otherwise the perimeter
                        intersection of the circles (C{bool}).
           @kwarg eps: The required I{minimal overlap} for C{B{area}=True} or the
                       I{intersection margin} if C{B{area}=False} (C{meter}, same
                       units as B{C{radius}}).
           @kwarg radius: Mean earth radius (C{meter}, conventionally).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{point2}} and
                        B{C{point3}} (C{bool}).

           @return: A L{Trilaterate5Tuple}C{(min, minPoint, max, maxPoint, n)} with
                    C{min} and C{max} in C{meter}, same units as B{C{eps}}, the
                    corresponding trilaterated points C{minPoint} and C{maxPoint}
                    as I{spherical} C{LatLon} and C{n}, the number of trilatered
                    points found for the given B{C{eps}}.

                    If only a single trilaterated point is found, C{min I{is} max},
                    C{minPoint I{is} maxPoint} and C{n = 1}.

                    For C{B{area}=True}, C{min} and C{max} are the smallest respectively
                    largest I{radial} overlap found.

                    For C{B{area}=False}, C{min} and C{max} represent the nearest
                    respectively farthest intersection margin.

                    If C{B{area}=True} and all 3 circles are concentric, C{n=0} and
                    C{minPoint} and C{maxPoint} are both the B{C{point#}} with the
                    smallest B{C{distance#}} C{min} and C{max} the largest B{C{distance#}}.

           @raise IntersectionError: Trilateration failed for the given B{C{eps}},
                                     insufficient overlap for C{B{area}=True} or
                                     no intersection or all (near-)concentric if
                                     C{B{area}=False}.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @raise ValueError: Coincident B{C{point2}} or B{C{point3}} or invalid
                              B{C{distance1}}, B{C{distance2}}, B{C{distance3}}
                              or B{C{radius}}.
        '''
        return _trilaterate5(self, distance1,
                             self.others(point2=point2), distance2,
                             self.others(point3=point3), distance3,
                             area=area, radius=radius, eps=eps, wrap=wrap)


_T00 = LatLon(0, 0, name='T00')  # reference instance (L{LatLon})


def areaOf(points, radius=R_M, wrap=False, polar=False):  # was wrap=True
    '''Calculate the area of a (spherical) polygon or composite (with the points joined by
       great circle arcs).

       @arg points: The polygon points or clips (L{LatLon}[], L{BooleanFHP} or L{BooleanGH}).
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter}, L{Ellipsoid},
                      L{Ellipsoid2}, L{Datum} or L{a_f2Tuple}) or C{None}.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the B{C{points}} (C{bool}).
       @kwarg polar: Use C{B{polar}=True} if the polygon encloses a pole (C{bool}), see
                     function L{ispolar<pygeodesy.points.ispolar>} and U{area of a polygon
                     enclosing a pole<https://GeographicLib.SourceForge.io/C++/doc/
                     classGeographicLib_1_1GeodesicExact.html#a3d7a9155e838a09a48dc14d0c3fac525>}.

       @return: Polygon area (C{meter} I{quared}, same units as B{C{radius}} or C{radians} if
                C{B{radius} is None}).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: Invalid B{C{radius}} or semi-circular polygon edge.

       @note: The area is based on I{Karney}'s U{'Area of a spherical polygon'
              <https://MathOverflow.net/questions/97711/ the-area-of-spherical-polygons>}, 3rd Answer.

       @see: Functions L{pygeodesy.areaOf}, L{sphericalNvector.areaOf}, L{pygeodesy.excessKarney},
             L{ellipsoidalExact.areaOf} and L{ellipsoidalKarney.areaOf}.
    '''
    if _MODS.booleans.isBoolean(points):
        return points._sum2(LatLon, areaOf, radius=radius, wrap=wrap)

    _at2, _t_2 = atan2, tan_2
    _un, _w180 = unrollPI, wrap180

    Ps = _T00.PointsIter(points, loop=1, wrap=wrap)
    p1 = p2 =  Ps[0]
    a1,  b1 =  p1.philam
    ta1, z1 = _t_2(a1), None

    A = Fsum()  # mean phi
    R = Fsum()  # see L{pygeodesy.excessKarney_}
    # ispolar: Summation of course deltas around pole is 0° rather than normally ±360°
    # <https://blog.Element84.com/determining-if-a-spherical-polygon-contains-a-pole.html>
    # XXX duplicate of function C{points.ispolar} to avoid copying all iterated points
    D = Fsum()
    for i, p2 in Ps.enumerate(closed=True):
        a2, b2 =  p2.philam
        db, b2 = _un(b1, b2, wrap=wrap and not Ps.looped)
        A  +=  a2
        ta2 = _t_2(a2)
        tdb = _t_2(db, points=i)
        R  += _at2(tdb * (ta1 + ta2),
                   _1_0 + ta1 * ta2)
        ta1, b1 = ta2, b2

        if not p2.isequalTo(p1, eps=EPS):
            z, z2 = _bearingTo2(p1, p2, wrap=wrap)
            if z1 is not None:
                D += _w180(z - z1)  # (z - z1 + 540) ...
            D += _w180(z2 - z)  # (z2 - z + 540) % 360 - 180
            p1, z1 = p2, z2

    R = abs(R * _2_0)
    if  abs(D) < _90_0 or polar:  # ispolar(points)
        R = abs(R - PI2)
    if radius:
        a  =  degrees(A.fover(len(A)))  # mean lat
        R *= _mean_radius(radius, a)**2
    return float(R)


def _destination2(a, b, r, t):
    '''(INTERNAL) Destination lat- and longitude in C{radians}.

       @arg a: Latitude (C{radians}).
       @arg b: Longitude (C{radians}).
       @arg r: Angular distance (C{radians}).
       @arg t: Bearing (compass C{radians}).

       @return: 2-Tuple (phi, lam) of (C{radians}, C{radiansPI}).
    '''
    # see <https://www.EdWilliams.org/avform.htm#LL>
    sa, ca, sr, cr, st, ct = sincos2_(a, r, t)
    ca *= sr

    a = asin1(ct * ca + cr * sa)
    d = atan2(st * ca,  cr - sa * sin(a))
    # note, in EdWilliams.org/avform.htm W is + and E is -
    return a, (b + d)  # (mod(b + d + PI, PI2) - PI)


def _int3d2(s, end, wrap, _i_, Vector, hs):
    # see <https://www.EdWilliams.org/intersect.htm> (5) ff
    # and similar logic in .ellipsoidalBaseDI._intersect3
    a1, b1 = s.philam

    if _isDegrees(end):  # bearing, get pseudo-end point
        a2, b2 = _destination2(a1, b1, PI_4, radians(end))
    else:  # must be a point
        s.others(end, name=_end_ + _i_)
        hs.append(end.height)
        a2, b2 = end.philam
        if wrap:
            a2, b2 = _Wrap.philam(a2, b2)

    db, b2 = unrollPI(b1, b2, wrap=wrap)
    if max(fabs(db), fabs(a2 - a1)) < EPS:
        raise _ValueError(_SPACE_(_line_ + _i_, _null_))
    # note, in EdWilliams.org/avform.htm W is + and E is -
    sb21, cb21, sb12, cb12 = sincos2_(db * _0_5,
                                    -(b1 + b2) * _0_5)
    cb21 *= sin(a1 - a2)  # sa21
    sb21 *= sin(a1 + a2)  # sa12
    x = Vector(fdot_(sb12, cb21, -cb12, sb21),
               fdot_(cb12, cb21,  sb12, sb21),
               cos(a1) * cos(a2) * sin(db))  # ll=start
    return x.unit(), (db, (a2 - a1))  # negated d


def _intdot(ds, a1, b1, a, b, wrap):
    # compute dot product ds . (-b + b1, a - a1)
    db, _ = unrollPI(b1, b, wrap=wrap)
    return fdot(ds, db, a - a1)


def intersecant2(center, circle, point, other, **radius_exact_height_wrap):
    '''Compute the intersections of a circle and a (great circle) line given as
       two points or as a point and a bearing from North.

       @arg center: Center of the circle (L{LatLon}).
       @arg circle: Radius of the circle (C{meter}, same units as the earth
                    B{C{radius}}) or a point on the circle (L{LatLon}).
       @arg point: A point on the (great circle) line (L{LatLon}).
       @arg other: An other point on the (great circle) line (L{LatLon}) or
                   the bearing at the B{C{point}} (compass C{degrees360}).
       @kwarg radius_exact_height_wrap: Optional keyword arguments, see method
                     L{intersecant2<pygeodesy.sphericalBase.LatLonSphericalBase.
                     intersecant2>} for further details.

       @return: 2-Tuple of the intersection points (representing a chord), each
                an instance of the B{C{point}} class.  Both points are the same
                instance if the (great circle) line is tangent to the circle.

       @raise IntersectionError: The circle and line do not intersect.

       @raise TypeError: If B{C{center}}, B{C{point}}, B{C{circle}} or B{C{other}}
                         not L{LatLon}.

       @raise UnitError: Invalid B{C{circle}}, B{C{other}}, B{C{radius}},
                         B{C{exact}}, B{C{height}} or B{C{napieradius}}.
    '''
    c = _T00.others(center=center)
    p = _T00.others(point=point)
    try:
        return _intersecant2(c, circle, p, other, **radius_exact_height_wrap)
    except (TypeError, ValueError) as x:
        raise _xError(x, center=center, circle=circle, point=point, other=other,
                                                     **radius_exact_height_wrap)


def _intersect(start1, end1, start2, end2, height=None, wrap=False,  # in.ellipsoidalBaseDI._intersect3
                                           LatLon=LatLon, **LatLon_kwds):
    # (INTERNAL) Intersect two (spherical) lines, see L{intersection}
    # above, separated to allow callers to embellish any exceptions

    s1, s2 = start1, start2
    if wrap:
        s2 = _Wrap.point(s2)
    hs = [s1.height, s2.height]

    a1, b1 = s1.philam
    a2, b2 = s2.philam
    db, b2 = unrollPI(b1, b2, wrap=wrap)
    r12 = vincentys_(a2, a1, db)
    if fabs(r12) < EPS:  # [nearly] coincident points
        a, b = favg(a1, a2), favg(b1, b2)

    # see <https://www.EdWilliams.org/avform.htm#Intersection>
    elif _isDegrees(end1) and _isDegrees(end2):  # both bearings
        sa1, ca1, sa2, ca2, sr12, cr12 = sincos2_(a1, a2, r12)

        x1, x2 = (sr12 * ca1), (sr12 * ca2)
        if isnear0(x1) or isnear0(x2):
            raise IntersectionError(_parallel_)
        # handle domain error for equivalent longitudes,
        # see also functions asin_safe and acos_safe at
        # <https://www.EdWilliams.org/avform.htm#Math>
        t12, t13 = acos1((sa2 - sa1 * cr12) / x1), radiansPI2(end1)
        t21, t23 = acos1((sa1 - sa2 * cr12) / x2), radiansPI2(end2)
        if sin(db) > 0:
            t21 = PI2 - t21
        else:
            t12 = PI2 - t12
        sx1, cx1, sx2, cx2 = sincos2_(wrapPI(t13 - t12),  # angle 2-1-3
                                      wrapPI(t21 - t23))  # angle 1-2-3)
        if isnear0(sx1) and isnear0(sx2):
            raise IntersectionError(_infinite_)
        sx3 = sx1 * sx2
# XXX   if sx3 < 0:
# XXX       raise ValueError(_ambiguous_)
        x3  = acos1(cr12 * sx3 - cx2 * cx1)
        r13 = atan2(sr12 * sx3,  cx2 + cx1 * cos(x3))

        a, b = _destination2(a1, b1, r13, t13)
        # like .ellipsoidalBaseDI,_intersect3, if this intersection
        # is "before" the first point, use the antipodal intersection
        if opposing_(t13, bearing_(a1, b1, a, b, wrap=wrap)):
            a, b = antipode_(a, b)  # PYCHOK PhiLam2Tuple

    else:  # end point(s) or bearing(s)
        _N_vector_ = _MODS.nvectorBase._N_vector_

        x1, d1 = _int3d2(s1, end1, wrap, _1_, _N_vector_, hs)
        x2, d2 = _int3d2(s2, end2, wrap, _2_, _N_vector_, hs)
        x = x1.cross(x2)
        if x.length < EPS:  # [nearly] colinear or parallel lines
            raise IntersectionError(_colinear_)
        a, b = x.philam
        # choose intersection similar to sphericalNvector
        if not (_intdot(d1, a1, b1, a, b, wrap) *
                _intdot(d2, a2, b2, a, b, wrap)) > 0:
            a, b = antipode_(a, b)  # PYCHOK PhiLam2Tuple

    h = fmean(hs) if height is None else Height(height)
    return _LL3Tuple(degrees90(a), degrees180(b), h,
                     intersection, LatLon, LatLon_kwds)


def intersection(start1, end1, start2, end2, height=None, wrap=False,
                                             **LatLon_and_kwds):
    '''Compute the intersection point of two lines, each defined by
       two points or by a start point and a bearing from North.

       @arg start1: Start point of the first line (L{LatLon}).
       @arg end1: End point of the first line (L{LatLon}) or the bearing
                  at the first start point (compass C{degrees360}).
       @arg start2: Start point of the second line (L{LatLon}).
       @arg end2: End point of the second line (L{LatLon}) or the bearing
                  at the second start point (compass C{degrees360}).
       @kwarg height: Optional height for the intersection point,
                      overriding the mean height (C{meter}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{start2}}
                    and both B{C{end*}} points (C{bool}).
       @kwarg LatLon_and_kwds: Optional class C{B{LatLon}=}L{LatLon} to use
                     for the intersection point and optionally additional
                     B{C{LatLon}} keyword arguments, ignored if C{B{LatLon}
                     is None}.

       @return: The intersection point as a (B{C{LatLon}}) or if C{B{LatLon}
                is None} a L{LatLon3Tuple}C{(lat, lon, height)}.  An alternate
                intersection point might be the L{antipode} to the returned result.

       @raise IntersectionError: Ambiguous or infinite intersection or colinear,
                                 parallel or otherwise non-intersecting lines.

       @raise TypeError: A B{C{start1}}, B{C{end1}}, B{C{start2}} or B{C{end2}}
                         point not L{LatLon}.

       @raise ValueError: Invalid B{C{height}} or C{null} line.
    '''
    s1 = _T00.others(start1=start1)
    s2 = _T00.others(start2=start2)
    try:
        return _intersect(s1, end1, s2, end2, height=height, wrap=wrap, **LatLon_and_kwds)
    except (TypeError, ValueError) as x:
        raise _xError(x, start1=start1, end1=end1, start2=start2, end2=end2)


def intersections2(center1, rad1, center2, rad2, radius=R_M, eps=_0_0,
                                                 height=None, wrap=False,  # was=True
                                               **LatLon_and_kwds):
    '''Compute the intersection points of two circles each defined by a
       center point and a radius.

       @arg center1: Center of the first circle (L{LatLon}).
       @arg rad1: Radius of the first circle (C{meter} or C{radians}, see
                  B{C{radius}}).
       @arg center2: Center of the second circle (L{LatLon}).
       @arg rad2: Radius of the second circle (C{meter} or C{radians}, see
                  B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter} or C{None} if B{C{rad1}},
                      B{C{rad2}} and B{C{eps}} are given in C{radians}).
       @kwarg eps: Required overlap (C{meter} or C{radians}, see B{C{radius}}).
       @kwarg height: Optional height for the intersection points (C{meter},
                      conventionally) or C{None} for the I{"radical height"}
                      at the I{radical line} between both centers.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll B{C{center2}}
                    (C{bool}).
       @kwarg LatLon_and_kwds: Optional class C{B{LatLon}=}L{LatLon} to use for
                     the intersection points and optionally additional B{C{LatLon}}
                     keyword arguments, ignored if C{B{LatLon} is None}.

       @return: 2-Tuple of the intersection points, each a B{C{LatLon}}
                instance or if C{B{LatLon} is None} a L{LatLon3Tuple}C{(lat,
                lon, height)}.  For abutting circles, both intersection
                points are the same instance, aka the I{radical center}.

       @raise IntersectionError: Concentric, antipodal, invalid or
                                 non-intersecting circles.

       @raise TypeError: If B{C{center1}} or B{C{center2}} not L{LatLon}.

       @raise ValueError: Invalid B{C{rad1}}, B{C{rad2}}, B{C{radius}},
                          B{C{eps}} or B{C{height}}.

       @note: Courtesy of U{Samuel Čavoj<https://GitHub.com/mrJean1/PyGeodesy/issues/41>}.

       @see: This U{Answer<https://StackOverflow.com/questions/53324667/
             find-intersection-coordinates-of-two-circles-on-earth/53331953>}.
    '''
    c1 = _T00.others(center1=center1)
    c2 = _T00.others(center2=center2)
    try:
        return _intersects2(c1, rad1, c2, rad2, radius=radius, eps=eps,
                                                height=height, wrap=wrap,
                                              **LatLon_and_kwds)
    except (TypeError, ValueError) as x:
        raise _xError(x, center1=center1, rad1=rad1,
                         center2=center2, rad2=rad2, wrap=wrap)


def _intersects2(c1, rad1, c2, rad2, radius=R_M, eps=_0_0,  # in .ellipsoidalBaseDI._intersects2
                                     height=None, too_d=None, wrap=False,  # was=True
                                     LatLon=LatLon, **LatLon_kwds):
    # (INTERNAL) Intersect two spherical circles, see L{intersections2}
    # above, separated to allow callers to embellish any exceptions

    def _dest3(bearing, h):
        a, b = _destination2(a1, b1, r1, bearing)
        return _LL3Tuple(degrees90(a), degrees180(b), h,
                         intersections2, LatLon, LatLon_kwds)

    a1, b1 = c1.philam
    a2, b2 = c2.philam
    if wrap:
        a2, b2 = _Wrap.philam(a2, b2)

    r1, r2, f = _rads3(rad1, rad2, radius)
    if f:  # swapped radii, swap centers
        a1, a2 = a2, a1  # PYCHOK swap!
        b1, b2 = b2, b1  # PYCHOK swap!

    db, b2 = unrollPI(b1, b2, wrap=wrap)
    d = vincentys_(a2, a1, db)  # radians
    if d < max(r1 - r2, EPS):
        raise IntersectionError(_near_(_concentric_))  # XXX ConcentricError?

    r = eps if radius is None else (m2radians(
        eps, radius=radius) if eps else _0_0)
    if r < _0_0:
        raise _ValueError(eps=r)

    x = fsumf_(r1, r2, -d)  # overlap
    if x > max(r, EPS):
        sd, cd, sr1, cr1, _, cr2 = sincos2_(d, r1, r2)
        x = sd * sr1
        if isnear0(x):
            raise _ValueError(_invalid_)
        x = acos1((cr2 - cd * cr1) / x)  # 0 <= x <= PI

    elif x < r:  # PYCHOK no cover
        t = (d * radius) if too_d is None else too_d
        raise IntersectionError(_too_(_Fmt.distant(t)))

    if height is None:  # "radical height"
        f = _radical2(d, r1, r2).ratio
        h = Height(favg(c1.height, c2.height, f=f))
    else:
        h = Height(height)

    b = bearing_(a1, b1, a2, b2, final=False, wrap=wrap)
    if x < EPS4:  # externally ...
        r = _dest3(b, h)
    elif x > _PI_EPS4:  # internally ...
        r = _dest3(b + PI, h)
    else:
        return _dest3(b + x, h), _dest3(b - x, h)
    return r, r  # ... abutting circles


@deprecated_function
def isPoleEnclosedBy(points, wrap=False):  # PYCHOK no cover
    '''DEPRECATED, use function L{pygeodesy.ispolar}.
    '''
    return ispolar(points, wrap=wrap)


def _LL3Tuple(lat, lon, height, where, LatLon, LatLon_kwds):
    '''(INTERNAL) Helper for L{intersection}, L{intersections2} and L{meanOf}.
    '''
    n = typename(where)
    if LatLon is None:
        r = LatLon3Tuple(lat, lon, height, name=n)
    else:
        kwds = _xkwds(LatLon_kwds, height=height, name=n)
        r = LatLon(lat, lon, **kwds)
    return r


def meanOf(points, height=None, wrap=False, LatLon=LatLon, **LatLon_kwds):
    '''Compute the I{geographic} mean of several points.

       @arg points: Points to be averaged (L{LatLon}[]).
       @kwarg height: Optional height at mean point, overriding the mean height
                      (C{meter}).
       @kwarg wrap: If C{True}, wrap or I{normalize} the B{C{points}} (C{bool}).
       @kwarg LatLon: Optional class to return the mean point (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword arguments,
                           ignored if C{B{LatLon} is None}.

       @return: The geographic mean and height (B{C{LatLon}}) or if C{B{LatLon}
                is None}, a L{LatLon3Tuple}C{(lat, lon, height)}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: No B{C{points}} or invalid B{C{height}}.
    '''
    def _N_vs(ps, w):
        Ps = _T00.PointsIter(ps, wrap=w)
        for p in Ps.iterate(closed=False):
            yield p._N_vector

    m = _MODS.nvectorBase
    # geographic, vectorial mean
    n =  m.sumOf(_N_vs(points, wrap), h=height, Vector=m.NvectorBase)
    lat, lon, h = n.latlonheight
    return _LL3Tuple(lat, lon, h, meanOf, LatLon, LatLon_kwds)


@deprecated_function
def nearestOn2(point, points, **closed_radius_LatLon_options):  # PYCHOK no cover
    '''DEPRECATED, use function L{sphericalTrigonometry.nearestOn3}.

       @return: ... 2-tuple C{(closest, distance)} of the C{closest}
                point (L{LatLon}) on the polygon and the C{distance}
                between the C{closest} and the given B{C{point}}.  The
                C{closest} is a B{C{LatLon}} or a L{LatLon2Tuple}C{(lat,
                lon)} if C{B{LatLon} is None} ...
    '''
    ll, d, _ = nearestOn3(point, points, **closed_radius_LatLon_options)  # PYCHOK 3-tuple
    if _xkwds_get(closed_radius_LatLon_options, LatLon=LatLon) is None:
        ll = LatLon2Tuple(ll.lat, ll.lon)
    return ll, d


def nearestOn3(point, points, closed=False, radius=R_M, wrap=False, adjust=True,
                                                        limit=9, **LatLon_and_kwds):
    '''Locate the point on a path or polygon closest to a reference point.

       Distances are I{approximated} using function L{equirectangular4
       <pygeodesy.equirectangular4>}, subject to the supplied B{C{options}}.

       @arg point: The reference point (L{LatLon}).
       @arg points: The path or polygon points (L{LatLon}[]).
       @kwarg closed: Optionally, close the polygon (C{bool}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                    B{C{points}} (C{bool}).
       @kwarg adjust: See function L{equirectangular4<pygeodesy.equirectangular4>}
                      (C{bool}).
       @kwarg limit: See function L{equirectangular4<pygeodesy.equirectangular4>}
                     (C{degrees}), default C{9 degrees} is about C{1,000 Km} (for
                     (mean spherical earth radius L{R_KM}).
       @kwarg LatLon_and_kwds: Optional class C{B{LatLon}=L{LatLon}} to return the
                     closest point and optionally additional C{B{LatLon}} keyword
                     arguments or specify C{B{LatLon}=None}.

       @return: A L{NearestOn3Tuple}C{(closest, distance, angle)} with the
                C{closest} point as B{C{LatLon}} or L{LatLon3Tuple}C{(lat,
                lon, height)} if C{B{LatLon} is None}.  The C{distance} is
                the L{equirectangular4<pygeodesy.equirectangular4>} distance
                between the C{closest} and the given B{C{point}} converted to
                C{meter}, same units as B{C{radius}}.  The C{angle} from the
                given B{C{point}} to the C{closest} is in compass C{degrees360},
                like function L{compassAngle<pygeodesy.compassAngle>}.  The
                C{height} is the (interpolated) height at the C{closest} point.

       @raise LimitError: Lat- and/or longitudinal delta exceeds the B{C{limit}},
                          see function L{equirectangular4<pygeodesy.equirectangular4>}.

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise ValueError: Invalid B{C{radius}}.

       @see: Functions L{equirectangular4<pygeodesy.equirectangular4>} and
             L{nearestOn5<pygeodesy.nearestOn5>}.
    '''
    t = _nearestOn5(point, points, closed=closed, wrap=wrap,
                                   adjust=adjust, limit=limit)
    d = degrees2m(t.distance, radius=radius)
    h = t.height
    n = typename(nearestOn3)

    LL, kwds = _xkwds_pop2(LatLon_and_kwds, LatLon=LatLon)
    r = LatLon3Tuple(t.lat, t.lon, h, name=n) if LL is None else \
                  LL(t.lat, t.lon, **_xkwds(kwds, height=h, name=n))
    return NearestOn3Tuple(r, d, t.angle, name=n)


def perimeterOf(points, closed=False, radius=R_M, wrap=True):
    '''Compute the perimeter of a (spherical) polygon or composite
       (with great circle arcs joining the points).

       @arg points: The polygon points or clips (L{LatLon}[], L{BooleanFHP}
                    or L{BooleanGH}).
       @kwarg closed: Optionally, close the polygon (C{bool}).
       @kwarg radius: Mean earth radius (C{meter}) or C{None}.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                    B{C{points}} (C{bool}).

       @return: Polygon perimeter (C{meter}, same units as B{C{radius}}
                or C{radians} if C{B{radius} is None}).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: Invalid B{C{radius}} or C{B{closed}=False} with
                          C{B{points}} a composite.

       @note: Distances are based on function L{vincentys_<pygeodesy.vincentys_>}.

       @see: Functions L{perimeterOf<pygeodesy.perimeterOf>},
             L{sphericalNvector.perimeterOf} and L{ellipsoidalKarney.perimeterOf}.
    '''
    def _rads(ps, c, w):  # angular edge lengths in radians
        Ps = _T00.PointsIter(ps, loop=1, wrap=w)
        a1, b1 = Ps[0].philam
        for p in Ps.iterate(closed=c):
            a2, b2 = p.philam
            db, b2 = unrollPI(b1, b2, wrap=w and not (c and Ps.looped))
            yield vincentys_(a2, a1, db)
            a1, b1 = a2, b2

    if _MODS.booleans.isBoolean(points):
        if not closed:
            raise _ValueError(closed=closed, points=_composite_)
        r = points._sum2(LatLon, perimeterOf, closed=True, radius=radius, wrap=wrap)
    else:
        r = fsum(_rads(points, closed, wrap))
    return _radians2m(r, radius)


def triangle7(latA, lonA, latB, lonB, latC, lonC, radius=R_M,
                                                  excess=excessAbc_,
                                                    wrap=False):
    '''Compute the angles, sides, and area of a (spherical) triangle.

       @arg latA: First corner latitude (C{degrees}).
       @arg lonA: First corner longitude (C{degrees}).
       @arg latB: Second corner latitude (C{degrees}).
       @arg lonB: Second corner longitude (C{degrees}).
       @arg latC: Third corner latitude (C{degrees}).
       @arg lonC: Third corner longitude (C{degrees}).
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple})
                      or C{None}.
       @kwarg excess: I{Spherical excess} callable (L{excessAbc_},
                      L{excessGirard_} or L{excessLHuilier_}).
       @kwarg wrap: If C{True}, wrap and L{unroll180<pygeodesy.unroll180>}
                    longitudes (C{bool}).

       @return: A L{Triangle7Tuple}C{(A, a, B, b, C, c, area)} with
                spherical angles C{A}, C{B} and C{C}, angular sides
                C{a}, C{b} and C{c} all in C{degrees} and C{area}
                in I{square} C{meter} or same units as B{C{radius}}
                I{squared} or if C{B{radius}=0} or C{None}, a
                L{Triangle8Tuple}C{(A, a, B, b, C, c, D, E)} with
                I{spherical deficit} C{D} and I{spherical excess}
                C{E} as the C{unit area}, all in C{radians}.
    '''
    t = triangle8_(Phid(latA=latA), Lamd(lonA=lonA),
                   Phid(latB=latB), Lamd(lonB=lonB),
                   Phid(latC=latC), Lamd(lonC=lonC),
                   excess=excess, wrap=wrap)
    return _t7Tuple(t, radius)


def triangle8_(phiA, lamA, phiB, lamB, phiC, lamC, excess=excessAbc_,
                                                     wrap=False):
    '''Compute the angles, sides, I{spherical deficit} and I{spherical
       excess} of a (spherical) triangle.

       @arg phiA: First corner latitude (C{radians}).
       @arg lamA: First corner longitude (C{radians}).
       @arg phiB: Second corner latitude (C{radians}).
       @arg lamB: Second corner longitude (C{radians}).
       @arg phiC: Third corner latitude (C{radians}).
       @arg lamC: Third corner longitude (C{radians}).
       @kwarg excess: I{Spherical excess} callable (L{excessAbc_},
                      L{excessGirard_} or L{excessLHuilier_}).
       @kwarg wrap: If C{True}, L{unrollPI<pygeodesy.unrollPI>} the
                    longitudinal deltas (C{bool}).

       @return: A L{Triangle8Tuple}C{(A, a, B, b, C, c, D, E)} with
                spherical angles C{A}, C{B} and C{C}, angular sides
                C{a}, C{b} and C{c}, I{spherical deficit} C{D} and
                I{spherical excess} C{E}, all in C{radians}.
    '''
    def _a_r(w, phiA, lamA, phiB, lamB, phiC, lamC):
        d, _ = unrollPI(lamB, lamC, wrap=w)
        a = vincentys_(phiC, phiB, d)
        return a, (phiB, lamB, phiC, lamC, phiA, lamA)  # rotate A, B, C

    def _A_r(a, sa, ca, sb, cb, sc, cc):
        s = sb * sc
        A = acos1((ca - cb * cc) / s) if isnon0(s) else a
        return A, (sb, cb, sc, cc, sa, ca)  # rotate sincos2_'s

    # notation: side C{a} is oposite to corner C{A}, etc.
    a, r = _a_r(wrap, phiA, lamA, phiB, lamB, phiC, lamC)
    b, r = _a_r(wrap, *r)
    c, _ = _a_r(wrap, *r)

    A, r = _A_r(a, *sincos2_(a, b, c))
    B, r = _A_r(b, *r)
    C, _ = _A_r(c, *r)

    D = fsumf_(PI2, -a, -b, -c)  # deficit aka defect
    E = excessGirard_(A, B, C)   if _isin(excess, excessGirard_,   True)  else (
        excessLHuilier_(a, b, c) if _isin(excess, excessLHuilier_, False) else
        excessAbc_(*max((A, b, c), (B, c, a), (C, a, b))))

    return Triangle8Tuple(A, a, B, b, C, c, D, E)


def _t7Tuple(t, radius):
    '''(INTERNAL) Convert a L{Triangle8Tuple} to L{Triangle7Tuple}.
    '''
    if radius:  # not _isin(radius, None, _0_0)
        r = radius if _isRadius(radius) else \
            _ellipsoidal_datum(radius).ellipsoid.Rmean
        A, B, C = map1(degrees, t.A, t.B, t.C)
        t = Triangle7Tuple(A, (r * t.a),
                           B, (r * t.b),
                           C, (r * t.c), t.E * r**2)
    return t


__all__ += _ALL_OTHER(Cartesian, LatLon,  # classes
                      areaOf,  # functions
                      intersecant2, intersection, intersections2, ispolar,
                      isPoleEnclosedBy,  # DEPRECATED, use ispolar
                      meanOf,
                      nearestOn2, nearestOn3,
                      perimeterOf,
                      sumOf,  # XXX == vector3d.sumOf
                      triangle7, triangle8_)

# **) MIT License
#
# Copyright (C) 2016-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
