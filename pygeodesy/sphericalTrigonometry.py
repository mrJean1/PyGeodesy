
# -*- coding: utf-8 -*-

u'''Trigonometric classes geodetic (lat-/longitude) L{LatLon} and
geocentric (ECEF) L{Cartesian} and functions L{areaOf}, L{intersection},
L{isPoleEnclosedBy}, L{meanOf}, L{nearestOn2} and L{perimeterOf},
I{all spherical}.

Pure Python implementation of geodetic (lat-/longitude) methods using
spherical trigonometry, transcribed from JavaScript originals by
I{(C) Chris Veness 2011-2016} published under the same MIT Licence**, see
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>}.

@newfield example: Example, Examples
'''

from pygeodesy.basics import EPS, PI, PI2, PI_2, PI_4, R_M, \
                             isscalar, map1, _xkwds
from pygeodesy.errors import _AssertionError, CrossError, crosserrors, \
                              IntersectionError, _ValueError, _xkwds_get
from pygeodesy.fmath import favg, fdot, fmean, fsum, fsum_
from pygeodesy.formy import antipode_, bearing_, vincentys_
from pygeodesy.interns import _1_, _2_, _coincident_, _colinear_, _end_, \
                              _fraction_, _invalid_, _item_sq, \
                              _near_concentric_, _not_convex_, _points_, \
                              _start_, _start1_, _start2_, _too_distant_fmt_
from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER
from pygeodesy.named import LatLon2Tuple, LatLon3Tuple, NearestOn3Tuple, \
                           _xnamed
from pygeodesy.nvectorBase import NvectorBase as _Nvector
from pygeodesy.points import _imdex2, ispolar, nearestOn5 as _nearestOn5
from pygeodesy.sphericalBase import _angular, CartesianSphericalBase, \
                                     LatLonSphericalBase, _rads3
from pygeodesy.units import Bearing_, Height, Radius, Radius_, Scalar
from pygeodesy.utily import acos1, asin1, degrees90, degrees180, degrees2m, \
                            iterNumpy2, radiansPI2, sincos2, tan_2, \
                            unrollPI, wrapPI
from pygeodesy.vector3d import _radical2, sumOf, Vector3d

from math import asin, atan2, copysign, cos, degrees, hypot, radians, sin

__all__ = _ALL_LAZY.sphericalTrigonometry
__version__ = '20.08.04'

_EPS_I2    = 4.0 * EPS
_PI_EPS_I2 = PI - _EPS_I2
if _PI_EPS_I2 >= PI:
    raise _AssertionError(_PI_EPS_I2=_PI_EPS_I2)


def _destination2(a, b, r, t):
    '''(INTERNAL) Destination phi- and longitude in C{radians}.

       @arg a: Latitude (C{radians}).
       @arg b: Longitude (C{radians}).
       @arg r: Angular distance (C{radians}).
       @arg t: Bearing (compass C{radians}).

       @return: 2-Tuple (phi, lam) of (C{radians}, C{radiansPI}).
    '''
    # see <https://www.EdWilliams.org/avform.htm#LL>
    sa, ca, sr, cr, st, ct = sincos2(a, r, t)

    a = asin1(ct * sr * ca + cr * sa)
    d = atan2(st * sr * ca, cr - sa * sin(a))
    # note, in EdWilliams.org/avform.htm W is + and E is -
    return a, b + d


class Cartesian(CartesianSphericalBase):
    '''Extended to convert geocentric, L{Cartesian} points to
       spherical, geodetic L{LatLon}.
    '''

    def toLatLon(self, **LatLon_datum_kwds):  # PYCHOK LatLon=LatLon
        '''Convert this cartesian point to an C{Nvector}-based
           geodetic point.

           @kwarg LatLon_datum_kwds: Optional L{LatLon}, B{C{datum}} and
                                     other keyword arguments, ignored if
                                     B{C{LatLon=None}}.  Use
                                     B{C{LatLon=...}} to override this
                                     L{LatLon} class or specify
                                     B{C{LatLon=None}}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}}
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon,
                    height, C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon_datum_kwds}}.
        '''
        kwds = _xkwds(LatLon_datum_kwds, LatLon=LatLon, datum=self.datum)
        return CartesianSphericalBase.toLatLon(self, **kwds)


class LatLon(LatLonSphericalBase):
    '''New point on spherical model earth model.

       @example:

       >>> p = LatLon(52.205, 0.119)  # height=0
    '''

    def _trackDistanceTo3(self, start, end, radius, wrap):
        '''(INTERNAL) Helper for along-/crossTrackDistanceTo.
        '''
        self.others(start, name=_start_)
        self.others(end, name=_end_)

        r = Radius_(radius)
        r = start.distanceTo(self, r, wrap=wrap) / r

        b = radians(start.initialBearingTo(self, wrap=wrap))
        e = radians(start.initialBearingTo(end, wrap=wrap))
        x = asin(sin(r) * sin(b - e))
        return r, x, (e - b)

    def alongTrackDistanceTo(self, start, end, radius=R_M, wrap=False):
        '''Compute the (signed) distance from the start to the closest
           point on the great circle path defined by a start and an
           end point.

           That is, if a perpendicular is drawn from this point to the
           great circle path, the along-track distance is the distance
           from the start point to the point where the perpendicular
           crosses the path.

           @arg start: Start point of great circle path (L{LatLon}).
           @arg end: End point of great circle path (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance along the great circle path (C{meter},
                    same units as B{C{radius}}), positive if after the
                    B{C{start}} toward the B{C{end}} point of the path or
                    negative if before the B{C{start}} point.

           @raise TypeError: The B{C{start}} or B{C{end}} point is not L{LatLon}.

           @raise ValueError: Invalid B{C{radius}}.

           @example:

           >>> p = LatLon(53.2611, -0.7972)

           >>> s = LatLon(53.3206, -1.7297)
           >>> e = LatLon(53.1887, 0.1334)
           >>> d = p.alongTrackDistanceTo(s, e)  # 62331.58
        '''
        r, x, b = self._trackDistanceTo3(start, end, radius, wrap)
        cx = cos(x)
        if abs(cx) > EPS:
            return copysign(acos1(cos(r) / cx), cos(b)) * radius
        else:
            return 0.0

    def bearingTo(self, other, wrap=False, raiser=False):  # PYCHOK no cover
        '''DEPRECATED, use method C{initialBearingTo}.
        '''
        return self.initialBearingTo(other, wrap=wrap, raiser=raiser)

    def crossingParallels(self, other, lat, wrap=False):
        '''Return the pair of meridians at which a great circle defined
           by this and an other point crosses the given latitude.

           @arg other: The other point defining great circle (L{LatLon}).
           @arg lat: Latitude at the crossing (C{degrees}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: 2-Tuple C{(lon1, lon2)}, both in C{degrees180} or
                    C{None} if the great circle doesn't reach B{C{lat}}.
        '''
        self.others(other)

        a1, b1 = self.philam
        a2, b2 = other.philam

        a = radians(lat)
        db, b2 = unrollPI(b1, b2, wrap=wrap)

        sa,  ca,  sa1, ca1, \
        sa2, ca2, sdb, cdb = sincos2(a, a1, a2, db)

        x = sa1 * ca2 * ca * sdb
        y = sa1 * ca2 * ca * cdb - ca1 * sa2 * ca
        z = ca1 * ca2 * sa * sdb

        h = hypot(x, y)
        if h < EPS or abs(z) > h:
            return None  # great circle doesn't reach latitude

        m = atan2(-y, x) + b1  # longitude at max latitude
        d = acos1(z / h)  # delta longitude to intersections
        return degrees180(m - d), degrees180(m + d)

    def crossTrackDistanceTo(self, start, end, radius=R_M, wrap=False):
        '''Compute the (signed) distance from this point to the great
           circle defined by a start and an end point.

           @arg start: Start point of great circle path (L{LatLon}).
           @arg end: End point of great circle path (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance to great circle (negative if to the
                    left or positive if to the right of the path).

           @raise TypeError: The B{C{start}} or B{C{end}} point is not L{LatLon}.

           @raise ValueError: Invalid B{C{radius}}.

           @example:

           >>> p = LatLon(53.2611, -0.7972)

           >>> s = LatLon(53.3206, -1.7297)
           >>> e = LatLon(53.1887, 0.1334)
           >>> d = p.crossTrackDistanceTo(s, e)  # -307.5
        '''
        _, x, _ = self._trackDistanceTo3(start, end, radius, wrap)
        return x * radius

    def destination(self, distance, bearing, radius=R_M, height=None):
        '''Locate the destination from this point after having
           travelled the given distance on the given initial bearing.

           @arg distance: Distance travelled (C{meter}, same units as
                          B{C{radius}}).
           @arg bearing: Bearing from this point (compass C{degrees360}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg height: Optional height at destination (C{meter}, same
                          units a B{C{radius}}).

           @return: Destination point (L{LatLon}).

           @raise ValueError: Invalid B{C{distance}}, B{C{bearing}},
                              B{C{radius}} or B{C{height}}.

           @example:

           >>> p1 = LatLon(51.4778, -0.0015)
           >>> p2 = p1.destination(7794, 300.7)
           >>> p2.toStr()  # '51.5135°N, 000.0983°W'

           @JSname: I{destinationPoint}.
        '''
        a, b = self.philam

        r, t = _angular(distance, radius), Bearing_(bearing)

        a, b = _destination2(a, b, r, t)
        h = self.height if height is None else Height(height)
        return self.classof(degrees90(a), degrees180(b), height=h)

    def distanceTo(self, other, radius=R_M, wrap=False):
        '''Compute the distance from this to an other point.

           @arg other: The other point (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance between this and the B{C{other}} point
                    (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: Invalid B{C{radius}}.

           @example:

           >>> p1 = LatLon(52.205, 0.119)
           >>> p2 = LatLon(48.857, 2.351);
           >>> d = p1.distanceTo(p2)  # 404300
        '''
        self.others(other)

        a1, b1 = self.philam
        a2, b2 = other.philam

        db, _ = unrollPI(b1, b2, wrap=wrap)
        r = vincentys_(a2, a1, db)
        return r * Radius(radius)

    def greatCircle(self, bearing):
        '''Compute the vector normal to great circle obtained by heading
           on the given initial bearing from this point.

           Direction of vector is such that initial bearing vector
           b = c × n, where n is an n-vector representing this point.

           @arg bearing: Bearing from this point (compass C{degrees360}).

           @return: Vector representing great circle (L{Vector3d}).

           @raise ValueError: Invalid B{C{bearing}}.

           @example:

           >>> p = LatLon(53.3206, -1.7297)
           >>> g = p.greatCircle(96.0)
           >>> g.toStr()  # (-0.794, 0.129, 0.594)
        '''
        a, b = self.philam

        t = Bearing_(bearing)

        sa, ca, sb, cb, st, ct = sincos2(a, b, t)

        return Vector3d(sb * ct - cb * sa * st,
                       -cb * ct - sb * sa * st,
                        ca * st)  # XXX .unit()?

    def initialBearingTo(self, other, wrap=False, raiser=False):
        '''Compute the initial bearing (forward azimuth) from this
           to an other point.

           @arg other: The other point (spherical L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).
           @kwarg raiser: Optionally, raise L{CrossError} (C{bool}),
                          use B{C{raiser}}=C{True} for behavior like
                          C{sphericalNvector.LatLon.initialBearingTo}.

           @return: Initial bearing (compass C{degrees360}).

           @raise CrossError: If this and the B{C{other}} point coincide,
                              provided B{C{raiser}} is C{True} and
                              L{crosserrors} is C{True}.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @example:

           >>> p1 = LatLon(52.205, 0.119)
           >>> p2 = LatLon(48.857, 2.351)
           >>> b = p1.initialBearingTo(p2)  # 156.2

           @JSname: I{bearingTo}.
        '''
        self.others(other)

        a1, b1 = self.philam
        a2, b2 = other.philam

        # XXX behavior like sphericalNvector.LatLon.initialBearingTo
        if raiser and crosserrors() and max(abs(a2 - a1), abs(b2 - b1)) < EPS:
            raise CrossError(_points_, self, txt=_coincident_)

        return degrees(bearing_(a1, b1, a2, b2, final=False, wrap=wrap))

    def intermediateTo(self, other, fraction, height=None, wrap=False):
        '''Locate the point at given fraction between this and an
           other point.

           @arg other: The other point (L{LatLon}).
           @arg fraction: Fraction between both points (float, between
                          0.0 for this and 1.0 for the other point).
           @kwarg height: Optional height, overriding the fractional
                          height (C{meter}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Intermediate point (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: Invalid B{C{fraction}} or B{C{height}}.

           @example:

           >>> p1 = LatLon(52.205, 0.119)
           >>> p2 = LatLon(48.857, 2.351)
           >>> p = p1.intermediateTo(p2, 0.25)  # 51.3721°N, 000.7073°E

           @JSname: I{intermediatePointTo}.
        '''
        self.others(other)

        f = Scalar(fraction, name=_fraction_)

        a1, b1 = self.philam
        a2, b2 = other.philam

        db, b2 = unrollPI(b1, b2, wrap=wrap)
        r = vincentys_(a2, a1, db)
        sr = sin(r)
        if abs(sr) > EPS:
            sa1, ca1, sa2, ca2, \
            sb1, cb1, sb2, cb2 = sincos2(a1, a2, b1, b2)

            A = sin((1 - f) * r) / sr
            B = sin(     f  * r) / sr

            x = A * ca1 * cb1 + B * ca2 * cb2
            y = A * ca1 * sb1 + B * ca2 * sb2
            z = A * sa1       + B * sa2

            a = atan2(z, hypot(x, y))
            b = atan2(y, x)

        else:  # points too close
            a = favg(a1, a2, f=f)
            b = favg(b1, b2, f=f)

        if height is None:
            h = self._havg(other, f=f)
        else:
            h = Height(height)
        return self.classof(degrees90(a), degrees180(b), height=h)

    def intersection(self, end1, other, end2, height=None, wrap=False):
        '''Locate the intersection point of two paths both defined
           by two points or a start point and bearing from North.

           @arg end1: End point of this path (L{LatLon}) or the
                      initial bearing at this point (compass
                      C{degrees360}).
           @arg other: Start point of the other path (L{LatLon}).
           @arg end2: End point of the other path (L{LatLon}) or
                      the initial bearing at the other start point
                      (compass C{degrees360}).
           @kwarg height: Optional height for intersection point,
                          overriding the mean height (C{meter}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: The intersection point (L{LatLon}).  An alternate
                    intersection point might be the L{antipode} to
                    the returned result.

           @raise IntersectionError: Intersection is ambiguous or infinite
                                     or the paths are coincident, colinear
                                     or parallel.

           @raise TypeError: If B{C{end1}}, B{C{other}} or B{C{end2}} is
                             not L{LatLon}.

           @raise ValueError: Invalid B{C{height}}.

           @example:

           >>> p = LatLon(51.8853, 0.2545)
           >>> s = LatLon(49.0034, 2.5735)
           >>> i = p.intersection(108.547, s, 32.435)  # '50.9078°N, 004.5084°E'
        '''
        return intersection(self, end1, other, end2,
                                  height=height, wrap=wrap,
                                  LatLon=self.classof)

    def intersections2(self, rad1, other, rad2, radius=R_M,
                             height=None, wrap=False):
        '''Compute the intersection points of two circles each defined
           by a center point and radius.

           @arg rad1: Radius of the this circle (C{meter} or C{radians},
                      see B{C{radius}}).
           @arg other: Center of the other circle (L{LatLon}).
           @arg rad2: Radius of the other circle (C{meter} or C{radians},
                      see B{C{radius}}).
           @kwarg radius: Mean earth radius (C{meter} or C{None} if both
                          B{C{rad1}} and B{C{rad2}} are given in C{radians}).
           @kwarg height: Optional height for the intersection points,
                          overriding the mean height (C{meter}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: 2-Tuple of the intersection points, each a L{LatLon}
                    instance.  For abutting circles, the intersection
                    points are the same instance.

           @raise IntersectionError: Concentric, antipodal, invalid or
                                     non-intersecting circles.

           @raise TypeError: If B{C{other}} is not L{LatLon}.

           @raise ValueError: Invalid B{C{rad1}}, B{C{rad2}}, B{C{radius}}
                              or B{C{height}}.
        '''
        c2 = self.others(other)

        return intersections2(self, rad1, c2, rad2, radius=radius,
                                                    height=height, wrap=wrap,
                                                    LatLon=self.classof)

    def isenclosedBy(self, points):
        '''Check whether a (convex) polygon encloses this point.

           @arg points: The polygon points (L{LatLon}[]).

           @return: C{True} if the polygon encloses this point,
                    C{False} otherwise.

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} are not L{LatLon}.

           @raise ValueError: Invalid B{C{points}}, non-convex polygon.

           @example:

           >>> b = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
           >>> p = LatLon(45,1, 1.1)
           >>> inside = p.isEnclosedBy(b)  # True
        '''
        n, points = self.points2(points, closed=True)

        n0 = self._N_vector

        if iterNumpy2(points):

            v1 = points[n-1]._N_vector
            v2 = points[n-2]._N_vector
            gc1 = v2.cross(v1)
            t0 = gc1.angleTo(n0) > PI_2
            for i in range(n):
                v2 = points[i]._N_vector
                gc = v1.cross(v2)
                v1 = v2

                ti = gc.angleTo(n0) > PI_2
                if ti != t0:
                    return False  # outside

                if gc1.angleTo(gc, vSign=n0) < 0:
                    raise _ValueError(_item_sq(points=i), points[i], txt=_not_convex_)
                gc1 = gc

        else:
            # get great-circle vector for each edge
            gc, v1 = [], points[n-1]._N_vector
            for i in range(n):
                v2 = points[i]._N_vector
                gc.append(v1.cross(v2))
                v1 = v2

            # check whether this point on same side of all
            # polygon edges (to the left or right depending
            # on anti-/clockwise polygon direction)
            t0 = gc[0].angleTo(n0) > PI_2  # True if on the right
            for i in range(1, n):
                ti = gc[i].angleTo(n0) > PI_2
                if ti != t0:  # different sides of edge i
                    return False  # outside

            # check for convex polygon (otherwise
            # the test above is not reliable)
            gc1 = gc[n-1]
            for i, gc2 in enumerate(gc):
                # angle between gc vectors, signed by direction of n0
                if gc1.angleTo(gc2, vSign=n0) < 0:
                    raise _ValueError(_item_sq(points=i), points[i], txt=_not_convex_)
                gc1 = gc2

        return True  # inside

    def isEnclosedBy(self, points):  # PYCHOK no cover
        '''DEPRECATED, use method C{isenclosedBy}.
        '''
        return self.isenclosedBy(points)

    def midpointTo(self, other, height=None, wrap=False):
        '''Find the midpoint between this and an other point.

           @arg other: The other point (L{LatLon}).
           @kwarg height: Optional height for midpoint, overriding
                          the mean height (C{meter}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Midpoint (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: Invalid B{C{height}}.

           @example:

           >>> p1 = LatLon(52.205, 0.119)
           >>> p2 = LatLon(48.857, 2.351)
           >>> m = p1.midpointTo(p2)  # '50.5363°N, 001.2746°E'
        '''
        self.others(other)

        # see <https://MathForum.org/library/drmath/view/51822.html>
        a1, b1 = self.philam
        a2, b2 = other.philam

        db, b2 = unrollPI(b1, b2, wrap=wrap)

        sa1, ca1, sa2, ca2, sdb, cdb = sincos2(a1, a2, db)

        x = ca2 * cdb + ca1
        y = ca2 * sdb

        a = atan2(sa1 + sa2, hypot(x, y))
        b = atan2(y, x) + b1

        if height is None:
            h = self._havg(other)
        else:
            h = Height(height)
        return self.classof(degrees90(a), degrees180(b), height=h)

    def nearestOn(self, point1, point2, radius=R_M, **options):
        '''Locate the point between two points closest and this point.

           Distances are approximated by function L{equirectangular_},
           subject to the supplied B{C{options}}.

           @arg point1: Start point (L{LatLon}).
           @arg point2: End point (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg options: Optional keyword arguments for function
                           L{equirectangular_}.

           @return: Closest point on the arc (L{LatLon}).

           @raise LimitError: Lat- and/or longitudinal delta exceeds
                              B{C{limit}}, see function L{equirectangular_}.

           @raise TypeError: If B{C{point1}} or B{C{point2}} is not L{LatLon}.

           @raise ValueError: Invalid B{C{radius}}.

           @see: Functions L{equirectangular_} and L{nearestOn5} and
                 method L{sphericalTrigonometry.LatLon.nearestOn3}.
        '''
        return self.nearestOn3([point1, point2], closed=False, radius=radius,
                                               **options)[0]

    def nearestOn2(self, points, closed=False, radius=R_M, **options):  # PYCHOK no cover
        '''DEPRECATED, use method L{sphericalTrigonometry.LatLon.nearestOn3}.

           @return: ... 2-Tuple C{(closest, distance)} of the closest
                    point (L{LatLon}) on the polygon and the distance
                    to that point from this point in C{meter}, same
                    units of B{C{radius}}.
        '''
        r = self.nearestOn3(points, closed=closed, radius=radius, **options)
        return tuple(r[:2])

    def nearestOn3(self, points, closed=False, radius=R_M, **options):
        '''Locate the point on a polygon closest to this point.

           Distances are approximated by function L{equirectangular_},
           subject to the supplied B{C{options}}.

           @arg points: The polygon points (L{LatLon}[]).
           @kwarg closed: Optionally, close the polygon (C{bool}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg options: Optional keyword arguments for function
                           L{equirectangular_}.

           @return: A L{NearestOn3Tuple}C{(closest, distance, angle)} of
                    the C{closest} point (L{LatLon}), the L{equirectangular_}
                    C{distance} between this and the C{closest} point in
                    C{meter}, same units as B{C{radius}}.  The C{angle}
                    from this to the C{closest} point is in compass
                    C{degrees360}, like function L{compassAngle}.

           @raise LimitError: Lat- and/or longitudinal delta exceeds
                              B{C{limit}}, see function L{equirectangular_}.

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} are not C{LatLon}.

           @raise ValueError: Invalid B{C{radius}}.

           @see: Functions L{compassAngle}, L{equirectangular_} and
                 L{nearestOn5}.
        '''
        lat, lon, d, c, h = _nearestOn5(self, points, closed=closed, **options)
        return NearestOn3Tuple(self.classof(lat, lon, height=h),
                               degrees2m(d, radius=radius), c)

    def toCartesian(self, **Cartesian_datum_kwds):  # PYCHOK Cartesian=Cartesian, datum=None
        '''Convert this point to C{Karney}-based cartesian (ECEF)
           coordinates.

           @kwarg Cartesian_datum_kwds: Optional L{Cartesian}, B{C{datum}}
                                        and other keyword arguments, ignored
                                        if B{C{Cartesian=None}}.  Use
                                        B{C{Cartesian=...}} to override
                                        this L{Cartesian} class or specify
                                        B{C{Cartesian=None}}.

           @return: The cartesian point (L{Cartesian}) or if B{C{Cartesian}}
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian_datum_kwds}}.
        '''
        kwds = _xkwds(Cartesian_datum_kwds, Cartesian=Cartesian, datum=self.datum)
        return LatLonSphericalBase.toCartesian(self, **kwds)


_T00 = LatLon(0, 0)  #: (INTERNAL) Reference instance (L{LatLon}).


def areaOf(points, radius=R_M, wrap=True):
    '''Calculate the area of a (spherical) polygon (with great circle
       arcs joining the points).

       @arg points: The polygon points (L{LatLon}[]).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: Polygon area (C{meter}, same units as B{C{radius}}, squared).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: Invalid B{C{radius}}.

       @note: The area is based on I{Karney}'s U{'Area of a spherical polygon'
              <https://OSGeo-org.1560.x6.nabble.com/
              Area-of-a-spherical-polygon-td3841625.html>}.

       @see: L{pygeodesy.areaOf}, L{sphericalNvector.areaOf} and
             L{ellipsoidalKarney.areaOf}.

       @example:

       >>> b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
       >>> areaOf(b)  # 8666058750.718977

       >>> c = LatLon(0, 0), LatLon(1, 0), LatLon(0, 1)
       >>> areaOf(c)  # 6.18e9
    '''
    n, points = _T00.points2(points, closed=True)

    # Area method due to Karney: for each edge of the polygon,
    #
    #                tan(Δλ/2) · (tan(φ1/2) + tan(φ2/2))
    #     tan(E/2) = ------------------------------------
    #                     1 + tan(φ1/2) · tan(φ2/2)
    #
    # where E is the spherical excess of the trapezium obtained by
    # extending the edge to the equator-circle vector for each edge

    def _exs(n, points):  # iterate over spherical edge excess
        a1, b1 = points[n-1].philam
        ta1 = tan_2(a1)
        for i in range(n):
            a2, b2 = points[i].philam
            db, b2 = unrollPI(b1, b2, wrap=wrap)
            ta2, tdb = map1(tan_2, a2, db)
            yield atan2(tdb * (ta1 + ta2), 1 + ta1 * ta2)
            ta1, b1 = ta2, b2

    s = fsum(_exs(n, points)) * 2

    if isPoleEnclosedBy(points):
        s = abs(s) - PI2

    return abs(s * Radius(radius)**2)


def _xb(a1, b1, end, a, b, wrap):
    # difference between the bearing to (a, b) and the given
    # bearing is negative if both are in opposite directions
    r = bearing_(a1, b1, a, b, wrap=wrap)
    return PI_2 - abs(wrapPI(r - radians(end)))


def _xdot(d, a1, b1, a, b, wrap):
    # compute dot product d . (-b + b1, a - a1)
    db, _ = unrollPI(b1, b, wrap=wrap)
    return fdot(d, db, a - a1)


def _x3d2(start, end, wrap, n, hs):
    # see <https://www.EdWilliams.org/intersect.htm> (5) ff
    a1, b1 = start.philam

    if isscalar(end):  # bearing, make a point
        a2, b2 = _destination2(a1, b1, PI_4, radians(end))
    else:  # must be a point
        _T00.others(end, name=_end_ + n)
        hs.append(end.height)
        a2, b2 = end.philam

    db, b2 = unrollPI(b1, b2, wrap=wrap)
    if max(abs(db), abs(a2 - a1)) < EPS:
        raise IntersectionError(start=start, end=end, txt='null path' + n)

    # note, in EdWilliams.org/avform.htm W is + and E is -
    b21, b12 = db * 0.5, -(b1 + b2) * 0.5

    sb21, cb21, sb12, cb12, \
    sa21,    _, sa12,    _ = sincos2(b21, b12, a1 - a2, a1 + a2)

    x = _Nvector(sa21 * sb12 * cb21 - sa12 * cb12 * sb21,
                 sa21 * cb12 * cb21 + sa12 * sb12 * sb21,
                 cos(a1) * cos(a2) * sin(db))  # ll=start
    return x.unit(), (db, (a2 - a1))  # negated d


def intersection(start1, end1, start2, end2, height=None, wrap=False,
                                             LatLon=LatLon, **LatLon_kwds):
    '''Compute the intersection point of two paths both defined
       by two points or a start point and bearing from North.

       @arg start1: Start point of the first path (L{LatLon}).
       @arg end1: End point ofthe first path (L{LatLon}) or
                  the initial bearing at the first start point
                  (compass C{degrees360}).
       @arg start2: Start point of the second path (L{LatLon}).
       @arg end2: End point of the second path (L{LatLon}) or
                  the initial bearing at the second start point
                  (compass C{degrees360}).
       @kwarg height: Optional height for the intersection point,
                      overriding the mean height (C{meter}).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).
       @kwarg LatLon: Optional class to return the intersection
                      point (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if B{C{LatLon=None}}.

       @return: The intersection point (B{C{LatLon}}) or a
                L{LatLon3Tuple}C{(lat, lon, height)} if B{C{LatLon}}
                is C{None}.  An alternate intersection point might
                be the L{antipode} to the returned result.

       @raise IntersectionError: Intersection is ambiguous or infinite
                                 or the paths are coincident, colinear
                                 or parallel.

       @raise TypeError: A B{C{start}} or B{C{end}} point not L{LatLon}.

       @raise ValueError: Invalid B{C{height}}.

       @example:

       >>> p = LatLon(51.8853, 0.2545)
       >>> s = LatLon(49.0034, 2.5735)
       >>> i = intersection(p, 108.547, s, 32.435)  # '50.9078°N, 004.5084°E'
    '''
    _T00.others(start1, name=_start1_)
    _T00.others(start2, name=_start2_)

    hs = [start1.height, start2.height]

    a1, b1 = start1.philam
    a2, b2 = start2.philam

    db, b2 = unrollPI(b1, b2, wrap=wrap)
    r12 = vincentys_(a2, a1, db)
    if abs(r12) < EPS:  # [nearly] coincident points
        a, b = favg(a1, a2), favg(b1, b2)

    # see <https://www.EdWilliams.org/avform.htm#Intersection>
    elif isscalar(end1) and isscalar(end2):  # both bearings
        sa1, ca1, sa2, ca2, sr12, cr12 = sincos2(a1, a2, r12)

        x1, x2 = (sr12 * ca1), (sr12 * ca2)
        if abs(x1) < EPS or abs(x2) < EPS:
            raise IntersectionError(start1=start1, end1=end1,
                                    start2=start2, end2=end2, txt='parallel')

        # handle domain error for equivalent longitudes,
        # see also functions asin_safe and acos_safe at
        # <https://www.EdWilliams.org/avform.htm#Math>
        t1, t2 = map1(acos1, (sa2 - sa1 * cr12) / x1,
                             (sa1 - sa2 * cr12) / x2)
        if sin(db) > 0:
            t12, t21 = t1, PI2 - t2
        else:
            t12, t21 = PI2 - t1, t2

        t13, t23 = map1(radiansPI2, end1, end2)
        x1, x2 = map1(wrapPI, t13 - t12,  # angle 2-1-3
                              t21 - t23)  # angle 1-2-3
        sx1, cx1, sx2, cx2 = sincos2(x1, x2)
        if sx1 == 0 and sx2 == 0:  # max(abs(sx1), abs(sx2)) < EPS
            raise IntersectionError(start1=start1, end1=end1,
                                    start2=start2, end2=end2, txt='infinite')
        sx3 = sx1 * sx2
#       if sx3 < 0:
#           raise IntersectionError(start1=start1, end1=end1,
#                                   start2=start2, end2=end2, txt=_ambiguous_)
        x3 = acos1(cr12 * sx3 - cx2 * cx1)
        r13 = atan2(sr12 * sx3, cx2 + cx1 * cos(x3))

        a, b = _destination2(a1, b1, r13, t13)
        # choose antipode for opposing bearings
        if _xb(a1, b1, end1, a, b, wrap) < 0 or \
           _xb(a2, b2, end2, a, b, wrap) < 0:
            a, b = antipode_(a, b)  # PYCHOK PhiLam2Tuple

    else:  # end point(s) or bearing(s)
        x1, d1 = _x3d2(start1, end1, wrap, _1_, hs)
        x2, d2 = _x3d2(start2, end2, wrap, _2_, hs)
        x = x1.cross(x2)
        if x.length < EPS:  # [nearly] colinear or parallel paths
            raise IntersectionError(start1=start1, end1=end1,
                                    start2=start2, end2=end2, txt=_colinear_)
        a, b = x.philam
        # choose intersection similar to sphericalNvector
        d1 = _xdot(d1, a1, b1, a, b, wrap)
        if d1:
            d2 = _xdot(d2, a2, b2, a, b, wrap)
            if (d2 < 0 and d1 > 0) or (d2 > 0 and d1 < 0):
                a, b = antipode_(a, b)  # PYCHOK PhiLam2Tuple

    h = fmean(hs) if height is None else Height(height)
    return _latlon3(degrees90(a), degrees180(b), h,
                    intersection, LatLon, **LatLon_kwds)


def intersections2(center1, rad1, center2, rad2, radius=R_M,
                                                 height=None, wrap=False,
                                                 LatLon=LatLon, **LatLon_kwds):
    '''Compute the intersection points of two circles each defined by
       a center point and radius.

       @arg center1: Center of the first circle (L{LatLon}).
       @arg rad1: Radius of the second circle (C{meter} or C{radians},
                  see B{C{radius}}).
       @arg center2: Center of the second circle (L{LatLon}).
       @arg rad2: Radius of the second circle (C{meter} or C{radians},
                  see B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter} or C{None} if both
                      B{C{rad1}} and B{C{rad2}} are given in C{radians}).
       @kwarg height: Optional height for the intersection points,
                      overriding the "radical height" at the "radical
                      line" between both centers (C{meter}).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).
       @kwarg LatLon: Optional class to return the intersection
                      points (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if B{C{LatLon=None}}.

       @return: 2-Tuple of the intersection points, each a B{C{LatLon}}
                instance or L{LatLon3Tuple}C{(lat, lon, height)} if
                B{C{LatLon}} is C{None}.  For abutting circles, both
                intersection points are the same instance.

       @raise IntersectionError: Concentric, antipodal, invalid or
                                 non-intersecting circles.

       @raise TypeError: If B{C{center1}} or B{C{center2}} not L{LatLon}.

       @raise ValueError: Invalid B{C{rad1}}, B{C{rad2}}, B{C{radius}} or
                          B{C{height}}.

       @note: Courtesy U{Samuel Čavoj<https://GitHub.com/mrJean1/PyGeodesy/issues/41>}.

       @see: This U{Answer<https://StackOverflow.com/questions/53324667/
             find-intersection-coordinates-of-two-circles-on-earth/53331953>}.
    '''

    c1 = _T00.others(center1, name='center1')
    c2 = _T00.others(center2, name='center2')

    try:
        return _intersect2(c1, rad1, c2, rad2, radius=radius,
                                               height=height, wrap=wrap,
                                               LatLon=LatLon, **LatLon_kwds)

    except (TypeError, ValueError) as x:
        raise IntersectionError(center1=center1, rad1=rad1,
                                center2=center2, rad2=rad2, txt=str(x))


def _intersect2(c1, rad1, c2, rad2, radius=R_M,  # in .ellipsoidalBase._intersect2
                                    height=None, wrap=False, too_d=None,
                                    LatLon=LatLon, **LatLon_kwds):
    # (INTERNAL) Intersect of two spherical circles, see L{intersections2}
    # above, separated to allow callers to embellish any exceptions

    def _dest1(bearing, h):
        a, b = _destination2(a1, b1, r1, bearing)
        return _latlon3(degrees90(a), degrees180(b), h,
                        intersections2, LatLon, **LatLon_kwds)

    r1, r2, f = _rads3(rad1, rad2, radius)
    if f:  # swapped
        c1, c2 = c2, c1  # PYCHOK swap

    a1, b1 = c1.philam
    a2, b2 = c2.philam

    db, b2 = unrollPI(b1, b2, wrap=wrap)
    d = vincentys_(a2, a1, db)  # radians
    if d < max(r1 - r2, EPS):
        raise ValueError(_near_concentric_)

    x = fsum_(r1, r2, -d)  # overlap
    if x > EPS:
        sd, cd, sr1, cr1, _, cr2 = sincos2(d, r1, r2)
        x = sd * sr1
        if abs(x) < EPS:
            raise ValueError(_invalid_)
        x = acos1((cr2 - cd * cr1) / x)  # 0 <= x <= PI

    elif x < 0:
        t = (d * radius) if too_d is None else too_d
        raise ValueError(_too_distant_fmt_ % (t,))

    if height is None:  # "radical height"
        f, _ = _radical2(d, r1, r2)  # "radical ratio"
        h = Height(favg(c1.height, c2.height, f=f))
    else:
        h = Height(height)

    b = bearing_(a1, b1, a2, b2, final=False, wrap=wrap)
    if x < _EPS_I2:  # externally ...
        t = _dest1(b, h)
    elif x > _PI_EPS_I2:  # internally ...
        t = _dest1(b + PI, h)
    else:
        return _dest1(b + x, h), _dest1(b - x, h)
    return t, t  # ... abutting circles


def isPoleEnclosedBy(points, wrap=False):  # PYCHOK no cover
    '''DEPRECATED, use function L{ispolar}.
    '''
    return ispolar(points, wrap=wrap)


def _latlon3(lat, lon, height, func, LatLon, **LatLon_kwds):
    '''(INTERNAL) Helper for L{intersection}, L{intersections2} and L{meanof}.
    '''
    if LatLon is None:
        r = LatLon3Tuple(lat, lon, height)
    else:
        kwds = _xkwds(LatLon_kwds, height=height)
        r = LatLon(lat, lon, **kwds)
    return _xnamed(r, func.__name__)


def meanOf(points, height=None, LatLon=LatLon, **LatLon_kwds):
    '''Compute the geographic mean of several points.

       @arg points: Points to be averaged (L{LatLon}[]).
       @kwarg height: Optional height at mean point, overriding
                      the mean height (C{meter}).
       @kwarg LatLon: Optional class to return the mean point
                      (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}}
                           keyword arguments, ignored if
                           B{C{LatLon=None}}.

       @return: Point at geographic mean and height (B{C{LatLon}})
                or a L{LatLon3Tuple}C{(lat, lon, height)} if
                B{C{LatLon}} is C{None}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: No B{C{points}} or invalid B{C{height}}.
    '''
    # geographic mean
    n, points = _T00.points2(points, closed=False)

    m = sumOf(points[i]._N_vector for i in range(n))
    lat, lon = m._N_vector.latlon

    if height is None:
        h = fmean(points[i].height for i in range(n))
    else:
        h = Height(height)
    return _latlon3(lat, lon, h, meanOf, LatLon, **LatLon_kwds)


def nearestOn2(point, points, **closed_radius_LatLon_options):  # PYCHOK no cover
    '''DEPRECATED, use function L{sphericalTrigonometry.nearestOn3}.

       @return: ... 2-tuple C{(closest, distance)} of the C{closest}
                point (L{LatLon}) on the polygon and the C{distance}
                between the C{closest} and the given B{C{point}}.  The
                C{closest} is a B{C{LatLon}} or a L{LatLon2Tuple}C{(lat,
                lon)} if B{C{LatLon}} is C{None} ...
    '''
    ll, d, _ = nearestOn3(point, points, **closed_radius_LatLon_options)
    if _xkwds_get(closed_radius_LatLon_options, LatLon=LatLon) is None:
        ll = LatLon2Tuple(ll.lat, ll.lon)
    return ll, d


def nearestOn3(point, points, closed=False, radius=R_M,
                              LatLon=LatLon, **options):
    '''Locate the point on a polygon closest to an other, reference point.

       Distances are approximated by function L{equirectangular_},
       subject to the supplied B{C{options}}.

       @arg point: The other, reference point (L{LatLon}).
       @arg points: The polygon points (L{LatLon}[]).
       @kwarg closed: Optionally, close the polygon (C{bool}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg LatLon: Optional class to return the closest point
                      (L{LatLon}) or C{None}.
       @kwarg options: Optional keyword arguments for function
                       L{equirectangular_}.

       @return: A L{NearestOn3Tuple}C{(closest, distance, angle)} with the
                C{closest} point as B{L{LatLon}} or L{LatLon3Tuple}C{(lat,
                lon, height)} if B{C{LatLon}} is C{None}.  The C{distance}
                is the L{equirectangular_} distance between the C{closest}
                and the given B{C{point}} in C{meter}, same units as
                B{C{radius}}.  The C{angle} from the given B{C{point}}
                to the C{closest} is in compass C{degrees360}, like function
                L{compassAngle}.  The C{height} is the (interpolated) height
                at the C{closest} point.

       @raise LimitError: Lat- and/or longitudinal delta exceeds the
                          B{C{limit}}, see function L{equirectangular_}.

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise ValueError: Invalid B{C{radius}}.

       @see: Functions L{equirectangular_} and L{nearestOn5}.
    '''
    lat, lon, d, c, h = _nearestOn5(point, points, closed=closed,
                                                   LatLon=None, **options)
    r = LatLon3Tuple(lat, lon, h) if LatLon is None else \
              LatLon(lat, lon, height=h)
    r = NearestOn3Tuple(r, degrees2m(d, radius=radius), c)
    return _xnamed(r, nearestOn3.__name__)


def perimeterOf(points, closed=False, radius=R_M, wrap=True):
    '''Compute the perimeter of a (spherical) polygon (with great circle
       arcs joining the points).

       @arg points: The polygon points (L{LatLon}[]).
       @kwarg closed: Optionally, close the polygon (C{bool}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: Polygon perimeter (C{meter}, same units as B{C{radius}}).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: Invalid B{C{radius}}.

       @note: This perimeter is based on the L{haversine} formula.

       @see: L{pygeodesy.perimeterOf}, L{sphericalNvector.perimeterOf}
             and L{ellipsoidalKarney.perimeterOf}.
    '''
    n, points = _T00.points2(points, closed=closed)

    def _rads(n, points, closed):  # angular edge lengths in radians
        i, m = _imdex2(closed, n)
        a1, b1 = points[i].philam
        for i in range(m, n):
            a2, b2 = points[i].philam
            db, b2 = unrollPI(b1, b2, wrap=wrap)
            yield vincentys_(a2, a1, db)
            a1, b1 = a2, b2

    r = fsum(_rads(n, points, closed))
    return r * Radius(radius)


__all__ += _ALL_OTHER(Cartesian, LatLon,  # classes
                      areaOf,  # functions
                      intersection, intersections2, ispolar,
                      isPoleEnclosedBy,  # DEPRECATED, use ispolar
                      meanOf,
                      nearestOn2, nearestOn3,
                      perimeterOf,
                      sumOf)  # == vector3d.sumOf

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
