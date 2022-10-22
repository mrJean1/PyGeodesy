
# -*- coding: utf-8 -*-

u'''Spherical, C{trigonometry}-based geodesy.

Trigonometric classes geodetic (lat-/longitude) L{LatLon} and
geocentric (ECEF) L{Cartesian} and functions L{areaOf}, L{intersection},
L{intersections2}, L{isPoleEnclosedBy}, L{meanOf}, L{nearestOn3} and
L{perimeterOf}, I{all spherical}.

Pure Python implementation of geodetic (lat-/longitude) methods using
spherical trigonometry, transcoded from JavaScript originals by
I{(C) Chris Veness 2011-2016} published under the same MIT Licence**, see
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import copysign0, isscalar, map1, signOf
from pygeodesy.constants import EPS, EPS1, EPS4, PI, PI2, PI_2, PI_4, R_M, \
                                isnear0, isnear1, isnon0, _0_0, _0_5, \
                                _1_0, _2_0, _90_0
from pygeodesy.datums import _ellipsoidal_datum, _mean_radius
from pygeodesy.errors import _AssertionError, CrossError, crosserrors, \
                             _ValueError, IntersectionError, _xError, \
                             _xkwds, _xkwds_get, _xkwds_pop
from pygeodesy.fmath import favg, fdot, fmean, hypot
from pygeodesy.fsums import Fsum, fsum, fsum_
from pygeodesy.formy import antipode_, bearing_, _bearingTo2, excessAbc, \
                            excessGirard, excessLHuilier, opposing_, _radical2, \
                            vincentys_
from pygeodesy.interns import _1_, _2_, _coincident_, _colinear_, _concentric_, \
                              _convex_, _end_, _infinite_, _invalid_, _LatLon_, \
                              _near_, _not_, _null_, _points_, _SPACE_, _too_
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS, _ALL_OTHER
# from pygeodesy.named import notImplemented  # from .points
from pygeodesy.namedTuples import LatLon2Tuple, LatLon3Tuple, \
                                  NearestOn3Tuple, Triangle7Tuple, \
                                  Triangle8Tuple
from pygeodesy.points import ispolar, nearestOn5 as _nearestOn5, \
                             notImplemented, Fmt as _Fmt  # XXX shadowed
from pygeodesy.props import deprecated_function, deprecated_method
from pygeodesy.sphericalBase import _angular, CartesianSphericalBase, \
                                     LatLonSphericalBase, _rads3, _trilaterate5
# from pygeodesy.streprs import Fmt as _Fmt  # from .points XXX shadowed
from pygeodesy.units import Bearing_, Height, Lam_, Phi_, Radius, \
                            Radius_, Scalar
from pygeodesy.utily import acos1, asin1, degrees90, degrees180, degrees2m, \
                            m2radians, radiansPI2, sincos2_, tan_2, unrollPI, \
                            wrap180, wrapPI
from pygeodesy.vector3d import sumOf, Vector3d

from math import asin, atan2, cos, degrees, radians, sin

__all__ = _ALL_LAZY.sphericalTrigonometry
__version__ = '22.10.12'

_parallel_ = 'parallel'
_path_     = 'path'

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

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}} is C{None},
                    an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon_and_kwds}} argument.
        '''
        kwds = _xkwds(LatLon_and_kwds, LatLon=LatLon, datum=self.datum)
        return CartesianSphericalBase.toLatLon(self, **kwds)


class LatLon(LatLonSphericalBase):
    '''New point on spherical model earth model.

       @example:

        >>> p = LatLon(52.205, 0.119)  # height=0
    '''

    def alongTrackDistanceTo(self, start, end, radius=R_M, wrap=False):
        '''Compute the (angular) distance (signed) from the start to
           the closest point on the great circle path defined by a
           start and an end point.

           That is, if a perpendicular is drawn from this point to the
           great circle path, the along-track distance is the distance
           from the start point to the point where the perpendicular
           crosses the path.

           @arg start: Start point of great circle path (L{LatLon}).
           @arg end: End point of great circle path (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None}.
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance along the great circle path (C{meter},
                    same units as B{C{radius}}) or C{radians} if
                    C{B{radius} is None}, positive if after the B{C{start}}
                    toward the B{C{end}} point of the path or negative
                    if before the B{C{start}} point.

           @raise TypeError: Invalid B{C{start}} or B{C{end}} point.

           @raise ValueError: Invalid B{C{radius}}.

           @example:

            >>> p = LatLon(53.2611, -0.7972)

            >>> s = LatLon(53.3206, -1.7297)
            >>> e = LatLon(53.1887, 0.1334)
            >>> d = p.alongTrackDistanceTo(s, e)  # 62331.58
        '''
        r, x, b = self._trackDistanceTo3(start, end, radius, wrap)
        cx = cos(x)
        return _0_0 if isnear0(cx) else \
               _r2m(copysign0(acos1(cos(r) / cx), cos(b)), radius)

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
        sa2, ca2, sdb, cdb = sincos2_(a, a1, a2, db)
        sa1 *= ca2 * ca

        x = sa1 * sdb
        y = sa1 * cdb - ca1 * sa2 * ca
        z = ca1 * sdb * ca2 * sa

        h = hypot(x, y)
        if h < EPS or abs(z) > h:  # PYCHOK no cover
            return None  # great circle doesn't reach latitude

        m = atan2(-y, x) + b1  # longitude at max latitude
        d = acos1(z / h)  # delta longitude to intersections
        return degrees180(m - d), degrees180(m + d)

    def crossTrackDistanceTo(self, start, end, radius=R_M, wrap=False):
        '''Compute the (angular) distance (signed) from this point to
           the great circle defined by a start and an end point.

           @arg start: Start point of great circle path (L{LatLon}).
           @arg end: End point of great circle path (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None}.
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance to great circle (negative if to the left
                    or positive if to the right of the path) (C{meter},
                    same units as B{C{radius}} or C{radians} if
                    B{C{radius}} is C{None}).

           @raise TypeError: The B{C{start}} or B{C{end}} point is not L{LatLon}.

           @raise ValueError: Invalid B{C{radius}}.

           @example:

            >>> p = LatLon(53.2611, -0.7972)

            >>> s = LatLon(53.3206, -1.7297)
            >>> e = LatLon(53.1887, 0.1334)
            >>> d = p.crossTrackDistanceTo(s, e)  # -307.5
        '''
        _, x, _ = self._trackDistanceTo3(start, end, radius, wrap)
        return _r2m(x, radius)

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
        '''
        a, b =  self.philam
        r, t = _angular(distance, radius), Bearing_(bearing)

        a, b = _destination2(a, b, r, t)
        h = self.height if height is None else Height(height)
        return self.classof(degrees90(a), degrees180(b), height=h)

    def distanceTo(self, other, radius=R_M, wrap=False):
        '''Compute the (angular) distance from this to an other point.

           @arg other: The other point (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None}.
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance between this and the B{C{other}} point
                    (C{meter}, same units as B{C{radius}} or
                    C{radians} if B{C{radius}} is C{None}).

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
        return _r2m(vincentys_(a2, a1, db), radius)

#   @Property_RO
#   def Ecef(self):
#       '''Get the ECEF I{class} (L{EcefVeness}), I{lazily}.
#       '''
#       return _MODS.ecef.EcefKarney

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
        sa, ca, sb, cb, st, ct = sincos2_(a, b, t)

        return Vector3d(sb * ct - cb * sa * st,
                       -cb * ct - sb * sa * st,
                        ca * st)  # XXX .unit()?

    def initialBearingTo(self, other, wrap=False, raiser=False):
        '''Compute the initial bearing (forward azimuth) from this
           to an other point.

           @arg other: The other point (spherical L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).
           @kwarg raiser: Optionally, raise L{CrossError} (C{bool}),
                          use C{B{raiser}=True} for behavior like
                          C{sphericalNvector.LatLon.initialBearingTo}.

           @return: Initial bearing (compass C{degrees360}).

           @raise CrossError: If this and the B{C{other}} point coincide,
                              provided both B{C{raiser}} is C{True} and
                              L{pygeodesy.crosserrors} is C{True}.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @example:

            >>> p1 = LatLon(52.205, 0.119)
            >>> p2 = LatLon(48.857, 2.351)
            >>> b = p1.initialBearingTo(p2)  # 156.2
        '''
        self.others(other)

        a1, b1 = self.philam
        a2, b2 = other.philam

        # XXX behavior like sphericalNvector.LatLon.initialBearingTo
        if raiser and crosserrors() and max(abs(a2 - a1), abs(b2 - b1)) < EPS:
            raise CrossError(_points_, self, txt=_coincident_)

        return degrees(bearing_(a1, b1, a2, b2, final=False, wrap=wrap))

    def intermediateTo(self, other, fraction, height=None, wrap=False):
        '''Locate the point at given fraction between (or along) this
           and an other point.

           @arg other: The other point (L{LatLon}).
           @arg fraction: Fraction between both points (C{scalar},
                          0.0 at this and 1.0 at the other point).
           @kwarg height: Optional height, overriding the intermediate
                          height (C{meter}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Intermediate point (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: Invalid B{C{fraction}} or B{C{height}}.

           @see: Methods C{midpointTo} and C{rhumbMidpointTo}.

           @example:

            >>> p1 = LatLon(52.205, 0.119)
            >>> p2 = LatLon(48.857, 2.351)
            >>> p = p1.intermediateTo(p2, 0.25)  # 51.3721°N, 000.7073°E
        '''
        f = Scalar(fraction=fraction)  # not high=_1_0
        if isnear0(f):  # PYCHOK no cover
            r = self
        elif isnear1(f) and not wrap:  # PYCHOK no cover
            r = self.others(other)
        else:
            a1, b1 = self.philam
            a2, b2 = self.others(other).philam

            db, b2 = unrollPI(b1, b2, wrap=wrap)
            r  = vincentys_(a2, a1, db)
            sr = sin(r)
            if isnon0(sr):
                sa1, ca1, sa2, ca2, \
                sb1, cb1, sb2, cb2 = sincos2_(a1, a2, b1, b2)

                t = f * r
                a = sin(r - t)  # / sr  superflous
                b = sin(    t)  # / sr  superflous

                x = a * ca1 * cb1 + b * ca2 * cb2
                y = a * ca1 * sb1 + b * ca2 * sb2
                z = a * sa1       + b * sa2

                a = atan2(z, hypot(x, y))
                b = atan2(y, x)

            else:  # PYCHOK no cover
                a = favg(a1, a2, f=f)  # coincident
                b = favg(b1, b2, f=f)

            h = self._havg(other, f=f) if height is None else Height(height)
            r = self.classof(degrees90(a), degrees180(b), height=h)
        return r

    def intersection(self, end1, other, end2, height=None, wrap=False):
        '''Compute the intersection point of two paths, each defined
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

           @raise IntersectionError: Ambiguous or infinite intersection
                                     or colinear, parallel or otherwise
                                     non-intersecting paths.

           @raise TypeError: If B{C{other}} is not L{LatLon} or B{C{end1}}
                             or B{C{end2}} not C{scalar} nor L{LatLon}.

           @raise ValueError: Invalid B{C{height}} or C{null} path.

           @example:

            >>> p = LatLon(51.8853, 0.2545)
            >>> s = LatLon(49.0034, 2.5735)
            >>> i = p.intersection(108.547, s, 32.435)  # '50.9078°N, 004.5084°E'
        '''
        try:
            s2 = self.others(other)
            return _intersect(self, end1, s2, end2, height=height, wrap=wrap,
                                                    LatLon=self.classof)
        except (TypeError, ValueError) as x:
            raise _xError(x, start1=self, end1=end1, other=other, end2=end2)

    def intersections2(self, rad1, other, rad2, radius=R_M, eps=_0_0,
                                                height=None, wrap=True):
        '''Compute the intersection points of two circles, each defined
           by a center point and radius.

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
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: 2-Tuple of the intersection points, each a L{LatLon}
                    instance.  For abutting circles, both intersection
                    points are the same instance, aka I{radical center}.

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
            raise _xError(x, center=self, rad1=rad1, other=other, rad2=rad2)

    def isenclosedBy(self, points):
        '''Check whether a (convex) polygon encloses this point.

           @arg points: The polygon points (L{LatLon}[]).

           @return: C{True} if the polygon encloses this point,
                    C{False} otherwise.

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} are not L{LatLon}.

           @raise ValueError: Invalid B{C{points}}, non-convex polygon.

           @see: Functions L{pygeodesy.isconvex}, L{pygeodesy.isenclosedBy}
                 and L{pygeodesy.ispolar} especially if the B{C{points}} may
                 enclose a pole or wrap around the earth longitudinally.

           @example:

            >>> b = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
            >>> p = LatLon(45,1, 1.1)
            >>> inside = p.isEnclosedBy(b)  # True
        '''
        Ps = self.PointsIter(points, loop=2, dedup=True)
        n0 = self._N_vector

        v2 = Ps[0]._N_vector
        v1 = Ps[1]._N_vector
        gc1 = v2.cross(v1)
        # check whether this point on same side of all
        # polygon edges (to the left or right depending
        # on anti-/clockwise polygon direction)
        t0 = gc1.angleTo(n0) > PI_2
        s0 = None
        # get great-circle vector for each edge
        for i, p in Ps.enumerate(closed=True):
            v2 = p._N_vector
            gc = v1.cross(v2)
            v1 = v2

            t = gc.angleTo(n0) > PI_2
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
                    raise _ValueError(t, p, txt=_not_(_convex_))
            gc1 = gc

        return True  # inside

    @deprecated_method
    def isEnclosedBy(self, points):  # PYCHOK no cover
        '''DEPRECATED, use method C{isenclosedBy}.'''
        return self.isenclosedBy(points)

    def midpointTo(self, other, height=None, fraction=_0_5, wrap=False):
        '''Find the midpoint between this and an other point.

           @arg other: The other point (L{LatLon}).
           @kwarg height: Optional height for midpoint, overriding
                          the mean height (C{meter}).
           @kwarg fraction: Midpoint location from this point (C{scalar}),
                            may be negative or greater than 1.0.
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Midpoint (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: Invalid B{C{height}}.

           @see: Methods C{intermediateTo} and C{rhumbMidpointTo}.

           @example:

            >>> p1 = LatLon(52.205, 0.119)
            >>> p2 = LatLon(48.857, 2.351)
            >>> m = p1.midpointTo(p2)  # '50.5363°N, 001.2746°E'
        '''
        if fraction is _0_5:
            self.others(other)

            # see <https://MathForum.org/library/drmath/view/51822.html>
            a1, b1 = self.philam
            a2, b2 = other.philam

            db, b2 = unrollPI(b1, b2, wrap=wrap)

            sa1, ca1, sa2, ca2, sdb, cdb = sincos2_(a1, a2, db)

            x = ca2 * cdb + ca1
            y = ca2 * sdb

            a = atan2(sa1 + sa2, hypot(x, y))
            b = atan2(y, x) + b1

            h = self._havg(other) if height is None else Height(height)
            r = self.classof(degrees90(a), degrees180(b), height=h)
        else:
            r = self.intermediateTo(other, fraction, height=height, wrap=wrap)
        return r

    def nearestOn(self, point1, point2, radius=R_M, **options):
        '''Locate the point between two points closest to this point.

           Distances are approximated by function L{pygeodesy.equirectangular_},
           subject to the supplied B{C{options}}.

           @arg point1: Start point (L{LatLon}).
           @arg point2: End point (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg options: Optional keyword arguments for function
                           L{pygeodesy.equirectangular_}.

           @return: Closest point on the great circle path (L{LatLon}).

           @raise LimitError: Lat- and/or longitudinal delta exceeds B{C{limit}},
                              see function L{pygeodesy.equirectangular_}.

           @raise NotImplementedError: Keyword argument C{B{within}=False}
                                       is not (yet) supported.

           @raise TypeError: Invalid B{C{point1}} or B{C{point2}}.

           @raise ValueError: Invalid B{C{radius}} or B{C{options}}.

           @see: Functions L{pygeodesy.equirectangular_} and L{pygeodesy.nearestOn5}
                 and method L{sphericalTrigonometry.LatLon.nearestOn3}.
        '''
        # remove kwarg B{C{within}} if present
        within = _xkwds_pop(options, within=True)
        if not within:
            notImplemented(self, within=within)

#       # UNTESTED - handle C{B{within}=False} and C{B{within}=True}
#       wrap = _xkwds_get(options, wrap=False)
#       a = self.alongTrackDistanceTo(point1, point2, radius=radius, wrap=wrap)
#       if abs(a) < EPS or (within and a < EPS):
#           return point1
#       d = point1.distanceTo(point2, radius=radius, wrap=wrap)
#       if isnear0(d):
#           return point1  # or point2
#       elif abs(d - a) < EPS or (a + EPS) > d:
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
                                               **options)[0]

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

    def nearestOn3(self, points, closed=False, radius=R_M, **options):
        '''Locate the point on a polygon closest to this point.

           Distances are approximated by function L{pygeodesy.equirectangular_},
           subject to the supplied B{C{options}}.

           @arg points: The polygon points (L{LatLon}[]).
           @kwarg closed: Optionally, close the polygon (C{bool}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg options: Optional keyword arguments for function
                           L{pygeodesy.equirectangular_}.

           @return: A L{NearestOn3Tuple}C{(closest, distance, angle)} of the
                    C{closest} point (L{LatLon}), the L{pygeodesy.equirectangular_}
                    C{distance} between this and the C{closest} point converted to
                    C{meter}, same units as B{C{radius}}.  The C{angle} from this
                    to the C{closest} point is in compass C{degrees360}, like
                    function L{pygeodesy.compassAngle}.

           @raise LimitError: Lat- and/or longitudinal delta exceeds B{C{limit}},
                              see function L{pygeodesy.equirectangular_}.

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} are not C{LatLon}.

           @raise ValueError: Invalid B{C{radius}} or B{C{options}}.

           @see: Functions L{pygeodesy.compassAngle}, L{pygeodesy.equirectangular_}
                 and L{pygeodesy.nearestOn5}.
        '''
        if _LatLon_ in options:
            raise _ValueError(**options)

        lat, lon, d, c, h = _nearestOn5(self, points, closed=closed, **options)
        return NearestOn3Tuple(self.classof(lat, lon, height=h),
                               degrees2m(d, radius=radius), c)

    def toCartesian(self, **Cartesian_datum_kwds):  # PYCHOK Cartesian=Cartesian, datum=None
        '''Convert this point to C{Karney}-based cartesian (ECEF)
           coordinates.

           @kwarg Cartesian_datum_kwds: Optional L{Cartesian}, B{C{datum}}
                                        and other keyword arguments, ignored
                                        if C{B{Cartesian} is None}.  Use
                                        C{B{Cartesian}=...} to override
                                        this L{Cartesian} class or specify
                                        C{B{Cartesian}=None}.

           @return: The cartesian point (L{Cartesian}) or if B{C{Cartesian}}
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian_datum_kwds}} argument.
        '''
        kwds = _xkwds(Cartesian_datum_kwds, Cartesian=Cartesian, datum=self.datum)
        return LatLonSphericalBase.toCartesian(self, **kwds)

    def _trackDistanceTo3(self, start, end, radius, wrap):
        '''(INTERNAL) Helper for .along-/crossTrackDistanceTo.
        '''
        self.others(start=start)
        self.others(end=end)

        r = Radius_(radius)
        r = start.distanceTo(self, r, wrap=wrap) / r

        b = radians(start.initialBearingTo(self, wrap=wrap))
        e = radians(start.initialBearingTo(end,  wrap=wrap))
        x = asin(sin(r) * sin(b - e))
        return r, x, (e - b)

    def triangle7(self, otherB, otherC, radius=R_M, wrap=False):
        '''Compute the angles, sides and area of a spherical triangle.

           @arg otherB: Second triangle point (C{LatLon}).
           @arg otherC: Third triangle point (C{LatLon}).
           @kwarg radius: Mean earth radius, ellipsoid or datum
                          (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                          L{Datum} or L{a_f2Tuple}) or C{None}.
           @kwarg wrap: Wrap/unroll angular distances (C{bool}).

           @return: L{Triangle7Tuple}C{(A, a, B, b, C, c, area)} or if
                    B{C{radius}} is C{None}, a L{Triangle8Tuple}C{(A,
                    a, B, b, C, c, D, E)}.

           @see: Function L{triangle7} and U{Spherical trigonometry
                 <https://WikiPedia.org/wiki/Spherical_trigonometry>}.
        '''
        self.others(otherB=otherB)
        self.others(otherC=otherC)

        r = self.philam + otherB.philam + otherC.philam
        t = triangle8_(*r, wrap=wrap)
        return self._xnamed(_t7Tuple(t, radius))

    def trilaterate5(self, distance1, point2, distance2, point3, distance3,
                           area=True, eps=EPS1, radius=R_M, wrap=False):
        '''Trilaterate three points by area overlap or perimeter intersection
           of three corresponding circles.

           @arg distance1: Distance to this point (C{meter}, same units
                           as B{C{radius}}).
           @arg point2: Second center point (C{LatLon}).
           @arg distance2: Distance to point2 (C{meter}, same units as
                           B{C{radius}}).
           @arg point3: Third center point (C{LatLon}).
           @arg distance3: Distance to point3 (C{meter}, same units as
                           B{C{radius}}).
           @kwarg area: If C{True} compute the area overlap, otherwise the
                        perimeter intersection of the circles (C{bool}).
           @kwarg eps: The required I{minimal overlap} for C{B{area}=True}
                       or the I{intersection margin} for C{B{area}=False}
                       (C{meter}, same units as B{C{radius}}).
           @kwarg radius: Mean earth radius (C{meter}, conventionally).
           @kwarg wrap: Wrap/unroll angular distances (C{bool}).

           @return: A L{Trilaterate5Tuple}C{(min, minPoint, max, maxPoint, n)}
                    with C{min} and C{max} in C{meter}, same units as B{C{eps}},
                    the corresponding trilaterated points C{minPoint} and
                    C{maxPoint} as I{spherical} C{LatLon} and C{n}, the number
                    of trilatered points found for the given B{C{eps}}.

                    If only a single trilaterated point is found, C{min I{is}
                    max}, C{minPoint I{is} maxPoint} and C{n = 1}.

                    For C{B{area}=True}, C{min} and C{max} are the smallest
                    respectively largest I{radial} overlap found.

                    For C{B{area}=False}, C{min} and C{max} represent the
                    nearest respectively farthest intersection margin.

                    If C{B{area}=True} and all 3 circles are concentric, C{n =
                    0} and C{minPoint} and C{maxPoint} are both the B{C{point#}}
                    with the smallest B{C{distance#}} C{min} and C{max} the
                    largest B{C{distance#}}.

           @raise IntersectionError: Trilateration failed for the given B{C{eps}},
                                     insufficient overlap for C{B{area}=True} or
                                     no intersection or all (near-)concentric for
                                     C{B{area}=False}.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @raise ValueError: Coincident B{C{points}} or invalid B{C{distance1}},
                              B{C{distance2}}, B{C{distance3}} or B{C{radius}}.
        '''
        return _trilaterate5(self, distance1,
                             self.others(point2=point2), distance2,
                             self.others(point3=point3), distance3,
                             area=area, radius=radius, eps=eps, wrap=wrap)


_T00 = LatLon(0, 0, name='T00')  # reference instance (L{LatLon})


def areaOf(points, radius=R_M, wrap=True):
    '''Calculate the area of a (spherical) polygon (with the points
       joined by great circle arcs).

       @arg points: The polygon points (L{LatLon}[]).
       @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                      L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple})
                      or C{None}.
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: Polygon area (C{meter} I{quared}, same units as
                B{C{radius}} or C{radians} if B{C{radius}} is C{None}).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: Invalid B{C{radius}} or semi-circular
                          polygon edge.

       @note: The area is based on I{Karney}'s U{'Area of a spherical
              polygon'<https://MathOverflow.net/questions/97711/
              the-area-of-spherical-polygons>}, 3rd Answer.

       @see: Functions L{pygeodesy.areaOf}, L{sphericalNvector.areaOf},
             L{pygeodesy.excessKarney}, L{ellipsoidalExact.areaOf} and
             L{ellipsoidalKarney.areaOf}.

       @example:

        >>> b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
        >>> areaOf(b)  # 8666058750.718977

        >>> c = LatLon(0, 0), LatLon(1, 0), LatLon(0, 1)
        >>> areaOf(c)  # 6.18e9
    '''
    Ps = _T00.PointsIter(points, loop=1)
    p1 = p2 = Ps[0]
    a1,  b1 = p1.philam
    ta1, z1 = tan_2(a1), None

    A = Fsum()  # mean phi
    E = Fsum()  # see L{pygeodesy.excessKarney_}
    # ispolar: Summation of course deltas around pole is 0° rather than normally ±360°
    # <https://blog.Element84.com/determining-if-a-spherical-polygon-contains-a-pole.html>
    # XXX duplicate of function C{points.ispolar} to avoid copying all iterated points
    D = Fsum()
    for i, p2 in Ps.enumerate(closed=True):
        a2, b2 = p2.philam
        db, b2 = unrollPI(b1, b2, wrap=wrap if i else False)
        ta2 = tan_2(a2)
        A += a2
        E += atan2(tan_2(db, points=i) * (ta1 + ta2),
                                   _1_0 + ta1 * ta2)
        ta1, b1 = ta2, b2

        if not p2.isequalTo(p1, eps=EPS):
            z, z2 = _bearingTo2(p1, p2, wrap=wrap)
            if z1 is not None:
                D += wrap180(z - z1)  # (z - z1 + 540) ...
            D += wrap180(z2 - z)  # (z2 - z + 540) % 360 - 180
            p1, z1 = p2, z2

    R = abs(E * _2_0)
    if abs(D) < _90_0:  # ispolar(points)
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


def _int3d2(start, end, wrap, _i_, Vector, hs):
    # see <https://www.EdWilliams.org/intersect.htm> (5) ff
    # and similar logic in .ellipsoidalBaseDI._intersect3
    a1, b1 = start.philam

    if isscalar(end):  # bearing, get pseudo-end point
        a2, b2 = _destination2(a1, b1, PI_4, radians(end))
    else:  # must be a point
        start.others(end, name=_end_ + _i_)
        hs.append(end.height)
        a2, b2 = end.philam

    db, b2 = unrollPI(b1, b2, wrap=wrap)
    if max(abs(db), abs(a2 - a1)) < EPS:
        raise _ValueError(_SPACE_(_path_ + _i_, _null_))
    # note, in EdWilliams.org/avform.htm W is + and E is -
    b21, b12 = db * _0_5, -(b1 + b2) * _0_5

    sb21, cb21, sb12, cb12, \
    sa21,    _, sa12,    _ = sincos2_(b21, b12, a1 - a2, a1 + a2)

    x = Vector(sa21 * sb12 * cb21 - sa12 * cb12 * sb21,
               sa21 * cb12 * cb21 + sa12 * sb12 * sb21,
               cos(a1) * cos(a2) * sin(db))  # ll=start
    return x.unit(), (db, (a2 - a1))  # negated d


def _intdot(ds, a1, b1, a, b, wrap):
    # compute dot product ds . (-b + b1, a - a1)
    db, _ = unrollPI(b1, b, wrap=wrap)
    return fdot(ds, db, a - a1)


def _intersect(start1, end1, start2, end2, height=None, wrap=False,  # in.ellipsoidalBaseDI._intersect3
                                           LatLon=None, **LatLon_kwds):
    # (INTERNAL) Intersect two (spherical) path, see L{intersection}
    # above, separated to allow callers to embellish any exceptions

    hs = [start1.height, start2.height]

    a1, b1 = start1.philam
    a2, b2 = start2.philam

    db, b2 = unrollPI(b1, b2, wrap=wrap)
    r12 = vincentys_(a2, a1, db)
    if abs(r12) < EPS:  # [nearly] coincident points
        a, b = favg(a1, a2), favg(b1, b2)

    # see <https://www.EdWilliams.org/avform.htm#Intersection>
    elif isscalar(end1) and isscalar(end2):  # both bearings
        sa1, ca1, sa2, ca2, sr12, cr12 = sincos2_(a1, a2, r12)

        x1, x2 = (sr12 * ca1), (sr12 * ca2)
        if isnear0(x1) or isnear0(x2):
            raise IntersectionError(_parallel_)
        # handle domain error for equivalent longitudes,
        # see also functions asin_safe and acos_safe at
        # <https://www.EdWilliams.org/avform.htm#Math>
        t1, t2 = acos1((sa2 - sa1 * cr12) / x1), \
                 acos1((sa1 - sa2 * cr12) / x2)
        if sin(db) > 0:
            t12, t21 = t1, PI2 - t2
        else:
            t12, t21 = PI2 - t1, t2
        t13, t23 = radiansPI2(end1), radiansPI2(end2)
        sx1, cx1, sx2, cx2 = sincos2_(wrapPI(t13 - t12),  # angle 2-1-3
                                      wrapPI(t21 - t23))  # angle 1-2-3)
        if isnear0(sx1) and isnear0(sx2):
            raise IntersectionError(_infinite_)
        sx3 = sx1 * sx2
# XXX   if sx3 < 0:
# XXX       raise ValueError(_ambiguous_)
        x3  = acos1(cr12 * sx3 - cx2 * cx1)
        r13 = atan2(sr12 * sx3, cx2 + cx1 * cos(x3))

        a, b = _destination2(a1, b1, r13, t13)
        # like .ellipsoidalBaseDI,_intersect3, if this intersection
        # is "before" the first point, use the antipodal intersection
        if opposing_(t13, bearing_(a1, b1, a, b, wrap=wrap)):
            a, b = antipode_(a, b)  # PYCHOK PhiLam2Tuple

    else:  # end point(s) or bearing(s)
        _N_vector_ = _MODS.nvectorBase._N_vector_

        x1, d1 = _int3d2(start1, end1, wrap, _1_, _N_vector_, hs)
        x2, d2 = _int3d2(start2, end2, wrap, _2_, _N_vector_, hs)
        x = x1.cross(x2)
        if x.length < EPS:  # [nearly] colinear or parallel paths
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
                                             LatLon=LatLon, **LatLon_kwds):
    '''Compute the intersection point of two paths, each defined
       by two points or a start point and bearing from North.

       @arg start1: Start point of the first path (L{LatLon}).
       @arg end1: End point of the first path (L{LatLon}) or
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
                           arguments, ignored if C{B{LatLon} is None}.

       @return: The intersection point as a (B{C{LatLon}}) or if
                C{B{LatLon} is None} a L{LatLon3Tuple}C{(lat, lon,
                height)}.  An alternate intersection point might
                be the L{antipode} to the returned result.

       @raise IntersectionError: Ambiguous or infinite intersection
                                 or colinear, parallel or otherwise
                                 non-intersecting paths.

       @raise TypeError: A B{C{start}} or B{C{end}} point not L{LatLon}.

       @raise ValueError: Invalid B{C{height}} or C{null} path.

       @example:

        >>> p = LatLon(51.8853, 0.2545)
        >>> s = LatLon(49.0034, 2.5735)
        >>> i = intersection(p, 108.547, s, 32.435)  # '50.9078°N, 004.5084°E'
    '''
    s1 = _T00.others(start1=start1)
    s2 = _T00.others(start2=start2)
    try:
        return _intersect(s1, end1, s2, end2, height=height, wrap=wrap,
                                              LatLon=LatLon, **LatLon_kwds)
    except (TypeError, ValueError) as x:
        raise _xError(x, start1=start1, end1=end1, start2=start2, end2=end2)


def intersections2(center1, rad1, center2, rad2, radius=R_M, eps=_0_0,
                                                 height=None, wrap=True,
                                                 LatLon=LatLon, **LatLon_kwds):
    '''Compute the intersection points of two circles each defined
       by a center point and a radius.

       @arg center1: Center of the first circle (L{LatLon}).
       @arg rad1: Radius of the first circle (C{meter} or C{radians},
                  see B{C{radius}}).
       @arg center2: Center of the second circle (L{LatLon}).
       @arg rad2: Radius of the second circle (C{meter} or C{radians},
                  see B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter} or C{None} if B{C{rad1}},
                      B{C{rad2}} and B{C{eps}} are given in C{radians}).
       @kwarg eps: Required overlap (C{meter} or C{radians}, see
                   B{C{radius}}).
       @kwarg height: Optional height for the intersection points (C{meter},
                      conventionally) or C{None} for the I{"radical height"}
                      at the I{radical line} between both centers.
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).
       @kwarg LatLon: Optional class to return the intersection
                      points (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon} is None}.

       @return: 2-Tuple of the intersection points, each a B{C{LatLon}}
                instance or if C{B{LatLon} is None} a L{LatLon3Tuple}C{(lat,
                lon, height)}.  For abutting circles, both intersection
                points are the same instance, aka I{radical center}.

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
                                                LatLon=LatLon, **LatLon_kwds)
    except (TypeError, ValueError) as x:
        raise _xError(x, center1=center1, rad1=rad1, center2=center2, rad2=rad2)


def _intersects2(c1, rad1, c2, rad2, radius=R_M, eps=_0_0,  # in .ellipsoidalBaseDI._intersects2
                                     height=None, too_d=None, wrap=True,
                                     LatLon=LatLon, **LatLon_kwds):
    # (INTERNAL) Intersect two spherical circles, see L{intersections2}
    # above, separated to allow callers to embellish any exceptions

    def _dest3(bearing, h):
        a, b = _destination2(a1, b1, r1, bearing)
        return _LL3Tuple(degrees90(a), degrees180(b), h,
                         intersections2, LatLon, LatLon_kwds)

    r1, r2, f = _rads3(rad1, rad2, radius)
    if f:  # swapped
        c1, c2 = c2, c1  # PYCHOK swap

    a1, b1 = c1.philam
    a2, b2 = c2.philam

    db, b2 = unrollPI(b1, b2, wrap=wrap)
    d = vincentys_(a2, a1, db)  # radians
    if d < max(r1 - r2, EPS):
        raise IntersectionError(_near_(_concentric_))  # XXX ConcentricError?

    r = eps if radius is None else (m2radians(
        eps, radius=radius) if eps else _0_0)
    if r < _0_0:
        raise _ValueError(eps=r)

    x = fsum_(r1, r2, -d)  # overlap
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


def _LL3Tuple(lat, lon, height, func, LatLon, LatLon_kwds):
    '''(INTERNAL) Helper for L{intersection}, L{intersections2} and L{meanOf}.
    '''
    n = func.__name__
    if LatLon is None:
        r = LatLon3Tuple(lat, lon, height, name=n)
    else:
        kwds = _xkwds(LatLon_kwds, height=height, name=n)
        r = LatLon(lat, lon, **kwds)
    return r


def meanOf(points, height=None, LatLon=LatLon, **LatLon_kwds):
    '''Compute the geographic mean of several points.

       @arg points: Points to be averaged (L{LatLon}[]).
       @kwarg height: Optional height at mean point, overriding the mean
                      height (C{meter}).
       @kwarg LatLon: Optional class to return the mean point (L{LatLon})
                      or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon} is None}.

       @return: The geographic mean and height (B{C{LatLon}}) or a
                L{LatLon3Tuple}C{(lat, lon, height)} if B{C{LatLon}}
                is C{None}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: No B{C{points}} or invalid B{C{height}}.
    '''
    # geographic mean
    m = sumOf(p._N_vector for p in _T00.PointsIter(points).iterate(closed=False))
    lat, lon = m._N_vector.latlon

    if height is None:
        h = fmean(p.height for p in _T00.PointsIter(points).iterate(closed=False))
    else:  # PYCHOK no cover
        h = Height(height)
    return _LL3Tuple(lat, lon, h, meanOf, LatLon, LatLon_kwds)


@deprecated_function
def nearestOn2(point, points, **closed_radius_LatLon_options):  # PYCHOK no cover
    '''DEPRECATED, use function L{sphericalTrigonometry.nearestOn3}.

       @return: ... 2-tuple C{(closest, distance)} of the C{closest}
                point (L{LatLon}) on the polygon and the C{distance}
                between the C{closest} and the given B{C{point}}.  The
                C{closest} is a B{C{LatLon}} or a L{LatLon2Tuple}C{(lat,
                lon)} if B{C{LatLon}} is C{None} ...
    '''
    ll, d, _ = nearestOn3(point, points, **closed_radius_LatLon_options)  # PYCHOK 3-tuple
    if _xkwds_get(closed_radius_LatLon_options, LatLon=LatLon) is None:
        ll = LatLon2Tuple(ll.lat, ll.lon)
    return ll, d


def nearestOn3(point, points, closed=False, radius=R_M,
                              LatLon=LatLon, **options):
    '''Locate the point on a path or polygon closest to a reference point.

       Distances are I{approximated} using function L{pygeodesy.equirectangular_},
       subject to the supplied B{C{options}}.

       @arg point: The reference point (L{LatLon}).
       @arg points: The path or polygon points (L{LatLon}[]).
       @kwarg closed: Optionally, close the polygon (C{bool}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg LatLon: Optional class to return the closest point (L{LatLon})
                      or C{None}.
       @kwarg options: Optional keyword arguments for function
                       L{pygeodesy.equirectangular_}.

       @return: A L{NearestOn3Tuple}C{(closest, distance, angle)} with the
                C{closest} point as B{C{LatLon}} or L{LatLon3Tuple}C{(lat,
                lon, height)} if B{C{LatLon}} is C{None}.  The C{distance}
                is the L{pygeodesy.equirectangular_} distance between the
                C{closest} and the given B{C{point}} converted to C{meter},
                same units as B{C{radius}}.  The C{angle} from the given
                B{C{point}} to the C{closest} is in compass C{degrees360},
                like function L{pygeodesy.compassAngle}.  The C{height} is
                the (interpolated) height at the C{closest} point.

       @raise LimitError: Lat- and/or longitudinal delta exceeds the B{C{limit}},
                          see function L{pygeodesy.equirectangular_}.

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise ValueError: Invalid B{C{radius}}.

       @see: Functions L{pygeodesy.equirectangular_} and L{pygeodesy.nearestOn5}.
    '''
    lat, lon, d, c, h = _nearestOn5(point, points, closed=closed,
                                                   LatLon=None, **options)
    d = degrees2m(d, radius=radius)
    n = nearestOn3.__name__
    r = LatLon3Tuple(lat, lon, h, name=n) if LatLon is None else \
              LatLon(lat, lon, height=h, name=n)
    return NearestOn3Tuple(r, d, c, name=n)


def perimeterOf(points, closed=False, radius=R_M, wrap=True):
    '''Compute the perimeter of a (spherical) polygon (with great circle
       arcs joining the points).

       @arg points: The polygon points (L{LatLon}[]).
       @kwarg closed: Optionally, close the polygon (C{bool}).
       @kwarg radius: Mean earth radius (C{meter}) or C{None}.
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: Polygon perimeter (C{meter}, same units as B{C{radius}}
                or C{radians} if B{C{radius}} is C{None}).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @raise ValueError: Invalid B{C{radius}}.

       @note: Distances are based on function L{pygeodesy.vincentys_}.

       @see: Functions L{pygeodesy.perimeterOf}, L{sphericalNvector.perimeterOf}
             and L{ellipsoidalKarney.perimeterOf}.
    '''
    def _rads(points, closed, wrap):  # angular edge lengths in radians
        Ps = _T00.PointsIter(points, loop=1)
        a1, b1 = Ps[0].philam
        for i, p in Ps.enumerate(closed=closed):
            a2, b2 = p.philam
            db, b2 = unrollPI(b1, b2, wrap=wrap if i else False)
            yield vincentys_(a2, a1, db)
            a1, b1 = a2, b2

    r = fsum(_rads(points, closed, wrap), floats=True)
    return _r2m(r, radius)


def _r2m(r, radius):
    '''(INTERNAL) Angular distance in C{radians} to C{meter}.
    '''
    if radius is not None:  # not in (None, _0_0)
        r *= R_M if radius is R_M else Radius(radius)
    return r


def triangle7(latA, lonA, latB, lonB, latC, lonC, radius=R_M,
                                                  excess=excessAbc,
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
       @kwarg excess: I{Spherical excess} callable (L{excessAbc},
                      L{excessGirard} or L{excessLHuilier}).
       @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

       @return: A L{Triangle7Tuple}C{(A, a, B, b, C, c, area)} with
                spherical angles C{A}, C{B} and C{C}, angular sides
                C{a}, C{b} and C{c} all in C{degrees} and C{area}
                in I{square} C{meter} or units of B{C{radius}}
                I{squared} or if B{C{radius}} is C{None}, a
                L{Triangle8Tuple}C{(A, a, B, b, C, c, D, E)} all in
                C{radians} with I{spherical excess} C{E} as the
                C{unit area} in C{radians}.
    '''
    t = triangle8_(Phi_(latA=latA), Lam_(lonA=lonA),
                   Phi_(latB=latB), Lam_(lonB=lonB),
                   Phi_(latC=latC), Lam_(lonC=lonC),
                   excess=excess, wrap=wrap)
    return _t7Tuple(t, radius)


def triangle8_(phiA, lamA, phiB, lamB, phiC, lamC, excess=excessAbc,
                                                     wrap=False):
    '''Compute the angles, sides, I{spherical deficit} and I{spherical
       excess} of a (spherical) triangle.

       @arg phiA: First corner latitude (C{radians}).
       @arg lamA: First corner longitude (C{radians}).
       @arg phiB: Second corner latitude (C{radians}).
       @arg lamB: Second corner longitude (C{radians}).
       @arg phiC: Third corner latitude (C{radians}).
       @arg lamC: Third corner longitude (C{radians}).
       @kwarg excess: I{Spherical excess} callable (L{excessAbc},
                      L{excessGirard} or L{excessLHuilier}).
       @kwarg wrap: Wrap and L{pygeodesy.unrollPI} longitudes (C{bool}).

       @return: A L{Triangle8Tuple}C{(A, a, B, b, C, c, D, E)} with
                spherical angles C{A}, C{B} and C{C}, angular sides
                C{a}, C{b} and C{c}, I{spherical deficit} C{D} and
                I{spherical excess} C{E}, all in C{radians}.
    '''
    def _a_r(w, phiA, lamA, phiB, lamB, phiC, lamC):
        d, _ = unrollPI(lamB, lamC, wrap=w)
        a = vincentys_(phiC, phiB, d)
        return a, (phiB, lamB, phiC, lamC, phiA, lamA)

    def _A_r(a, sa, ca, sb, cb, sc, cc):
        s = sb * sc
        A = acos1((ca - cb * cc) / s) if isnon0(s) else a
        return A, (sb, cb, sc, cc, sa, ca)  # rotate sincos2's

    # notation: side C{a} is oposite to corner C{A}, etc.
    a, r = _a_r(wrap, phiA, lamA, phiB, lamB, phiC, lamC)
    b, r = _a_r(wrap, *r)
    c, _ = _a_r(wrap, *r)

    A, r = _A_r(a, *sincos2_(a, b, c))
    B, r = _A_r(b, *r)
    C, _ = _A_r(c, *r)

    D = fsum_(PI2, -a, -b, -c, floats=True)  # deficit aka defect
    E = excessGirard(A, B, C)   if excess in (excessGirard, True) else (
        excessLHuilier(a, b, c) if excess in (excessLHuilier, False) else
        excessAbc(*max((A, b, c), (B, c, a), (C, a, b))))

    return Triangle8Tuple(A, a, B, b, C, c, D, E)


def _t7Tuple(t, radius):
    '''(INTERNAL) Convert a L{Triangle8Tuple} to L{Triangle7Tuple}.
    '''
    if radius:  # not in (None, _0_0)
        r = radius if isscalar(radius) else \
            _ellipsoidal_datum(radius).ellipsoid.Rmean
        A, B, C = map1(degrees, t.A, t.B, t.C)
        t = Triangle7Tuple(A, (r * t.a),
                           B, (r * t.b),
                           C, (r * t.c), t.E * r**2)
    return t


__all__ += _ALL_OTHER(Cartesian, LatLon,  # classes
                      areaOf,  # functions
                      intersection, intersections2, ispolar,
                      isPoleEnclosedBy,  # DEPRECATED, use ispolar
                      meanOf,
                      nearestOn2, nearestOn3,
                      perimeterOf,
                      sumOf,  # == vector3d.sumOf
                      triangle7, triangle8_)

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
