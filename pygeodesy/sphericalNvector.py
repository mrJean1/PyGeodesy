
# -*- coding: utf-8 -*-

u'''Spherical, C{N-vector}-based geodesy.

N-vector-based classes geodetic (lat-/longitude) L{LatLon}, geocentric
(ECEF) L{Cartesian} and L{Nvector} and functions L{areaOf}, L{intersection},
L{meanOf}, L{nearestOn3}, L{perimeterOf}, L{sumOf}, L{triangulate} and
L{trilaterate}, I{all spherical}.

Pure Python implementation of n-vector-based spherical geodetic (lat-/longitude)
methods, transcoded from JavaScript originals by I{(C) Chris Veness 2011-2016},
published under the same MIT Licence**.  See U{Vector-based geodesy
<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>} and
U{Module latlon-nvector-spherical
<https://www.Movable-Type.co.UK/scripts/geodesy/docs/module-latlon-nvector-spherical.html>}.

Tools for working with points and paths on (a spherical model of) the
earth’s surface using using n-vectors rather than the more common
spherical trigonometry.  N-vectors make many calculations much simpler,
and easier to follow, compared with the trigonometric equivalents.

Based on Kenneth Gade’s U{‘Non-singular Horizontal Position Representation’
<https://www.NavLab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf>},
The Journal of Navigation (2010), vol 63, nr 3, pp 395-417.

Note that the formulations below take x => 0°N,0°E, y => 0°N,90°E and
z => 90°N while Gade uses x => 90°N, y => 0°N,90°E, z => 0°N,0°E.

Also note that on a spherical earth model, an n-vector is equivalent
to a normalised version of an (ECEF) cartesian coordinate.
'''
# make sure int/int division yields float quosient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isscalar, _xinstanceof
from pygeodesy.datums import Datums
from pygeodesy.errors import _xkwds
from pygeodesy.fmath import fmean, fsum
from pygeodesy.interns import EPS, EPS0, PI, PI2, PI_2, R_M, _end_, \
                             _Nv00_, _other_, _point_, _points_, \
                             _pole_, _0_0, _0_5, _1_0
from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER
from pygeodesy.named import notImplemented
from pygeodesy.namedTuples import NearestOn3Tuple
from pygeodesy.nvectorBase import NvectorBase, NorthPole, LatLonNvectorBase, \
                                  sumOf as _sumOf, _triangulate, _trilaterate
from pygeodesy.points import ispolar  # PYCHOK exported
from pygeodesy.props import deprecated_function, deprecated_method
from pygeodesy.sphericalBase import _angular, CartesianSphericalBase, \
                                     LatLonSphericalBase
from pygeodesy.units import Bearing, Bearing_, Height, Radius, Scalar
from pygeodesy.utily import degrees360, sincos2, sincos2_, sincos2d

from math import atan2

__all__ = _ALL_LAZY.sphericalNvector
__version__ = '21.11.30'

_paths_ = 'paths'


class Cartesian(CartesianSphericalBase):
    '''Extended to convert geocentric, L{Cartesian} points to
       L{Nvector} and n-vector-based, spherical L{LatLon}.
    '''

    def toLatLon(self, **LatLon_and_kwds):  # PYCHOK LatLon=LatLon
        '''Convert this cartesian to an C{Nvector}-based geodetic point.

           @kwarg LatLon_and_kwds: Optional L{LatLon} and L{LatLon} keyword
                                   arguments, like C{datum}.  Use C{B{LatLon}=...}
                                   to override this L{LatLon} class or specify
                                   C{B{LatLon}=None}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}} is set
                    to C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon_and_kwds}} argument.
        '''
        kwds = _xkwds(LatLon_and_kwds, LatLon=LatLon, datum=self.datum)
        return CartesianSphericalBase.toLatLon(self, **kwds)

    def toNvector(self, **Nvector_and_kwds):  # PYCHOK Datums.WGS84
        '''Convert this cartesian to L{Nvector} components, I{including height}.

           @kwarg Nvector_and_kwds: Optional L{Nvector} and L{Nvector} keyword
                                    arguments, like C{datum}.  Use C{B{Nvector}=...}
                                    to override this L{Nvector} class or specify
                                    C{B{Nvector}=None}.

           @return: The C{n-vector} components (L{Nvector}) or if B{C{Nvector}}
                    is set to C{None}, a L{Vector4Tuple}C{(x, y, z, h)}

           @raise TypeError: Invalid B{C{Nvector_and_kwds}} argument.
        '''
        # ll = CartesianBase.toLatLon(self, LatLon=LatLon,
        #                                    datum=datum or self.datum)
        # kwds = _xkwds(kwds, Nvector=Nvector)
        # return ll.toNvector(**kwds)
        kwds = _xkwds(Nvector_and_kwds, Nvector=Nvector, datum=self.datum)
        return CartesianSphericalBase.toNvector(self, **kwds)


class LatLon(LatLonNvectorBase, LatLonSphericalBase):
    '''New n-vector based point on a spherical earth model.

       Tools for working with points and paths on (a spherical
       model of) the earth's surface using vector-based methods.

       @example:

        >>> from sphericalNvector import LatLon
        >>> p = LatLon(52.205, 0.119)
    '''
    _Nv = None  # cached_toNvector L{Nvector})

    def _update(self, updated, *attrs):  # PYCHOK args
        '''(INTERNAL) Zap cached attributes if updated.
        '''
        if updated:  # reset caches
            LatLonNvectorBase._update(self, updated, _Nv=self._Nv)  # special case
            LatLonSphericalBase._update(self, updated, *attrs)

    def alongTrackDistanceTo(self, start, end, radius=R_M):
        '''Compute the (signed) distance from the start to the closest
           point on the great circle path defined by a start and an
           end point.

           That is, if a perpendicular is drawn from this point to the
           great circle path, the along-track distance is the distance
           from the start point to the point where the perpendicular
           crosses the path.

           @arg start: Start point of great circle path (L{LatLon}).
           @arg end: End point of great circle path (L{LatLon}) or
                     initial bearing from start point (compass
                     C{degrees360}).
           @kwarg radius: Mean earth radius (C{meter}).

           @return: Distance along the great circle path (positive if
                    after the start toward the end point of the path
                    or negative if before the start point).

           @raise TypeError: If B{C{start}} or B{C{end}} point is not L{LatLon}.

           @raise Valuerror: Some points coincide.

           @example:

            >>> p = LatLon(53.2611, -0.7972)

            >>> s = LatLon(53.3206, -1.7297)
            >>> e = LatLon(53.1887, 0.1334)
            >>> d = p.alongTrackDistanceTo(s, e)  # 62331.58
        '''
        self.others(start=start)
        gc, _, _ = self._gc3(start, end, _end_)

        p = self.toNvector()
        a = gc.cross(p).cross(gc)  # along-track point gc × p × gc
        return start.toNvector().angleTo(a, vSign=gc) * radius

    @deprecated_method
    def bearingTo(self, other, **unused):  # PYCHOK no cover
        '''DEPRECATED, use method L{initialBearingTo}.
        '''
        return self.initialBearingTo(other)

    def crossTrackDistanceTo(self, start, end, radius=R_M):
        '''Compute the (signed) distance from this point to great circle
           defined by a start and end point.

           @arg start: Start point of great circle path (L{LatLon}).
           @arg end: End point of great circle path (L{LatLon}) or
                     initial bearing from start point (compass
                     C{degrees360}).
           @kwarg radius: Mean earth radius (C{meter}).

           @return: Distance to great circle (negative if to the
                    left or positive if to the right of the path).

           @raise TypeError: If B{C{start}} or B{C{end}} point is not L{LatLon}.

           @raise Valuerror: Some points coincide.

           @example:

            >>> p = LatLon(53.2611, -0.7972)

            >>> s = LatLon(53.3206, -1.7297)
            >>> d = p.crossTrackDistanceTo(s, 96)  # -305.7

            >>> e = LatLon(53.1887, 0.1334)
            >>> d = p.crossTrackDistanceTo(s, e)  # -307.5
        '''
        self.others(start=start)
        gc, _, _ = self._gc3(start, end, _end_)

        p = self.toNvector()
        return (gc.angleTo(p) - PI_2) * radius

    def destination(self, distance, bearing, radius=R_M, height=None):
        '''Locate the destination from this point after having
           travelled the given distance on the given bearing.

           @arg distance: Distance travelled (C{meter}, same units
                          as B{C{radius}}).
           @arg bearing: Bearing from this point (compass C{degrees360}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg height: Optional height at destination, overriding the
                          default height (C{meter}, same units as B{C{radius}}).

           @return: Destination point (L{LatLon}).

           @raise Valuerror: Polar coincidence or invalid B{C{distance}},
                             B{C{bearing}}, B{C{radius}} or B{C{height}}.

           @example:

            >>> p = LatLon(51.4778, -0.0015)
            >>> q = p.destination(7794, 300.7)
            >>> q.toStr()  # 51.513546°N, 000.098345°W

           @JSname: I{destinationPoint}.
        '''
        a = _angular(distance, radius)
        sa, ca, sb, cb = sincos2_(a, Bearing_(bearing))

        p = self.toNvector()
        e = NorthPole.cross(p, raiser=_pole_).unit()  # east vector at p
        n = p.cross(e)  # north vector at p
        q = n.times(cb).plus(e.times(sb))  # direction vector @ p
        n = p.times(ca).plus(q.times(sa))
        return n.toLatLon(height=height, LatLon=self.classof)  # Nvector(n.x, n.y, n.z).toLatLon(...)

    def distanceTo(self, other, radius=R_M, wrap=False):
        '''Compute the distance from this to an other point.

           @arg other: The other point (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None}.
           @kwarg wrap: Wrap/unroll the angular distance (C{bool}).

           @return: Distance between this and the B{C{other}} point
                    (C{meter}, same units as B{C{radius}} or C{radians}
                    if B{C{radius}} is C{None}).

           @raise TypeError: Invalid B{C{other}} point.

           @example:

            >>> p = LatLon(52.205, 0.119)
            >>> q = LatLon(48.857, 2.351);
            >>> d = p.distanceTo(q)  # 404.3 km
        '''
        self.others(other)

        r = abs(self.toNvector().angleTo(other.toNvector(), wrap=wrap))
        return r if radius is None else (Radius(radius) * r)

#   @Property_RO
#   def Ecef(self):
#       '''Get the ECEF I{class} (L{EcefVeness}), I{lazily}.
#       '''
#       from pygeodesy.ecef import EcefKarney
#       return EcefKarney

    def _gc3(self, start, end, namend, raiser=_points_):
        '''(INTERNAL) Return great circle, start and end Nvectors.
        '''
        s = start.toNvector()
        if isscalar(end):  # bearing
            gc = s.greatCircle(end)
            e = None
        else:
            self.others(end, name=namend)
            e = end.toNvector()
            gc = s.cross(e, raiser=raiser)  # XXX .unit()?
        return gc, s, e

    def greatCircle(self, bearing):
        '''Compute the vector normal to great circle obtained by
           heading on the given bearing from this point.

           Direction of vector is such that initial bearing vector
           b = c × n, where n is an n-vector representing this point.

           @arg bearing: Bearing from this point (compass C{degrees360}).

           @return: N-vector representing the great circle (L{Nvector}).
        '''
        a, b = self.philam
        t = Bearing_(bearing)

        sa, ca, sb, cb, st, ct = sincos2_(a, b, t)
        return Nvector(sb * ct - sa * cb * st,
                      -cb * ct - sa * sb * st,
                       ca * st, name=self.name)  # XXX .unit()

    def greatCircleTo(self, other):
        '''Compute the vector normal to great circle obtained by
           heading from this to an other point or on a given bearing.

           Direction of vector is such that initial bearing vector
           b = c × n, where n is an n-vector representing this point.

           @arg other: The other point (L{LatLon}) or the bearing from
                       this point (compass C{degrees360}).

           @return: N-vector representing the great circle (L{Nvector}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise Valuerror: Points coincide.

           @example:

            >>> p = LatLon(53.3206, -1.7297)
            >>> gc = p.greatCircle(96.0)
            >>> gc.toStr()  # (-0.79408, 0.12856, 0.59406)

            >>> q = LatLon(53.1887, 0.1334)
            >>> g = p.greatCircleTo(q)
            >>> g.toStr()  # (-0.79408, 0.12859, 0.59406)
        '''
        gc, _, _ = self._gc3(self, other, _other_)
        return gc.unit()

    def initialBearingTo(self, other, **unused):
        '''Compute the initial bearing (forward azimuth) from this
           to an other point.

           @arg other: The other point (L{LatLon}).
           @arg unused: Optional keyword argument B{C{wrap}} ignored.

           @return: Initial bearing (compass C{degrees360}).

           @raise Crosserror: This point coincides with the B{C{other}}
                              point or the C{NorthPole}, provided
                              L{pygeodesy.crosserrors} is C{True}.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @example:

            >>> p1 = LatLon(52.205, 0.119)
            >>> p2 = LatLon(48.857, 2.351)
            >>> b = p1.initialBearingTo(p2)  # 156.2

           @JSname: I{bearingTo}.
        '''
        self.others(other)
        # see <https://MathForum.org/library/drmath/view/55417.html>
        n = self.toNvector()
#       gc1 = self.greatCircleTo(other)
        gc1 = n.cross(other.toNvector(), raiser=_points_)  # .unit()
#       gc2 = self.greatCircleTo(NorthPole)
        gc2 = n.cross(NorthPole, raiser=_pole_)  # .unit()
        return degrees360(gc1.angleTo(gc2, vSign=n))

    def intermediateChordTo(self, other, fraction, height=None):
        '''Locate the point projected from the point at given fraction
           on a straight line (chord) between this and an other point.

           @arg other: The other point (L{LatLon}).
           @arg fraction: Fraction between both points (float, between
                          0.0 for this and 1.0 for the other point).
           @kwarg height: Optional height at the intermediate point,
                          overriding the fractional height (C{meter}).

           @return: Intermediate point (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @example:

            >>> p = LatLon(52.205, 0.119)
            >>> q = LatLon(48.857, 2.351)
            >>> i = p.intermediateChordTo(q, 0.25)  # 51.3723°N, 000.7072°E

           @JSname: I{intermediatePointOnChordTo}, I{intermediatePointDirectlyTo}.
        '''
        self.others(other)

        f = Scalar(fraction=fraction)
        i = other.toNvector().times(f).plus(
             self.toNvector().times(1 - f))
#       i = other.toNvector() * f + \
#            self.toNvector() * (1 - f))

        h = self._havg(other, f=f) if height is None else Height(height)
        return i.toLatLon(height=h, LatLon=self.classof)  # Nvector(i.x, i.y, i.z).toLatLon(...)

    def intermediateTo(self, other, fraction, height=None, **unused):  # wrap=False
        '''Locate the point at a given fraction between this and an
           other point.

           @arg other: The other point (L{LatLon}).
           @arg fraction: Fraction between both points (C{float}, between
                          0.0 for this and 1.0 for the other point).
           @kwarg height: Optional height at the intermediate point,
                          overriding the fractional height (C{meter}).

           @return: Intermediate point (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise Valuerror: Points coincide or invalid B{C{height}}.

           @example:

            >>> p = LatLon(52.205, 0.119)
            >>> q = LatLon(48.857, 2.351)
            >>> i = p.intermediateTo(q, 0.25)  # 51.3721°N, 000.7074°E

           @JSname: I{intermediatePointTo}.
        '''
        q = self.others(other).toNvector()
        p = self.toNvector()
        f = Scalar(fraction=fraction)

        x = p.cross(q, raiser=_points_)
        d = x.unit().cross(p)  # unit(p × q) × p
        # angular distance α, tan(α) = |p × q| / p ⋅ q
        s, c = sincos2(atan2(x.length, p.dot(q)) * f)  # interpolated
        i = p.times(c).plus(d.times(s))  # p * cosα + d * sinα

        h = self._havg(other, f=f) if height is None else Height(height)
        return i.toLatLon(height=h, LatLon=self.classof)  # Nvector(i.x, i.y, i.z).toLatLon(...)

    def intersection(self, end1, start2, end2, height=None):
        '''Locate the intersection point of two paths each defined
           by two points or a start point and bearing from North.

           @arg end1: End point of the first path (L{LatLon}) or the
                      initial bearing at this point (compass C{degrees360}).
           @arg start2: Start point of the second path (L{LatLon}).
           @arg end2: End point of the second path (L{LatLon}) or the
                      initial bearing at the second point (compass
                      C{degrees}).
           @kwarg height: Optional height at the intersection point,
                          overriding the mean height (C{meter}).

           @return: The intersection point (L{LatLon}) or C{None}
                    if no unique intersection exists.

           @raise TypeError: If B{C{start2}}, B{C{end1}} or B{C{end2}} point
                             is not L{LatLon}.

           @raise ValueError: Intersection is ambiguous or infinite or
                              the paths are parallel, coincident or null.

           @example:

            >>> s = LatLon(51.8853, 0.2545)
            >>> e = LatLon(49.0034, 2.5735)
            >>> i = s.intersection(108.55, e, 32.44)  # 50.9076°N, 004.5086°E
        '''
        return intersection(self, end1, start2, end2,
                            height=height, LatLon=self.classof)

    def isenclosedBy(self, points):
        '''Check whether this point is enclosed by a (convex) polygon.

           @arg points: The polygon points (L{LatLon}[]).

           @return: C{True} if the polygon encloses this point,
                    C{False} otherwise.

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} are not L{LatLon}.

           @example:

            >>> b = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
            >>> p = LatLon(45.1, 1.1)
            >>> inside = p.isenclosedBy(b)  # True

           @JSname: I{enclosedBy}.
        '''
        # sum subtended angles of each edge (using n0, the
        # normal vector to this point for sign of α)
        def _subtangles(Ps, n0):
            vs1 = n0.minus(Ps[0].toNvector())
            for p in Ps.iterate(closed=True):
                vs2 = n0.minus(p.toNvector())
                yield vs1.angleTo(vs2, vSign=n0)  # PYCHOK false
                vs1 = vs2

        # Note, this method uses angle summation test: on a plane,
        # angles for an enclosed point will sum to 360°, angles for
        # an exterior point will sum to 0°.  On a sphere, enclosed
        # point angles will sum to less than 360° (due to spherical
        # excess), exterior point angles will be small but non-zero.
        s = fsum(_subtangles(self.PointsIter(points, loop=1),
                             self.toNvector()))  # normal vector
        # XXX are winding number optimisations equally applicable to
        # spherical surface?
        return abs(s) > PI

    @deprecated_method
    def isEnclosedBy(self, points):  # PYCHOK no cover
        '''DEPRECATED, use method C{isenclosedBy}.'''
        return self.isenclosedBy(points)

    def iswithin(self, point1, point2):
        '''Check whether this point is between two other points.

           If this point is not on the great circle arc defined by
           both points, return whether it is within the area bound
           by perpendiculars to the great circle at each point (in
           the same hemispere).

           @arg point1: Start point of the arc (L{LatLon}).
           @arg point2: End point of the arc (L{LatLon}).

           @return: C{True} if this point is within the arc,
                    C{False} otherwise.

           @raise TypeError: If B{C{point1}} or B{C{point2}} is not L{LatLon}.

           @JSname: I{isBetween}.
        '''
        n0 = self.toNvector()
        n1 = self.others(point1=point1).toNvector()
        n2 = self.others(point2=point2).toNvector()

        # corner case, null arc
        if n1.isequalTo(n2):
            return n0.isequalTo(n1) or n0.isequalTo(n2)  # PYCHOK returns

        if n0.dot(n1) < 0 or n0.dot(n2) < 0:  # different hemisphere
            return False  # PYCHOK returns

        # get vectors representing d0=p0->p1 and d2=p2->p1 and the
        # dot product d0⋅d2 tells us if p0 is on the p2 side of p1 or
        # on the other side (similarly for d0=p0->p2 and d1=p1->p2
        # and dot product d0⋅d1 and p0 on the p1 side of p2 or not)
        return n0.minus(n1).dot(n2.minus(n1)) >= 0 and \
               n0.minus(n2).dot(n1.minus(n2)) >= 0

    @deprecated_method
    def isWithin(self, point1, point2):  # PYCHOK no cover
        '''DEPRECATED, use method C{iswithin}.'''
        return self.iswithin(point1, point2)

    def midpointTo(self, other, height=None):
        '''Find the midpoint between this and an other point.

           @arg other: The other point (L{LatLon}).
           @kwarg height: Optional height at the midpoint, overriding
                          the mean height (C{meter}).

           @return: Midpoint (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @example:

            >>> p1 = LatLon(52.205, 0.119)
            >>> p2 = LatLon(48.857, 2.351)
            >>> m = p1.midpointTo(p2)  # '50.5363°N, 001.2746°E'
        '''
        self.others(other)

        m = self.toNvector().plus(other.toNvector())
        h = self._havg(other) if height is None else height
        return m.toLatLon(height=h, LatLon=self.classof)

    def nearestOn(self, point1, point2, height=None, within=True, wrap=False):
        '''Locate the point on the great circle arc between two
           points closest to this point.

           @arg point1: Start point of the arc (L{LatLon}).
           @arg point2: End point of the arc (L{LatLon}).
           @kwarg height: Optional height, overriding the mean height
                          for the point within the arc (C{meter}), or
                          C{None} to interpolate the height.
           @kwarg within: If C{True} return the closest point between
                          both given points, otherwise the closest
                          point elsewhere on the arc (C{bool}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Closest point on the arc (L{LatLon}).

           @raise NotImplementedError: Keyword argument C{B{wrap}=True}
                                       not supported.

           @raise TypeError: Invalid B{C{point1}} or B{C{point2}}.

           @example:

            >>> s1 = LatLon(51.0, 1.0)
            >>> s2 = LatLon(51.0, 2.0)

            >>> s = LatLon(51.0, 1.9)
            >>> p = s.nearestOn(s1, s2)  # 51.0004°N, 001.9000°E

            >>> d = p.distanceTo(s)  # 42.71 m

            >>> s = LatLon(51.0, 2.1)
            >>> p = s.nearestOn(s1, s2)  # 51.0000°N, 002.0000°E

           @JSname: I{nearestPointOnSegment}.
        '''
        if wrap:  # wrap=True throws C{NotImplementedError} always.
            notImplemented(self, wrap=wrap)

        if self.iswithin(point1, point2) and not point1.isequalTo(point2, EPS):
            # closer to arc than to its endpoints,
            # find the closest point on the arc
            gc1 = point1.toNvector().cross(point2.toNvector())
            gc2 = self.toNvector().cross(gc1)
            n   = gc1.cross(gc2)

        elif within:  # for backward compatibility
            return point1 if self.distanceTo(point1) < self.distanceTo(point2) else point2

        else:  # handle beyond arc extent by .vector3d.nearestOn
            n1 = point1.toNvector()
            n2 = point2.toNvector()
            n  = self.toNvector().nearestOn(n1, n2, within=False)
            if n is n1:
                return point1
            elif n is n2:
                return point2

        p = n.toLatLon(height=height or 0, LatLon=self.classof)
        if height in (None, False):  # interpolate height within extent
            d =  point1.distanceTo(point2)
            f = (point1.distanceTo(p) / d) if d > EPS0 else _0_5
            p.height = point1._havg(point2, f=max(_0_0, min(f, _1_0)))
        return p

    # @deprecated_method
    def nearestOn2(self, points, **closed_radius_height):  # PYCHOK no cover
        '''DEPRECATED, use method L{sphericalNvector.LatLon.nearestOn3}.

           @return: ... 2-Tuple C{(closest, distance)} of the C{closest}
                    point (L{LatLon}) on the polygon and the C{distance}
                    to that point from this point ...
        '''
        r = self.nearestOn3(points, **closed_radius_height)
        return r.closest, r.distance

    def nearestOn3(self, points, closed=False, radius=R_M, height=None):
        '''Locate the point on a path or polygon (with great circle
           arcs joining consecutive points) closest to this point.

           The closest point is either on within the extent of any great
           circle arc or the nearest of the arc's end points.

           @arg points: The path or polygon points (L{LatLon}[]).
           @kwarg closed: Optionally, close the polygon (C{bool}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None}.
           @kwarg height: Optional height, overriding the mean height
                          for a point within the arc (C{meter}).

           @return: A L{NearestOn3Tuple}C{(closest, distance, angle)} of
                    the C{closest} point (L{LatLon}), the C{distance}
                    between this and the C{closest} point in C{meter},
                    same units as B{C{radius}} or in C{radians} if
                    B{C{radius}} is C{None} and the C{angle} from this
                    to the C{closest} point in compass C{degrees360}.

           @raise TypeError: Some B{C{points}} are not C{LatLon}.

           @raise ValueError: No B{C{points}}.
        '''
        Ps = self.PointsIter(points, loop=1)

        R = self.distanceTo
        N = self.nearestOn

        c = p1 = Ps[0]
        r = R(c, radius=None)  # radians
        for p2 in Ps.iterate(closed=closed):
            p = N(p1, p2, height=height)
            d = R(p, radius=None)  # radians
            if d < r:
                c, r = p, d
            p1 = p2
        d = r if radius is None else (Radius(radius) * r)
        return NearestOn3Tuple(c, d, degrees360(r))

    def toCartesian(self, **Cartesian_and_kwds):  # PYCHOK Cartesian=Cartesian, datum=None
        '''Convert this point to C{Nvector}-based cartesian (ECEF) coordinates.

           @kwarg Cartesian_and_kwds: Optional L{Cartesian} and L{Cartesian} keyword
                                      arguments, like C{datum}.  Use C{B{Cartesian}=...}
                                      to override this L{Cartesian} class or specify
                                      C{B{Cartesian}=None}.

           @return: The cartesian point (L{Cartesian}) or if B{C{Cartesian}} is
                    set to C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian_and_kwds}} argument.
        '''
        kwds = _xkwds(Cartesian_and_kwds, Cartesian=Cartesian, datum=self.datum)
        return LatLonSphericalBase.toCartesian(self, **kwds)

    def toNvector(self, **Nvector_and_kwds):  # PYCHOK signature
        '''Convert this point to L{Nvector} components, I{including height}.

           @kwarg Nvector_and_kwds: Optional L{Nvector} and L{Nvector} keyword
                                    arguments.  Use C{B{Nvector}=...} to override
                                    this L{Nvector} class or specify
                                    C{B{Nvector}=None}.

           @return: The C{n-vector} components (L{Nvector}) or if B{C{Nvector}} is
                    set to C{None}, a L{Vector4Tuple}C{(x, y, z, h)} if B{C{Nvector}}.

           @raise TypeError: Invalid B{C{Nvector_and_kwds}} argument.

           @example:

            >>> p = LatLon(45, 45)
            >>> n = p.toNvector()
            >>> n.toStr()  # [0.50, 0.50, 0.70710]

           @JSname: I{toVector}.
        '''
        kwds = _xkwds(Nvector_and_kwds, Nvector=Nvector)
        return LatLonNvectorBase.toNvector(self, **kwds)


class Nvector(NvectorBase):
    '''An n-vector is a position representation using a (unit) vector
       normal to the earth's surface.  Unlike lat-/longitude points,
       n-vectors have no singularities or discontinuities.

       For many applications, n-vectors are more convenient to work
       with than other position representations like lat-/longitude,
       earth-centred earth-fixed (ECEF) vectors, UTM coordinates, etc.

       On a spherical model earth, an n-vector is equivalent to an
       earth-centred earth-fixed (ECEF) vector.

       Note commonality with L{ellipsoidalNvector.Nvector}.
    '''
    _datum = Datums.Sphere  # default datum (L{Datum})

    def toCartesian(self, **Cartesian_and_kwds):  # PYCHOK Cartesian=Cartesian
        '''Convert this n-vector to C{Nvector}-based cartesian
           (ECEF) coordinates.

           @kwarg Cartesian_and_kwds: Optional L{Cartesian} and L{Cartesian} keyword
                                      arguments, like C{h}.  Use C{B{Cartesian}=...}
                                      to override this L{Cartesian} class or specify
                                      C{B{Cartesian}=None}.

           @return: The cartesian point (L{Cartesian}) or if B{C{Cartesian}} is
                    set to C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian_and_kwds}} argument.
        '''
        kwds = _xkwds(Cartesian_and_kwds, h=self.h, Cartesian=Cartesian)
        return NvectorBase.toCartesian(self, **kwds)  # class or .classof

    def toLatLon(self, **LatLon_and_kwds):  # PYCHOK height=None, LatLon=LatLon
        '''Convert this n-vector to an C{Nvector}-based geodetic point.

           @kwarg LatLon_and_kwds: Optional L{LatLon} and L{LatLon} keyword
                                   arguments, like C{height}.  Use C{B{LatLon}=...}
                                   to override this L{LatLon} class or specify
                                   C{B{LatLon}=None}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}} is set
                    to C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon_and_kwds}} argument.

           @raise ValueError: Invalid B{C{height}}.
        '''
        kwds = _xkwds(LatLon_and_kwds, height=self.h, LatLon=LatLon)
        return NvectorBase.toLatLon(self, **kwds)  # class or .classof

    def greatCircle(self, bearing):
        '''Compute the n-vector normal to great circle obtained by
           heading on given compass bearing from this point as its
           n-vector.

           Direction of vector is such that initial bearing vector
           b = c × p.

           @arg bearing: Initial compass bearing (C{degrees}).

           @return: N-vector representing great circle (L{Nvector}).

           @raise Valuerror: Polar coincidence.

           @example:

            >>> n = LatLon(53.3206, -1.7297).toNvector()
            >>> gc = n.greatCircle(96.0)  # [-0.794, 0.129, 0.594]
        '''
        s, c = sincos2d(Bearing(bearing))

        e = NorthPole.cross(self, raiser=_pole_)  # easting
        n = self.cross(e, raiser=_point_)  # northing

        e = e.times(c / e.length)
        n = n.times(s / n.length)
        return n.minus(e)


_Nvll = LatLon(_0_0, _0_0, name=_Nv00_)  # reference instance (L{LatLon})


def areaOf(points, radius=R_M):
    '''Calculate the area of a (spherical) polygon (with great circle
       arcs joining consecutive points).

       @arg points: The polygon points (L{LatLon}[]).
       @kwarg radius: Mean earth radius (C{meter}) or C{None}.

       @return: Polygon area (C{meter} I{squared} , same units as
                B{C{radius}}, or C{radians} if B{C{radius}} is C{None}).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @see: Functions L{pygeodesy.areaOf}, L{sphericalTrigonometry.areaOf}
             and L{ellipsoidalKarney.areaOf}.

       @example:

        >>> b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
        >>> areaOf(b)  # 8666058750.718977
    '''
    def _interangles(Ps):
        # use vector to 1st point as plane normal for sign of α
        n0 = Ps[0].toNvector()

        v2 = Ps[0]._N_vector  # XXX v2 == no?
        v1 = Ps[1]._N_vector
        gc = v2.cross(v1)
        for p in Ps.iterate(closed=True):
            v2 = p._N_vector
            gc1 = v1.cross(v2)
            v1 = v2
            yield gc.angleTo(gc1, vSign=n0)
            gc = gc1

    # sum interior angles: depending on whether polygon is cw or ccw,
    # angle between edges is π−α or π+α, where α is angle between
    # great-circle vectors; so sum α, then take n·π − |Σα| (cannot
    # use Σ(π−|α|) as concave polygons would fail)
    s = fsum(_interangles(_Nvll.PointsIter(points, loop=2)))
    # using Girard’s theorem: A = [Σθᵢ − (n−2)·π]·R²
    # (PI2 - abs(s) == (n*PI - abs(s)) - (n-2)*PI)
    r = abs(PI2 - abs(s))
    return r if radius is None else (r * Radius(radius)**2)


def intersection(start1, end1, start2, end2,
                 height=None, LatLon=LatLon, **LatLon_kwds):
    '''Locate the intersection of two paths each defined by two
       points or by a start point and an initial bearing.

       @arg start1: Start point of the first path (L{LatLon}).
       @arg end1: End point of the first path (L{LatLon}) or the
                  initial bearing at the first start point
                  (compass C{degrees360}).
       @arg start2: Start point of the second path (L{LatLon}).
       @arg end2: End point of the second path (L{LatLon}) or the
                  initial bearing at the second start point
                  (compass C{degrees360}).
       @kwarg height: Optional height at the intersection point,
                      overriding the mean height (C{meter}).
       @kwarg LatLon: Optional class to return the intersection
                      point (L{LatLon}).
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon} is None}.

       @return: The intersection point (B{C{LatLon}}) or 3-tuple
                (C{degrees90}, C{degrees180}, height) if B{C{LatLon}}
                is C{None} or C{None} if no unique intersection
                exists.

       @raise TypeError: If B{C{start*}} or B{C{end*}} is not L{LatLon}.

       @raise ValueError: Intersection is ambiguous or infinite or
                          the paths are parallel, coincident or null.

       @example:

        >>> p = LatLon(51.8853, 0.2545)
        >>> q = LatLon(49.0034, 2.5735)
        >>> i = intersection(p, 108.55, q, 32.44)  # 50.9076°N, 004.5086°E
    '''
    _Nvll.others(start1=start1)
    _Nvll.others(start2=start2)

    # If gc1 and gc2 are great circles through start and end points
    # (or defined by start point and bearing), then the candidate
    # intersections are simply gc1 × gc2 and gc2 × gc1.  Most of the
    # work is deciding the correct intersection point to select!  If
    # bearing is given, that determines the intersection, but if both
    # paths are defined by start/end points, take closer intersection.
    gc1, s1, e1 = _Nvll._gc3(start1, end1, 'end1')
    gc2, s2, e2 = _Nvll._gc3(start2, end2, 'end2')

    hs = start1.height, start2.height
    # there are two (antipodal) candidate intersection
    # points ... we have to choose the one to return
    i1 = gc1.cross(gc2, raiser=_paths_)
    # postpone computing i2 until needed
    # i2 = gc2.cross(gc1, raiser=_paths_)

    # selection of intersection point depends on how
    # paths are defined (by bearings or endpoints)
    if e1 and e2:  # endpoint+endpoint
        d = sumOf((s1, s2, e1, e2)).dot(i1)
        hs += end1.height, end2.height
    elif e1 and not e2:  # endpoint+bearing
        # gc2 x v2 . i1 +ve means v2 bearing points to i1
        d = gc2.cross(s2).dot(i1)
        hs += end1.height,
    elif e2 and not e1:  # bearing+endpoint
        # gc1 x v1 . i1 +ve means v1 bearing points to i1
        d = gc1.cross(s1).dot(i1)
        hs += end2.height,
    else:  # bearing+bearing
        # if gc x v . i1 is +ve, initial bearing is
        # towards i1, otherwise towards antipodal i2
        d1 = gc1.cross(s1).dot(i1)  # +ve means p1 bearing points to i1
        d2 = gc2.cross(s2).dot(i1)  # +ve means p2 bearing points to i1
        if d1 > 0 and d2 > 0:
            d = 1  # both point to i1
        elif d1 < 0 and d2 < 0:
            d = -1  # both point to i2
        else:  # d1, d2 opposite signs
            # intersection is at further-away intersection point,
            # take opposite intersection from mid- point of v1
            # and v2 [is this always true?] XXX changed to always
            # get intersection p1 bearing points to, aka being
            # located "after" p1 along the bearing at p1, like
            # function .sphericalTrigonometry._intersect and
            # .ellipsoidalBaseDI._intersect3
            d = d1  # neg(s1.plus(s2).dot(i1))

    i = i1 if d > 0 else gc2.cross(gc1, raiser=_paths_)

    h = fmean(hs) if height is None else height
    kwds = _xkwds(LatLon_kwds, height=h, LatLon=LatLon)
    return i.toLatLon(**kwds)  # Nvector(i.x, i.y, i.z).toLatLon(...)


def meanOf(points, height=None, LatLon=LatLon, **LatLon_kwds):
    '''Compute the geographic mean of the supplied points.

       @arg points: Array of points to be averaged (L{LatLon}[]).
       @kwarg height: Optional height, overriding the mean height
                      (C{meter}).
       @kwarg LatLon: Optional class to return the mean point
                      (L{LatLon}).
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon} is None}.

       @return: Point at geographic mean and mean height (B{C{LatLon}}).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not C{LatLon}.
    '''
    Ps = _Nvll.PointsIter(points)
    # geographic mean
    m = sumOf(p._N_vector for p in Ps.iterate(closed=False))
    kwds = _xkwds(LatLon_kwds, height=height, LatLon=LatLon,
                               name=meanOf.__name__)
    return m.toLatLon(**kwds)


@deprecated_function
def nearestOn2(point, points, **closed_radius_height):  # PYCHOK no cover
    '''DEPRECATED, use method L{sphericalNvector.nearestOn3}.

       @return: ... 2-Tuple C{(closest, distance)} of the C{closest}
                point (L{LatLon}) on the polygon and the C{distance}
                between the C{closest} and the given B{C{point}} ...
    '''
    r = nearestOn3(point, points, **closed_radius_height)
    return r.closest, r.distance


def nearestOn3(point, points, closed=False, radius=R_M, height=None):
    '''Locate the point on a polygon (with great circle arcs
       joining consecutive points) closest to an other point.

       If the given point is within the extent of any great circle
       arc, the closest point is on that arc.  Otherwise, the
       closest is the nearest of the arc's end points.

       @arg point: The other, reference point (L{LatLon}).
       @arg points: The polygon points (L{LatLon}[]).
       @kwarg closed: Optionally, close the polygon (C{bool}).
       @kwarg radius: Mean earth radius (C{meter}) or C{None}.
       @kwarg height: Optional height, overriding the mean height
                      for a point within the arc (C{meter}).

       @return: A L{NearestOn3Tuple}C{(closest, distance, angle)} of
                the C{closest} point (L{LatLon}) on the polygon, the
                C{distance} and the C{angle} between the C{closest}
                and the given B{C{point}}.  The C{distance} is in
                C{meter}, same units as B{C{radius}} or in C{radians}
                if B{C{radius}} is C{None}, the C{angle} is in compass
                C{degrees360}.

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} or B{C{point}} not C{LatLon}.
    '''
    _xinstanceof(LatLon, point=point)

    return point.nearestOn3(points, closed=closed, radius=radius, height=height)


def perimeterOf(points, closed=False, radius=R_M):
    '''Compute the perimeter of a (spherical) polygon (with great circle
       arcs joining consecutive points).

       @arg points: The polygon points (L{LatLon}[]).
       @kwarg closed: Optionally, close the polygon (C{bool}).
       @kwarg radius: Mean earth radius (C{meter}) or C{None}.

       @return: Polygon perimeter (C{meter}, same units as B{C{radius}}
                or C{radians} if B{C{radius}} is C{None}).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @see: Functions L{pygeodesy.perimeterOf}, L{sphericalTrigonometry.perimeterOf}
             and L{ellipsoidalKarney.perimeterOf}.
    '''
    def _rads(Ps, closed):  # angular edge lengths in radians
        v1 = Ps[0]._N_vector
        for p in Ps.iterate(closed=closed):
            v2 = p._N_vector
            yield v1.angleTo(v2)
            v1 = v2

    r = fsum(_rads(_Nvll.PointsIter(points, loop=1), closed))
    return r if radius is None else (Radius(radius) * r)


def sumOf(nvectors, Vector=Nvector, h=None, **Vector_kwds):
    '''Return the vectorial sum of two or more n-vectors.

       @arg nvectors: Vectors to be added (L{Nvector}[]).
       @kwarg Vector: Optional class for the vectorial sum (L{Nvector}).
       @kwarg h: Optional height, overriding the mean height (C{meter}).
       @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword arguments.

       @return: Vectorial sum (B{C{Vector}}).

       @raise VectorError: No B{C{nvectors}}.
    '''
    return _sumOf(nvectors, Vector=Vector, h=h, **Vector_kwds)


def triangulate(point1, bearing1, point2, bearing2,
                height=None, LatLon=LatLon, **LatLon_kwds):
    '''Locate a point given two known points and the initial bearings
       from those points.

       @arg point1: First reference point (L{LatLon}).
       @arg bearing1: Bearing at the first point (compass C{degrees360}).
       @arg point2: Second reference point (L{LatLon}).
       @arg bearing2: Bearing at the second point (compass C{degrees360}).
       @kwarg height: Optional height at the triangulated point, overriding
                      the mean height (C{meter}).
       @kwarg LatLon: Optional class to return the triangulated point
                      (L{LatLon}).
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon} is None}.

       @return: Triangulated point (B{C{LatLon}}).

       @raise TypeError: If B{C{point1}} or B{C{point2}} is not L{LatLon}.

       @raise Valuerror: Points coincide.

       @example:

        >>> p = LatLon("47°18.228'N","002°34.326'W")  # Basse Castouillet
        >>> q = LatLon("47°18.664'N","002°31.717'W")  # Basse Hergo
        >>> t = triangulate(p, 7, q, 295)  # 47.323667°N, 002.568501°W'
    '''
    return _triangulate(_Nvll.others(point1=point1), bearing1,
                        _Nvll.others(point2=point2), bearing2,
                         height=height, LatLon=LatLon, **LatLon_kwds)


def trilaterate(point1, distance1, point2, distance2, point3, distance3,  # PYCHOK args
                                   radius=R_M, height=None, useZ=False,
                                   LatLon=LatLon, **LatLon_kwds):
    '''Locate a point at given distances from three other points.

       @arg point1: First point (L{LatLon}).
       @arg distance1: Distance to the first point (C{meter}, same units
                       as B{C{radius}}).
       @arg point2: Second point (L{LatLon}).
       @arg distance2: Distance to the second point (C{meter}, same units
                       as B{C{radius}}).
       @arg point3: Third point (L{LatLon}).
       @arg distance3: Distance to the third point (C{meter}, same units
                       as B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg height: Optional height at the trilaterated point, overriding
                      the IDW height (C{meter}, same units as B{C{radius}}).
       @kwarg useZ: Include Z component iff non-NaN, non-zero (C{bool}).
       @kwarg LatLon: Optional class to return the trilaterated
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword arguments,
                           ignored if C{B{LatLon} is None}.

       @return: Trilaterated point (B{C{LatLon}}).

       @raise IntersectionError: No intersection, trilateration failed.

       @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{point3}}.

       @raise ValueError: Coincident B{C{points}} or invalid B{C{distance1}},
                          B{C{distance2}}, B{C{distance3}} or B{C{radius}}.

       @see: U{Trilateration<https://WikiPedia.org/wiki/Trilateration>}.
    '''
    return _trilaterate(_Nvll.others(point1=point1), distance1,
                        _Nvll.others(point2=point2), distance2,
                        _Nvll.others(point3=point3), distance3,
                         radius=radius, height=height, useZ=useZ,
                         LatLon=LatLon, **LatLon_kwds)


__all__ += _ALL_OTHER(Cartesian, LatLon, Nvector,  # classes
                      areaOf,  # functions
                      intersection, ispolar,
                      meanOf,
                      nearestOn2, nearestOn3,
                      perimeterOf,
                      sumOf,
                      triangulate, trilaterate)

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
