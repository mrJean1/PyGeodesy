
# -*- coding: utf-8 -*-

u'''N-vector-based classes geodetic (lat-/longitude) L{LatLon}, geocentric
(ECEF) L{Cartesian} and L{Nvector} and functions L{areaOf}, L{intersection},
L{meanOf}, L{nearestOn2}, L{triangulate} and L{trilaterate}, I{all spherical}.

Pure Python implementation of n-vector-based spherical geodetic (lat-/longitude)
methods, transcribed from JavaScript originals by I{(C) Chris Veness 2011-2016},
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

@newfield example: Example, Examples
'''

from pygeodesy.basics import EPS, EPS_2, PI, PI2, PI_2, R_M, \
                             isscalar, map1, _xinstanceof, _xkwds
from pygeodesy.datum import Datums
from pygeodesy.ecef import EcefKarney
from pygeodesy.errors import _ValueError
from pygeodesy.fmath import fidw, fmean, fsum, fsum_
from pygeodesy.interns import _1_, _2_, _bearing_, _coincident_, \
                              _distance_, _end_, _fraction_, \
                              _other_, _point_, _points_, _pole_, \
                              _start_, _start1_, _start2_
from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER
from pygeodesy.named import NearestOn3Tuple
from pygeodesy.nvectorBase import NvectorBase, NorthPole, LatLonNvectorBase, \
                                  sumOf as _sumOf
from pygeodesy.points import _imdex2, ispolar  # PYCHOK exported
from pygeodesy.sphericalBase import _angular, CartesianSphericalBase, \
                                     LatLonSphericalBase
from pygeodesy.streprs import unstr
from pygeodesy.units import Bearing, Bearing_, Height, Radius, Radius_, Scalar
from pygeodesy.utily import degrees360, iterNumpy2, sincos2, sincos2d

from math import atan2, fabs, sqrt

__all__ = _ALL_LAZY.sphericalNvector
__version__ = '20.07.24'

_paths_ = 'paths'


class Cartesian(CartesianSphericalBase):
    '''Extended to convert geocentric, L{Cartesian} points to
       L{Nvector} and n-vector-based, spherical L{LatLon}.
    '''

    def toLatLon(self, **LatLon_datum_kwds):  # PYCHOK LatLon=LatLon
        '''Convert this cartesian to an C{Nvector}-based geodetic point.

           @kwarg LatLon_datum_kwds: Optional L{LatLon}, B{C{datum}} and
                  other keyword arguments, ignored if B{C{LatLon=None}}.
                  Use B{C{LatLon=...}} to override the L{LatLon} class
                  or specify B{C{LatLon=None}}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}} is
                    C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon_datum_kwds}}.
        '''
        kwds = _xkwds(LatLon_datum_kwds, LatLon=LatLon, datum=self.datum)
        return CartesianSphericalBase.toLatLon(self, **kwds)

    def toNvector(self, **Nvector_datum_kwds):  # PYCHOK Datums.WGS84
        '''Convert this cartesian to L{Nvector} components,
           I{including height}.

           @kwarg Nvector_datum_kwds: Optional L{Nvector}, B{C{datum}} and
                  other keyword arguments, ignored if B{C{Nvector=None}}.
                  Use B{C{Nvector=...}} to override the L{Nvector} class
                  or specify B{C{Nvector=None}}.

           @return: The C{n-vector} components (L{Nvector}) or a
                    L{Vector4Tuple}C{(x, y, z, h)} if B{C{Nvector}}
                    is C{None}.

           @raise TypeError: Invalid B{C{Nvector_datum_kwds}}.
        '''
        # ll = CartesianBase.toLatLon(self, LatLon=LatLon,
        #                                    datum=datum or self.datum)
        # kwds = _xkwds(kwds, Nvector=Nvector)
        # return ll.toNvector(**kwds)
        kwds = _xkwds(Nvector_datum_kwds, Nvector=Nvector, datum=self.datum)
        return CartesianSphericalBase.toNvector(self, **kwds)


class LatLon(LatLonNvectorBase, LatLonSphericalBase):
    '''New n-vector based point on a spherical earth model.

       Tools for working with points and paths on (a spherical
       model of) the earth's surface using vector-based methods.

       @example:

       >>> from sphericalNvector import LatLon
       >>> p = LatLon(52.205, 0.119)
    '''
    _Nv = None  #: (INTERNAL) cache _toNvector L{Nvector}).

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
        self.others(start, name=_start_)
        gc, _, _ = self._gc3(start, end, _end_)

        p = self.toNvector()
        a = gc.cross(p).cross(gc)  # along-track point gc × p × gc
        return start.toNvector().angleTo(a, vSign=gc) * radius

    def bearingTo(self, other, **unused):  # PYCHOK no cover
        '''DEPRECATED, use method C{initialBearingTo}.
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
        self.others(start, name=_start_)
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

           @raise Valuerror: Polar coincidence ior invalid B{C{distance}},
                             B{C{bearing}}, B{C{radius}} or B{C{height}}.

           @example:

           >>> p = LatLon(51.4778, -0.0015)
           >>> q = p.destination(7794, 300.7)
           >>> q.toStr()  # 51.513546°N, 000.098345°W

           @JSname: I{destinationPoint}.
        '''
        a = _angular(distance, radius)
        sa, ca, sb, cb = sincos2(a, Bearing_(bearing))

        p = self.toNvector()
        e = NorthPole.cross(p, raiser=_pole_).unit()  # east vector at p
        n = p.cross(e)  # north vector at p
        q = n.times(cb).plus(e.times(sb))  # direction vector @ p
        n = p.times(ca).plus(q.times(sa))
        return n.toLatLon(height=height, LatLon=self.classof)  # Nvector(n.x, n.y, n.z).toLatLon(...)

    def distanceTo(self, other, radius=R_M, **unused):  # for -DistanceTo
        '''Compute the distance from this to an other point.

           @arg other: The other point (L{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}).

           @return: Distance between this and the B{C{other}} point
                    (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(48.857, 2.351);
           >>> d = p.distanceTo(q)  # 404.3 km
        '''
        self.others(other)

        return self.toNvector().angleTo(other.toNvector()) * Radius(radius)

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

        sa, ca, sb, cb, st, ct = sincos2(a, b, t)
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
                              L{crosserrors} is C{True}.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @example:

           >>> p1 = LatLon(52.205, 0.119)
           >>> p2 = LatLon(48.857, 2.351)
           >>> b = p1.initialBearingTo(p2)  # 156.2

           @JSname: I{bearingTo}.
        '''
        self.others(other, name=_other_)
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

        f = Scalar(fraction, name=_fraction_)
        i = other.toNvector().times(f).plus(
             self.toNvector().times(1 - f))
#       i = other.toNvector() * f + \
#            self.toNvector() * (1 - f))

        h = self._havg(other, f=f) if height is None else Height(height)
        return i.toLatLon(height=h, LatLon=self.classof)  # Nvector(i.x, i.y, i.z).toLatLon(...)

    def intermediateTo(self, other, fraction, height=None):
        '''Locate the point at a given fraction between this and an
           other point.

           @arg other: The other point (L{LatLon}).
           @arg fraction: Fraction between both points (float, between
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
        f = Scalar(fraction, name=_fraction_)

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
           >>> inside = p.isEnclosedBy(b)  # True

           @JSname: I{enclosedBy}.
        '''
        n, points = self.points2(points, closed=True)

        # use normal vector to this point for sign of α
        n0 = self.toNvector()

        if iterNumpy2(points):

            def _subtangles(n, points, n0):  # iterate
                vs = n0.minus(points[n-1].toNvector())
                for i in range(n):
                    vs1 = n0.minus(points[i].toNvector())
                    yield vs.angleTo(vs1, vSign=n0)  # PYCHOK false
                    vs = vs1

            # sum subtended angles
            s = fsum(_subtangles(n, points, n0))

        else:
            # get vectors from this to each point
            vs = [n0.minus(points[i].toNvector()) for i in range(n)]
            # close the polygon
            vs.append(vs[0])

            # sum subtended angles of each edge (using n0 to determine sign)
            s = fsum(vs[i].angleTo(vs[i+1], vSign=n0) for i in range(n))

        # Note, this method uses angle summation test: on a plane,
        # angles for an enclosed point will sum to 360°, angles for
        # an exterior point will sum to 0°.  On a sphere, enclosed
        # point angles will sum to less than 360° (due to spherical
        # excess), exterior point angles will be small but non-zero.

        # XXX are winding number optimisations equally applicable to
        # spherical surface?
        return abs(s) > PI

    def isEnclosedBy(self, points):  # PYCHOK no cover
        '''DEPRECATED, use method C{isenclosedBy}.
        '''
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
        n1 = self.others(point1, name='point1').toNvector()
        n2 = self.others(point2, name='point2').toNvector()

        # corner case, null arc
        if n1.isequalTo(n2):
            return n0.isequalTo(n1) or n0.isequalTo(n2)  # PYCHOK returns

        if n0.dot(n1) < 0 or n0.dot(n2) < 0:  # different hemisphere
            return False  # PYCHOK returns

        # get vectors representing d0=p0->p1 and d2=p2->p1
        # and dot product d0⋅d2 tells us if p0 is on the
        # p2 side of p1 or on the other side (similarly
        # for d0=p0->p2 and d1=p1->p2 and dot product
        # d0⋅d1 and p0 on the p1 side of p2 or not)
        return n0.minus(n1).dot(n2.minus(n1)) >= 0 and \
               n0.minus(n2).dot(n1.minus(n2)) >= 0

    def isWithin(self, point1, point2):  # PYCHOK no cover
        '''DEPRECATED, use method C{iswithin}.
        '''
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

    def nearestOn(self, point1, point2, height=None):
        '''Locate the point on the great circle arc between two
           points closest to this point.

           If this point is within the extent of the arc between both
           end points, return the closest point on the arc.  Otherwise,
           return the closest of the arc's end points.

           @arg point1: Start point of the arc (L{LatLon}).
           @arg point2: End point of the arc (L{LatLon}).
           @kwarg height: Optional height, overriding the mean height
                          for the point within the arc (C{meter}).

           @return: Closest point on the arc (L{LatLon}).

           @raise TypeError: If B{C{point1}} or B{C{point2}} is not L{LatLon}.

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
        if self.isWithin(point1, point2) and not point1.isequalTo(point2, EPS):
            # closer to arc than to its endpoints,
            # find the closest point on the arc
            gc1 = point1.toNvector().cross(point2.toNvector())
            gc2 = self.toNvector().cross(gc1)
            p = gc1.cross(gc2).toLatLon(height=height or 0, LatLon=self.classof)
            if height is None:  # interpolate height
                d = point1.distanceTo(point2)
                f = 0.5 if d < EPS else point1.distanceTo(p) / d
                p.height = point1._havg(point2, f=f)

        # beyond arc extent, take closer endpoint
        elif self.distanceTo(point1) < self.distanceTo(point2):
            p = point1
        else:
            p = point2

        return p

    def nearestOn2(self, points, **closed_radius_height):  # PYCHOK no cover
        '''DEPRECATED, use method L{sphericalNvector.LatLon.nearestOn3}.

           @return: ... 2-Tuple C{(closest, distance)} of the C{closest}
                    point (L{LatLon}) on the polygon and the C{distance}
                    to that point from this point ...
        '''
        r = self.nearestOn3(points, **closed_radius_height)
        return r.closest, r.distance

    def nearestOn3(self, points, closed=False, radius=R_M, height=None):
        '''Locate the point on a polygon (with great circle arcs
           joining consecutive points) closest to this point.

           If this point is within the extent of any great circle
           arc, the closest point is on that arc.  Otherwise,
           the closest is the nearest of the arc's end points.

           @arg points: The polygon points (L{LatLon}[]).
           @kwarg closed: Optionally, close the polygon (C{bool}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg height: Optional height, overriding the mean height
                          for a point within the arc (C{meter}).

           @return: A L{NearestOn3Tuple}C{(closest, distance, angle)} of
                    the C{closest} point (L{LatLon}), the C{distance}
                    between this and the C{closest} point in C{meter},
                    same units as B{C{radius}} and the C{angle} from
                    this to the C{closest} point in compass C{degrees360}.

           @raise TypeError: Some B{C{points}} are not C{LatLon}.

           @raise ValueError: No B{C{points}}.
        '''
        n, points = self.points2(points, closed=closed)

        i, m = _imdex2(closed, n)
        c = p2 = points[i]
        r = self.distanceTo(c, radius=1)  # force radians
        for i in range(m, n):
            p1, p2 = p2, points[i]
            p = self.nearestOn(p1, p2, height=height)
            t = self.distanceTo(p, radius=1)  # force radians
            if t < r:
                c, r = p, t

        return NearestOn3Tuple(c, r * Radius(radius), degrees360(r))

    def toCartesian(self, **Cartesian_datum_kwds):  # PYCHOK Cartesian=Cartesian, datum=None
        '''Convert this point to C{Nvector}-based cartesian (ECEF)
           coordinates.

           @kwarg Cartesian_datum_kwds: Optional L{Cartesian}, B{C{datum}} or
                                        other keyword arguments, ignored if
                                        B{C{Cartesian=None}}.  Use
                                        B{C{Cartesian=...}} to override this
                                        L{Cartesian} class or specify
                                        B{C{Cartesian=None}}.

           @return: The cartesian point (L{Cartesian}) or if B{C{Cartesian}}
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian_datum_kwds}}.
        '''
        kwds = _xkwds(Cartesian_datum_kwds, Cartesian=Cartesian,
                                                datum=self.datum)
        return LatLonSphericalBase.toCartesian(self, **kwds)

    def toNvector(self, **Nvector_kwds):  # PYCHOK signature
        '''Convert this point to L{Nvector} components, I{including
           height}.

           @kwarg Nvector_kwds: Optional L{Nvector} keyword arguments,
                                ignored if B{C{Nvector=None}}.  Use
                                B{C{Nvector=...}} to override this
                                L{Nvector} class or specify
                                B{C{Nvector=None}}

           @return: The C{n-vector} components (L{Nvector}) or a
                    L{Vector4Tuple}C{(x, y, z, h)} if B{C{Nvector}}
                    is C{None}.

           @raise TypeError: Invalid B{C{Nvector_kwds}}.

           @example:

           >>> p = LatLon(45, 45)
           >>> n = p.toNvector()
           >>> n.toStr()  # [0.50, 0.50, 0.70710]

           @JSname: I{toVector}.
        '''
        kwds = _xkwds(Nvector_kwds, Nvector=Nvector)
        return LatLonNvectorBase.toNvector(self, **kwds)

    def triangulate(self, bearing1, other, bearing2, height=None):
        '''Locate a point given this and an other point and a bearing
           at this and the other point.

           @arg bearing1: Bearing at this point (compass C{degrees360}).
           @arg other: The other point (L{LatLon}).
           @arg bearing2: Bearing at the other point (compass C{degrees360}).
           @kwarg height: Optional height at the triangulated point,
                          overriding the mean height (C{meter}).

           @return: Triangulated point (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise Valuerror: Points coincide.

           @example:

           >>> p = LatLon("47°18.228'N","002°34.326'W")  # Basse Castouillet
           >>> q = LatLon("47°18.664'N","002°31.717'W")  # Basse Hergo
           >>> t = p.triangulate(7, q, 295)  # 47.323667°N, 002.568501°W'
        '''
        return triangulate(self, bearing1, other, bearing2,
                                 height=height, LatLon=self.classof)

    def trilaterate(self, distance1, point2, distance2, point3, distance3,
                          radius=R_M, height=None, useZ=False):
        '''Locate a point at given distances from this and two other points.
           See also U{Trilateration<https://WikiPedia.org/wiki/Trilateration>}.

           @arg distance1: Distance to this point (C{meter}, same units
                           as B{C{radius}}).
           @arg point2: Second reference point (L{LatLon}).
           @arg distance2: Distance to point2 (C{meter}, same units as
                           B{C{radius}}).
           @arg point3: Third reference point (L{LatLon}).
           @arg distance3: Distance to point3 (C{meter}, same units as
                           B{C{radius}}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg height: Optional height at trilaterated point, overriding
                          the mean height (C{meter}, same units as B{C{radius}}).
           @kwarg useZ: Include Z component iff non-NaN, non-zero (C{bool}).

           @return: Trilaterated point (L{LatLon}).

           @raise TypeError: Some B{C{points}} are not L{LatLon}.

           @raise ValueError: Distance(s) exceeds trilateration or
                              some B{C{points}} coincide.
        '''
        return trilaterate(self, distance1, point2, distance2,
                                            point3, distance3,
                                 radius=radius, height=height,
                                 LatLon=self.classof, useZ=useZ)


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
    _datum = Datums.Sphere  #: (INTERNAL) Default datum (L{Datum}).
    _Ecef  = EcefKarney     #: (INTERNAL) Preferred C{Ecef...} class.

    def toCartesian(self, **Cartesian_h_kwds):  # PYCHOK Cartesian=Cartesian
        '''Convert this n-vector to C{Nvector}-based cartesian
           (ECEF) coordinates.

           @kwarg Cartesian_h_kwds: Optional L{Cartesian}, B{C{h}} and
                                    other keyword arguments, ignored if
                                    B{C{Cartesian=None}}.  Use
                                    B{C{Cartesian=...}} to override this
                                    L{Cartesian} class or specify
                                    B{C{Cartesian=None}}.

           @return: The cartesian point (L{Cartesian}) or if B{C{Cartesian}}
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian}}, B{C{h}} or other
                             B{C{Cartesian_h_kwds}}.
        '''
        kwds = _xkwds(Cartesian_h_kwds, h=self.h, Cartesian=Cartesian)
        return NvectorBase.toCartesian(self, **kwds)  # class or .classof

    def toLatLon(self, **LatLon_height_kwds):  # PYCHOK height=None, LatLon=LatLon
        '''Convert this n-vector to an C{Nvector}-based geodetic point.

           @kwarg LatLon_height_kwds: Optional L{LatLon}, B{C{height}} and
                                      other keyword arguments, ignored if
                                      B{C{LatLon=None}}.  Use
                                      B{C{LatLon=...}} to override this
                                      L{LatLon} class or specify
                                      B{C{LatLon=None}}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}}
                    is C{None}, a L{LatLon3Tuple}C{(lat, lon, height)}.

           @raise TypeError: Invalid B{C{LatLon}}, B{C{height}} or
                             other B{C{LatLon_height_kwds}}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        kwds = _xkwds(LatLon_height_kwds, height=self.h, LatLon=LatLon)
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


_Nvll = LatLon(0, 0, name='Nv00')  #: (INTERNAL) Reference instance (L{LatLon}).


def areaOf(points, radius=R_M):
    '''Calculate the area of a (spherical) polygon (with great circle
       arcs joining consecutive points).

       @arg points: The polygon points (L{LatLon}[]).
       @kwarg radius: Mean earth radius (C{meter}).

       @return: Polygon area (C{meter}, same units as B{C{radius}}, squared).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @see: L{pygeodesy.areaOf}, L{sphericalTrigonometry.areaOf} and
             L{ellipsoidalKarney.areaOf}.

       @example:

       >>> b = LatLon(45, 1), LatLon(45, 2), LatLon(46, 2), LatLon(46, 1)
       >>> areaOf(b)  # 8666058750.718977
    '''
    n, points = _Nvll.points2(points, closed=True)

    # use vector to 1st point as plane normal for sign of α
    n0 = points[0].toNvector()

    if iterNumpy2(points):

        def _interangles(n, points, n0):  # iterate
            v2 = points[n-2]._N_vector
            v1 = points[n-1]._N_vector
            gc = v2.cross(v1)
            for i in range(n):
                v2 = points[i]._N_vector
                gc1 = v1.cross(v2)
                v1 = v2

                yield gc.angleTo(gc1, vSign=n0)
                gc = gc1

        # sum interior angles
        s = fsum(_interangles(n, points, n0))

    else:
        # get great-circle vector for each edge
        gc, v1 = [], points[n-1]._N_vector
        for i in range(n):
            v2 = points[i]._N_vector
            gc.append(v1.cross(v2))  # PYCHOK false, does have .cross
            v1 = v2
        gc.append(gc[0])  # XXX needed?

        # sum interior angles: depending on whether polygon is cw or ccw,
        # angle between edges is π−α or π+α, where α is angle between
        # great-circle vectors; so sum α, then take n·π − |Σα| (cannot
        # use Σ(π−|α|) as concave polygons would fail)
        s = fsum(gc[i].angleTo(gc[i+1], vSign=n0) for i in range(n))

    # using Girard’s theorem: A = [Σθᵢ − (n−2)·π]·R²
    # (PI2 - abs(s) == (n*PI - abs(s)) - (n-2)*PI)
    return abs(PI2 - abs(s)) * Radius(radius)**2


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
                           arguments, ignored if B{C{LatLon=None}}.

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
    _Nvll.others(start1, name=_start1_)
    _Nvll.others(start2, name=_start2_)

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
            # intersection is at further-away intersection
            # point, take opposite intersection from mid-
            # point of v1 and v2 [is this always true?]
            d = -s1.plus(s2).dot(i1)

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
                           arguments, ignored if B{C{LatLon=None}}.

       @return: Point at geographic mean and mean height (B{C{LatLon}}).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not C{LatLon}.
    '''
    n, points = _Nvll.points2(points, closed=False)
    # geographic mean
    m = sumOf(points[i]._N_vector for i in range(n))
    kwds = _xkwds(LatLon_kwds, height=height, LatLon=LatLon)
    return m.toLatLon(**kwds)


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
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg height: Optional height, overriding the mean height
                      for a point within the arc (C{meter}).

       @return: A L{NearestOn3Tuple}C{(closest, distance, angle)} of
                the C{closest} point (L{LatLon}) on the polygon, the
                C{distance} and the C{angle} between the C{closest}
                and the given B{C{point}}.  The C{distance} is in
                C{meter}, same units as B{C{radius}}, the C{angle}
                is in compass C{degrees360}.

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
       @kwarg radius: Mean earth radius (C{meter}).

       @return: Polygon perimeter (C{meter}, same units as B{C{radius}}).

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not L{LatLon}.

       @see: L{pygeodesy.perimeterOf}, L{sphericalTrigonometry.perimeterOf}
             and L{ellipsoidalKarney.perimeterOf}.
    '''
    n, points = _Nvll.points2(points, closed=closed)

    def _rads(n, points, closed):  # angular edge lengths in radians
        i, m = _imdex2(closed, n)
        v1 = points[i]._N_vector
        for i in range(m, n):
            v2 = points[i]._N_vector
            yield v1.angleTo(v2)
            v1 = v2

    r = fsum(_rads(n, points, closed))
    return r * Radius(radius)


def sumOf(nvectors, Vector=Nvector, h=None, **Vector_kwds):
    '''Return the vectorial sum of two or more n-vectors.

       @arg nvectors: Vectors to be added (L{Nvector}[]).
       @kwarg Vector: Optional class for the vectorial sum (L{Nvector}).
       @kwarg h: Optional height, overriding the mean height (C{meter}).
       @kwarg kwds: Optional, additional B{C{Vector}} keyword arguments.

       @return: Vectorial sum (B{C{Vector}}).

       @raise VectorError: No B{C{nvectors}}.
    '''
    return _sumOf(nvectors, Vector=Vector, h=h, **Vector_kwds)


def triangulate(point1, bearing1, point2, bearing2,
                height=None, LatLon=LatLon, **LatLon_kwds):
    '''Locate a point given two known points and initial bearings
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
                           arguments, ignored if B{C{LatLon=None}}.

       @return: Triangulated point (B{C{LatLon}}).

       @raise TypeError: If B{C{point1}} or B{C{point2}} is not L{LatLon}.

       @raise Valuerror: Points coincide.

       @example:

       >>> p = LatLon("47°18.228'N","002°34.326'W")  # Basse Castouillet
       >>> q = LatLon("47°18.664'N","002°31.717'W")  # Basse Hergo
       >>> t = triangulate(p, 7, q, 295)  # 47.323667°N, 002.568501°W'
    '''
    _Nvll.others(point1, name='point1')
    _Nvll.others(point2, name='point2')

    if point1.isequalTo(point2, EPS):
        raise _ValueError(points=point2, txt=_coincident_)

    def _gc(p, b, i):
        n = p.toNvector()
        de = NorthPole.cross(n, raiser=_pole_).unit()  # east vector @ n
        dn = n.cross(de)  # north vector @ n
        s, c = sincos2d(Bearing(b, name=_bearing_ + i))
        dest = de.times(s)
        dnct = dn.times(c)
        d = dnct.plus(dest)  # direction vector @ n
        return n.cross(d)  # great circle point + bearing

    gc1 = _gc(point1, bearing1, _1_)  # great circle p1 + b1
    gc2 = _gc(point2, bearing2, _2_)  # great circle p2 + b2

    n = gc1.cross(gc2, raiser=_points_)  # n-vector of intersection point

    h = point1._havg(point2) if height is None else Height(height)
    kwds = _xkwds(LatLon_kwds, height=h, LatLon=LatLon)
    return n.toLatLon(**kwds)  # Nvector(n.x, n.y, n.z).toLatLon(...)


def trilaterate(point1, distance1, point2, distance2, point3, distance3,
                radius=R_M, height=None, useZ=False, LatLon=LatLon, **LatLon_kwds):
    '''Locate a point at given distances from three other points.
       See also U{Trilateration<https://WikiPedia.org/wiki/Trilateration>}.

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
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if B{C{LatLon=None}}.

       @return: Trilaterated point (B{C{LatLon}}).

       @raise TypeError: If B{C{point1}}, B{C{point2}} or B{C{point3}}
                         is not L{LatLon}.

       @raise ValueError: Invalid B{C{distance1}}, B{C{distance2}},
                          B{C{distance3}} or B{C{radius}}, or some
                          B{C{distances}} exceed trilateration or
                          some B{C{points}} coincide.
    '''
    def _nd2(p, d, r, i, *qs):
        # return Nvector and angular distance squared
        _Nvll.others(p, name=_point_ + i)
        for q in qs:
            if p.isequalTo(q, EPS):
                raise _ValueError(points=p, txt=_coincident_)
        return p.toNvector(), (Scalar(d, name=_distance_ + i) / r)**2

    r = Radius_(radius)

    n1, d12 = _nd2(point1, distance1, r, _1_)
    n2, d22 = _nd2(point2, distance2, r, _2_, point1)
    n3, d32 = _nd2(point3, distance3, r, '3', point1, point2)

    # the following uses x,y coordinate system with origin at n1, x axis n1->n2
    y = n3.minus(n1)
    x = n2.minus(n1)

    d = x.length  # distance n1->n2
    if d > EPS_2:  # and y.length > EPS_2:
        X = x.unit()  # unit vector in x direction n1->n2
        i = X.dot(y)  # signed magnitude of x component of n1->n3
        Y = y.minus(X.times(i)).unit()  # unit vector in y direction
        j = Y.dot(y)  # signed magnitude of y component of n1->n3
        if abs(j) > EPS_2:
            # courtesy Carlos Freitas <https://GitHub.com/mrJean1/PyGeodesy/issues/33>
            x = fsum_(d12, -d22, d**2) / (2 * d)  # n1->intersection x- and ...
            y = fsum_(d12, -d32, i**2, j**2) / (2 * j) - (x * i / j)  # ... y-component

            n = n1.plus(X.times(x)).plus(Y.times(y))  # .plus(Z.times(z))
            if useZ:  # include non-NaN, non-zero Z component
                z = fsum_(d12, -(x**2), -(y**2))
                if z > EPS:
                    Z = X.cross(Y)  # unit vector perpendicular to plane
                    n = n.plus(Z.times(sqrt(z)))

            if height is None:
                h = fidw((point1.height, point2.height, point3.height),
                         map1(fabs, distance1, distance2, distance3))
            else:
                h = Height(height)
            kwds = _xkwds(LatLon_kwds, height=h, LatLon=LatLon)
            return n.toLatLon(**kwds)  # Nvector(n.x, n.y, n.z).toLatLon(...)

    # no intersection, d < EPS_2 or abs(j) < EPS_2
    t = unstr(trilaterate.__name__, point1, distance1,
                                    point2, distance2,
                                    point3, distance3, useZ=useZ, d=d)
    raise _ValueError('no intersection', txt=t)


__all__ += _ALL_OTHER(Cartesian, LatLon, Nvector,  # classes
                      areaOf,  # functions
                      intersection, ispolar,
                      meanOf,
                      nearestOn2,
                      perimeterOf,
                      sumOf,
                      triangulate, trilaterate)

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
