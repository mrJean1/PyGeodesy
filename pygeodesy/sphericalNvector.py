
# -*- coding: utf-8 -*-

u'''N-vector-based spherical geodetic (lat-/longitude) classes L{LatLon}
and L{Nvector} and functions L{areaOf}, L{intersection}, L{meanOf},
L{nearestOn2}, L{triangulate} and L{trilaterate}.

Pure Python implementation of n-vector-based spherical geodetic (lat-/longitude)
methods, transcribed from JavaScript originals by I{(C) Chris Veness 2011-2016},
published under the same MIT Licence**.  See U{Vector-based geodesy
<http://www.Movable-Type.co.UK/scripts/latlong-vectors.html>} and
U{Module latlon-nvector-spherical
<http://www.Movable-Type.co.UK/scripts/geodesy/docs/module-latlon-nvector-spherical.html>}.

Tools for working with points and paths on (a spherical model of) the
earth’s surface using using n-vectors rather than the more common
spherical trigonometry.  N-vectors make many calculations much simpler,
and easier to follow, compared with the trigonometric equivalents.

Based on Kenneth Gade’s U{‘Non-singular Horizontal Position Representation’
<http://www.NavLab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf>},
The Journal of Navigation (2010), vol 63, nr 3, pp 395-417.

Note that the formulations below take x => 0°N,0°E, y => 0°N,90°E and
z => 90°N while Gade uses x => 90°N, y => 0°N,90°E, z => 0°N,0°E.

Also note that on a spherical earth model, an n-vector is equivalent
to a normalised version of an (ECEF) cartesian coordinate.

@newfield example: Example, Examples
'''

from datum import R_M
from fmath import EPS, fmean, fsum, fsum_, isscalar
from nvector import NorthPole, LatLonNvectorBase, \
                    Nvector as NvectorBase, sumOf
from points import _imdex2, ispolar  # PYCHOK ispolar
from sphericalBase import LatLonSphericalBase
from utily import PI, PI2, PI_2, degrees360, iterNumpy2

from math import atan2, cos, radians, sin

# all public contants, classes and functions
__all__ = ('LatLon', 'Nvector',  # classes
           'areaOf',  # functions
           'intersection', 'ispolar',
           'meanOf',
           'nearestOn2',
           'triangulate', 'trilaterate')
__version__ = '18.10.26'


class LatLon(LatLonNvectorBase, LatLonSphericalBase):
    '''New n-vector based point on spherical earth model.

       Tools for working with points and paths on (a spherical
       model of) the earth's surface using vector-based methods.

       @example:

       >>> from sphericalNvector import LatLon
       >>> p = LatLon(52.205, 0.119)
    '''
    _Nv = None  #: (INTERNAL) cache _toNvector L{Nvector}).

    def _gc3(self, start, end, namend, raiser='points'):
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

    def _update(self, updated):
        '''(INTERNAL) Clear caches id updated.
        '''
        if updated:  # reset caches
            if self._Nv:
                self._Nv._fromll = None
                self._Nv = None
            LatLonNvectorBase._update(self, updated)
            LatLonSphericalBase._update(self, updated)

    def alongTrackDistanceTo(self, start, end, radius=R_M):
        '''Compute the (signed) distance from the start to the closest
           point on the great circle path defined by a start and an
           end point.

           That is, if a perpendicular is drawn from this point to the
           great circle path, the along-track distance is the distance
           from the start point to the point where the perpendicular
           crosses the path.

           @param start: Start point of great circle path (L{LatLon}).
           @param end: End point of great circle path (L{LatLon}) or
                       initial bearing from start point (compass
                       C{degrees360}).
           @keyword radius: Optional, mean earth radius (C{meter}).

           @return: Distance along the great circle path (positive if
                    after the start toward the end point of the path
                    or negative if before the start point).

           @raise TypeError: The I{start} or I{end} point is not L{LatLon}.

           @raise Valuerror: Some points coincide.

           @example:

           >>> p = LatLon(53.2611, -0.7972)

           >>> s = LatLon(53.3206, -1.7297)
           >>> e = LatLon(53.1887, 0.1334)
           >>> d = p.alongTrackDistanceTo(s, e)  # 62331.58
        '''
        self.others(start, name='start')
        gc, _, _ = self._gc3(start, end, 'end')

        p = self.toNvector()
        a = gc.cross(p).cross(gc)  # along-track point gc × p × gc
        return start.toNvector().angleTo(a, vSign=gc) * radius

    def bearingTo(self, other, **unused):
        '''DEPRECATED, use method C{initialBearingTo}.
        '''
        return self.initialBearingTo(other)

    def copy(self):
        '''Copy this point.

           @return: The copy (L{LatLon} or subclass thereof).
        '''
        return LatLonNvectorBase.copy(self)

    def crossTrackDistanceTo(self, start, end, radius=R_M):
        '''Compute the (signed) distance from this point to great circle
           defined by a start and end point.

           @param start: Start point of great circle path (L{LatLon}).
           @param end: End point of great circle path (L{LatLon}) or
                       initial bearing from start point (compass
                       C{degrees360}).
           @keyword radius: Optional, mean earth radius (C{meter}).

           @return: Distance to great circle (negative if to the
                    left or positive if to the right of the path).

           @raise TypeError: The I{start} or I{end} point is not L{LatLon}.

           @raise Valuerror: Some points coincide.

           @example:

           >>> p = LatLon(53.2611, -0.7972)

           >>> s = LatLon(53.3206, -1.7297)
           >>> d = p.crossTrackDistanceTo(s, 96)  # -305.7

           >>> e = LatLon(53.1887, 0.1334)
           >>> d = p.crossTrackDistanceTo(s, e)  # -307.5
        '''
        self.others(start, name='start')
        gc, _, _ = self._gc3(start, end, 'end')

        p = self.toNvector()
        return (gc.angleTo(p) - PI_2) * radius

    def destination(self, distance, bearing, radius=R_M, height=None):
        '''Locate the destination from this point after having
           travelled the given distance on the given bearing.

           @param distance: Distance travelled (C{meter}, same units
                            as I{radius}).
           @param bearing: Bearing from this point (compass C{degrees360}).
           @keyword radius: Optional, mean earth radius (C{meter}).
           @keyword height: Optional height at destination, overriding
                            the default height (C{meter}, same units as
                            I{radius}).

           @return: Destination point (L{LatLon}).

           @raise Valuerror: Polar coincidence.

           @example:

           >>> p = LatLon(51.4778, -0.0015)
           >>> q = p.destination(7794, 300.7)
           >>> q.toStr()  # 51.513546°N, 000.098345°W

           @JSname: I{destinationPoint}.
        '''
        p = self.toNvector()

        e = NorthPole.cross(p, raiser='pole').unit()  # east vector at p
        n = p.cross(e)  # north vector at p

        t = radians(bearing)
        q = n.times(cos(t)).plus(e.times(sin(t)))  # direction vector @ p

        r = float(distance) / float(radius)  # angular distance in radians
        n = p.times(cos(r)).plus(q.times(sin(r)))
        return n.toLatLon(height=height, LatLon=self.classof)  # Nvector(n.x, n.y, n.z).toLatLon(...)

    def distanceTo(self, other, radius=R_M):
        '''Compute the distance from this to an other point.

           @param other: The other point (L{LatLon}).
           @keyword radius: Optional, mean earth radius (C{meter}).

           @return: Distance between this and the I{other} point
                    (C{meter}, same units as I{radius}).

           @raise TypeError: The I{other} point is not L{LatLon}.

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(48.857, 2.351);
           >>> d = p.distanceTo(q)  # 404.3 km
        '''
        self.others(other)

        return self.toNvector().angleTo(other.toNvector()) * radius

    def greatCircle(self, bearing):
        '''Compute the vector normal to great circle obtained by
           heading on the given bearing from this point.

           Direction of vector is such that initial bearing vector
           b = c × n, where n is an n-vector representing this point.

           @param bearing: Bearing from this point (compass C{degrees360}).

           @return: N-vector representing great circle (L{Nvector}).

           @example:

           >>> p = LatLon(53.3206, -1.7297)
           >>> gc = p.greatCircle(96.0)
           >>> gc.toStr()  # [-0.794, 0.129, 0.594]
        '''
        a, b = self.to2ab()
        t = radians(bearing)

        ca, sa = cos(a), sin(a)
        cb, sb = cos(b), sin(b)
        ct, st = cos(t), sin(t)

        return Nvector(sb * ct - sa * cb * st,
                      -cb * ct - sa * sb * st,
                       ca * st, name=self.name)  # XXX .unit()

    def greatCircleTo(self, other):
        '''Compute the vector normal to great circle obtained by
           heading from this to an other point or on a given bearing.

           Direction of vector is such that initial bearing vector
           b = c × n, where n is an n-vector representing this point.

           @param other: The other point (L{LatLon}) or the bearing
                         from this point (compass C{degrees360}).

           @return: N-vector representing great circle (L{Nvector}).

           @raise TypeError: The I{other} point is not L{LatLon}.

           @raise Valuerror: Points coincide.

           @example:

           >>> p = LatLon(53.3206, -1.7297)
           >>> gc = p.greatCircle(96.0)
           >>> gc.toStr()  # (-0.79408, 0.12856, 0.59406)

           >>> q = LatLon(53.1887, 0.1334)
           >>> g = p.greatCircleTo(q)
           >>> g.toStr()  # (-0.79408, 0.12859, 0.59406)
        '''
        gc, _, _ = self._gc3(self, other, 'other')
        return gc.unit()

    def initialBearingTo(self, other, **unused):
        '''Compute the initial bearing (forward azimuth) from this
           to an other point.

           @param other: The other point (L{LatLon}).
           @param unused: Optional keyword argument I{wrap} ignored.

           @return: Initial bearing (compass C{degrees360}).

           @raise Crosserror: This point coincides with the I{other}
                              point or the C{NorthPole}, provided
                              L{crosserrors} is C{True}.

           @raise TypeError: The I{other} point is not L{LatLon}.

           @example:

           >>> p1 = LatLon(52.205, 0.119)
           >>> p2 = LatLon(48.857, 2.351)
           >>> b = p1.initialBearingTo(p2)  # 156.2

           @JSname: I{bearingTo}.
        '''
        self.others(other, name='other')
        # see <http://MathForum.org/library/drmath/view/55417.html>
        n = self.toNvector()
#       gc1 = self.greatCircleTo(other)
        gc1 = n.cross(other.toNvector(), raiser='points')  # .unit()
#       gc2 = self.greatCircleTo(NorthPole)
        gc2 = n.cross(NorthPole, raiser='pole')  # .unit()
        return degrees360(gc1.angleTo(gc2, vSign=n))

    def intermediateChordTo(self, other, fraction, height=None):
        '''Locate the point projected from the point at given fraction
           on a straight line (chord) between this and an other point.

           @param other: The other point (L{LatLon}).
           @param fraction: Fraction between both points (float, 0.0 =
                            this point, 1.0 = other point).
           @keyword height: Optional height at the intermediate point,
                            overriding the fractional height (C{meter}).

           @return: Intermediate point (L{LatLon}).

           @raise TypeError: The I{other} point is not L{LatLon}.

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(48.857, 2.351)
           >>> i = p.intermediateChordTo(q, 0.25)  # 51.3723°N, 000.7072°E

           @JSname: I{intermediatePointOnChordTo}, I{intermediatePointDirectlyTo}.
        '''
        self.others(other)

        i = other.toNvector().times(fraction).plus(
             self.toNvector().times(1 - fraction))
#       i = other.toNvector() * fraction + \
#            self.toNvector() * (1 - fraction))

        return i.toLatLon(height=height, LatLon=self.classof)  # Nvector(i.x, i.y, i.z).toLatLon(...)

    def intermediateTo(self, other, fraction, height=None):
        '''Locate the point at a given fraction between this and an
           other point.

           @param other: The other point (L{LatLon}).
           @param fraction: Fraction between both points (float, 0.0 =
                            this point, 1.0 = the other point).
           @keyword height: Optional height at the intermediate point,
                            overriding the fractional height (C{meter}).

           @return: Intermediate point (L{LatLon}).

           @raise TypeError: The I{other} point is not L{LatLon}.

           @raise Valuerror: Points coincide.

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(48.857, 2.351)
           >>> i = p.intermediateTo(q, 0.25)  # 51.3721°N, 000.7074°E

           @JSname: I{intermediatePointTo}.
        '''
        self.others(other)

        p = self.toNvector()
        q = other.toNvector()

        x = p.cross(q, raiser='points')
        d = x.unit().cross(p)  # unit(p × q) × p
        # angular distance tan(a) = |p × q| / p ⋅ q
        a = atan2(x.length, p.dot(q)) * fraction  # interpolated
        i = p.times(cos(a)).plus(d.times(sin(a)))  # p * cosα + d * sinα

        return i.toLatLon(height=height, LatLon=self.classof)  # Nvector(i.x, i.y, i.z).toLatLon(...)

    def intersection(self, end1, start2, end2, height=None):
        '''Locates the point of intersection of two paths each defined
           by two points or a start point and bearing from North.

           @param end1: End point of first path (L{LatLon}) or the
                        initial bearing at this point (compass
                        C{degrees360}).
           @param start2: Start point of the second path (L{LatLon}).
           @param end2: End point of second path (L{LatLon}) or the
                        initial bearing at the second point (compass
                        C{degrees}).
           @keyword height: Optional height at the intersection point,
                            overriding the mean height (C{meter}).

           @return: Intersection point (L{LatLon}).

           @raise TypeError: The I{start2}, I{end1} or I{end2} is
                             not L{LatLon}.

           @raise Valuerror: Points coincide.

           @example:

           >>> s = LatLon(51.8853, 0.2545)
           >>> e = LatLon(49.0034, 2.5735)
           >>> i = s.intersection(108.55, e, 32.44)  # 50.9076°N, 004.5086°E
        '''
        return intersection(self, end1, start2, end2,
                            height=height, LatLon=self.classof)

    def isenclosedBy(self, points):
        '''Check whether this point is enclosed by a (convex) polygon.

           @param points: The polygon points (L{LatLon}[]).

           @return: C{True} if the polygon encloses this point,
                    C{False} otherwise.

           @raise TypeError: Some I{points} are not L{LatLon}.

           @raise ValueError: Insufficient number of I{points}.

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

    def isEnclosedBy(self, points):
        '''DEPRECATED, use method C{isenclosedBy}.
        '''
        return self.isenclosedBy(points)

    def iswithin(self, point1, point2):
        '''Check whether this point is between two other points.

           If this point is not on the great circle arc defined by
           both points, return whether it is within the area bound
           by perpendiculars to the great circle at each point (in
           the same hemispere).

           @param point1: Start point of the arc (L{LatLon}).
           @param point2: End point of the arc (L{LatLon}).

           @return: C{True} if this point is within the arc,
                    C{False} otherwise.

           @raise TypeError: If I{point1} or I{point2} is not L{LatLon}.

           @JSname: I{isBetween}.
        '''
        self.others(point1, name='point1')
        self.others(point2, name='point2')

        n0 = self.toNvector()
        n1 = point1.toNvector()
        n2 = point2.toNvector()

        # corner case, null arc
        if n1.isequalTo(n2):
            return n0.isequalTo(n1) or n0.isequalTo(n2)

        if n0.dot(n1) < 0 or n0.dot(n2) < 0:
            return False  # different hemisphere

        # get vectors representing d0=p0->p1 and d2=p2->p1
        # and dot product d0⋅d2 tells us if p0 is on the
        # p2 side of p1 or on the other side (similarly
        # for d0=p0->p2 and d1=p1->p2 and dot product
        # d0⋅d1 and p0 on the p1 side of p2 or not)
        return n0.minus(n1).dot(n2.minus(n1)) >= 0 and \
               n0.minus(n2).dot(n1.minus(n2)) >= 0

    def isWithin(self, point1, point2):
        '''DEPRECATED, use method C{iswithin}.
        '''
        return self.iswithin(point1, point2)

    def midpointTo(self, other, height=None):
        '''Find the midpoint between this and an other point.

           @param other: The other point (L{LatLon}).
           @keyword height: Optional height at the midpoint,
                            overriding the mean height (C{meter}).

           @return: Midpoint (L{LatLon}).

           @raise TypeError: The I{other} point is not L{LatLon}.

           @example:

           >>> p1 = LatLon(52.205, 0.119)
           >>> p2 = LatLon(48.857, 2.351)
           >>> m = p1.midpointTo(p2)  # '50.5363°N, 001.2746°E'
        '''
        self.others(other)

        m = self.toNvector().plus(other.toNvector())
        return m.toLatLon(height=height, LatLon=self.classof)

    def nearestOn(self, point1, point2, height=None):
        '''Locate the point on the great circle arc between two
           points closest to this point.

           If this point is within the extent of the arc between both
           end points, return the closest point on the arc.  Otherwise,
           return the closest of the arc's end points.

           @param point1: Start point of the arc (L{LatLon}).
           @param point2: End point of the arc (L{LatLon}).
           @keyword height: Optional height, overriding the mean height
                            for the point within the arc (C{meter}).

           @return: Closest point on the arc (L{LatLon}).

           @raise TypeError: If I{point1} or I{point2} is not L{LatLon}.

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
            p = gc1.cross(gc2).toLatLon(height=height, LatLon=self.classof)

        # beyond arc extent, take closer endpoint
        elif self.distanceTo(point1) < self.distanceTo(point2):
            p = point1
        else:
            p = point2

        return p

    def nearestOn2(self, points, closed=False, radius=R_M, height=None):
        '''Locate the point on a polygon (with great circle arcs
           joining consecutive points) closest to this point.

           If this point is within the extent of any great circle
           arc, the closest point is on that arc.  Otherwise,
           the closest is the nearest of the arc's end points.

           @param points: The polygon points (L{LatLon}[]).
           @keyword closed: Optionally, close the polygon (C{bool}).
           @keyword radius: Optional, mean earth radius (C{meter}).
           @keyword height: Optional height, overriding the mean height
                            for a point within the arc (C{meter}).

           @return: 2-Tuple (closest, distance) of the closest point
                    (L{LatLon}) on the polygon and the distance to
                    that point from the given point in C{meter}, same
                    units of I{radius}.

           @raise TypeError: Some I{points} are not C{LatLon}.

           @raise ValueError: No I{points}.
        '''
        n, points = self.points2(points, closed=closed)

        i, m = _imdex2(closed, n)
        c = p2 = points[i]
        d = self.distanceTo(c, radius=radius)
        for i in range(m, n):
            p1, p2 = p2, points[i]
            p = self.nearestOn(p1, p2, height=height)
            t = self.distanceTo(p, radius=radius)
            if t < d:
                c, d = p, t
        return c, d

    def toNvector(self):
        '''Convert this (geodetic) point to a (spherical) n-vector
           normal to the earth's surface.

           @return: N-vector representing this point (L{Nvector}).

           @example:

           >>> p = LatLon(45, 45)
           >>> n = p.toNvector()
           >>> n.toStr()  # [0.50, 0.50, 0.70710]

           @JSname: I{toVector}.
        '''
        if self._Nv is None:
            x, y, z, h = self.to4xyzh()
            self._Nv = Nvector(x, y, z, h=h, ll=self, name=self.name)
        return self._Nv

    def triangulate(self, bearing1, other, bearing2, height=None):
        '''Locate a point given this and an other point and a bearing
           at this and the other point.

           @param bearing1: Bearing at this point (compass C{degrees360}).
           @param other: The other point (L{LatLon}).
           @param bearing2: Bearing at the other point (compass
                            C{degrees360}).
           @keyword height: Optional height at the triangulated point,
                            overriding the mean height (C{meter}).

           @return: Triangulated point (L{LatLon}).

           @raise TypeError: The I{other} point is not L{LatLon}.

           @raise Valuerror: Points coincide.

           @example:

           >>> p = LatLon("47°18.228'N","002°34.326'W")  # Basse Castouillet
           >>> q = LatLon("47°18.664'N","002°31.717'W")  # Basse Hergo
           >>> t = p.triangulate(7, q, 295)  # 47.323667°N, 002.568501°W'
        '''
        return triangulate(self, bearing1, other, bearing2,
                                 height=height, LatLon=self.classof)

    def trilaterate(self, distance1, point2, distance2, point3, distance3,
                          radius=R_M, height=None):
        '''Locate a point at given distances from this and two other points.
           See also U{Trilateration<http://WikiPedia.org/wiki/Trilateration>}.

           @param distance1: Distance to this point (C{meter}, same units
                             as I{radius}).
           @param point2: Second reference point (L{LatLon}).
           @param distance2: Distance to point2 (C{meter}, same units as
                             I{radius}).
           @param point3: Third reference point (L{LatLon}).
           @param distance3: Distance to point3 (C{meter}, same units as
                             I{radius}).
           @keyword radius: Optional, mean earth radius (C{meter}).
           @keyword height: Optional height at trilaterated point,
                            overriding the mean height (C{meter}, same
                            units as I{radius}).

           @return: Trilaterated point (L{LatLon}).

           @raise TypeError: One of the I{points} is not L{LatLon}.

           @raise ValueError: Distance(s) exceeds trilateration or
                              some I{points} coincide.
        '''
        return trilaterate(self, distance1, point2, distance2,
                                            point3, distance3,
                                 radius=radius, height=height,
                                 LatLon=self.classof)


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

    def toLatLon(self, height=None, LatLon=LatLon):
        '''Convert this n-vector to a (spherical) geodetic point.

           @keyword height: Optional height above earth radius,
                            overriding the default height (C{meter}).
           @keyword LatLon: Optional (spherical) LatLon (sub-)class
                            to use for the point (L{LatLon}) or C{None}.

           @return: Point equivalent to this n-vector (L{LatLon}) or
                    3-tuple (C{degrees90}, C{degrees180}, height) if
                    I{LatLon} is C{None}.

           @example:

           >>> n = Nvector(0.5, 0.5, 0.7071)
           >>> p = n.toLatLon()  # 45.0°N, 45.0°E
        '''
        return self._toLLh(LatLon, height)

    def greatCircle(self, bearing):
        '''Compute the n-vector normal to great circle obtained by
           heading on given compass bearing from this point as its
           n-vector.

           Direction of vector is such that initial bearing vector
           b = c × p.

           @param bearing: Initial compass bearing (C{degrees}).

           @return: N-vector representing great circle (L{Nvector}).

           @raise Valuerror: Polar coincidence.

           @example:

           >>> n = LatLon(53.3206, -1.7297).toNvector()
           >>> gc = n.greatCircle(96.0)  # [-0.794, 0.129, 0.594]
        '''
        t = radians(bearing)

        e = NorthPole.cross(self, raiser='pole')  # easting
        n = self.cross(e, raiser='point')  # northing

        e = e.times(cos(t) / e.length)
        n = n.times(sin(t) / n.length)
        return n.minus(e)


_Nvll = LatLon(0, 0, name='Nv00')  #: (INTERNAL) Reference instance (L{LatLon}).


def areaOf(points, radius=R_M):
    '''Calculate the area of a (spherical) polygon (with great circle
       arcs joining consecutive points).

       @param points: The polygon points (L{LatLon}[]).
       @keyword radius: Optional, mean earth radius (C{meter}).

       @return: Polygon area (C{meter}, same units as I{radius}, squared).

       @raise TypeError: Some I{points} are not L{LatLon}.

       @raise ValueError: Insufficient number of I{points}.

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
            v2 = points[n-2].toNvector()
            v1 = points[n-1].toNvector()
            gc = v2.cross(v1)
            for i in range(n):
                v2 = points[i].toNvector()
                gc1 = v1.cross(v2)
                v1 = v2

                yield gc.angleTo(gc1, vSign=n0)
                gc = gc1

        # sum interior angles
        s = fsum(_interangles(n, points, n0))

    else:
        # get great-circle vector for each edge
        gc, v1 = [], points[n-1].toNvector()
        for i in range(n):
            v2 = points[i].toNvector()
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
    return abs(PI2 - abs(s)) * radius**2


def intersection(start1, end1, start2, end2,
                 height=None, LatLon=LatLon):
    '''Locate the intersection of two paths each defined by two
       points or by a start point and an initial bearing.

       @param start1: Start point of the first path (L{LatLon}).
       @param end1: End point of first path (L{LatLon}) or the
                    initial bearing at the first start point
                    (compass C{degrees360}).
       @param start2: Start point of the second path (L{LatLon}).
       @param end2: End point of second path (L{LatLon}) or the
                    initial bearing at the second start point
                    (compass C{degrees360}).
       @keyword height: Optional height at the intersection point,
                        overriding the default height (C{meter}).
       @keyword LatLon: Optional (sub-)class for the intersection
                        point (L{LatLon}).

       @return: Intersection point (L{LatLon}) or C{None} if no
                unique intersection exists.

       @raise TypeError: Start or end point(s) not L{LatLon}.

       @raise Valuerror: Paths coincide.

       @example:

       >>> p = LatLon(51.8853, 0.2545)
       >>> q = LatLon(49.0034, 2.5735)
       >>> i = intersection(p, 108.55, q, 32.44)  # 50.9076°N, 004.5086°E
    '''
    _Nvll.others(start1, name='start1')
    _Nvll.others(start2, name='start2')

    # If gc1 and gc2 are great circles through start and end points
    # (or defined by start point and bearing), then the candidate
    # intersections are simply gc1 × gc2 and gc2 × gc1.  Most of the
    # work is deciding the correct intersection point to select!  If
    # bearing is given, that determines the intersection, but if both
    # paths are defined by start/end points, take closer intersection.
    gc1, s1, e1 = _Nvll._gc3(start1, end1, 'end1')
    gc2, s2, e2 = _Nvll._gc3(start2, end2, 'end2')

    # there are two (antipodal) candidate intersection
    # points ... we have to choose the one to return
    i1 = gc1.cross(gc2, raiser='paths')
    i2 = gc2.cross(gc1, raiser='paths')

    # selection of intersection point depends on how
    # paths are defined (by bearings or endpoints)
    if e1 and e2:  # endpoint+endpoint
        d = sumOf((s1, s2, e1, e2)).dot(i1)
    elif e1 and not e2:  # endpoint+bearing
        # gc2 x v2 . i1 +ve means v2 bearing points to i1
        d = gc2.cross(s2).dot(i1)
    elif e2 and not e1:  # bearing+endpoint
        # gc1 x v1 . i1 +ve means v1 bearing points to i1
        d = gc1.cross(s1).dot(i1)
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

    i = i1 if d > 0 else i2
    return i.toLatLon(height=height, LatLon=LatLon)  # Nvector(i.x, i.y, i.z).toLatLon(...)


def meanOf(points, height=None, LatLon=LatLon):
    '''Compute the geographic mean of the supplied points.

       @param points: Array of points to be averaged (L{LatLon}[]).
       @keyword height: Optional height, overriding the mean height (C{meter}).
       @keyword LatLon: Optional (sub-)class for the mean point (L{LatLon}).

       @return: Point at geographic mean and mean height (L{LatLon}).

       @raise ValueError: Insufficient number of I{points}.
    '''
    n, points = _Nvll.points2(points, closed=False)
    # geographic mean
    m = sumOf(points[i].toNvector() for i in range(n))
    a, b, h = m.to3llh()

    return LatLon(a, b, height=h if height is None else height)


def nearestOn2(point, points, closed=False, radius=R_M, height=None):
    '''Locate the point on a polygon (with great circle arcs
       joining consecutive points) closest to an other point.

       If the given point is within the extent of any great circle
       arc, the closest point is on that arc.  Otherwise, the
       closest is the nearest of the arc's end points.

       @param point: The other, reference point (L{LatLon}).
       @param points: The polygon points (L{LatLon}[]).
       @keyword closed: Optionally, close the polygon (C{bool}).
       @keyword radius: Optional, mean earth radius (C{meter}).
       @keyword height: Optional height, overriding the mean height
                        for a point within the arc (C{meter}).

       @return: 2-Tuple (I{closest}, I{distance}) of the closest point
                (L{LatLon}) on the polygon and the I{distance} from the
                I{closest} to the other I{point} in C{meter}, same
                units as I{radius}.

       @raise TypeError: Some I{points} or the point not C{LatLon}.

       @raise ValueError: No I{points}.
    '''
    if not isinstance(point, LatLon):
        raise TypeError('%s not %r: %r' % ('point', LatLon, point))

    return point.nearestOn2(points, closed=closed, radius=radius, height=height)


def triangulate(point1, bearing1, point2, bearing2,
                height=None, LatLon=LatLon):
    '''Locate a point given two known points and initial bearings
       from those points.

       @param point1: First reference point (L{LatLon}).
       @param bearing1: Bearing at the first point (compass C{degrees360}).
       @param point2: Second reference point (L{LatLon}).
       @param bearing2: Bearing at the second point (compass C{degrees360}).
       @keyword height: Optional height at the triangulated point,
                        overriding the mean height (C{meter}).
       @keyword LatLon: Optional (sub-)class for the triangulated point
                        (L{LatLon}).

       @return: Triangulated point (L{LatLon}).

       @raise TypeError: One of the I{points} is not L{LatLon}.

       @raise Valuerror: Points coincide.

       @example:

       >>> p = LatLon("47°18.228'N","002°34.326'W")  # Basse Castouillet
       >>> q = LatLon("47°18.664'N","002°31.717'W")  # Basse Hergo
       >>> t = triangulate(p, 7, q, 295)  # 47.323667°N, 002.568501°W'
    '''
    _Nvll.others(point1, name='point1')
    _Nvll.others(point2, name='point2')

    if point1.isequalTo(point2, EPS):
        raise ValueError('%s %s: %r' % ('coincident', 'points', point2))

    def _gc(p, b):
        n, t = p.toNvector(), radians(b)
        de = NorthPole.cross(n, raiser='pole').unit()  # east vector @ n
        dn = n.cross(de)  # north vector @ n
        dest = de.times(sin(t))
        dnct = dn.times(cos(t))
        d = dnct.plus(dest)  # direction vector @ n
        return n.cross(d)  # great circle point + bearing

    gc1 = _gc(point1, bearing1)  # great circle p1 + b1
    gc2 = _gc(point2, bearing2)  # great circle p2 + b2

    n = gc1.cross(gc2, raiser='points')  # n-vector of intersection point

    if height is None:
        h = point1._havg(point2)
    else:
        h = height
    return n.toLatLon(height=h, LatLon=LatLon)  # Nvector(n.x, n.y, n.z).toLatLon(...)


def trilaterate(point1, distance1, point2, distance2, point3, distance3,
                radius=R_M, height=None, LatLon=LatLon):
    '''Locate a point at given distances from three other points.
       See also U{Trilateration<http://WikiPedia.org/wiki/Trilateration>}.

       @param point1: First point (L{LatLon}).
       @param distance1: Distance to the first point (C{meter}, same units
                         as I{radius}).
       @param point2: Second point (L{LatLon}).
       @param distance2: Distance to the second point (C{meter}, same units
                         as I{radius}).
       @param point3: Third point (L{LatLon}).
       @param distance3: Distance to the third point (C{meter}, same units
                         as I{radius}).
       @keyword radius: Optional, mean earth radius (C{meter}).
       @keyword height: Optional height at the trilaterated point, overriding
                        the mean height (C{meter}, same units as I{radius}).
       @keyword LatLon: Optional (sub-)class for the trilaterated point
                        (L{LatLon}).

       @return: Trilaterated point (L{LatLon}).

       @raise TypeError: One of the I{points} is not L{LatLon}.

       @raise ValueError: Invalid I{radius}, some I{distances} exceed
                          trilateration or some I{points} coincide.
    '''
    def _nd2(p, d, name, *qs):
        # return Nvector and radial distance squared
        _Nvll.others(p, name=name)
        for q in qs:
            if p.isequalTo(q, EPS):
                raise ValueError('%s %s: %r' % ('coincident', 'points', p))
        return p.toNvector(), (float(d) / radius)**2

    if float(radius or 0) < EPS:
        raise ValueError('%s %s: %r' % ('radius', 'invalid', radius))

    n1, d12 = _nd2(point1, distance1, 'point1')
    n2, d22 = _nd2(point2, distance2, 'point2', point1)
    n3, d32 = _nd2(point3, distance3, 'point3', point1, point2)

    # the following uses x,y coordinate system with origin at n1, x axis n1->n2
    x = n2.minus(n1)
    y = n3.minus(n1)
    z = 0

    d = x.length  # distance n1->n2
    if d > EPS:  # and (y.length * 2) > EPS:
        X = x.unit()  # unit vector in x direction n1->n2
        i = X.dot(y)  # signed magnitude of x component of n1->n3
        Y = y.minus(X.times(i)).unit()  # unit vector in y direction
        j = Y.dot(y)  # signed magnitude of y component of n1->n3
        if abs(j) > EPS:
            x = fsum_(d12, -d22, d**2) / d  # n1->intersection x- and ...
            y = fsum_(d12, -d32, i**2, j**2, -2 * x * i) / j  # ... y-component
            z = (x**2 + y**2) * 0.25

#   z = sqrt(d12 - z)  # z will be NaN for no intersections
    if not 0 < z < d12:
        raise ValueError('no %s for %r, %r, %r at %r, %r, %r' %
                        ('trilaterate', point1, point2, point3,
                          distance1, distance2, distance3))
#   Z = X.cross(Y)  # unit vector perpendicular to plane
    # note don't use Z component; assume points at same height
    n = n1.plus(X.times(x)).plus(Y.times(y))  # .plus(Z.times(z))

    if height is None:
        h = fmean((point1.height, point2.height, point3.height))
    else:
        h = height
    return n.toLatLon(height=h, LatLon=LatLon)  # Nvector(n.x, n.y, n.z).toLatLon(...)

# **) MIT License
#
# Copyright (C) 2016-2018 -- mrJean1 at Gmail dot com
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
