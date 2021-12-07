
# -*- coding: utf-8 -*-

u'''Extended 3-D vector class L{Vector3d} and functions.

Function L{intersection3d3}, L{intersections2}, L{parse3d}, L{sumOf},
L{trilaterate2d2} and L{trilaterate3d2}.
'''

from pygeodesy.basics import len2, map2
from pygeodesy.errors import IntersectionError, _TypeError, _ValueError, \
                             VectorError, _xError, _xkwds, _xkwds_popitem
from pygeodesy.fmath import fdot, fsum, fsum1_
from pygeodesy.formy import _radical2
from pygeodesy.interns import EPS, EPS0, EPS1, EPS4, MISSING, NN, \
                             _COMMA_, _concentric_, _datum_, _h_, _height_, \
                             _intersection_, _name_, _near_, _negative_, \
                             _no_, _too_, _xyz_, _y_, _z_, _0_0, _1_0
from pygeodesy.iters import Fmt, PointsIter
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import _xnamed, _xotherError
from pygeodesy.namedTuples import Intersection3Tuple, NearestOn2Tuple, \
                                  NearestOn6Tuple, Vector3Tuple  # Vector4Tuple
# from pygeodesy.streprs import Fmt  # from .iters
from pygeodesy.units import _fi_j2, Radius, Radius_
from pygeodesy.vector3dBase import Vector3dBase

from math import sqrt

__all__ = _ALL_LAZY.vector3d
__version__ = '21.11.28'


class Vector3d(Vector3dBase):
    '''Extended 3-D vector.

       In a geodesy context, these may be used to represent:
        - n-vector representing a normal to a point on earth's surface
        - earth-centered, earth-fixed cartesian (ECEF)
        - great circle normal to vector
        - motion vector on earth's surface
        - etc.
    '''

    def circin6(self, point2, point3, eps=EPS4):
        '''Return the radius and center of the I{inscribed} aka I{In- circle}
           of a (3-D) triangle formed by this and two other points.

           @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                        C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
           @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                        C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
           @kwarg eps: Tolerance for function L{pygeodesy.trilaterate3d2} if
                       C{B{useZ} is True} otherwise L{pygeodesy.trilaterate2d2}.

           @return: L{Circin6Tuple}C{(radius, center, deltas, cA, cB, cC)}.  The
                    C{center} and contact points C{cA}, C{cB} and C{cC}, each an
                    instance of this (sub-)class, are co-planar with this and the
                    two given points.

           @raise ImportError: Package C{numpy} not found, not installed or older
                               than version 1.10.

           @raise IntersectionError: Near-coincident or -colinear points or
                                     a trilateration or C{numpy} issue.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @see: Function L{pygeodesy.circin6}, U{Incircle
                 <https://MathWorld.Wolfram.com/Incircle.html>} and U{Contact
                 Triangle<https://MathWorld.Wolfram.com/ContactTriangle.html>}.
        '''
        from pygeodesy.vector2d import _circin6
        try:
            return _circin6(self, point2, point3, eps=eps, useZ=True)
        except (AssertionError, TypeError, ValueError) as x:
            raise _xError(x, point=self, point2=point2, point3=point3)

    def circum3(self, point2, point3, circum=True, eps=EPS4):
        '''Return the radius and center of the smallest circle I{through} or
           I{containing} this and two other (3-D) points.

           @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}
                        or C{Vector4Tuple}).
           @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}
                        or C{Vector4Tuple}).
           @kwarg circum: If C{True} return the C{circumradius} and C{circumcenter},
                          always, ignoring the I{Meeus}' Type I case (C{bool}).
           @kwarg eps: Tolerance passed to function L{pygeodesy.trilaterate3d2}.

           @return: A L{Circum3Tuple}C{(radius, center, deltas)}.  The C{center}, an
                    instance of this (sub-)class, is co-planar with this and the two
                    given points.

           @raise ImportError: Package C{numpy} not found, not installed or older than
                               version 1.10.

           @raise IntersectionError: Near-concentric, -coincident or -colinear points
                                     or a trilateration or C{numpy} issue.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @see: Function L{pygeodesy.circum3} and methods L{circum4_} and L{meeus2}.
        '''
        from pygeodesy.vector2d import _circum3
        try:
            return _circum3(self, point2, point3, circum=circum, eps=eps, useZ=True,
                                                  clas=self.classof)
        except (AssertionError, TypeError, ValueError) as x:
            raise _xError(x, point=self, point2=point2, point3=point3, circum=circum)

    def circum4_(self, *points):
        '''Best-fit a sphere through this and two or more other (3-D) points.

           @arg points: Other points (each a C{Cartesian}, L{Vector3d}, C{Vector3Tuple}
                        or C{Vector4Tuple}).

           @return: L{Circum4Tuple}C{(radius, center, rank, residuals)} with C{center}
                    an instance if this (sub-)class.

           @raise ImportError: Package C{numpy} not found, not installed or
                               older than version 1.10.

           @raise NumPyError: Some C{numpy} issue.

           @raise PointsError: Too few B{C{points}}.

           @raise TypeError: One of the B{C{points}} invalid.

           @see: Function L{pygeodesy.circum4_} and methods L{circum3} and L{meeus2}.
        '''
        from pygeodesy.vector2d import circum4_
        return circum4_(self, *points, Vector=self.classof)

    def iscolinearWith(self, point1, point2, eps=EPS):
        '''Check whether this and two other (3-D) points are colinear.

           @arg point1: One point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}
                        or C{Vector4Tuple}).
           @arg point2: An other point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}
                        or C{Vector4Tuple}).
           @kwarg eps: Tolerance (C{scalar}), same units as C{x},
                       C{y}, and C{z}.

           @return: C{True} if this point is colinear with B{C{point1}} and
                    B{C{point2}}, C{False} otherwise.

           @raise TypeError: Invalid B{C{point1}} or B{C{point2}}.

           @see: Method L{nearestOn}.
        '''
        from pygeodesy.vector2d import _iscolinearWith
        v = self if self.name else _otherV3d(NN_OK=False, this=self)
        return _iscolinearWith(v, point1, point2, eps=eps)

    def meeus2(self, point2, point3, circum=False):
        '''Return the radius and I{Meeus}' Type of the smallest circle
           I{through} or I{containing} this and two other (3-D) points.

           @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}
                        or C{Vector4Tuple}).
           @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}
                        or C{Vector4Tuple}).
           @kwarg circum: If C{True} return the C{circumradius} and C{circumcenter}
                          always, ignoring the I{Meeus}' Type I case (C{bool}).

           @return: L{Meeus2Tuple}C{(radius, Type)}.

           @raise IntersectionError: Coincident or colinear points, iff C{B{circum}=True}.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @see: Function L{pygeodesy.meeus2} and methods L{circum3} and L{circum4_}.
        '''
        from pygeodesy.vector2d import _meeus4, Meeus2Tuple
        try:
            r, t, _, _ = _meeus4(self, point2, point3, circum=circum, clas=self.classof)
        except (TypeError, ValueError) as x:
            raise _xError(x, point=self, point2=point2, point3=point3, circum=circum)
        return Meeus2Tuple(r, t)

    def nearestOn(self, point1, point2, within=True):
        '''Locate the point between two points closest to this point.

           @arg point1: Start point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}
                        or C{Vector4Tuple}).
           @arg point2: End point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}
                        or C{Vector4Tuple}).
           @kwarg within: If C{True} return the closest point between
                          the given points, otherwise the closest
                          point on the extended line through both
                          points (C{bool}).

           @return: Closest point, either B{C{point1}} or B{C{point2}} or an
                    instance of this (sub-)class.

           @raise TypeError: Invalid B{C{point1}} or B{C{point2}}.

           @see: Method L{sphericalTrigonometry.LatLon.nearestOn3} and
                 U{3-D Point-Line Distance<https://MathWorld.Wolfram.com/
                 Point-LineDistance3-Dimensional.html>}.
        '''
        return _nearestOn2(self, point1, point2, within=within).closest

    def nearestOn6(self, points, closed=False, useZ=True):  # eps=EPS
        '''Locate the point on a path or polygon closest to this point.

           The closest point is either on and within the extent of a polygon
           edge or the nearest of that edge's end points.

           @arg points: The path or polygon points (C{Cartesian}, L{Vector3d},
                        C{Vector3Tuple} or C{Vector4Tuple}[]).
           @kwarg closed: Optionally, close the path or polygon (C{bool}).
           @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=0} (C{bool}).

           @return: A L{NearestOn6Tuple}C{(closest, distance, fi, j, start, end)}
                    with the C{closest}, the C{start} and the C{end} point each
                    an instance of this point's (sub-)class.

           @raise PointsError: Insufficient number of B{C{points}}

           @raise TypeError: Non-cartesian B{C{points}}.

           @note: Distances measured with method L{Vector3d.equirectangular}.

           @see: Function L{nearestOn6}.
        '''
        return nearestOn6(self, points, closed=closed, useZ=useZ)  # Vector=self.classof

    def parse(self, str3d, sep=_COMMA_, name=NN):
        '''Parse an C{"x, y, z"} string to a L{Vector3d} instance.

           @arg str3d: X, y and z string (C{str}), see function L{parse3d}.
           @kwarg sep: Optional separator (C{str}).
           @kwarg name: Optional instance name (C{str}), overriding this name.

           @return: The instance (L{Vector3d}).

           @raise VectorError: Invalid B{C{str3d}}.
        '''
        return parse3d(str3d, sep=sep, Vector=self.classof, name=name or self.name)

    def radii11(self, point2, point3):
        '''Return the radii of the C{Circum-}, C{In-}, I{Soddy} and C{Tangent}
           circles of a (3-D) triangle.

           @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                        C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
           @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                        C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).

           @return: L{Radii11Tuple}C{(rA, rB, rC, cR, rIn, riS, roS, a, b, c, s)}.

           @raise TriangleError: Near-coincident or -colinear points.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @see: Function L{pygeodesy.radii11}, U{Incircle
                 <https://MathWorld.Wolfram.com/Incircle.html>}, U{Soddy Circles
                 <https://MathWorld.Wolfram.com/SoddyCircles.html>} and U{Tangent
                 Circles<https://MathWorld.Wolfram.com/TangentCircles.html>}.
        '''
        from pygeodesy.vector2d import _radii11ABC
        try:
            return _radii11ABC(self, point2, point3, useZ=True)[0]
        except (TypeError, ValueError) as x:
            raise _xError(x, point=self, point2=point2, point3=point3)

    def soddy4(self, point2, point3, eps=EPS4):
        '''Return the radius and center of the C{inner} I{Soddy} circle of a
           (3-D) triangle.

           @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                        C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
           @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                        C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
           @kwarg eps: Tolerance for function L{pygeodesy.trilaterate3d2} if
                       C{B{useZ} is True} otherwise L{pygeodesy.trilaterate2d2}.

           @return: L{Soddy4Tuple}C{(radius, center, deltas, outer)}.  The C{center},
                    an instance of B{C{point1}}'s (sub-)class, is co-planar with the
                    three given points.

           @raise ImportError: Package C{numpy} not found, not installed or older
                               than version 1.10.

           @raise IntersectionError: Near-coincident or -colinear points or
                                     a trilateration or C{numpy} issue.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @see: Function L{pygeodesy.soddy4}.
        '''
        from pygeodesy.vector2d import soddy4
        return soddy4(self, point2, point3, eps=eps, useZ=True)

    def trilaterate2d2(self, radius, center2, radius2, center3, radius3, eps=EPS, z=0):
        '''Trilaterate this and two other circles, each given as a (2-D) center
           and a radius.

           @arg radius: Radius of this circle (same C{units} as this C{x} and C{y}.
           @arg center2: Center of the 2nd circle (C{Cartesian}, L{Vector3d},
                         C{Vector2Tuple}, C{Vector3Tuple} or C{Vector4Tuple}).
           @arg radius2: Radius of this circle (same C{units} as this C{x} and C{y}.
           @arg center3: Center of the 3rd circle (C{Cartesian}, L{Vector3d},
                         C{Vector2Tuple}, C{Vector3Tuple} or C{Vector4Tuple}).
           @arg radius3: Radius of the 3rd circle (same C{units} as this C{x} and C{y}.
           @kwarg eps: Check the trilaterated point I{delta} on all 3 circles (C{scalar})
                       or C{None}.
           @kwarg z: Optional Z component of the trilaterated point (C{scalar}).

           @return: Trilaterated point, an instance of this (sub-)class with C{z=B{z}}.

           @raise IntersectionError: No intersection, near-concentric or -colinear
                                     centers, trilateration failed some other way
                                     or the trilaterated point is off one circle
                                     by more than B{C{eps}}.

           @raise TypeError: Invalid B{C{center2}} or B{C{center3}}.

           @raise UnitError: Invalid B{C{radius1}}, B{C{radius2}} or B{C{radius3}}.

           @see: Function L{pygeodesy.trilaterate2d2}.
        '''
        from pygeodesy.vector2d import _trilaterate2d2

        def _xyr3(r, **name_v):
            v = _otherV3d(useZ=False, **name_v)
            return v.x, v.y, r

        try:
            return _trilaterate2d2(*(_xyr3(radius,  center=self) +
                                     _xyr3(radius2, center2=center2) +
                                     _xyr3(radius3, center3=center3)),
                                      eps=eps, Vector=self.classof, z=z)
        except (AssertionError, TypeError, ValueError) as x:
            raise _xError(x, center=self,     radius=radius,
                             center2=center2, radius2=radius2,
                             center3=center3, radius3=radius3)

    def trilaterate3d2(self, radius, center2, radius2, center3, radius3, eps=EPS):
        '''Trilaterate this and two other spheres, each given as a (3-D) center
           and a radius.

           @arg radius: Radius of this sphere (same C{units} as this C{x}, C{y}
                        and C{z}).
           @arg center2: Center of the 2nd sphere (C{Cartesian}, L{Vector3d},
                         C{Vector3Tuple} or C{Vector4Tuple}).
           @arg radius2: Radius of this sphere (same C{units} as this C{x}, C{y}
                         and C{z}).
           @arg center3: Center of the 3rd sphere (C{Cartesian}, , L{Vector3d},
                         C{Vector3Tuple} or C{Vector4Tuple}).
           @arg radius3: Radius of the 3rd sphere (same C{units} as this C{x}, C{y}
                         and C{z}).
           @kwarg eps: Tolerance (C{scalar}), same units as C{x}, C{y}, and C{z}.

           @return: 2-Tuple with two trilaterated points, each an instance of this
                    (sub-)class.  Both points are the same instance if all three
                    spheres intersect or abut in a single point.

           @raise ImportError: Package C{numpy} not found, not installed or
                               older than version 1.10.

           @raise IntersectionError: Near-concentric, -colinear, too distant or
                                     non-intersecting spheres or C{numpy} issue.

           @raise TypeError: Invalid B{C{center2}} or B{C{center3}}.

           @raise UnitError: Invalid B{C{radius}}, B{C{radius2}} or B{C{radius3}}.

           @note: Package U{numpy<https://pypi.org/project/numpy>} is required,
                  version 1.10 or later.

           @see: Norrdine, A. U{I{An Algebraic Solution to the Multilateration
                 Problem}<https://www.ResearchGate.net/publication/
                 275027725_An_Algebraic_Solution_to_the_Multilateration_Problem>}
                 and U{I{implementation}<https://www.ResearchGate.net/publication/
                 288825016_Trilateration_Matlab_Code>}.
        '''
        from pygeodesy.vector2d import _trilaterate3d2
        try:
            return _trilaterate3d2(_otherV3d(center=self, NN_OK=False),
                                    Radius_(radius, low=eps),
                                    center2, radius2, center3, radius3,
                                    eps=eps, clas=self.classof)
        except (AssertionError, TypeError, ValueError) as x:
            raise _xError(x, center=self,     radius=radius,
                             center2=center2, radius2=radius2,
                             center3=center3, radius3=radius3)


def _intersect3d3(start1, end1, start2, end2, eps=EPS, useZ=False):
    # (INTERNAL) Intersect two lines, see L{intersection} above,
    # separated to allow callers to embellish any exceptions

    def _outside(t, d2):  # -1 before start#, +1 after end#
        return -1 if t < 0 else (+1 if t > d2 else 0)  # XXX d2 + eps?

    s1 = _otherV3d(useZ=useZ, start1=start1)
    e1 = _otherV3d(useZ=useZ, end1=end1)
    s2 = _otherV3d(useZ=useZ, start2=start2)
    e2 = _otherV3d(useZ=useZ, end2=end2)

    a  = e1.minus(s1)
    b  = e2.minus(s2)
    c  = s2.minus(s1)

    ab = a.cross(b)
    d  = abs(c.dot(ab))
    e  = max(EPS0, eps or _0_0)
    if d > EPS0 and ab.length > e:
        d = d / ab.length
        if d > e:  # argonic, skew lines distance
            raise IntersectionError(skew_d=d, txt=_no_(_intersection_))

    # co-planar, non-skew lines
    ab2 = ab.length2
    if ab2 < e:  # colinear, parallel or null line(s)
        x = a.length2 > b.length2
        if x:  # make a shortest
            a,  b  = b,  a
            s1, s2 = s2, s1
            e1, e2 = e2, e1
        if b.length2 < e:  # both null
            if c.length < e:
                return s1, 0, 0
            elif e2.minus(e1).length < e:
                return e1, 0, 0
        elif a.length2 < e:  # null (s1, e1), non-null (s2, e2)
            # like _nearestOn2(s1, s2, e2, within=False, eps=e)
            t = s1.minus(s2).dot(b)
            v = s2.plus(b.times(t / b.length2))
            if s1.minus(v).length < e:
                o = _outside(t, b.length2)
                return (v, o, 0) if x else (v, 0, (o * 2))
        raise IntersectionError(length2=ab2, txt=_no_(_intersection_))

    cb =  c.cross(b)
    t  =  cb.dot(ab)
    o1 = _outside(t, ab2)
    v  =  s1.plus(a.times(t / ab2))
    o2 = _outside(v.minus(s2).dot(b), b.length2) * 2
    return v, o1, o2


def intersection3d3(start1, end1, start2, end2, eps=EPS, useZ=True,
                                              **Vector_and_kwds):
    '''Compute the intersection point of two lines, each defined by or
       through a start and end point (3-D).

       @arg start1: Start point of the first line (C{Cartesian}, L{Vector3d},
                    C{Vector3Tuple} or C{Vector4Tuple}).
       @arg end1: End point of the first line (C{Cartesian}, L{Vector3d},
                  C{Vector3Tuple} or C{Vector4Tuple}).
       @arg start2: Start point of the second line (C{Cartesian}, L{Vector3d},
                    C{Vector3Tuple} or C{Vector4Tuple}).
       @arg end2: End point of the second line (C{Cartesian}, L{Vector3d},
                  C{Vector3Tuple} or C{Vector4Tuple}).
       @kwarg eps: Tolerance for skew line distance and length (C{EPS}).
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=0} (C{bool}).
       @kwarg Vector_and_kwds: Optional class C{B{Vector}=None} to return the
                               intersection points and optional, additional B{C{Vector}}
                               keyword arguments, otherwise B{C{start1}}'s (sub-)class.

       @return: An L{Intersection3Tuple}C{(point, outside1, outside2)} with
                C{point} an instance of B{C{Vector}} or B{C{start1}}'s (sub-)class.

       @raise IntersectionError: Invalid, skew, non-co-planar or otherwise
                                 non-intersecting lines.

       @see: U{Line-line intersection<https://MathWorld.Wolfram.com/Line-LineIntersection.html>}
             and U{line-line distance<https://MathWorld.Wolfram.com/Line-LineDistance.html>},
             U{skew lines<https://MathWorld.Wolfram.com/SkewLines.html>} and U{point-line
             distance<https://MathWorld.Wolfram.com/Point-LineDistance3-Dimensional.html>}.
    '''
    try:
        v, o1, o2 = _intersect3d3(start1, end1, start2, end2, eps=eps, useZ=useZ)
    except (TypeError, ValueError) as x:
        raise _xError(x, start1=start1, end1=end1, start2=start2, end2=end2)
    v = _nVc(v, **_xkwds(Vector_and_kwds, clas=start1.classof,
                                          name=intersection3d3.__name__))
    return Intersection3Tuple(v, o1, o2)


def intersections2(center1, radius1, center2, radius2, sphere=True, **Vector_and_kwds):
    '''Compute the intersection of two spheres or circles, each defined by a
       (3-D) center point and a radius.

       @arg center1: Center of the first sphere or circle (C{Cartesian}, L{Vector3d},
                     C{Vector3Tuple} or C{Vector4Tuple}).
       @arg radius1: Radius of the first sphere or circle (same units as the
                     B{C{center1}} coordinates).
       @arg center2: Center of the second sphere or circle (C{Cartesian}, L{Vector3d},
                     C{Vector3Tuple} or C{Vector4Tuple}).
       @arg radius2: Radius of the second sphere or circle (same units as the
                     B{C{center1}} and B{C{center2}} coordinates).
       @kwarg sphere: If C{True} compute the center and radius of the intersection of
                      two spheres.  If C{False}, ignore the C{z}-component and compute
                      the intersection of two circles (C{bool}).
       @kwarg Vector_and_kwds: Optional class C{B{Vector}=None} to return the
                               intersection points and optional, additional B{C{Vector}}
                               keyword arguments, otherwise B{C{center1}}'s (sub-)class.

       @return: If B{C{sphere}} is C{True}, a 2-tuple of the C{center} and C{radius}
                of the intersection of the I{spheres}.  The C{radius} is C{0.0} for
                abutting spheres (and the C{center} is aka I{radical center}).

                If B{C{sphere}} is C{False}, a 2-tuple with the two intersection
                points of the I{circles}.  For abutting circles, both points are
                the same instance, aka I{radical center}.

       @raise IntersectionError: Concentric, invalid or non-intersecting spheres
                                 or circles.

       @raise TypeError: Invalid B{C{center1}} or B{C{center2}}.

       @raise UnitError: Invalid B{C{radius1}} or B{C{radius2}}.

       @see: U{Sphere-Sphere<https://MathWorld.Wolfram.com/Sphere-SphereIntersection.html>} and
             U{Circle-Circle<https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>}
             Intersection.
    '''
    try:
        return _intersects2(center1, Radius_(radius1=radius1),
                            center2, Radius_(radius2=radius2), sphere=sphere,
                            clas=center1.classof, **Vector_and_kwds)
    except (TypeError, ValueError) as x:
        raise _xError(x, center1=center1, radius1=radius1, center2=center2, radius2=radius2)


def _intersects2(center1, r1, center2, r2, sphere=True, too_d=None,  # in CartesianEllipsoidalBase.intersections2,
                                         **clas_Vector_and_kwds):  # .ellipsoidalBaseDI._intersections2
    # (INTERNAL) Intersect two spheres or circles, see L{intersections2}
    # above, separated to allow callers to embellish any exceptions

    def _nV3(x, y, z):
        v = Vector3d(x, y, z)
        n = intersections2.__name__
        return _nVc(v, **_xkwds(clas_Vector_and_kwds, name=n))

    def _xV3(c1, u, x, y):
        xy1 = x, y, _1_0  # transform to original space
        return _nV3(fdot(xy1, u.x, -u.y, c1.x),
                    fdot(xy1, u.y,  u.x, c1.y), _0_0)

    c1 = _otherV3d(useZ=sphere, center1=center1)
    c2 = _otherV3d(useZ=sphere, center2=center2)

    if r1 < r2:  # r1, r2 == R, r
        c1, c2 = c2, c1
        r1, r2 = r2, r1

    m = c2.minus(c1)
    d = m.length
    if d < max(r2 - r1, EPS):
        raise IntersectionError(_near_(_concentric_))

    o = fsum1_(-d, r1, r2)  # overlap == -(d - (r1 + r2))
    # compute intersections with c1 at (0, 0) and c2 at (d, 0), like
    # <https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>
    if o > EPS:  # overlapping, r1, r2 == R, r
        x = _radical2(d, r1, r2).xline
        y = _1_0 - (x / r1)**2
        if y > EPS:
            y = r1 * sqrt(y)  # y == a / 2
        elif y < 0:
            raise IntersectionError(_negative_)
        else:  # abutting
            y = _0_0
    elif o < 0:
        t = d if too_d is None else too_d
        raise IntersectionError(_too_(Fmt.distant(t)))
    else:  # abutting
        x, y = r1, _0_0

    u = m.unit()
    if sphere:  # sphere center and radius
        c = c1 if x < EPS  else (
            c2 if x > EPS1 else c1.plus(u.times(x)))
        t = _nV3(c.x, c.y, c.z), Radius(y)

    elif y > 0:  # intersecting circles
        t = _xV3(c1, u, x, y), _xV3(c1, u, x, -y)
    else:  # abutting circles
        t = _xV3(c1, u, x, 0)
        t = t, t
    return t


def iscolinearWith(point, point1, point2, eps=EPS, useZ=True):
    '''Check whether a point is colinear with two other (2- or 3-D) points.

       @arg point: The point (L{Vector3d}, C{Vector3Tuple} or C{Vector4Tuple}).
       @arg point1: First point (L{Vector3d}, C{Vector3Tuple} or C{Vector4Tuple}).
       @arg point2: Second point (L{Vector3d}, C{Vector3Tuple} or C{Vector4Tuple}).
       @kwarg eps: Tolerance (C{scalar}), same units as C{x}, C{y} and C{z}.
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=0} (C{bool}).

       @return: C{True} if B{C{point}} is colinear B{C{point1}} and B{C{point2}},
                        C{False} otherwise.

       @raise TypeError: Invalid B{C{point}}, B{C{point1}} or B{C{point2}}.

       @see: Function L{nearestOn}.
    '''
    p = _otherV3d(useZ=useZ, point=point)

    from pygeodesy.vector2d import _iscolinearWith
    return _iscolinearWith(p, point1, point2, eps=eps, useZ=useZ)


def nearestOn(point, point1, point2, within=True, useZ=True, Vector=None, **Vector_kwds):
    '''Locate the point between two points closest to a reference (2- or 3-D).

       @arg point: Reference point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}
                                    or C{Vector4Tuple}).
       @arg point1: Start point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple} or
                                 C{Vector4Tuple}).
       @arg point2: End point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple} or
                               C{Vector4Tuple}).
       @kwarg within: If C{True} return the closest point between both given
                      points, otherwise the closest point on the extended line
                      through both points (C{bool}).
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=0} (C{bool}).
       @kwarg Vector: Class to return closest point (C{Cartesian}, L{Vector3d}
                      or C{Vector3Tuple}) or C{None}.
       @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword arguments,
                           ignored if C{B{Vector} is None}.

       @return: Closest point, either B{C{point1}} or B{C{point2}} or an instance
                of the B{C{point}}'s (sub-)class or B{C{Vector}} if not C{None}.

       @raise TypeError: Invalid B{C{point}}, B{C{point1}} or B{C{point2}}.

       @see: U{3-D Point-Line Distance<https://MathWorld.Wolfram.com/Point-LineDistance3-Dimensional.html>},
             C{Cartesian} and C{LatLon} methods C{nearestOn}, method L{sphericalTrigonometry.LatLon.nearestOn3}
             and function L{sphericalTrigonometry.nearestOn3}.
    '''
    p0 = _otherV3d(useZ=useZ, point =point)
    p1 = _otherV3d(useZ=useZ, point1=point1)
    p2 = _otherV3d(useZ=useZ, point2=point2)

    p, _ = _nearestOn2(p0, p1, p2, within=within)
    if Vector is not None:
        p = Vector(p.x, p.y, **_xkwds(Vector_kwds, z=p.z, name=nearestOn.__name__))
    elif p is p1:
        p = point1
    elif p is p2:
        p = point2
    else:  # ignore Vector_kwds
        p = point.classof(p.x, p.y, Vector_kwds.get(_z_, p.z), name=nearestOn.__name__)
    return p


def _nearestOn2(p0, p1, p2, within=True, eps=EPS):
    # (INTERNAL) Closest point and fraction, see L{nearestOn} above,
    # separated to allow callers to embellish any exceptions
    p21 = p2.minus(p1)
    d2 = p21.length2
    if d2 < eps:  # coincident
        p = p1  # ~= p2
        t = 0
    else:  # see comments in .points.nearestOn5
        t = p0.minus(p1).dot(p21) / d2
        if within and t < eps:
            p = p1
            t = 0
        elif within and t > (_1_0 - eps):
            p = p2
            t = 1
        else:
            p = p1.plus(p21.times(t))
    return NearestOn2Tuple(p, t)


def nearestOn6(point, points, closed=False, useZ=True, **Vector_and_kwds):  # eps=EPS
    '''Locate the point on a path or polygon closest to a reference point.

       The closest point is either on and within the extent of a polygon edge or
       the nearest of that edge's end points.

       @arg point: Reference point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}
                                    or C{Vector4Tuple}).
       @arg points: The path or polygon points (C{Cartesian}, L{Vector3d},
                    C{Vector3Tuple} or C{Vector4Tuple}[]).
       @kwarg closed: Optionally, close the path or polygon (C{bool}).
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=0} (C{bool}).
       @kwarg Vector_and_kwds: Optional class C{B{Vector}=None} to return the
                               closest point and optional, additional B{C{Vector}}
                               keyword arguments, otherwise B{C{point}}'s (sub-)class.

       @return: A L{NearestOn6Tuple}C{(closest, distance, fi, j, start, end)} with
                the C{closest}, the C{start} and the C{end} point each an instance
                of the B{C{Vector}} keyword argument of if {B{Vector}=None} or not
                specified, an instance of the B{C{point}}'s (sub-)class.

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Non-cartesian B{C{point}} and B{C{points}}.

       @note: Distances measured with method L{Vector3d.equirectangular}.

       @see: Method C{LatLon.nearestOn6} or function L{nearestOn5} for geodetic points.
    '''
    r  = _otherV3d(useZ=useZ, point=point)
    D2 = r.equirectangular  # distance squared

    Ps = PointsIter(points, loop=1, name=nearestOn6.__name__)
    p1 = c = s = e = _otherV3d(useZ=useZ, i=0, points=Ps[0])
    c2 = D2(c)  # == r.minus(c).length2

    f = i = 0  # p1..p2 == points[i]..[j]
    for j, p2 in Ps.enumerate(closed=closed):
        p2 = _otherV3d(useZ=useZ, i=j, points=p2)
        p, t = _nearestOn2(r, p1, p2)  # within=True, eps=EPS
        d2 = D2(p)  # == r.minus(p).length2
        if d2 < c2:
            c2, c, s, e, f = d2, p, p1, p2, (i + t)
        p1, i = p2, j

    f, j = _fi_j2(f, len(Ps))  # like .ellipsoidalBaseDI._nearestOn2_

    kwds = _xkwds(Vector_and_kwds, clas=point.classof, name=Ps.name)
    v = _nVc(c, **kwds)
    s = _nVc(s, **kwds) if s is not c else v
    e = _nVc(e, **kwds) if e is not c else v
    return NearestOn6Tuple(v, sqrt(c2), f, j, s, e)


def _nVc(v, clas=None, name=NN, Vector=None, **Vector_kwds):  # in .vector2d
    # return a named C{Vector} or C{clas} instance
    if Vector is not None:
        v = Vector(v.x, v.y, v.z, **Vector_kwds)
    elif clas is not None:
        v = clas(v.x, v.y, v.z)  # ignore Vector_kwds
    return _xnamed(v, name) if name else v


def _otherV3d(useZ=True, NN_OK=True, i=None, **name_v):  # in .CartesianEllipsoidalBase.intersections2,
    # check named vector instance, return Vector3d            .Ellipsoid.height4, .formy.hartzell, .vector2d
    def _name_i(name, i):
        return name if i is None else Fmt.SQUARE(name, i)

    name, v = _xkwds_popitem(name_v)
    if useZ and isinstance(v, Vector3dBase):
        return v if NN_OK or v.name else v.copy(name=_name_i(name, i))
    try:
        return Vector3d(v.x, v.y, (v.z if useZ else _0_0), name=_name_i(name, i))
    except AttributeError:  # no _x_ or _y_ attr
        pass
    raise _xotherError(Vector3d(0, 0, 0), v, name=_name_i(name, i), up=2)


def parse3d(str3d, sep=_COMMA_, Vector=Vector3d, **Vector_kwds):
    '''Parse an C{"x, y, z"} string.

       @arg str3d: X, y and z values (C{str}).
       @kwarg sep: Optional separator (C{str}).
       @kwarg Vector: Optional class (L{Vector3d}).
       @kwarg Vector_kwds: Optional B{C{Vector}} keyword arguments,
                           ignored if C{B{Vector} is None}.

       @return: A B{C{Vector}} instance or if B{C{Vector}} is C{None},
                a named L{Vector3Tuple}C{(x, y, z)}.

       @raise VectorError: Invalid B{C{str3d}}.
    '''
    try:
        v = [float(v.strip()) for v in str3d.split(sep)]
        n = len(v)
        if n != 3:
            raise _ValueError(len=n)
    except (TypeError, ValueError) as x:
        raise VectorError(str3d=str3d, txt=str(x))
    return _xnamed((Vector3Tuple(v) if Vector is None else  # *v
                    Vector(*v, **Vector_kwds)), parse3d.__name__)


def sumOf(vectors, Vector=Vector3d, **Vector_kwds):
    '''Compute the vectorial sum of several vectors.

       @arg vectors: Vectors to be added (L{Vector3d}[]).
       @kwarg Vector: Optional class for the vectorial sum (L{Vector3d}).
       @kwarg Vector_kwds: Optional B{C{Vector}} keyword arguments,
                           ignored if C{B{Vector} is None}.

       @return: Vectorial sum as B{C{Vector}} or if B{C{Vector}} is
                C{None}, a named L{Vector3Tuple}C{(x, y, z)}.

       @raise VectorError: No B{C{vectors}}.
    '''
    n, vectors = len2(vectors)
    if n < 1:
        raise VectorError(vectors=n, txt=MISSING)

    v = Vector3Tuple(fsum(v.x for v in vectors),
                     fsum(v.y for v in vectors),
                     fsum(v.z for v in vectors))
    return _xnamed((v if Vector is None else
                         Vector(*v, **Vector_kwds)), sumOf.__name__)


def trilaterate2d2(x1, y1, radius1, x2, y2, radius2, x3, y3, radius3,
                                    eps=None, **Vector_and_kwds):
    '''Trilaterate three circles, each given as a (2-D) center and a radius.

       @arg x1: Center C{x} coordinate of the 1st circle (C{scalar}).
       @arg y1: Center C{y} coordinate of the 1st circle (C{scalar}).
       @arg radius1: Radius of the 1st circle (C{scalar}).
       @arg x2: Center C{x} coordinate of the 2nd circle (C{scalar}).
       @arg y2: Center C{y} coordinate of the 2nd circle (C{scalar}).
       @arg radius2: Radius of the 2nd circle (C{scalar}).
       @arg x3: Center C{x} coordinate of the 3rd circle (C{scalar}).
       @arg y3: Center C{y} coordinate of the 3rd circle (C{scalar}).
       @arg radius3: Radius of the 3rd circle (C{scalar}).
       @kwarg eps: Check the trilaterated point I{delta} on all 3
                   circles (C{scalar}) or C{None}.
       @kwarg Vector_and_kwds: Optional class C{B{Vector}=None} to return the
                               trilateration and optional, additional B{C{Vector}}
                               keyword arguments, otherwise (L{Vector3d}).

       @return: Trilaterated point as C{B{Vector}(x, y, **B{Vector_kwds})}
                or L{Vector2Tuple}C{(x, y)} if C{B{Vector} is None}..

       @raise IntersectionError: No intersection, near-concentric or -colinear
                                 centers, trilateration failed some other way
                                 or the trilaterated point is off one circle
                                 by more than B{C{eps}}.

       @raise UnitError: Invalid B{C{radius1}}, B{C{radius2}} or B{C{radius3}}.

       @see: U{Issue #49<https://GitHub.com/mrJean1/PyGeodesy/issues/49>},
             U{Find X location using 3 known (X,Y) location using trilateration
             <https://math.StackExchange.com/questions/884807>} and function
             L{pygeodesy.trilaterate3d2}.
    '''
    from pygeodesy.vector2d import _trilaterate2d2
    return _trilaterate2d2(x1, y1, radius1, x2, y2, radius2, x3, y3, radius3,
                                            eps=eps, **Vector_and_kwds)


def trilaterate3d2(center1, radius1, center2, radius2, center3, radius3,
                                     eps=EPS, **Vector_and_kwds):
    '''Trilaterate three spheres, each given as a (3-D) center and a radius.

       @arg center1: Center of the 1st sphere (C{Cartesian}, L{Vector3d},
                     C{Vector3Tuple} or C{Vector4Tuple}).
       @arg radius1: Radius of the 1st sphere (same C{units} as C{x}, C{y}
                     and C{z}).
       @arg center2: Center of the 2nd sphere (C{Cartesian}, L{Vector3d},
                     C{Vector3Tuple} or C{Vector4Tuple}).
       @arg radius2: Radius of this sphere (same C{units} as C{x}, C{y}
                     and C{z}).
       @arg center3: Center of the 3rd sphere (C{Cartesian}, L{Vector3d},
                     C{Vector3Tuple} or C{Vector4Tuple}).
       @arg radius3: Radius of the 3rd sphere (same C{units} as C{x}, C{y}
                     and C{z}).
       @kwarg eps: Tolerance (C{scalar}), same units as C{x}, C{y}, and C{z}.
       @kwarg Vector_and_kwds: Optional class C{B{Vector}=None} to return the
                               trilateration and optional, additional B{C{Vector}}
                               keyword arguments, otherwise B{C{center1}}'s
                               (sub-)class.

       @return: 2-Tuple with two trilaterated points, each a B{C{Vector}}
                instance.  Both points are the same instance if all three
                spheres abut/intersect in a single point.

       @raise ImportError: Package C{numpy} not found, not installed or
                           older than version 1.10.

       @raise IntersectionError: Near-concentric, -colinear, too distant or
                                 non-intersecting spheres.

       @raise NumPyError: Some C{numpy} issue.

       @raise TypeError: Invalid B{C{center1}}, B{C{center2}} or B{C{center3}}.

       @raise UnitError: Invalid B{C{radius1}}, B{C{radius2}} or B{C{radius3}}.

       @see: Norrdine, A. U{I{An Algebraic Solution to the Multilateration
             Problem}<https://www.ResearchGate.net/publication/
             275027725_An_Algebraic_Solution_to_the_Multilateration_Problem>},
             the U{I{implementation}<https://www.ResearchGate.net/publication/
             288825016_Trilateration_Matlab_Code>} and function
             L{pygeodesy.trilaterate2d2}.
    '''
    from pygeodesy.vector2d import _trilaterate3d2
    try:
        return _trilaterate3d2(_otherV3d(center1=center1, NN_OK=False),
                                Radius_(radius1=radius1, low=eps),
                                center2, radius2, center3, radius3, eps=eps,
                                clas=center1.classof, **Vector_and_kwds)
    except (AssertionError, TypeError, ValueError) as x:
        raise _xError(x, center1=center1, radius1=radius1,
                         center2=center2, radius2=radius2,
                         center3=center3, radius3=radius3)


def _xyzn4(xyz, y, z, Error=_TypeError):  # see: .ecef
    '''(INTERNAL) Get an C{(x, y, z, name)} 4-tuple.
    '''
    try:
        t = xyz.x, xyz.y, xyz.z
        n = getattr(xyz, _name_, NN)
    except AttributeError:
        t = xyz, y, z
        n = NN
    try:
        x, y, z = map2(float, t)
    except (TypeError, ValueError) as x:
        d = dict(zip((_xyz_, _y_, _z_), t))
        raise Error(txt=str(x), **d)

    return x, y, z, n


def _xyzhdn6(xyz, y, z, height, datum, ll, Error=_TypeError):  # by .cartesianBase, .nvectorBase
    '''(INTERNAL) Get an C{(x, y, z, h, d, name)} 6-tuple.
    '''
    x, y, z, n = _xyzn4(xyz, y, z, Error=Error)

    h = height or getattr(xyz, _height_, None) \
               or getattr(xyz, _h_, None) \
               or getattr(ll,  _height_, None)

    d = datum or getattr(xyz, _datum_, None) \
              or getattr(ll,  _datum_, None)

    return x, y, z, h, d, n


__all__ += _ALL_DOCS(intersections2, sumOf, Vector3dBase)

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
