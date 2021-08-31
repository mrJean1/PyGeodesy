
# -*- coding: utf-8 -*-

u'''Extended 3-D vector class L{Vector3d} and functions.

Function L{intersection3d3}, L{intersections2}, L{iscolinearWith},
L{nearestOn}, L{parse3d}, L{sumOf}, L{trilaterate2d2} and
L{trilaterate3d2}.
'''

from pygeodesy.basics import isnear0, len2, map1, map2, _xnumpy
from pygeodesy.errors import _and, _AssertionError, IntersectionError, \
                              NumPyError, PointsError, _TypeError, _ValueError, \
                              VectorError, _xError, _xkwds, _xkwds_popitem
from pygeodesy.fmath import fdot, fsum, fsum_, fsum1_, hypot, hypot2_
from pygeodesy.formy import _radical2
from pygeodesy.interns import EPS, EPS0, EPS02, EPS1, EPS4, INF, MISSING, NN, \
                             _EPS4e8, _and_, _center_, _coincident_, _colinear_, \
                             _concentric_, _COMMA_, _COMMASPACE_, _datum_, _few_, \
                             _h_, _height_, _intersection_, _invalid_, _name_, \
                             _near_, _negative_, _no_, _radius_, _s_, _SPACE_, \
                             _too_, _xyz_, _y_, _z_, _0_0, _0_5, _1_0, \
                             _1_0_T, _2_0, _4_0
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import _NamedTuple, _Pass, _xnamed, _xotherError
from pygeodesy.namedTuples import Intersection3Tuple, LatLon3Tuple, Vector2Tuple, \
                                  Vector3Tuple  # Vector4Tuple
from pygeodesy.streprs import Fmt
from pygeodesy.units import Float, Int, Meter, Radius, Radius_
from pygeodesy.vector3dBase import Vector3dBase

from contextlib import contextmanager
from math import sqrt

__all__ = _ALL_LAZY.vector3d
__version__ = '21.08.31'

_cA_        = 'cA'
_cB_        = 'cB'
_cC_        = 'cC'
_deltas_    = 'deltas'
_of_        = 'of'
_outer_     = 'outer'
_raise_     = 'raise'  # PYCHOK used!
_rank_      = 'rank'
_residuals_ = 'residuals'
_Type_      = 'Type'
_with_      = 'with'


class Circin6Tuple(_NamedTuple):
    '''6-Tuple C{(radius, center, deltas, cA, cB, cC)} with the C{radius}, the
       trilaterated C{center} and contact points of the I{inscribed} aka I{In-
       circle} of a triangle.  The C{center} is I{un}ambiguous if C{deltas} is
       C{None}, otherwise C{center} is the mean and C{deltas} the differences of
       the L{trilaterate3d2} results.  Contact points C{cA}, C{cB} and C{cC} are
       the points of tangency, aka the corners of the U{Contact Triangle
       <https://MathWorld.Wolfram.com/ContactTriangle.html>}.
    '''
    _Names_ = (_radius_, _center_, _deltas_, _cA_,  _cB_,  _cC_)
    _Units_ = ( Radius,  _Pass,    _Pass,    _Pass, _Pass, _Pass)


class Circum3Tuple(_NamedTuple):  # in .latlonBase
    '''3-Tuple C{(radius, center, deltas)} with the C{circumradius} and trilaterated
       C{circumcenter} of the C{circumcircle} through 3 points (aka {Meeus}' Type II
       circle) or the C{radius} and C{center} of the smallest I{Meeus}' Type I circle.
       The C{center} is I{un}ambiguous if C{deltas} is C{None}, otherwise C{center}
       is the mean and C{deltas} the differences of the L{trilaterate3d2} results.
    '''
    _Names_ = (_radius_, _center_, _deltas_)
    _Units_ = ( Radius,  _Pass,    _Pass)


class Circum4Tuple(_NamedTuple):
    '''4-Tuple C{(radius, center, rank, residuals)} with C{radius} and C{center}
       of a sphere I{least-squares} fitted through given points and the C{rank}
       and C{residuals} -if any- from U{numpy.linalg.lstsq
       <https://NumPy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html>}.
    '''
    _Names_ = (_radius_, _center_, _rank_, _residuals_)
    _Units_ = ( Radius,  _Pass,     Int,   _Pass)


class Meeus2Tuple(_NamedTuple):
    '''2-Tuple C{(radius, Type)} with C{radius} and I{Meeus}' C{Type} of the smallest
       circle I{containing} 3 points.  C{Type} is C{None} for a I{Meeus}' Type II
       C{circumcircle} passing through all 3 points.  Otherwise C{Type} is the center
       of a I{Meeus}' Type I circle with 2 points on (a diameter of) and 1 point
       inside the circle.
    '''
    _Names_ = (_radius_, _Type_)
    _Units_ = ( Radius,  _Pass)


class Radii11Tuple(_NamedTuple):
    '''11-Tuple C{(rA, rB, rC, cR, rIn, riS, roS, a, b, c, s)} with the C{Tangent}
       circle radii C{rA}, C{rB} and C{rC}, the C{circumradius} C{cR}, the C{Incircle}
       radius C{rIn} aka C{inradius}, the inner and outer I{Soddy} circle radii C{riS}
       and C{roS} and the sides C{a}, C{b} and C{c} and semi-perimeter C{s} of a
       triangle, all in C{meter} conventionally.

       @note: C{Circumradius} C{cR} and outer I{Soddy} radius C{roS} may be C{INF}.
    '''
    _Names_ = ('rA', 'rB', 'rC', 'cR', 'rIn', 'riS', 'roS', 'a', 'b', 'c', _s_)
    _Units_ = ( Meter,) * len(_Names_)


class Soddy4Tuple(_NamedTuple):
    '''4-Tuple C{(radius, center, deltas, outer)} with C{radius} and trilaterated
       C{center} of the I{inner} I{Soddy} circle and the radius of the C{outer}
       I{Soddy} circle.  The C{center} is I{un}ambiguous if C{deltas} is C{None},
       otherwise C{center} is the mean and C{deltas} the differences of the
       L{trilaterate3d2} results.

       @note: The outer I{Soddy} radius C{outer} may be C{INF}.
    '''
    _Names_ = (_radius_, _center_, _deltas_, _outer_)
    _Units_ = ( Radius,  _Pass,    _Pass,     Radius)


class Vector3d(Vector3dBase):
    '''Extended 3-D vector.

       In a geodesy context, these may be used to represent:
        - n-vector representing a normal to a point on earth's surface
        - earth-centered, earth-fixed cartesian (ECEF)
        - great circle normal to vector
        - motion vector on earth's surface
        - etc.
    '''
    _numpy = None  # module numpy imported once by function _numpy below

    def circin6(self, point2, point3, eps=EPS4):
        '''Return the radius and center of the I{inscribed} aka I{In- circle}
           of a (3-D) triangle formed by this and two other points.

           @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                        C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
           @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                        C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
           @kwarg eps: Tolerance for function L{trilaterate3d2} if C{B{useZ} is True}
                       else L{trilaterate2d2}.

           @return: L{Circin6Tuple}C{(radius, center, deltas, cA, cB, cC)}.  The
                    C{center} and contact points C{cA}, C{cB} and C{cC}, each an
                    instance of this (sub-)class, are co-planar with this and the
                    two given points.

           @raise ImportError: Package C{numpy} not found, not installed or older
                               than version 1.10.

           @raise IntersectionError: Near-coincident or colinear points or
                                     a trilateration or C{numpy} issue.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @see: Function L{pygeodesy.circin6}, U{Incircle
                 <https://MathWorld.Wolfram.com/Incircle.html>} and U{Contact
                 Triangle<https://MathWorld.Wolfram.com/ContactTriangle.html>}.
        '''
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
           @kwarg eps: Tolerance passed to function L{trilaterate3d2}.

           @return: A L{Circum3Tuple}C{(radius, center, deltas)}.  The C{center}, an
                    instance of this (sub-)class, is co-planar with this and the two
                    given points.

           @raise ImportError: Package C{numpy} not found, not installed or older than
                               version 1.10.

           @raise IntersectionError: Near-concentric, coincident or colinear points or
                                     a trilateration or C{numpy} issue.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @see: Function L{pygeodesy.circum3} and methods L{circum4_} and L{meeus2}.
        '''
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
                 U{3-D Point-Line distance<https://MathWorld.Wolfram.com/
                 Point-LineDistance3-Dimensional.html>}.
        '''
        return nearestOn(self, point1, point2, within=within)

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

           @raise IntersectionError: Near-coincident or colinear points.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @see: Function L{pygeodesy.radii11}, U{Incircle
                 <https://MathWorld.Wolfram.com/Incircle.html>}, U{Soddy Circles
                 <https://MathWorld.Wolfram.com/SoddyCircles.html>} and U{Tangent
                 Circles<https://MathWorld.Wolfram.com/TangentCircles.html>}.
        '''
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
           @kwarg eps: Tolerance for function L{trilaterate3d2} if C{B{useZ} is True}
                       else L{trilaterate2d2}.

           @return: L{Soddy4Tuple}C{(radius, center, deltas, outer)}.  The C{center},
                    an instance of B{C{point1}}'s (sub-)class, is co-planar with the
                    three given points.

           @raise ImportError: Package C{numpy} not found, not installed or older
                               than version 1.10.

           @raise IntersectionError: Near-coincident or colinear points or
                                     a trilateration or C{numpy} issue.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @see: Function L{pygeodesy.soddy4}.
        '''
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

           @raise IntersectionError: No intersection, colinear or near-concentric
                                     centers, trilateration failed some other way
                                     or the trilaterated point is off one circle
                                     by more than B{C{eps}}.

           @raise TypeError: Invalid B{C{center2}} or B{C{center3}}.

           @raise UnitError: Invalid B{C{radius1}}, B{C{radius2}} or B{C{radius3}}.

           @see: Function L{pygeodesy.trilaterate2d2}.
        '''
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

           @raise IntersectionError: Near-concentric, colinear, too distant or
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
        try:
            return _trilaterate3d2(_otherV3d(center=self, NN_OK=False),
                                    Radius_(radius, low=eps),
                                    center2, radius2, center3, radius3,
                                    eps=eps, clas=self.classof)
        except (AssertionError, TypeError, ValueError) as x:
            raise _xError(x, center=self,     radius=radius,
                             center2=center2, radius2=radius2,
                             center3=center3, radius3=radius3)


def circin6(point1, point2, point3, eps=EPS4, useZ=True):
    '''Return the radius and center of the I{inscribed} aka I{In- circle}
       of a (2- or 3-D) triangle.

       @arg point1: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @kwarg eps: Tolerance for function L{trilaterate3d2} if C{B{useZ} is True}
                   else L{trilaterate2d2}.
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=0} (C{bool}).

       @return: L{Circin6Tuple}C{(radius, center, deltas, cA, cB, cC)}.  The
                C{center} and contact points C{cA}, C{cB} and C{cC}, each an
                instance of B{C{point1}}'s (sub-)class, are co-planar with
                the three given points.

       @raise ImportError: Package C{numpy} not found, not installed or older
                           than version 1.10 and C{B{useZ} is True}.

       @raise IntersectionError: Near-coincident or colinear points or
                                 a trilateration or C{numpy} issue.

       @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{point3}}.

       @see: Functions L{radii11} and L{circum3}, U{Incircle
             <https://MathWorld.Wolfram.com/Incircle.html>} and U{Contact Triangle
             <https://MathWorld.Wolfram.com/ContactTriangle.html>}.
    '''
    try:
        return _circin6(point1, point2, point3, eps=eps, useZ=useZ)
    except (AssertionError, TypeError, ValueError) as x:
        raise _xError(x, point1=point1, point2=point2, point3=point3)


def _circin6(point1, point2, point3, eps=EPS4, useZ=True, dLL3=False, **Vector_kwds):
    # (INTERNAL) Radius, center, deltas, 3 contact points

    def _fraction(r, a):
        return (r / a) if a > EPS0 else _0_5

    def _contact2(a, p2, r2, p3, r3, V, V_kwds):
        c = p2.intermediateTo(p3, _fraction(r2, a)) if r2 > r3 else \
            p3.intermediateTo(p2, _fraction(r3, a))
        C = V(c.x, c.y, c.z, **V_kwds)
        return c, C

    t, p1, p2, p3 = _radii11ABC(point1, point2, point3, useZ=useZ)
    V, r1, r2, r3 =  point1.classof, t.rA, t.rB, t.rC

    c1, cA = _contact2(t.a, p2, r2, p3, r3, V, _xkwds(Vector_kwds, name=_cA_))
    c2, cB = _contact2(t.b, p3, r3, p1, r1, V, _xkwds(Vector_kwds, name=_cB_))
    c3, cC = _contact2(t.c, p1, r1, p2, r2, V, _xkwds(Vector_kwds, name=_cC_))

    r    =  t.rIn
    c, d = _tricenter3d2(c1, r, c2, r, c3, r, eps=eps, useZ=useZ, dLL3=dLL3,
                                              Vector=V, name=circin6.__name__,
                                            **Vector_kwds)
    return Circin6Tuple(r, c, d, cA, cB, cC)


def circum3(point1, point2, point3, circum=True, eps=EPS4, useZ=True):
    '''Return the radius and center of the smallest circle I{through} or I{containing}
       three (2- or 3-D) points.

       @arg point1: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple} or
                    C{Vector4Tuple}).
       @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple} or
                    C{Vector4Tuple}).
       @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple} or
                    C{Vector4Tuple}).
       @kwarg circum: If C{True} return the C{circumradius} and C{circumcenter}
                      always, ignoring the I{Meeus}' Type I case (C{bool}).
       @kwarg eps: Tolerance for function L{trilaterate3d2} if C{B{useZ} is True}
                   else L{trilaterate2d2}.
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=0} (C{bool}).

       @return: A L{Circum3Tuple}C{(radius, center, deltas)}.  The C{center}, an
                instance of B{C{point1}}'s (sub-)class, is co-planar with the three
                given points.

       @raise ImportError: Package C{numpy} not found, not installed or older
                           than version 1.10 and C{B{useZ} is True}.

       @raise IntersectionError: Near-coincident or colinear points or
                                 a trilateration or C{numpy} issue.

       @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{point3}}.

       @see: U{Jean Meeus, "Astronomical Algorithms", 2nd Ed. 1998, page 127ff
             <http://www.Agopax.IT/Libri_astronomia/pdf/Astronomical%20Algorithms.pdf>},
             U{circumradius<https://MathWorld.Wolfram.com/Circumradius.html>},
             U{circumcircle<https://MathWorld.Wolfram.com/Circumcircle.html>} and
             functions L{pygeodesy.circum4_} and L{pygeodesy.meeus2}.
    '''
    try:
        p1 = _otherV3d(useZ=useZ, point1=point1)
        return _circum3(p1, point2, point3, circum=circum, eps=eps, useZ=useZ,
                                            clas=point1.classof)
    except (AssertionError, TypeError, ValueError) as x:
        raise _xError(x, point1=point1, point2=point2, point3=point3, circum=circum)


def _circum3(p1, point2, point3, circum=True, eps=EPS4, useZ=True, dLL3=False,
                                 clas=Vector3d, **clas_kwds):  # in .latlonBase
    # (INTERNAL) Radius, center, deltas
    r, d, p2, p3 = _meeus4(p1, point2, point3, circum=circum, useZ=useZ,
                                               clas=clas, **clas_kwds)
    if d is None:  # Meeus' Type II or circum=True
        kwds = _xkwds(clas_kwds, eps=eps, Vector=clas, name=circum3.__name__)
        c, d = _tricenter3d2(p1, r, p2, r, p3, r, useZ=useZ, dLL3=dLL3, **kwds)
    else:  # Meeus' Type I
        c, d = d, None
    return Circum3Tuple(r, c, d)


def circum4_(*points, **Vector_and_kwds):
    '''Best-fit a sphere through three or more (3-D) points.

       @arg points: The points (each a C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    or C{Vector4Tuple}).
       @kwarg Vector_and_kwds: Optional class C{B{Vector}=None} to return the center
                               and optional, additional B{C{Vector}} keyword arguments,
                               otherwise the first B{C{points}}' (sub-)class.

       @return: L{Circum4Tuple}C{(radius, center, rank, residuals)} with C{center} an
                instance of C{B{points}[0]}' (sub-)class or B{C{Vector}} if specified.

       @raise ImportError: Package C{numpy} not found, not installed or older than
                           version 1.10.

       @raise NumPyError: Some C{numpy} issue.

       @raise PointsError: Too few B{C{points}}.

       @raise TypeError: One of the B{C{points}} is invalid.

       @see: U{Charles F. Jekel, "Least Squares Sphere Fit", Sep 13, 2015
             <https://Jekel.me/2015/Least-Squares-Sphere-Fit/>} and U{Appendix A
             <https://hdl.handle.net/10019.1/98627>}, U{numpy.linalg.lstsq
             <https://NumPy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html>},
             U{Eberly 6<https://www.sci.Utah.EDU/~balling/FEtools/doc_files/LeastSquaresFitting.pdf>}
             and functions L{pygeodesy.circum3} and L{pygeodesy.meeus2}.
    '''
    n, ps = len2(points)
    if n < 3:
        raise PointsError(points=n, txt=_too_(_few_))

    A, b = [], []
    for i, p in enumerate(ps):
        v = _otherV3d(useZ=True, i=i, points=p)
        A.append(v.times(_2_0).xyz + _1_0_T)
        b.append(v.length2)

    with _numpy(n, circum4_) as np:
        A = np.array(A).reshape((n, 4))
        b = np.array(b).reshape((n, 1))
        C, R, rk, _ = np.linalg.lstsq(A, b, rcond=None)  # to silence warning
        C = map1(float, *C)
        R = map1(float, *R)  # empty if rk < 4 or n <= 4

    n = circum4_.__name__
    c = Vector3d(*C[:3], name=n)
    r = Radius(sqrt(fsum_(C[3], *c.x2y2z2)), name=n)

    c = _nVc(c, **_xkwds(Vector_and_kwds, clas=ps[0].classof, name=n))
    return Circum4Tuple(r, c, rk, R)


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
            # like _nearestOn(s1, s2, e2, within=False, eps=e)
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


def iscolinearWith(point, point1, point2, eps=EPS):
    '''Check whether a point is colinear with two other (3-D) points.

       @arg point: The point (L{Vector3d}, C{Vector3Tuple} or
                              C{Vector4Tuple}).
       @arg point1: First point (L{Vector3d}, C{Vector3Tuple}
                                 or C{Vector4Tuple}).
       @arg point2: Second point (L{Vector3d}, C{Vector3Tuple}
                                  or C{Vector4Tuple}).
       @kwarg eps: Tolerance (C{scalar}), same units as C{x},
                   C{y} and C{z}.

       @return: C{True} if B{C{point}} is colinear B{C{point1}}
                and B{C{point2}}, C{False} otherwise.

       @raise TypeError: Invalid B{C{point}}, B{C{point1}} or
                         B{C{point2}}.

       @see: Function L{vector3d.nearestOn}.
    '''
    return _iscolinearWith(_otherV3d(point=point),
                            point1, point2, eps=eps)


def _iscolinearWith(v, p1, p2, eps=EPS):
    # (INTERNAL) Check colinear, see L{isColinear} above,
    # separated to allow callers to embellish any exceptions
    v1 = _otherV3d(point1=p1)
    v2 = _otherV3d(point2=p2)
    n  = _nearestOn(v, v1, v2, within=False, eps=eps)
    return n is v1 or n.minus(v).length2 < eps


def meeus2(point1, point2, point3, circum=False, useZ=True):
    '''Return the radius and I{Meeus}' Type of the smallest circle I{through}
       or I{containing} three (3-D) points.

       @arg point1: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @kwarg circum: If C{True} return the non-zero C{circumradius} always,
                      ignoring the I{Meeus}' Type I case (C{bool}).
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=0} (C{bool}).

       @return: L{Meeus2Tuple}C{(radius, Type)}.

       @raise IntersectionError: Near-coincident or colinear points, iff C{B{circum}=True}.

       @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{point3}}.

       @see: U{Jean Meeus, "Astronomical Algorithms", 2nd Ed. 1998, page 127ff
             <http://www.Agopax.IT/Libri_astronomia/pdf/Astronomical%20Algorithms.pdf>},
             U{circumradius<https://MathWorld.Wolfram.com/Circumradius.html>},
             U{circumcircle<https://MathWorld.Wolfram.com/Circumcircle.html>} and
             functions L{pygeodesy.circum3} and L{pygeodesy.circum4_}.
    '''
    try:
        A = _otherV3d(useZ=useZ, point1=point1)
        r, t, _, _ = _meeus4(A, point2, point3, circum=circum, useZ=useZ, clas=point1.classof)
    except (TypeError, ValueError) as x:
        raise _xError(x, point1=point1, point2=point2, point3=point3, circum=circum)
    return Meeus2Tuple(r, t)


def _meeus4(A, point2, point3, circum=False, useZ=True, clas=None, **clas_kwds):
    # (INTERNAL) Radius and Meeus' Type
    B = p2 = _otherV3d(useZ=useZ, point2=point2)
    C = p3 = _otherV3d(useZ=useZ, point3=point3)

    a = B.minus(C).length2
    b = C.minus(A).length2
    c = A.minus(B).length2
    if a < b:
        a, b, A, B = b, a, B, A
    if a < c:
        a, c, A, C = c, a, C, A

    if a > EPS02 and (circum or a < (b + c)):  # circumradius
        b = sqrt(b / a)
        c = sqrt(c / a)
        r = fsum1_(_1_0, b, c) * fsum1_(_1_0, b, -c) * fsum1_(-_1_0, b, c) * fsum1_(_1_0, -b, c)
        if r < EPS02:
            raise IntersectionError(_coincident_ if b < EPS0 or c < EPS0 else (
                                    _colinear_ if _iscolinearWith(A, B, C) else _invalid_))
        r = sqrt(a / r) * b * c
        t = None  # Meeus' Type II
    else:  # obtuse or right angle
        r = 0 if a < EPS02 else sqrt(a) * _0_5
        t = B.plus(C).times(_0_5)  # Meeus' Type I
        if clas is not None:
            t = clas(t.x, t.y, t.z, **_xkwds(clas_kwds, name=meeus2.__name__))
    return r, t, p2, p3


def nearestOn(point, point1, point2, within=True, Vector=None, **Vector_kwds):
    '''Locate the point between two points closest to a reference (3-D).

       @arg point: Reference point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}
                                    or C{Vector4Tuple}).
       @arg point1: Start point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple} or
                                 C{Vector4Tuple}).
       @arg point2: End point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple} or
                               C{Vector4Tuple}).
       @kwarg within: If C{True} return the closest point between both given
                      points, otherwise the closest point on the extended line
                      through both points (C{bool}).
       @kwarg Vector: Class to return closest point (C{Cartesian}, L{Vector3d}
                      or C{Vector3Tuple}) or C{None}.
       @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword arguments,
                           ignored if C{B{Vector} is None}.

       @return: Closest point, either B{C{point1}} or B{C{point2}} or an instance
                of the B{C{point}}'s (sub-)class or B{C{Vector}} if not C{None}.

       @raise TypeError: Invalid B{C{point}}, B{C{point1}} or B{C{point2}}.

       @see: U{3-D Point-Linedistance<https://MathWorld.Wolfram.com/Point-LineDistance3-Dimensional.html>},
             C{Cartesian} and C{LatLon} methods C{nearestOn}, method L{sphericalTrigonometry.LatLon.nearestOn3}
             and function L{sphericalTrigonometry.nearestOn3}.
    '''
    p0 = _otherV3d(point =point)
    p1 = _otherV3d(point1=point1)
    p2 = _otherV3d(point2=point2)

    p = _nearestOn(p0, p1, p2, within=within)
    if Vector is not None:
        p = Vector(p.x, p.y, p.z, **_xkwds(Vector_kwds, name=nearestOn.__name__))
    elif p is p1:
        p = point1
    elif p is p2:
        p = point2
    else:  # ignore Vector_kwds
        p = point.classof(p.x, p.y, p.z, name=nearestOn.__name__)
    return p


def _nearestOn(p0, p1, p2, within=True, eps=EPS):
    # (INTERNAL) Get closest point, see L{nearestOn} above,
    # separated to allow callers to embellish any exceptions
    p21 = p2.minus(p1)
    d2 = p21.length2
    if d2 < eps:  # coincident
        p = p1  # ... or point2
    else:
        t = p0.minus(p1).dot(p21) / d2
        p = p1 if (within and t < eps) else (
            p2 if (within and t > (_1_0 - eps)) else
            p1.plus(p21.times(t)))
    return p


def _null_space2(numpy, A, eps):
    # (INTERNAL) Return the nullspace and rank of matrix A
    # @see: <https://SciPy-Cookbook.ReadTheDocs.io/items/RankNullspace.html>,
    # <https://NumPy.org/doc/stable/reference/generated/numpy.linalg.svd.html>,
    # <https://StackOverflow.com/questions/19820921>,
    # <https://StackOverflow.com/questions/2992947> and
    # <https://StackOverflow.com/questions/5889142>
    A = numpy.array(A)
    m = max(numpy.shape(A))
    if m != 4:  # for this usage
        raise _AssertionError(shape=m, txt=_null_space2.__name__)
    # if needed, square A, pad with zeros
    A = numpy.resize(A, m * m).reshape(m, m)
#   try:  # no numpy.linalg.null_space <https://docs.SciPy.org/doc/>
#       return scipy.linalg.null_space(A)  # XXX no scipy.linalg?
#   except AttributeError:
#       pass
    _, s, v = numpy.linalg.svd(A)
    t = max(eps, eps * s[0])  # tol, s[0] is largest singular
    r = numpy.sum(s > t)  # rank
    if r == 3:  # get null_space
        n = numpy.transpose(v[r:])
        s = numpy.shape(n)
        if s != (m, 1):  # bad null_space shape
            raise _AssertionError(shape=s, txt=_null_space2.__name__)
        e = float(numpy.max(numpy.abs(numpy.dot(A, n))))
        if e > t:  # residual not near-zero
            raise _AssertionError(res=e, tol=t, txt=_null_space2.__name__)
    else:  # coincident, colinear, concentric centers, ambiguous, etc.
        n = None
    # del A, s, vh  # release numpy
    return n, r


@contextmanager  # <https://www.python.org/dev/peps/pep-0343/> Examples
def _numpy(arg, where):
    # get numpy with any errors raised as NumPyError
    np = Vector3d._numpy
    if np is None:  # get numpy, once or ImportError
        Vector3d._numpy = np = _xnumpy(where, 1, 10)  # macOS' Python 2.7 numpy 1.8 OK
    try:  # <https://NumPy.org/doc/stable/reference/generated/numpy.seterr.html>
        e = np.seterr(all=_raise_)  # throw FloatingPointError for numpy errors
        yield np
    except Exception as x:  # mostly FloatingPointError?
        raise NumPyError(x.__class__.__name__, arg, txt=str(x))
    finally:  # restore numpy error handling
        np.seterr(**e)


def _nVc(v, clas=None, name=NN, Vector=None, **Vector_kwds):
    # return a named C{Vector} or C{clas} instance
    if Vector is not None:
        v = Vector(v.x, v.y, v.z, **Vector_kwds)
    elif clas is not None:
        v = clas(v.x, v.y, v.z)  # ignore Vector_kwds
    return _xnamed(v, name) if name else v


def _otherV3d(useZ=True, NN_OK=True, i=None, **name_v):  # in .CartesianEllipsoidalBase.intersections2, .Ellipsoid.height4, .formy.hartzell
    # check named vector instance, return Vector3d
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

       @return: New B{C{Vector}} or if B{C{Vector}} is C{None},
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
    return _xnamed((Vector3Tuple(*v) if Vector is None else
                    Vector(*v, **Vector_kwds)), parse3d.__name__)


def radii11(point1, point2, point3, useZ=True):
    '''Return the radii of the C{In-}, I{Soddy} and C{Tangent} circles of a
       (2- or 3-D) triangle.

       @arg point1: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=0} (C{bool}).

       @return: L{Radii11Tuple}C{(rA, rB, rC, cR, rIn, riS, roS, a, b, c, s)}.

       @raise IntersectionError: Near-coincident or colinear points.

       @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{point3}}.

       @see: U{Circumradius<https://MathWorld.Wolfram.com/Circumradius.html>},
             U{Incircle<https://MathWorld.Wolfram.com/Incircle.html>}, U{Soddy
             Circles<https://MathWorld.Wolfram.com/SoddyCircles.html>} and
             U{Tangent Circles<https://MathWorld.Wolfram.com/TangentCircles.html>}.
    '''
    try:
        return _radii11ABC(point1, point2, point3, useZ=useZ)[0]
    except (TypeError, ValueError) as x:
        raise _xError(x, point1=point1, point2=point2, point3=point3)


def _radii11ABC(point1, point2, point3, useZ=True):
    # (INTERNAL) Tangent, Circum, Incircle, Soddy radii, sides and semi-perimeter
    A = _otherV3d(useZ=useZ, point1=point1, NN_OK=False)
    B = _otherV3d(useZ=useZ, point2=point2, NN_OK=False)
    C = _otherV3d(useZ=useZ, point3=point3, NN_OK=False)

    a = B.minus(C).length
    b = C.minus(A).length
    c = A.minus(B).length

    s = fsum1_(a, b, c) * _0_5  # semi-perimeter
    if s > EPS0:
        rs = (s - a), (s - b), (s - c)
        r3, r2, r1 = sorted(rs)  # r3 <= r2 <= r1
        if r3 > EPS0:  # and r2 > EPS0 and r1 > EPS0
            r3_r1 = r3 / r1
            r3_r2 = r3 / r2
            # t = r1 * r2 * r3 * (r1 + r2 + r3)
            #   = r1**2 * r2 * r3 * (1 + r2 / r1 + r3 / r1)
            #   = (r1 * r2)**2 * (r3 / r2) * (1 + r2 / r1 + r3 / r1)
            t = r3_r2 * fsum1_(_1_0, r2 / r1, r3_r1)  # * (r1 * r2)**2
            if t > EPS02:
                t = sqrt(t) * _2_0  # * r1 * r2
                # d = r1 * r2 + r2 * r3 + r3 * r1
                #   = r1 * (r2 + r2 * r3 / r1 + r3)
                #   = r1 * r2 * (1 + r3 / r1 + r3 / r2)
                d = fsum1_(_1_0, r3_r1, r3_r2)  # * r1 * r2
                # si/o = r1 * r2 * r3 / (r1 * r2 * (d +/- t))
                #      = r3 / (d +/- t)
                si =  r3 / (d + t)
                so = (r3 / (d - t)) if d > t else INF
                # ci = sqrt(r1 * r2 * r3 / s)
                #    = r1 * sqrt(r2 * r3 / r1 / s)
                ci = r1 * sqrt(r2 * r3_r1 / s)
                # co = a * b * c / (4 * ci * s)
                t  = _4_0 * s * ci
                co = (a * b * c / t) if t > EPS0 else INF
                r1, r2, r3 = rs  # original order
                t = Radii11Tuple(r1, r2, r3, co, ci, si, so, a, b, c, s)
                return t, A, B, C

    raise IntersectionError(_near_(_coincident_) if min(a, b, c) < EPS0 else (
                            _colinear_ if _iscolinearWith(A, B, C) else _invalid_))


def soddy4(point1, point2, point3, eps=EPS4, useZ=True):
    '''Return the radius and center of the C{inner} I{Soddy} circle of a
       (2- or 3-D) triangle.

       @arg point1: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @kwarg eps: Tolerance for function L{trilaterate3d2} if C{B{useZ} is True}
                   else L{trilaterate2d2}.
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=0} (C{bool}).

       @return: L{Soddy4Tuple}C{(radius, center, deltas, outer)}.  The C{center},
                an instance of B{C{point1}}'s (sub-)class, is co-planar with the
                three given points.

       @raise ImportError: Package C{numpy} not found, not installed or older
                           than version 1.10 and C{B{useZ} is True}.

       @raise IntersectionError: Near-coincident or colinear points or
                                 a trilateration or C{numpy} issue.

       @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{point3}}.

       @see: Functions L{radii11} and L{circum3}.
    '''
    t, p1, p2, p3 = _radii11ABC(point1, point2, point3, useZ=useZ)

    r    =  t.riS
    c, d = _tricenter3d2(p1, t.rA + r,
                         p2, t.rB + r,
                         p3, t.rC + r, eps=eps, useZ=useZ,
                                       Vector=point1.classof, name=soddy4.__name__)
    return Soddy4Tuple(r, c, d, t.roS)


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


def _tricenter3d2(p1, r1, p2, r2, p3, r3, eps=EPS4, useZ=True, dLL3=False, **kwds):
    # (INTERNAL) Trilaterate and disambiguate the 3-D center
    d, kwds = None, _xkwds(kwds, eps=eps, coin=True)
    if useZ and p1.z != p2.z != p3.z:  # ignore z if all match
        a, b = _trilaterate3d2(p1, r1, p2, r2, p3, r3, **kwds)
        if a is b:  # no unambiguity
            c = a  # == b
        else:
            c = a.plus(b).times(_0_5)  # mean
            if not a.isconjugateTo(b, minum=0, eps=eps):
                if dLL3:  # deltas as (lat, lon, height)
                    a = a.toLatLon()
                    b = b.toLatLon()
                    d = LatLon3Tuple(b.lat    - a.lat,
                                     b.lon    - a.lon,
                                     b.height - a.height, name=_deltas_)
                else:
                    d = b.minus(a)  # vectorial deltas
    else:
        if useZ:  # pass z to Vector if given
            kwds = _xkwds(kwds, z=p1.z)
        c = _trilaterate2d2(p1.x, p1.y, r1,
                            p2.x, p2.y, r2,
                            p3.x, p3.y, r3, **kwds)
    return c, d


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

       @raise IntersectionError: No intersection, colinear or near-concentric
                                 centers, trilateration failed some other way
                                 or the trilaterated point is off one circle
                                 by more than B{C{eps}}.

       @raise UnitError: Invalid B{C{radius1}}, B{C{radius2}} or B{C{radius3}}.

       @see: U{Issue #49<https://GitHub.com/mrJean1/PyGeodesy/issues/49>},
             U{Find X location using 3 known (X,Y) location using trilateration
             <https://math.StackExchange.com/questions/884807>} and L{trilaterate3d2}.
    '''
    return _trilaterate2d2(x1, y1, radius1, x2, y2, radius2, x3, y3, radius3,
                                            eps=eps, **Vector_and_kwds)


def _trilaterate2d2(x1, y1, radius1, x2, y2, radius2, x3, y3, radius3,
                                             coin=False, eps=None,
                                             Vector=None, **Vector_kwds):
    # (INTERNAL) Trilaterate three circles, see L{trilaterate2d2} above

    def _abct4(x1, y1, r1, x2, y2, r2, coin):
        a =  x2 - x1
        b =  y2 - y1
        t = _trinear2far(r1, r2, hypot(a, b), coin)
        c = _0_0 if t else (hypot2_(r1, x2, y2) - hypot2_(r2, x1, y1))
        return a, b, c, t

    def _astr(**kwds):  # kwds as (name=value, ...) strings
        return Fmt.PAREN(_COMMASPACE_(*(Fmt.EQUAL(*t) for t in kwds.items())))

    r1 = Radius_(radius1=radius1)
    r2 = Radius_(radius2=radius2)
    r3 = Radius_(radius3=radius3)

    a, b, c, t = _abct4(x1, y1, r1, x2, y2, r2, coin)
    if t:
        raise IntersectionError(_and(_astr(x1=x1, y1=y1, radius1=r1),
                                     _astr(x2=x2, y2=y2, radius2=r2)), txt=t)

    d, e, f, t = _abct4(x2, y2, r2, x3, y3, r3, coin)
    if t:
        raise IntersectionError(_and(_astr(x2=x2, y2=y2, radius2=r2),
                                     _astr(x3=x3, y3=y3, radius3=r3)), txt=t)

    _, _, _, t = _abct4(x3, y3, r3, x1, y1, r1, coin)
    if t:
        raise IntersectionError(_and(_astr(x3=x3, y3=y3, radius3=r3),
                                     _astr(x1=x1, y1=y1, radius1=r1)), txt=t)

    q = (a * e - b * d) * _2_0
    if isnear0(q):
        t = _no_(_intersection_)
        raise IntersectionError(_and(_astr(x1=x1, y1=y1, radius1=r1),
                                     _astr(x2=x2, y2=y2, radius2=r2),
                                     _astr(x3=x3, y3=y3, radius3=r3)), txt=t)
    t = Vector2Tuple((c * e - b * f) / q,
                     (a * f - c * d) / q, name=trilaterate2d2.__name__)

    if eps and eps > 0:
        for x, y, r in ((x1, y1, r1), (x2, y2, r2), (x3, y3, r3)):
            d = hypot(x - t.x, y - t.y)
            e = abs(d - r)
            if e > eps:
                t = _and(Float(delta=e).toRepr(), r.toRepr(),
                         Float(distance=d).toRepr(), t.toRepr())
                raise IntersectionError(t, txt=Fmt.exceeds_eps(eps))

    if Vector is not None:
        t = Vector(t.x, t.y, **_xkwds(Vector_kwds, name=t.name))
    return t


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

       @raise IntersectionError: Near-concentric, colinear, too distant or
                                 non-intersecting spheres.

       @raise NumPyError: Some C{numpy} issue.

       @raise TypeError: Invalid B{C{center1}}, B{C{center2}} or B{C{center3}}.

       @raise UnitError: Invalid B{C{radius1}}, B{C{radius2}} or B{C{radius3}}.

       @see: Norrdine, A. U{I{An Algebraic Solution to the Multilateration
             Problem}<https://www.ResearchGate.net/publication/
             275027725_An_Algebraic_Solution_to_the_Multilateration_Problem>}
             and U{I{implementation}<https://www.ResearchGate.net/publication/
             288825016_Trilateration_Matlab_Code>} and L{trilaterate2d2}.
    '''
    try:
        return _trilaterate3d2(_otherV3d(center1=center1, NN_OK=False),
                                Radius_(radius1=radius1, low=eps),
                                center2, radius2, center3, radius3, eps=eps,
                                clas=center1.classof, **Vector_and_kwds)
    except (AssertionError, NumPyError, TypeError, ValueError) as x:
        raise _xError(x, center1=center1, radius1=radius1,
                         center2=center2, radius2=radius2,
                         center3=center3, radius3=radius3)


def _trilaterate3d2(c1, r1, c2, r2, c3, r3, eps=EPS, coin=False,  # MCCABE 14
                                          **clas_Vector_and_kwds):
    # (INTERNAL) Intersect three spheres or circles, see L{trilaterate3d2}
    # above, separated to allow callers to embellish any exceptions, like
    # C{FloatingPointError}s from C{numpy}

    def _F3d2(F):
        # map numpy 4-vector to floats tuple and Vector3d
        T = map2(float, F)
        return T, Vector3d(*T[1:])

    def _N3(t01, x, z):
        # compute x, y and z and return as Vector
        v = x.plus(z.times(t01))
        n = trilaterate3d2.__name__
        return _nVc(v, **_xkwds(clas_Vector_and_kwds, name=n))

    def _perturbe5(eps, r):
        # perturbe radii to handle corner cases like this
        # <https://GitHub.com/mrJean1/PyGeodesy/issues/49>
        yield _0_0
        if eps and eps > 0:
            p = max(eps,  EPS)
            yield  p
            yield -min(p, r)
            q = max(eps, _EPS4e8)
            if q > p:
                yield  q
                yield -min(q, r)

    def _roots(numpy, *coeffs):
        # only real, non-complex roots of a polynomial, if any
        rs = numpy.polynomial.polynomial.polyroots(coeffs)
        return tuple(float(r) for r in rs if not numpy.iscomplex(r))

    c2 = _otherV3d(center2=c2, NN_OK=False)
    c3 = _otherV3d(center3=c3, NN_OK=False)
    rs = [r1, Radius_(radius2=r2, low=eps),
              Radius_(radius3=r3, low=eps)]

    # get null_space Z, pseudo-inverse A and vector B, once
    A = [(_1_0_T + c.times(-_2_0).xyz) for c in (c1, c2, c3)]  # 3 x 4
    with _numpy(None, trilaterate3d2) as np:
        Z, _ = _null_space2(np, A, eps)
        A    =  np.linalg.pinv(A)  # Moore-Penrose pseudo-inverse
    if Z is None:  # coincident, colinear, concentric, etc.
        raise _trilaterror(c1, r1, c2, r2, c3, r3, eps, coin)
    Z, z = _F3d2(Z)
    z2 =  z.length2
    bs = [c.length2 for c in (c1, c2, c3)]

    for p in _perturbe5(eps, min(rs)):
        b = [((r + p)**2 - b) for r, b in zip(rs, bs)]  # 1 x 3 or 3 x 1
        with _numpy(p, trilaterate3d2) as np:
            X, x = _F3d2(np.dot(A, b))
            # quadratic polynomial coefficients, ordered (^0, ^1, ^2)
            t = _roots(np, fdot(X, -_1_0, *x.xyz),
                          (fdot(Z, -_0_5, *x.xyz) * _2_0), z2)
            if t:
                break
    else:  # coincident, concentric, colinear, too distant, no intersection, etc.
        raise _trilaterror(c1, r1, c2, r2, c3, r3, eps, coin)

    v = _N3(t[0], x, z)
    if len(t) < 2:  # one intersection
        t = v, v
    elif abs(t[0] - t[1]) < eps:  # abutting
        t = v, v
    else:  # "lowest" intersection first (to avoid test failures)
        u = _N3(t[1], x, z)
        t = (u, v) if u.x < v.x else (v, u)
    return t


def _trilaterror(c1, r1, c2, r2, c3, r3, eps, coin):
    # return IntersectionError with the cause of the error

    def _no_intersection():
        t = _no_(_intersection_)
        if coin:
            def _reprs(*crs):
                return _and(*map(repr, crs))

            r =  repr(r1) if r1 == r2 == r3 else _reprs(r1, r2, r3)
            t = _SPACE_(t, _of_, _reprs(c1, c2, c3), _with_, _radius_, r)
        return t

    def _txt(c1, r1, c2, r2):
        t = _trinear2far(r1, r2, c1.minus(c2).length, coin)
        return _SPACE_(c1.name, _and_, c2.name, t) if t else t

    t = _txt(c1, r1, c2, r2) or \
        _txt(c1, r1, c3, r3) or \
        _txt(c2, r2, c3, r3) or (
        _colinear_ if _iscolinearWith(c1, c2, c3, eps=eps) else
        _no_intersection())
    return IntersectionError(t, txt=None)


def _trinear2far(r1, r2, h, coin):
    # check for near-coincident/-concentric or too distant spheres/circles
    return _too_(Fmt.distant(h)) if h > (r1 + r2) else (_near_(
           _coincident_ if coin else _concentric_) if h < abs(r1 - r2) else NN)


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
# Copyright (C) 2016-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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
