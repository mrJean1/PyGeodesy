
# -*- coding: utf-8 -*-

u'''2- or 3-D vectorial functions L{circin6}, L{circum3}, L{circum4_},
L{iscolinearWith}, L{meeus2}, L{nearestOn}, L{radii11} and L{soddy4}.
'''

from pygeodesy.basics import len2, map2, _xnumpy
from pygeodesy.constants import EPS, EPS0, EPS02, EPS4, INF, INT0, \
                               _EPS4e8, isnear0, _0_0, _0_25, _0_5, _N_0_5, \
                               _1_0, _1_0_1T, _N_1_0, _2_0, _N_2_0, _4_0
from pygeodesy.errors import _and, _AssertionError, IntersectionError, NumPyError, \
                              PointsError, TriangleError, _xError, _xkwds
from pygeodesy.fmath import fabs, fdot, hypot, hypot2_, sqrt
from pygeodesy.fsums import Fsum, fsumf_, fsum1f_
from pygeodesy.interns import NN, _a_, _and_, _b_, _c_, _center_, _coincident_, \
                             _colinear_, _concentric_, _COMMASPACE_, _few_, \
                             _intersection_, _invalid_, _near_, _no_, _of_, \
                             _radius_, _rIn_, _s_, _SPACE_, _too_, _with_
# from pygeodesy.lazily import _ALL_LAZY  # from .named
from pygeodesy.named import _ALL_LAZY, _NamedTuple, _Pass, Property_RO
from pygeodesy.namedTuples import LatLon3Tuple, Vector2Tuple
# from pygeodesy.props import Property_RO  # from .named
from pygeodesy.streprs import Fmt, unstr
from pygeodesy.units import Float, Int, Meter, Radius, Radius_
from pygeodesy.vector3d import iscolinearWith, nearestOn, _nearestOn2, _nVc, _otherV3d, \
                               trilaterate2d2, trilaterate3d2, Vector3d  # PYCHOK unused

from contextlib import contextmanager
# from math import fabs, sqrt  # from .fmath

__all__ = _ALL_LAZY.vector2d
__version__ = '23.06.08'

_cA_        = 'cA'
_cB_        = 'cB'
_cC_        = 'cC'
_deltas_    = 'deltas'
_outer_     = 'outer'
_raise_     = 'raise'  # PYCHOK used!
_rank_      = 'rank'
_residuals_ = 'residuals'
_Type_      = 'Type'


class Circin6Tuple(_NamedTuple):
    '''6-Tuple C{(radius, center, deltas, cA, cB, cC)} with the C{radius}, the
       trilaterated C{center} and contact points of the I{inscribed} aka I{In-
       circle} of a triangle.  The C{center} is I{un}ambiguous if C{deltas} is
       C{None}, otherwise C{center} is the mean and C{deltas} the differences of
       the L{pygeodesy.trilaterate3d2} results.  Contact points C{cA}, C{cB} and
       C{cC} are the points of tangency, aka the corners of the U{Contact Triangle
       <https://MathWorld.Wolfram.com/ContactTriangle.html>}.
    '''
    _Names_ = (_radius_, _center_, _deltas_, _cA_,  _cB_,  _cC_)
    _Units_ = ( Radius,  _Pass,    _Pass,    _Pass, _Pass, _Pass)


class Circum3Tuple(_NamedTuple):  # in .latlonBase
    '''3-Tuple C{(radius, center, deltas)} with the C{circumradius} and trilaterated
       C{circumcenter} of the C{circumcircle} through 3 points (aka {Meeus}' Type II
       circle) or the C{radius} and C{center} of the smallest I{Meeus}' Type I circle.
       The C{center} is I{un}ambiguous if C{deltas} is C{None}, otherwise C{center}
       is the mean and C{deltas} the differences of the L{pygeodesy.trilaterate3d2}
       results.
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
       and C{roS}, the sides C{a}, C{b} and C{c} and semi-perimeter C{s} of a triangle,
       all in C{meter} conventionally.

       @note: C{Circumradius} C{cR} and outer I{Soddy} radius C{roS} may be C{INF}.
    '''
    _Names_ = ('rA', 'rB', 'rC', 'cR', _rIn_, 'riS', 'roS', _a_, _b_, _c_, _s_)
    _Units_ = ( Meter,) * len(_Names_)


class Soddy4Tuple(_NamedTuple):
    '''4-Tuple C{(radius, center, deltas, outer)} with C{radius} and (trilaterated)
       C{center} of the I{inner} I{Soddy} circle and the radius of the C{outer}
       I{Soddy} circle.  The C{center} is I{un}ambiguous if C{deltas} is C{None},
       otherwise C{center} is the mean and C{deltas} the differences of the
       L{pygeodesy.trilaterate3d2} results.

       @note: The outer I{Soddy} radius C{outer} may be C{INF}.
    '''
    _Names_ = (_radius_, _center_, _deltas_, _outer_)
    _Units_ = ( Radius,  _Pass,    _Pass,     Radius)


def circin6(point1, point2, point3, eps=EPS4, useZ=True):
    '''Return the radius and center of the I{inscribed} aka I{Incircle} of
       a (2- or 3-D) triangle.

       @arg point1: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @kwarg eps: Tolerance for function L{pygeodesy.trilaterate3d2} if
                   C{B{useZ} is True} else L{pygeodesy.trilaterate2d2}.
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=INT0} (C{bool}).

       @return: L{Circin6Tuple}C{(radius, center, deltas, cA, cB, cC)}.  The
                C{center} and contact points C{cA}, C{cB} and C{cC}, each an
                instance of B{C{point1}}'s (sub-)class, are co-planar with
                the three given points.

       @raise ImportError: Package C{numpy} not found, not installed or older
                           than version 1.10 and C{B{useZ} is True}.

       @raise IntersectionError: Near-coincident or -colinear points or
                                 a trilateration or C{numpy} issue.

       @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{point3}}.

       @see: Functions L{radii11} and L{circum3}, U{Contact Triangle
             <https://MathWorld.Wolfram.com/ContactTriangle.html>} and
             U{Incircle<https://MathWorld.Wolfram.com/Incircle.html>}.
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
                       **_xkwds(Vector_kwds, Vector=V, name=circin6.__name__))
    return Circin6Tuple(r, c, d, cA, cB, cC)


def circum3(point1, point2, point3, circum=True, eps=EPS4, useZ=True):
    '''Return the radius and center of the smallest circle I{through} or
       I{containing} three (2- or 3-D) points.

       @arg point1: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple} or
                    C{Vector4Tuple}).
       @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple} or
                    C{Vector4Tuple}).
       @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple} or
                    C{Vector4Tuple}).
       @kwarg circum: If C{True} return the C{circumradius} and C{circumcenter}
                      always, ignoring the I{Meeus}' Type I case (C{bool}).
       @kwarg eps: Tolerance for function L{pygeodesy.trilaterate3d2} if C{B{useZ}
                   is True} else L{pygeodesy.trilaterate2d2}.
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=INT0} (C{bool}).

       @return: A L{Circum3Tuple}C{(radius, center, deltas)}.  The C{center}, an
                instance of B{C{point1}}'s (sub-)class, is co-planar with the three
                given points.

       @raise ImportError: Package C{numpy} not found, not installed or older
                           than version 1.10 and C{B{useZ} is True}.

       @raise IntersectionError: Near-coincident or -colinear points or
                                 a trilateration or C{numpy} issue.

       @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{point3}}.

       @see: Functions L{pygeodesy.circum4_} and L{pygeodesy.meeus2} and Meeus, J.
             U{I{Astronomical Algorithms}<http://www.Agopax.IT/Libri_astronomia/pdf/
             Astronomical%20Algorithms.pdf>}, 2nd Ed. 1998, page 127ff, U{circumradius
             <https://MathWorld.Wolfram.com/Circumradius.html>} and U{circumcircle
             <https://MathWorld.Wolfram.com/Circumcircle.html>}.
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


def circum4_(*points, **useZ_Vector_and_kwds):
    '''Best-fit a sphere through three or more (3-D) points.

       @arg points: The points (each a C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    or C{Vector4Tuple}).
       @kwarg useZ_Vector_and_kwds: Keyword arguments C{B{useZ}=True} (C{bool})
                   to use the Z components, otherwise force all C{z=INT0}, class
                   C{B{Vector}=None} to return the center point with optionally,
                   additional nB{C{Vector}} keyword arguments, otherwise the
                   first B{C{points}}' (sub-)class is used.

       @return: L{Circum4Tuple}C{(radius, center, rank, residuals)} with C{center} an
                instance of C{B{points}[0]}' (sub-)class or B{C{Vector}} if specified.

       @raise ImportError: Package C{numpy} not found, not installed or older than
                           version 1.10.

       @raise NumPyError: Some C{numpy} issue.

       @raise PointsError: Too few B{C{points}}.

       @raise TypeError: One of the B{C{points}} is invalid.

       @see: Functions L{pygeodesy.circum3} and L{pygeodesy.meeus2}, Jekel, Charles F. U{I{Least
             Squares Sphere Fit}<https://Jekel.me/2015/Least-Squares-Sphere-Fit/>} Sep 13, 2015,
             U{Appendix A<https://hdl.handle.net/10019.1/98627>}, U{numpy.linalg.lstsq<https://
             NumPy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html>} and U{Eberly 6
             <https://www.sci.Utah.EDU/~balling/FEtools/doc_files/LeastSquaresFitting.pdf>}.
    '''
    def _useZ_kwds(useZ=True, **kwds):
        return useZ, kwds

    n, ps = len2(points)
    if n < 3:
        raise PointsError(points=n, txt=_too_(_few_))
    useZ, kwds = _useZ_kwds(**useZ_Vector_and_kwds)

    A, b = [], []
    for i, p in enumerate(ps):
        v = _otherV3d(useZ=useZ, i=i, points=p)
        A.append(v.times(_2_0).xyz + _1_0_1T)
        b.append(v.length2)

    with _numpy(circum4_, n=n) as _np:
        A = _np.array(A).reshape((n, 4))
        b = _np.array(b).reshape((n, 1))
        C, R, rk, _ = _np.least_squares4(A, b, rcond=None)  # to silence warning
        C = map2(float, C)
        R = map2(float, R)  # empty if rk < 4 or n <= 4

    n = circum4_.__name__
    c = Vector3d(*C[:3], name=n)
    r = Radius(sqrt(fsumf_(C[3], *c.x2y2z2)), name=n)

    c = _nVc(c, **_xkwds(kwds, clas=ps[0].classof, name=n))
    return Circum4Tuple(r, c, rk, R)


def _iscolinearWith(p, point1, point2, eps=EPS, useZ=True):
    # (INTERNAL) Check colinear, see L{iscolinearWith} above,
    # separated to allow callers to embellish any exceptions
    p1 = _otherV3d(useZ=useZ, point1=point1)
    p2 = _otherV3d(useZ=useZ, point2=point2)
    n, _ = _nearestOn2(p, p1, p2, within=False, eps=eps)
    return n is p1 or n.minus(p).length2 < eps


def meeus2(point1, point2, point3, circum=False, useZ=True):
    '''Return the radius and I{Meeus}' Type of the smallest circle I{through}
       or I{containing} three (2- or 3-D) points.

       @arg point1: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @kwarg circum: If C{True} return the C{circumradius} and C{circumcenter}
                      always, overriding I{Meeus}' Type II case (C{bool}).
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=INT0} (C{bool}).

       @return: L{Meeus2Tuple}C{(radius, Type)}, with C{Type} the C{circumcenter}
                iff C{B{circum}=True}.

       @raise IntersectionError: Near-coincident or -colinear points, iff C{B{circum}=True}.

       @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{point3}}.

       @see: Functions L{pygeodesy.circum3} and L{pygeodesy.circum4_} and Meeus, J.
             U{I{Astronomical Algorithms}<http://www.Agopax.IT/Libri_astronomia/pdf/
             Astronomical%20Algorithms.pdf>}, 2nd Ed. 1998, page 127ff, U{circumradius
             <https://MathWorld.Wolfram.com/Circumradius.html>} and U{circumcircle
             <https://MathWorld.Wolfram.com/Circumcircle.html>}.
    '''
    try:
        A = _otherV3d(useZ=useZ, point1=point1)
        return _meeus2(A, point2, point3, circum, useZ=useZ, clas=point1.classof)
    except (TypeError, ValueError) as x:
        raise _xError(x, point1=point1, point2=point2, point3=point3, circum=circum)


def _meeus2(A, point2, point3, circum, useZ=True, **clas_and_kwds):  # in .vector3d
    # (INTERNAL) Radius and center or Meeus' Type
    f = _circum3 if circum else _meeus4
    t =  f(A, point2, point3, circum=circum, useZ=useZ, **clas_and_kwds)[:2]
    return Meeus2Tuple(t)


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
        R = Fsum(_1_0, b, c) * Fsum(_1_0, b, -c) * Fsum(_1_0, -b, c) * Fsum(_N_1_0, b, c)
        r = R.fover(a)
        if r < EPS02:
            raise IntersectionError(_coincident_ if b < EPS0 or c < EPS0 else (
                                    _colinear_ if _iscolinearWith(A, B, C) else _invalid_))
        r = b * c / sqrt(r)
        t = None  # Meeus' Type II
    else:  # obtuse or right angle at A
        r = sqrt(a * _0_25) if a > EPS02 else INT0
        t = B.plus(C).times(_0_5)  # Meeus' Type I
        if clas is not None:
            t = clas(t.x, t.y, t.z, **_xkwds(clas_kwds, name=meeus2.__name__))
    return r, t, p2, p3


class _numpy(object):  # see also .geodesicw._wargs
    '''(INTERNAL) Partial C{NumPy} wrapper.
    '''
    @contextmanager  # <https://www.python.org/dev/peps/pep-0343/> Examples
    def __call__(self, where, *args, **kwds):
        '''(INTERNAL) Yield self with any errors raised as L{NumPyError}.
        '''
        np = self.np
        try:  # <https://NumPy.org/doc/stable/reference/generated/numpy.seterr.html>
            e = np.seterr(all=_raise_)  # throw FloatingPointError for numpy errors
            yield self
        except Exception as x:  # mostly FloatingPointError?
            t = unstr(where, *args, **kwds)
            raise NumPyError(t, cause=x)  # _xError2?
        finally:  # restore numpy error handling
            np.seterr(**e)

    @Property_RO
    def array(self):
        return self.np.array

    @Property_RO
    def least_squares4(self):
        '''Linear least-squares function.
        '''
        return self.np.linalg.lstsq

    @Property_RO
    def np(self):
        '''Import numpy 1.10+ once.
        '''
        return _xnumpy(self.__class__, 1, 10)

    def null_space2(self, A, rcond=None):
        '''Return the C{null_space} and C{rank} of matrix B{C{A}}.

           @see: U{Source<https://docs.SciPy.org/doc/scipy/reference/generated/scipy.linalg.null_space.html>}
                 U{SciPY Cookbook<https://SciPy-Cookbook.ReadTheDocs.io/items/RankNullspace.html>}, U{here
                 <https://NumPy.org/doc/stable/reference/generated/numpy.linalg.svd.html>}, U{here
                 <https://StackOverflow.com/questions/19820921>}, U{here
                 <https://StackOverflow.com/questions/2992947>} and U{here
                 <https://StackOverflow.com/questions/5889142>}.
        '''
        def _Error(txt=self.null_space2.__name__, **kwds):
            return _AssertionError(txt=txt, **kwds)

        np = self.np
        A  = np.array(A)
        m  = max(A.shape)
        if m != 4:  # for this usage
            raise _Error(shape=m)
        # if needed, square A, pad with zeros
        A = np.resize(A, m * m).reshape(m, m)
#       try:  # no np.linalg.null_space <https://docs.SciPy.org/doc/>
#           Z = scipy.linalg.null_space(A)  # XXX no scipy.linalg?
#           return Z, ...
#       except AttributeError:
#           pass
        U, S, V = np.linalg.svd(A)
        s = max(EPS, rcond) if rcond else (EPS * max(U.shape[0], V.shape[1]))
        t = max(EPS, float(np.max(S) * s))  # abs_tol, rel_tol * largest singular
        r = int(np.sum(S > t))  # rank
        if r == 3:  # get null_space
            Z = np.transpose(V[r:])
            s = map2(int, Z.shape)
            if s != (m, 1):  # bad null_space shape
                raise _Error(shape=s, m=m)
            D = A.dot(Z)  # near-zeros-vector
            n = float(np.linalg.norm(D, INF))  # INF = max(fabs(D)), 2 = hypot_(*D)
            if n > t:  # largest exceed tol
                raise _Error(dot=tuple(D.ravel()), norm=n, tol=t)
        else:  # coincident, colinear, concentric centers, ambiguous, etc.
            Z = None
        # del A, S, U, V  # release numpy
        return Z, r

    @Property_RO
    def pseudo_inverse(self):
        '''Moore-Penrose pseudo-inverse function.
        '''
        return self.np.linalg.pinv

    def real_roots(self, *coeffs):
        '''Compute the real, non-complex roots of a polynomial.
        '''
        np = self.np
        rs = np.polynomial.polynomial.polyroots(coeffs)
        return tuple(float(r) for r in rs if not np.iscomplex(r))

_numpy = _numpy()  # PYCHOK singleton


def radii11(point1, point2, point3, useZ=True):
    '''Return the radii of the C{In-}, I{Soddy} and C{Tangent} circles of a
       (2- or 3-D) triangle.

       @arg point1: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=INT0} (C{bool}).

       @return: L{Radii11Tuple}C{(rA, rB, rC, cR, rIn, riS, roS, a, b, c, s)}.

       @raise TriangleError: Near-coincident or -colinear points.

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

    S = Fsum(a, b, c) * _0_5
    s = float(S)  # semi-perimeter
    if s > EPS0:
        rs = float(S - a), float(S - b), float(S - c)
        r3, r2, r1 = sorted(rs)  # r3 <= r2 <= r1
        if r3 > EPS0:  # and r2 > EPS0 and r1 > EPS0
            r3_r1 = r3 / r1
            r3_r2 = r3 / r2
            # t = r1 * r2 * r3 * (r1 + r2 + r3)
            #   = r1**2 * r2 * r3 * (1 + r2 / r1 + r3 / r1)
            #   = (r1 * r2)**2 * (r3 / r2) * (1 + r2 / r1 + r3 / r1)
            t = r3_r2 * fsum1f_(_1_0, r2 / r1, r3_r1)  # * (r1 * r2)**2
            if t > EPS02:
                t = sqrt(t) * _2_0  # * r1 * r2
                # d = r1 * r2 + r2 * r3 + r3 * r1
                #   = r1 * (r2 + r2 * r3 / r1 + r3)
                #   = r1 * r2 * (1 + r3 / r1 + r3 / r2)
                d = fsum1f_(_1_0, r3_r1, r3_r2)  # * r1 * r2
                # si/o = r1 * r2 * r3 / (r1 * r2 * (d +/- t))
                #      = r3 / (d +/- t)
                si =  r3 / (d + t)
                so = (r3 / (d - t)) if d > t else INF
                # ci = sqrt(r1 * r2 * r3 / s)
                #    = r1 * sqrt(r2 * r3 / r1 / s)
                ci = r1 * sqrt(r2 * r3_r1 / s)
                # co = a * b * c / (4 * ci * s)
                t  = ci * s * _4_0
                co = (a * b * c / t) if t > EPS0 else INF
                r1, r2, r3 = rs  # original order
                t = Radii11Tuple(r1, r2, r3, co, ci, si, so, a, b, c, s)
                return t, A, B, C

    raise TriangleError(_near_(_coincident_) if min(a, b, c) < EPS0 else (
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
       @kwarg eps: Tolerance for function L{pygeodesy.trilaterate3d2} if
                   C{B{useZ} is True} otherwise L{pygeodesy.trilaterate2d2}.
       @kwarg useZ: If C{True}, use the Z components, otherwise force C{z=INT0} (C{bool}).

       @return: L{Soddy4Tuple}C{(radius, center, deltas, outer)}.  The C{center},
                an instance of B{C{point1}}'s (sub-)class, is co-planar with the
                three given points.  The C{outer} I{Soddy} radius may be C{INF}.

       @raise ImportError: Package C{numpy} not found, not installed or older
                           than version 1.10 and C{B{useZ} is True}.

       @raise IntersectionError: Near-coincident or -colinear points or
                                 a trilateration or C{numpy} issue.

       @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{point3}}.

       @see: Functions L{radii11} and L{circum3} and U{Soddy Circles
             <https://MathWorld.Wolfram.com/SoddyCircles.html>}.
    '''
    t, p1, p2, p3 = _radii11ABC(point1, point2, point3, useZ=useZ)

    r    =  t.riS
    c, d = _tricenter3d2(p1, t.rA + r,
                         p2, t.rB + r,
                         p3, t.rC + r, eps=eps, useZ=useZ,
                                       Vector=point1.classof, name=soddy4.__name__)
    return Soddy4Tuple(r, c, d, t.roS)


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


def _trilaterate2d2(x1, y1, radius1, x2, y2, radius2, x3, y3, radius3,
                                             coin=False, eps=None,
                                             Vector=None, **Vector_kwds):
    # (INTERNAL) Trilaterate three circles, see L{pygeodesy.trilaterate2d2}

    def _abct4(x1, y1, r1, x2, y2, r2):
        a =  x2 - x1
        b =  y2 - y1
        t = _tri3near2far(r1, r2, hypot(a, b), coin)
        c = _0_0 if t else (hypot2_(r1, x2, y2) - hypot2_(r2, x1, y1))
        return a, b, c, t

    def _astr(**kwds):  # kwds as (name=value, ...) strings
        return Fmt.PAREN(_COMMASPACE_(*(Fmt.EQUALg(*t) for t in kwds.items())))

    r1 = Radius_(radius1=radius1)
    r2 = Radius_(radius2=radius2)
    r3 = Radius_(radius3=radius3)

    a, b, c, t = _abct4(x1, y1, r1, x2, y2, r2)
    if t:
        raise IntersectionError(_and(_astr(x1=x1, y1=y1, radius1=r1),
                                     _astr(x2=x2, y2=y2, radius2=r2)), txt=t)

    d, e, f, t = _abct4(x2, y2, r2, x3, y3, r3)
    if t:
        raise IntersectionError(_and(_astr(x2=x2, y2=y2, radius2=r2),
                                     _astr(x3=x3, y3=y3, radius3=r3)), txt=t)

    _, _, _, t = _abct4(x3, y3, r3, x1, y1, r1)
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

    if eps and eps > 0:  # check distances to center vs radius
        for x, y, r in ((x1, y1, r1), (x2, y2, r2), (x3, y3, r3)):
            d = hypot(x - t.x, y - t.y)
            e = fabs(d - r)
            if e > eps:
                t = _and(Float(delta=e).toRepr(), r.toRepr(),
                         Float(distance=d).toRepr(), t.toRepr())
                raise IntersectionError(t, txt=Fmt.exceeds_eps(eps))

    if Vector is not None:
        t = Vector(t.x, t.y, **_xkwds(Vector_kwds, name=t.name))
    return t


def _trilaterate3d2(c1, r1, c2, r2, c3, r3, eps=EPS, coin=False,
                                          **clas_Vector_and_kwds):
    # (INTERNAL) Intersect three spheres or circles, see function
    # L{pygeodesy.trilaterate3d2}, separated to allow callers to
    # embellish exceptions, like C{FloatingPointError}s from C{numpy}

    def _F3d2(F):
        # map numpy 4-vector to floats tuple and Vector3d
        T = map2(float, F)
        return T, Vector3d(*T[1:])

    def _N3(t01, x, z):
        # compute x, y and z and return as B{C{clas}} or B{C{Vector}}
        v = x.plus(z.times(t01))
        n = trilaterate3d2.__name__
        return _nVc(v, **_xkwds(clas_Vector_and_kwds, name=n))

    c2 = _otherV3d(center2=c2, NN_OK=False)
    c3 = _otherV3d(center3=c3, NN_OK=False)
    rs = (r1, Radius_(radius2=r2, low=EPS),
              Radius_(radius3=r3, low=EPS))

    # get matrix A[3 x 4], its pseudo-inverse and null_space Z
    A = [(_1_0_1T + c.times(_N_2_0).xyz) for c in (c1, c2, c3)]
    with _numpy(trilaterate3d2, A=A, eps=eps) as _np:
        Z, _ = _np.null_space2(A, eps)
        if Z is not None:
            Z, z = _F3d2(Z)  # [4 x 1]
            z2 =  z.length2
            A  = _np.pseudo_inverse(A)  # [4 x 3]
            bs = [c.length2 for c in (c1, c2, c3)]
            # perturbe radii and vector b slightly by eps and eps * 4
            for p in _tri5perturb(eps, min(rs)):
                b = [((r + p)**2 - b) for r, b in zip(rs, bs)]  # [3 x 1]
                X, x = _F3d2(A.dot(b))
                # quadratic polynomial, coefficients ordered (^0, ^1, ^2)
                t = _np.real_roots(fdot(X, _N_1_0, *x.xyz),
                                   fdot(Z, _N_0_5, *x.xyz) * _2_0, z2)
                if t:
                    v = _N3(t[0], x, z)
                    if len(t) < 2:  # one intersection
                        t = v, v
                    elif fabs(t[0] - t[1]) < eps:  # abutting
                        t = v, v
                    else:  # "lowest" intersection first (to avoid test failures)
                        u = _N3(t[1], x, z)
                        t = (u, v) if u.x < v.x else (v, u)
                    return t

    # coincident, concentric, colinear, too distant, no intersection:
    # create the explanation and and throw an IntersectionError

    def _no_intersection(coin):
        t = _no_(_intersection_)
        if coin:
            def _reprs(*crs):
                return _and(*map(repr, crs))

            r =  repr(r1) if r1 == r2 == r3 else _reprs(r1, r2, r3)
            t = _SPACE_(t, _of_, _reprs(c1, c2, c3), _with_, _radius_, r)
        elif Z is None:
            t = _COMMASPACE_(t, _no_(_numpy.null_space2.__name__))
        return t

    t = _tri4near2far(c1, r1, c2, r2, coin) or \
        _tri4near2far(c1, r1, c3, r3, coin) or \
        _tri4near2far(c2, r2, c3, r3, coin) or (
        _colinear_ if _iscolinearWith(c1, c2, c3, eps=eps) else
        _no_intersection(coin))
    raise IntersectionError(t, txt=None)


def _tri3near2far(r1, r2, h, coin):
    # check for near-coincident/-concentric or too distant spheres/circles
    return _too_(Fmt.distant(h)) if h > (r1 + r2) else (_near_(
           _coincident_ if coin else _concentric_) if h < fabs(r1 - r2) else NN)


def _tri4near2far(c1, r1, c2, r2, coin):
    # check for near-coincident/-concentric or too distant spheres/circles
    t = _tri3near2far(r1, r2, c1.minus(c2).length, coin)
    return _SPACE_(c1.name, _and_, c2.name, t) if t else NN


def _tri5perturb(eps, r):
    # perturb the radii to handle this corner case
    # <https://GitHub.com/mrJean1/PyGeodesy/issues/49>
    yield _0_0
    if eps and eps > 0:
        p = max(eps, EPS)
        yield  p
        m = min(p, r)
        yield -m
        q = max(eps * _4_0, _EPS4e8)
        if q > p:
            yield q
            q = min(q, r)
            if q > m:
                yield -q

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
