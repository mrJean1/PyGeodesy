
# -*- coding: utf-8 -*-

u'''Class C{Ellipse} for 2-D ellipse attributes, like perimeter, area, etc.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

# from pygeodesy.basics import islistuple  # _MODS
from pygeodesy.constants import EPS, EPS_2, INT0, NEG0, PI, PI_2, PI3_2, PI2, \
                               _0_0, _1_0, _4_0, _isfinite, _over, _1_over  # _N_1_0
from pygeodesy.constants import _0_5, _3_0, _10_0, MANT_DIG as _DIG53  # PYCHOK used!
# from pygeodesy.ellipsoids import Ellipsoid  # _MODS
from pygeodesy.errors import _ConvergenceError, _ValueError, _xkwds, _xkwds_pop2
from pygeodesy.fmath import euclid, fhorner, fmean_, hypot, polar2d
from pygeodesy.fsums import _fsum  # PYCHOK used!
from pygeodesy.internals import typename,  _DOT_, _UNDER_
# from pygeodesy.interns import _DOT_, _UNDER_  # from .internals
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _NamedBase,  unstr
from pygeodesy.props import Property_RO, property_RO, property_ROnce
# from pygeodesy.streprs import unstr  # from .named
# from pygeodesy.triaxials import Triaxial_ , TriaxialError # _MODS
from pygeodesy.units import Degrees, Meter, Meter2, Radians, Radius, Scalar
from pygeodesy.utily import atan2, atan2p, sincos2, sincos2d
# from pygeodesy.vector3d import Vector3d  # _MODS

from math import degrees, fabs, radians, sqrt
# import operator as _operator  # from .fmath

__all__ = _ALL_LAZY.ellipses
__version__ = '26.03.27'

_TOL53    =  sqrt(EPS_2)     # sqrt(pow(_0_5, _DIG53))
_TOL53_53 = _TOL53 / _DIG53  # "flat" b/a tolerance, 1.9e-10
# assert _DIG53 == 53


class Ellipse(_NamedBase):
    '''Class to compute various attributes of a 2-D ellipse.
    '''
#   _ab3   = (a, b, a * b)  # unordered
    _flat  =  False
    _maxit = _DIG53
#   _Pab4  = (r, P, a, b)  # a >= b, ordered
    _4_PI_ = _4_0 / PI - (14 / 11)
    _tEPS  =  None

    def __init__(self, a, b, **name):
        '''New L{Ellipse} with semi-axes B{C{a}} and B{C{b}}.

           The ellipse is C{oblate} if C{B{a} > B{b}}, C{prolate} if
           C{B{a} < B{b}}, C{circular} if C{B{a} == B{b}} and C{"flat"}
           if C{min(B{a}, B{b}) <<< max(B{a}, B{b})}.

           @arg a: X semi-axis length (C{meter}, conventionally).
           @arg b: Y semi-axis length (C{meter}, conventionally).

           @raise ValueError: Invalid B{C{a}} or B{C{b}}.

           @see: U{Ellipse<https://MathWorld.Wolfram.com/Ellipse.html>}.
        '''
        if name:
            self.name = name
        self._ab3 = a, b, (a * b)  # unordered

        r = a < b
        if r:  # prolate
            a, b = b, a
        if b < 0 or not _isfinite(a):  # PYCHOK no cover
            raise self._Error(None)
        if a > b:
            if _isFlat(a, b):
                self._flat = True
                P = _4_0 * a
#               b = _0_0
            else:  # pro-/oblate
                P = None
        else:  # circular
            P = a * PI2
#           b = a
        self._Pab4 = r, P, a, b  # ordered

    @Property_RO
    def a(self):
        '''Get semi-axis C{B{a}} of this ellipse (C{meter}, conventionally).
        '''
        a, _, _ = self._ab3
        return Meter(a=a)

    @Property_RO
    def apses2(self):
        '''Get 2-tuple C{(apoapsis, periapsis)} with the U{apo-<https://MathWorld.Wolfram.com/Apoapsis.html>}
           and U{periapsis<https://MathWorld.Wolfram.com/Periapsis.html>} of this ellipse, both C{meter}.
        '''
        _, _, a, p = self._Pab4
        c = self.c
        if c:  # a != p
            p  = a - c
            a += c
        return a, p

    def arc(self, deg2, deg1=0):
        '''Compute the length of U{elliptic arc<https://www.JohnDCook.com/blog/2022/11/02/elliptic-arc-length/>}
           C{(B{deg2} - B{deg1})}, both counter-clockwise from semi-axis B{C{a}} to B{C{b}} of the ellipse.

           @arg deg2: End angle of the elliptic arc (C{degrees}).
           @kwarg deg1: Start angle of the elliptic arc (C{degrees}).

           @return: Arc length, signed (C{meter}, conventionally).
        '''
        return self.arc_(radians(deg2), (radians(deg1) if deg1 else _0_0))

    def arc_(self, rad2, rad1=0):
        '''Compute the length of U{elliptic arc<https://www.JohnDCook.com/blog/2022/11/02/elliptic-arc-length/>}
           C{(B{rad2} - B{rad1})}, both counter-clockwise from semi-axis B{C{a}} to B{C{b}} of the ellipse.

           @arg rad2: End angle of the elliptic arc (C{radians}).
           @kwarg rad1: Start angle of the elliptic arc (C{radians}).

           @return: Arc length, signed (C{meter}, conventionally).
        '''
        r, R, a, _ = self._Pab4
        if R is None:
            _e =  self._ellipe or self._ellipE
            k  =  self.e2
            r  =  PI_2 if r else _0_0
            R  = _arc(_e, k, r + rad2)
            r += rad1
            if r:
                R -= _arc(_e, k, r)
        else:
            a = (rad2 - rad1) / PI2
        return Meter(arc=R * a)

    @Property_RO
    def area(self):
        '''Get the area of this ellipse (C{meter**2}, conventionally).
        '''
        _, _, ab = self._ab3
        return Meter2(area=ab * PI)

    @Property_RO
    def b(self):
        '''Get semi-axis C{B{b}} of this ellipse (C{meter}, conventionally).
        '''
        _, b, _ = self._ab3
        return Meter(b=b)

    @Property_RO
    def c(self):
        '''Get the U{linear eccentricity<https://WikiPedia.org/wiki/Ellipse#Linear_eccentricity>}
           C{c}, I{unsigned} (C{meter}, conventionally).
        '''
        return Meter(c=fabs(self.foci))

    @Property_RO
    def e(self):
        '''Get the eccentricity (C{scalar, 0 <= B{e} <= 1}).
        '''
        e2 = self.e2
        return Scalar(e=sqrt(e2) if 0 < e2 < 1 else e2)

    @Property_RO
    def e2(self):
        '''Get the eccentricity I{squared} (C{scalar, 0 <= B{e2} <= 1}).
        '''
        # C{e2} is aka C{k}, Elliptic C{k2} and SciPy's C{m}
        _, _, a, b = self._Pab4
        e2 = ((_1_0 - (b / a)**2) if 0 < b else _1_0) if b < a else _0_0
        return Scalar(e2=e2)

    @Property_RO
    def _Ek(self):
        '''(INTERNAL) Get the C{Elliptic(k)} instance.
        '''
        return _MODS.elliptic._Ek(self.e2)

    def _ellipE(self, k, phi=None):  # PYCHOK k
        '''(INTERNAL) Get the in-/complete integral of the 2nd kind.
        '''
        # assert k == self._Ek.k2
        return self._Ek.cE if phi is None else self._Ek.fE(phi)

    @property_ROnce
    def _ellipe(self):
        '''(INTERNAL) Wrap functions C{scipy.special.ellipe} and C{-.ellipeinc}, I{once}.
        '''
        try:
            from scipy.special import ellipe, ellipeinc

            def _ellipe(k, phi=None):
                r = ellipe(k) if phi is None else ellipeinc(phi, k)
                return float(r)

        except (AttributeError, ImportError):
            _ellipe = None
        return _ellipe  # overwrite property_ROnce

    def _Error(self, where, **cause):  # PYCHOK no cover
        '''(INTERNAL) Build an L{EllipseError}.
        '''
        t = self.named3
        u = unstr(t, a=self.a, b=self.b)
        if where:
            t =  typename(where, where)
            u = _DOT_(u, t)
        return EllipseError(u, **cause)

    @Property_RO
    def foci(self):
        '''Get the U{linear eccentricity<https://WikiPedia.org/wiki/Ellipse#Linear_eccentricity>},
           I{signed} (C{meter}, conventionally), C{positive} if this ellipse is oblate, C{negative}
           if prolate or C{0} if circular.  See also property L{Ellipse.c}.
        '''
        c = float(self.e)
        if c:
            r, _, a, _ = self._Pab4
            c  *= a
            if r:  # prolate
                c = -c
        return Meter(foci=c)  # signed

    @property_ROnce
    def _GKs(self):
        '''(INTERNAL) Compute the coefficients for property C{.perimeterGK}, I{once}.
        '''
        # U{numerators<https://OEIS.org/A056981>}, U{denominators<https://OEIS.org/A056982>}
        return (1, 1 / 4, 1 / 64, 1 / 256, 25 / 16384, 49 / 65536,
                441 / 1048576, 1089 / 4194304)  # overwrite property_ROnce

    def hartzell4(self, x, y, los=False):
        '''Compute the intersection of this ellipse with a Line-Of-Sight from Point-Of-View
           C{(B{x}, B{y})} I{outside} this ellipse.

           @kwarg los: Line-Of-Sight, I{direction} to the ellipse (L{Los}, L{Vector3d},
                       L{Vector2Tuple} or 2-tuple C{(dx, dy)}) or C{True} for the I{normal,
                       perpendicular, plumb} to this ellipse or C{False} or C{None} to
                       point to its center.

           @return: L{Vector4Tuple}C{(x, y, z, h)} with coordinates C{x}, C{y} and C{z=0}
                    of the intersection and C{h} the distance to "Point-Of-View" C{(B{x},
                    B{y})} I{along the} B{C{los}}, all in C{meter}, conventionally.

           @raise EllipseError: Invalid B{C{x}}, B{C{y}} or B{C{los}} or B{C{los}} points
                                outside or away from this ellipse.

           @see: Function L{hartzell4<triaxials.triaxial5.hartzell4>} for further details.
        '''
        V3d = _MODS.vector3d.Vector3d
        if los not in (True, False, None):
            try:
                los = V3d(los.x, los.y, 0)
            except (AttributeError, TypeError):
                if _MODS.basics.islistuple(los, minum=2):
                    los = V3d(*map(float, los[:2]))
        return self._triaxialX(self.hartzell4, V3d(x, y, 0), los=los)

    def height4(self, x, y, **normal_eps):
        '''Compute the projection on and distance to this ellipse from a point C{(B{x}, B{y})}
           in- or outside this ellipse.

           @kwarg normal_eps: With default C{B{normal}=True} the projection is I{perpendicular,
                         plumb} to this ellipse, otherwise C{radially} to its center (C{bool}).
                         Tolerance C{B{eps}=EPS} for root finding and validation (C{scalar}),
                         use a negative value to skip validation.

           @return: L{Vector4Tuple}C{(x, y, z, h)} with coordinates C{x}, C{y} and C{z=0}
                    of the projection on or the intersection with the ellipse and C{h} the
                    I{signed, normal distance} to the ellipse in C{meter}, conventionally.
                    Positive C{h} indicates, C{x} and/or C{y} are outside the ellipse,
                    negative C{h} means inside.

           @raise EllipseError: Invalid B{C{x}}, B{C{y}} or B{C{eps}}, no convergence in
                  root finding or validation failed.

           @see: Methods L{Ellipse.normal3d}, L{Ellipse.normal4} and function L{height4
                 <triaxials.triaxial5.height4>}.
        '''
        return self._triaxialX(self.height4, x, y, 0, **normal_eps)

    def _HGKs(self, h, maxit):
        '''(INTERNAL) Yield the terms for property C{.perimeterHGK}.
        '''
        s = t = _1_0
        yield s
        for u in range(-1, maxit * 2, 2):
            t *= u / (u + 3) * h
            t2 = t**2
            yield t2
            p  = s
            s += t2
            if s == p:  # 44 trips
                break
        else:  # PYCHOK no cover
            raise _ConvergenceError(maxit, s, p)

    @property_RO
    def isCircular(self):
        '''Is this ellipse circular? (C{bool})
        '''
        return self.a == self.b

    @property_RO
    def isFlat(self):
        '''Is this ellipse "flat", too pro-/oblate? (C{bool})
        '''
        return self._flat

    @property_RO
    def isOblate(self):
        '''Is this ellipse oblate (foci on semi-axis C{a})? (C{bool})
        '''
        return self.a > self.b

    @property_RO
    def isProlate(self):
        '''Is this ellipse prolate (foci on semi-axis C{b})? (C{bool})
        '''
        return self.a < self.b

    @Property_RO
    def lati(self):
        '''Get the U{semi-latus rectum<https://WikiPedia.org/wiki/Ellipse#Semi-latus_rectum>},
           I{signed} (C{meter}, conventionally), C{positive} if this ellipse is oblate or
           circular, C{0} if "flat" and oblate, C{negative} if prolate or C{NEG0} if "flat"
           and prolate.  See also property L{Ellipse.p}.
        '''
        r, _, a, p = self._Pab4
        if 0 < p < a:
            p *= p / a
        if r:
            p = (-p) if p else NEG0
        return Meter(lati=p)  # signed

    def normal3d(self, deg_x, y=None, **length):
        '''Get a 3-D vector I{perpendicular to} this ellipse from point C{(B{x}, B{y})}
           C{on} this ellipse or at C{B{deg} degrees} along this ellipse.

           @kwarg length: Optional, signed C{B{length}=1} in out-/inward direction
                          (C{scalar}).

           @return: A C{Vector3d(x_, y_, z_=0)} normalized to B{C{length}}, pointing
                    out- or inward for postive respectively negative B{C{length}}.

           @raise EllipseError: Invalid B{C{x}} and/or B{C{y}}.

           @see: Methods L{Ellipse.height4}, L{Ellipse.normal4}, L{Ellipse.sideOf} and
                 C{Triaxial_.normal3d}.
        '''
        return self._triaxialX(self.normal3d, *self._xy03(deg_x, y), **length)

    def normal4(self, deg_x, y=None, **height_normal):
        '''Compute a point at B{C{height}} above or below this ellipse point C{(B{x},
           B{y})} C{on} this ellipse or at C{B{deg} degrees} along this ellipse.

           @kwarg height_normal: The desired distance C{B{height}=0} in- or outside this
                         ellipse (C{meter}, conventionally) and C{B{normal}=True},  If
                         C{B{normal}=True}, the B{C{height}} is I{perpendicular, plumb}
                         to this ellipse, otherwise C{radially} to its center (C{bool}).

           @return: L{Vector4Tuple}C{(x, y, z, h)} with coordinates C{x}, C{y} and C{z=0}
                    and C{h} the I{signed, normal distance} to the ellipse in C{meter},
                    conventionally.  Positive C{h} indicates, C{x} and/or C{y} are outside
                    the ellipse, negative C{h} means inside.

           @raise EllipseError: Invalid B{C{x}} and/or B{C{y}}.

           @see: Methods L{Ellipse.height4}, L{Ellipse.normal3d}, L{Ellipse.sideOf} and
                 C{Triaxial_.normal4}.
        '''
        return self._triaxialX(self.normal4, *self._xy03(deg_x, y), **height_normal)

    @Property_RO
    def p(self):
        '''Get the U{semi-latus rectum<https://WikiPedia.org/wiki/Ellipse#Semi-latus_rectum>}
           C{p (aka B{𝓁}, script-small-l)}, I{unsigned} (C{meter}, conventionally).
        '''
        return Meter(p=fabs(self.lati))

    @Property_RO
    def perimeterAGM(self):
        '''Compute the perimeter of this ellipse using the U{Arithmetic-Geometric Mean
           <https://PaulBourke.net/geometry/ellipsecirc>} formula (C{meter}, conventionally).
        '''
        _, P, a, b = self._Pab4
        if P is None:
            t  = _TOL53
            m  = -1
            c  = a + b
            ds = [c**2]
            _d = ds.append
            for _ in range(self._maxit):  # 4..5 trips
                b  = sqrt(a * b)
                a  = c * _0_5
                c  = a + b
                d  = a - b
                m *= 2
                _d(d**2 * m)
                if d <= (b * t):
                    break
            else:  # PYCHOK no cover
                raise _ConvergenceError(self._maxit, _over(d, b), t)
            P = _over(_fsum(ds) * PI, c)  # nonfinites=True
        return Meter(perimeterAGM=P)

    @Property_RO
    def perimeter4Arc3(self):
        '''Compute the perimeter (and arcs) of this ellipse using the U{4-Arc
           <https://PaulBourke.net/geometry/ellipsecirc>} (aka 4-Center)
           approximation as a 3-Tuple C{(P, Ra, Rb)} with perimeter C{P}, arc
           radii C{Ra} and C{Rb} at the respective semi-axes (all in C{meter},
           conventionally).
        '''
        r, P, a, b = self._Pab4
        if P is None:
            h = hypot(a, b)
            t = atan2(b, a)
            s, c = sincos2(t)
            L = (h - (a - b)) * _0_5
            a = _over(L, c)
            b = _over(h - L, s)
            P = (t * b + (PI_2 - t) * a) * _4_0
        elif a > b:  # flat
            a, b = _0_0, _1_over(b)  # INF
#       else:  # circular
#           pass
        if r:
            a, b = b, a
        return Meter(perimeter4Arc=P), Radius(Ra=a), Radius(Rb=b)

#   @Property_RO
#   def perimeterBPA(self):
#       '''Compute the perimeter of this ellipse using the U{Bronshtein Padé
#          Approximant <https://www.math.TTU.edu/~pearce/papers/schov.pdf>}
#          (C{meter}, conventionally).
#       '''
#       P, h = self._Ph2
#       if h:
#           h *=  h
#           P *= _over(h**2 * _3_0 - _64_0, h * _16_0 - _64_0) * PI
#       return Meter(perimeterBPA=P)

    @Property_RO
    def perimeterCR(self):
        '''Compute the perimeter of this ellipse using U{Rackauckas'
           <https://www.ChrisRackauckas.com/assets/Papers/ChrisRackauckas-The_Circumference_of_an_Ellipse.pdf>}
           approximation, also U{here<https://ExtremeLearning.com.AU/a-formula-for-the-perimeter-of-an-ellipse>} and
           U{here<http://www.EByte.IT/library/docs/math05a/EllipsePerimeterApprox05.html>} (C{meter}, conventionally).
        '''
        P, h = self._Ph2
        if h:
            h *=  h
            P *= _over(fhorner(h, 135168,  -85760, -5568, 3867),
                       fhorner(h, 135168, -119552, 22208, -345)) * PI
        return Meter(perimeterCR=P)

    @Property_RO
    def perimeterGK(self):
        '''Compute the perimeter of this ellipse using the U{Gauss-Kummer
           <https://www.JohnDCook.com/blog/2023/05/28/approximate-ellipse-perimeter>}
           series, C{B{b / a} > 0.75} (C{meter}, conventionally).
        '''
        P, h = self._Ph2
        if h:
            P *= fhorner(h**2, *self._GKs) * PI
        return Meter(perimeterGK=P)

    @Property_RO
    def perimeterHGK(self):
        '''Compute the perimeter of this ellipse using the U{Hypergeometric Gauss-Kummer
           <https://web.Tecnico.ULisboa.PT/~mcasquilho/compute/com/,ellips/PerimeterOfEllipse.pdf>}
           series (C{meter}, conventionally).
        '''
        P, h = self._Ph2
        if h:
            hs =  self._HGKs(h, self._maxit)
            P *= _fsum(hs) * PI  # nonfinites=True
        return Meter(perimeterHGK=P)

#   @Property_RO
#   def perimeterJWPA(self):
#       '''Compute the perimeter of this ellipse using the U{Jacobson-Waadeland
#          Padé Approximant <https://www.math.TTU.edu/~pearce/papers/schov.pdf>}
#          (C{meter}, conventionally).
#       '''
#       P, h = self._Ph2
#       if h:
#           h *=  h
#           P *= _over(fhorner(h, 256,  -48, -21),
#                      fhorner(h, 256, -112,   3)) * PI
#       return Meter(perimeterJWPA=P)

    @Property_RO
    def perimeter2k(self):
        '''Compute the perimeter of this ellipse using the complete integral
           of the 2nd kind, C{Elliptic.cE} (C{meter}, conventionally).
        '''
        return Meter(perimeter2k=self._perimeter2k(self._ellipE))

    @Property_RO
    def perimeter2k_(self):
        '''Compute the perimeter of this ellipse using U{SciPy's ellipe
           <https://www.JohnDCook.com/perimeter_ellipse.html>} function
           if available, otherwise use property C{perimeter2k} (C{meter},
           conventionally).
        '''
        return Meter(perimeter2k_=self._perimeter2k(self._ellipe or self._ellipE))

    def _perimeter2k(self, _ellip):
        '''(INTERNAL) Helper for methods C{.PE2k} and C{.Pe2k}.
        '''
        _, P, a, _ = self._Pab4
        if P is None:  # see .ellipsoids.Ellipsoid.L
            k =  self.e2
            P = _ellip(k) * a * _4_0
        return P

#   @Property_RO
#   def perimeterLS(self):
#       '''Compute the perimeter of this ellipse using the U{Linderholm-Segal
#          <https://www.JohnDCook.com/blog/2021/03/24/perimeter-of-an-ellipse>}
#          formula, aka C{3/2 norm} (C{meter}, conventionally).
#       '''
#       _, P, a, b = self._Pab4
#       if P is None:
#           n = pow(a, _1_5) + pow(b, _1_5)
#           P = pow(n * _0_5, _2_3rd) * PI2
#       return Meter(perimeterLS=P)

    @Property_RO
    def perimeter2R(self):
        '''Compute the perimeter of this ellipse using U{Ramanujan's 2nd
           <https://PaulBourke.net/geometry/ellipsecirc>} approximation,
           C{B{b / a} > 0.9} (C{meter}, conventionally).
        '''
        P, h = self._Ph2
        if h:
            P *= _2RC(h, _1_0)
        return Meter(perimeter2R=P)

    @Property_RO
    def perimeter2RC(self):
        '''Compute the perimeter of this ellipse using U{Cantrell Ramanujan's 2nd
           <http://www.EByte.IT/library/docs/math05a/EllipsePerimeterApprox05.html>}
           approximation, C{B{b / a} > 0.9} (C{meter}, conventionally).
        '''
        P, h = self._Ph2
        if h:
            P *= _2RC(h, pow(h, 24) * self._4_PI_ + _1_0)
        return Meter(perimeter2RC=P)

#   @Property_RO
#   def perimeterSPA(self):
#       '''Compute the perimeter of this ellipse using the U{Selmer Padé Approximant
#          <http://www.EByte.IT/library/docs/math05a/EllipsePerimeterApprox05.html>}
#          (C{meter}, conventionally).
#       '''
#       P, h = self._Ph2
#       if h:
#           h *=  h
#           P *= _over(h * _3_0 + _16_0, _16_0 - h) * PI
#       return Meter(perimeterSPA=P)

    @Property_RO
    def _Ph2(self):
        _, P, a, b = self._Pab4
        if P is None:
            if b:
                P =  a + b
                h = (a - b) / P
            else:
                P =  a
                h = _1_0
        else:
            h = None
        return P, h

    def point(self, deg_x, y=None):
        '''Return the point I{on} this ellipse at C{B{deg}} or C{atan2d(B{y}, B{x})
           degrees} along this ellipse.

           @return: A 2-tuple C{(x, y)}.
        '''
        s, c = sincos2d(deg_x) if y is None else self._sc2(deg_x, y, None)
        return (self.a * c), (self.b * s)

    def points(self, np, nq=4, ccw=False, ended=False, eps=EPS):  # MCCABE 13
        '''Yield up to C{np} points along this ellipse, each a 2-tuple C{(x, y)},
           starting at semi-axis C{+a}, in (counter-)clockwise order and distributed
           evenly along the minor semi-axis.

           @arg np: Number of points to generate (C{int}).
           @kwarg nq: Number of quarters to cover (C{int}, 1..4).
           @kwarg ccw: Use C{B{ccw}=True} for counter-clockwise order (C{bool}).
           @kwarg ended: If C{True}, include the last quadrant's end point (C{bool}).
           @kwarg eps: Tolerance for duplicate points (C{meter}, conventionally).

           @see: U{Directrix<https://MathWorld.Wolfram.com/ConicSectionDirectrix.html>}.
        '''
        a, b, _ = self._ab3
        if min(a, b) > eps and not self.isFlat:
            q  =  max(min(int(nq),  4), 1)
            n  =  max(    int(np) // q, 1)
            if not ccw:
                b = -b
            ps = list(_q1ps(a, b, n, eps))
            for p in ps:  # 1st quadrant
                yield p
            p0 =  ps.pop(0)  # E
            pq = _0_0, b  # N/S
            if q > 1:
                yield pq
                for x, y in reversed(ps):  # 2nd
                    yield (-x), y
                pq = (-a), _0_0  # W
                if q > 2:
                    yield pq
                    for x, y in ps:  # 3rd
                        yield (-x), (-y)
                    pq = _0_0, (-b)  # S/N
                    if q > 3:
                        yield pq
                        for x, y in reversed(ps):  # 4th
                            yield x, (-y)
                        pq = p0
            if ended:
                yield pq
        else:  # "flat"
            p0 = a, b
            yield p0
            if max(a, b) > eps:
                yield (-a), (-b)
                if ended:
                    yield p0

    def polar2d(self, deg_x, y=None):
        '''For a point at C{B{deg}} or C{atan2d(B{y}, B{x}) degrees} along this
           ellipse, return 2-tuple C{(radius, angle)} with the polar U{radius
           <https://WikiPedia.org/wiki/Ellipse#Polar_form_relative_to_center>}
           from the center (C{meter}, conventionally) and C{angle} in C{degrees}.
        '''
        return polar2d(*self.point(deg_x, y))

    @Property_RO
    def R1(self):
        '''Get this ellipse' I{arithmetic mean} radius, C{(2 * a + b) / 3} (C{meter}, conventionally).
        '''
        _, _, a, r = self._Pab4
        if r:
            r = fmean_(a, a, r) if a > r else a
        return Radius(R1=r or _0_0)

    @Property_RO
    def R2(self):
        '''Get this ellipse' I{authalic} radius, C{sqrt(B{a} * B{b})} (C{meter}, conventionally).
        '''
        a, b, ab = self._ab3
        return Radius(R2=(sqrt(ab) if a != b else float(a)) if ab else _0_0)

    Rauthalic = Rgeometric = R2

    def Roc(self, deg_x, y=None, eps=None):
        '''Compute the U{radius of curvature<https://WikiPedia.org/wiki/Radius_of_curvature>}
           at a point C{B{deg}} or C{atan2d(B{y}, B{x}) degrees} along this ellipse.

           @see: Method L{Roc_<Ellipse.Roc_>} for ruther details.
        '''
        x = radians(deg_x) if y is None else deg_x
        return self.Roc_(x, y, eps=eps)

    def Roc_(self, rad_x, y=None, eps=None):
        '''Compute the U{radius of curvature<https://WikiPedia.org/wiki/Radius_of_curvature>}
           at a point C{B{rad}} or C{atan2(B{y}, B{x}) radians} along this ellipse.

           @kwarg eps: See method C{sideOf}, use C{B{eps}=0} to permit any points.

           @return: Curvature (C{meter}, conventionally).

           @raise ValueError: Point C{(B{x}, B{y})} off this ellipse, unless C{B{eps}=0}.
        '''
        try:
            a, b, ab = self._ab3
            if b != a:
                s, c = sincos2(rad_x) if y is None else self._sc2(rad_x, y, eps)
                r = _over(hypot(a * s, b * c)**3, ab)
            else:  # circular
                r = float(a)
        except Exception as X:
            raise self._Error(self.Roc_, cause=X)
        return Radius(Roc=r)

    @Property_RO
    def Rrectifying(self):
        '''Get this ellipse' I{rectifying} radius, C{perimeter2k_ / PI2} (C{meter}, conventionally).
        '''
        return Radius(Rrectifying=self.perimeter2k_ / PI2)

    def _sc2(self, x, y, eps):
        '''(INTERNAL) Helper for methods C{.point}, C{.polar}, C{.Roc_} and C{.slope_}.
        '''
        if eps and eps > 0:
            s = self._sideOf(x, y, eps)
            if s:
                raise _ValueError(x=x, y=y, eps=eps, sideOf=s)
        h = hypot(x, y)
        s = _over(y, h)
        c = _over(x, h)
        return s, c

    def _sideOf(self, x, y, eps):
        '''(INTERNAL) Helper for methods C{._sc2} and C{.sideOf}.
        '''
        a, b, ab = self._ab3
        s = ab or max(a, b)
        if s:
            s = (hypot(x * b, y * a) - s) / s
#           s = max(_N_1_0, min(_1_0, s))
        else:  # dot
            s = _1_0 if x or y else _0_0
        return INT0 if fabs(s) < eps else s

    def sideOf(self, x, y, eps=EPS):
        '''Return a C{positive}, C{negative} or C{0} fraction if point C{(B{x}, B{y})}
           is C{outside}, C{inside} respectively C{on} this ellipse.
        '''
        try:
            return Scalar(sideOf=self._sideOf(x, y, eps))
        except Exception as X:
            raise self._Error(self.sideOf, x=x, y=y, cause=X)

    def slope(self, deg_x, y=None, eps=None):
        '''Compute the U{tangent slope<https://WikiPedia.org/wiki/Ellipse#Tangent_slope_as_parameter>}
           at a point C{B{deg}} or C{atan2d(B{y}, B{x}) degrees} along this ellipse.

           @return: Slope (C{degrees}), negative for C{0 <= B{deg} < 90}.

           @see: Method L{slope_<Ellipse.slope_>} for further details.
        '''
        x = radians(deg_x) if y is None else deg_x
        return Degrees(slope=degrees(self.slope_(x, y, eps=eps)))

    def slope_(self, rad_x, y=None, eps=None):
        '''Compute the U{tangent slope<https://WikiPedia.org/wiki/Ellipse#Tangent_slope_as_parameter>}
           at a point C{B{rad}} or C{atan2(B{y}, B{x}) radians} along this ellipse.

           @kwarg eps: See method C{sideOf}, use C{B{eps}=0} to permit any points.

           @return: Slope (C{radians}), negative for C{0 <= B{rad} < PI/2}.

           @raise ValueError: C{(B{x}, B{y})} off this ellipse, unless C{B{eps}=0}.
        '''
        # <https://UNacademy.com/content/jee/study-material/mathematics/equation-of-a-tangent-to-the-ellipse/>
        s, c = sincos2(rad_x) if y is None else self._sc2(rad_x, y, eps)
        r = atan2p(-self.b * c, self.a * s)
        if r >= PI3_2:
            r -= PI2
        return Radians(slope=r or _0_0)  # no -0.0

    def toEllipsoid(self, **Ellipsoid_and_kwds):
        '''Return an L{Ellipsoid<pygeodesy.Ellipsoid>} from this ellipse'
           C{a} and C{b} semi-axes.

           @kwarg Ellipsoid_and_kwds: Optional C{B{Ellipsoid}=Ellipsoid} class
                                and additional C{Ellipsoid} keyword arguments.
        '''
        E, kwds = _xkwds_pop2(Ellipsoid_and_kwds, Ellipsoid=
                                 _MODS.ellipsoids.Ellipsoid)
        return E(self.a, b=self.b, **_xkwds(kwds, name=self.name))

    def toStr(self, prec=8, terse=2, **sep_name):  # PYCHOK signature
        '''Return this ellipse as a text string.

           @kwarg prec: Number of decimal digits, unstripped (C{int}).
           @kwarg terse: Limit the number of items (C{int}, 0...9),
                         use C{B{terse}=0} or C{=None} for all.
           @kwarg sep_name: Optional C{B{name}=NN} (C{str}) or C{None}
                      to exclude this ellipse' name and separator
                      C{B{sep}=", "} to join the items (C{str}).

           @return: This C{Ellipse}' attributes (C{str}).
        '''
        E = Ellipse
        t = E.a, E.b
        if (terse or 0) != 2:
            t += E.c, E.e, E.e2, E.p, E.area, E.perimeter2k, E.R2
            if terse:
                t = t[:terse]
        return self._instr(prec=prec, props=t, **sep_name)

    def toTriaxial_(self, c=EPS, **Triaxial_and_kwds):  # like .Ellipse5Tuple.toTriaxial_
        '''Return a L{Triaxial_<pygeodesy.Triaxial_>} from this ellipse' semi-axes.

           @kwarg c: Near-zero, minor semi-axis (C{meter}, conventionally).
           @kwarg Triaxial_and_kwds: Optional C{B{Triaxial}=Triaxial_} class and
                               additional C{Triaxial} keyword arguments.
        '''
        T, kwds = _xkwds_pop2(Triaxial_and_kwds, Triaxial=_MODS.triaxials.Triaxial_)
        return T(self.a, b=self.b, c=c, **_xkwds(kwds, name=self.name or _UNDER_))  # 'NN'

    def _triaxialX(self, method, *args, **kwds):
        '''(INTERNAL) Invoke a triaxial method and map exceptions to L{EllipseError}s.
        '''
        try:
            t = self._tEPS
            if t is None:
                self._tEPS = t = self.toTriaxial_(EPS)
            _m = getattr(t, method.__name__)
            return _m(*args, **kwds)
        except Exception as x:
            raise self._Error(method, Triaxial_=t, cause=x)

    def _xy03(self, deg_x, y):
        if y is None:
            y, x = sincos2d(deg_x)
            y *= self.b
            x *= self.a
        else:
            x  = float(deg_x)
            y  = float(y)
        return x, y, 0


class EllipseError(_ValueError):
    '''Raised for any L{Ellipse} or C{ellipses} issue.
    '''
    pass  # ...


def _arc(_e, k, r):
    # in C{Ellipse.arc_}
    t, r = divmod(r, PI2)
    R = _e(k, r)  # phi=r
    if t:  # + t * perimeter
        t *= _e(k) * _4_0
        R +=  t
    return R


def _isFlat(a, b):  # in .triaxials.bases
    # is C{b <<< a}?
    return b < (a * _TOL53_53)


def _q1ps(a, b, n, eps):
    # yield the 1st quadrant C{Ellipse.points}
    if a > b:  # oblate
        def _yx2(i):
            y = i / n
            return y, sqrt(_1_0 - y**2)

    elif a < b:  # prolate
        def _yx2(i):  # PYCHOK redef
            x = (n - i) / n
            return sqrt(_1_0 - x**2), x

    else:  # circular
        r = PI_2 / n
        def _yx2(i):  # PYCHOK redef
            return sincos2(r * i)

    p = a, _0_0  # == p0
    yield p
    for i in range(1, n):
        y, x = _yx2(i)
        y *= b
        x *= a
        if euclid(x, y, *p) > eps:
            p = x, y
            yield p


def _2RC(h, r):  # in Ellipse.perimeter2R and .perimeter2RC
    h *= _3_0 * h
    r +=  h / (sqrt(_4_0 - h) + _10_0)
    return r * PI

# **) MIT License
#
# Copyright (C) 2026-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
