
# -*- coding: utf-8 -*-

u'''(INTERNAL) Private, 3-D vector base class C{Vector3dBase}.

A pure Python implementation of vector-based functions by I{(C) Chris Veness
2011-2024} published under the same MIT Licence**, see U{Vector-based geodesy
<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.
'''

from pygeodesy.basics import _copysign, islistuple, isscalar, map1, map2, \
                             _signOf, _zip
from pygeodesy.constants import EPS, EPS0, INT0, PI, PI2, _copysignINF, \
                               _float0, isnear0, isnear1, isneg0, \
                               _pos_self, _1_0
from pygeodesy.errors import CrossError, VectorError, _xcallable, _xError
from pygeodesy.fmath import euclid_, fdot, fdot_, hypot_, hypot2_  # _MODS.fmath.fma
from pygeodesy.interns import _coincident_, _colinear_, _COMMASPACE_, _xyz_
from pygeodesy.lazily import _ALL_LAZY, _ALL_DOCS, _ALL_MODS as _MODS
from pygeodesy.named import _NamedBase, _NotImplemented, _xother3
# from pygeodesy.namedTuples import Vector3Tuple  # _MODS
from pygeodesy.props import deprecated_method, Property, Property_RO, \
                            property_doc_, property_RO, _update_all
from pygeodesy.streprs import Fmt, strs, unstr
from pygeodesy.units import Float, Scalar
from pygeodesy.utily import atan2, sincos2,  fabs

from math import ceil as _ceil, floor as _floor, trunc as _trunc

__all__ = _ALL_LAZY.vector3dBase
__version__ = '24.11.24'


class Vector3dBase(_NamedBase):  # sync __methods__ with .fsums.Fsum
    '''(INTERNAL) Generic 3-D vector base class.
    '''
    _crosserrors = True  # un/set by .errors.crosserrors

    _ll = None  # original latlon, '_fromll'
#   _x  = INT0  # X component
#   _y  = INT0  # Y component
#   _z  = INT0  # Z component

    def __init__(self, x_xyz, y=INT0, z=INT0, ll=None, **name):
        '''New L{Vector3d} or C{Vector3dBase} instance.

           The vector may be normalised or use x, y, z for position and
           distance from earth centre or height relative to the surface
           of the earth' sphere or ellipsoid.

           @arg x_xyz: X component of vector (C{scalar}) or a (3-D) vector
                       (C{Cartesian}, L{Ecef9Tuple}, C{Nvector}, L{Vector3d},
                       L{Vector3Tuple}, L{Vector4Tuple} or a C{tuple} or
                       C{list} of 3+ C{scalar} items).
           @kwarg y: Y component of vector (C{scalar}), ignored if B{C{x_xyz}}
                     is not C{scalar}, otherwise same units as B{C{x_xyz}}.
           @kwarg z: Z component of vector (C{scalar}), ignored if B{C{x_xyz}}
                     is not C{scalar}, otherwise same units as B{C{x_xyz}}.
           @kwarg ll: Optional latlon reference (C{LatLon}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise VectorError: Invalid B{C{x_xyz}}, B{C{y}} or B{C{z}}.
        '''
        self._x, \
        self._y, \
        self._z = _xyz3(type(self), x_xyz, y, z) if isscalar(x_xyz) else \
                  _xyz3(type(self), x_xyz)
        if ll:
            self._ll = ll
        if name:
            self.name = name

    def __abs__(self):
        '''Return the norm of this vector.

           @return: Norm, unit length (C{float});
        '''
        return self.length

    def __add__(self, other):
        '''Add this to an other vector (L{Vector3d}).

           @return: Vectorial sum (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.plus(other)

    def __bool__(self):  # PYCHOK PyChecker
        '''Is this vector non-zero?

           @see: Method C{bools}.
        '''
        return bool(self.x or self.y or self.z)

    def __ceil__(self):  # PYCHOK no cover
        '''Return a vector with the C{ceil} of these components.

           @return: Ceil-ed (L{Vector3d}).
        '''
        return self._mapped(_ceil)

    def __cmp__(self, other):  # Python 2-
        '''Compare this and an other vector (L{Vector3d}).

           @return: -1, 0 or +1 (C{int}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return _signOf(self.length, self._other_cmp(other))

    cmp = __cmp__

    def __divmod__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __eq__(self, other):
        '''Is this vector equal to an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if equal, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.isequalTo(other, eps=EPS0)

    def __float__(self):  # PYCHOK no cover
        '''Not implemented, see method C{float}.'''
        return _NotImplemented(self)  # must return C{float}

    def __floor__(self):  # PYCHOK no cover
        '''Return a vector with the C{floor} of these components.

           @return: Floor-ed (L{Vector3d}).
        '''
        return self._mapped(_floor)

    def __floordiv__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __ge__(self, other):
        '''Is this vector longer than or equal to an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if so, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.length >= self._other_cmp(other)

#   def __getitem__(self, key):
#       '''Return C{item} at index or slice C{[B{key}]}.
#       '''
#       return self.xyz[key]

    def __gt__(self, other):
        '''Is this vector longer than an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if so, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.length > self._other_cmp(other)

    def __hash__(self):  # PYCHOK no cover
        '''Return this instance' C{hash}.
        '''
        return hash(self.xyz)  # XXX id(self)?

    def __iadd__(self, other):
        '''Add this and an other vector I{in-place}, C{this += B{other}}.

           @arg other: The other vector (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self._xyz(self.plus(other))

    def __ifloordiv__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __imatmul__(self, other):  # PYCHOK Python 3.5+
        '''Cross multiply this and an other vector I{in-place}, C{this @= B{other}}.

           @arg other: The other vector (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.

           @see: Luciano Ramalho, "Fluent Python", O'Reilly, 2016 p. 397+, 2022 p. 578+.
        '''
        return self._xyz(self.cross(other))

    def __imod__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __imul__(self, scalar):
        '''Multiply this vector by a scalar I{in-place}, C{this *= B{scalar}}.

           @arg scalar: Factor (C{scalar}).

           @raise TypeError: Non-scalar B{C{scalar}}.
        '''
        return self._xyz(self.times(scalar))

    def __int__(self):  # PYCHOK no cover
        '''Not implemented, see method C{ints}.'''
        return _NotImplemented(self)  # must return C{int}

    def __ipow__(self, other, *mod):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other, *mod)

    def __isub__(self, other):
        '''Subtract an other vector from this one I{in-place}, C{this -= B{other}}.

           @arg other: The other vector (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self._xyz(self.minus(other))

#   def __iter__(self):
#       '''Return an C{iter}ator over this vector's components.
#       '''
#       return iter(self.xyz3)

    def __itruediv__(self, scalar):
        '''Divide this vector by a scalar I{in-place}, C{this /= B{scalar}}.

           @arg scalar: The divisor (C{scalar}).

           @raise TypeError: Non-scalar B{C{scalar}}.
        '''
        return self._xyz(self.dividedBy(scalar))

    def __le__(self, other):  # Python 3+
        '''Is this vector shorter than or equal to an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if so, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.length <= self._other_cmp(other)

#   def __len__(self):
#       '''Return C{3}, always.
#       '''
#       return len(self.xyz)

    def __lt__(self, other):  # Python 3+
        '''Is this vector shorter than an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if so, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.length < self._other_cmp(other)

    def __matmul__(self, other):  # PYCHOK Python 3.5+
        '''Compute the cross product of this and an other vector, C{this @ B{other}}.

           @arg other: The other vector (L{Vector3d}).

           @return: Cross product (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.cross(other)

    def __mod__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __mul__(self, scalar):
        '''Multiply this vector by a scalar, C{this * B{scalar}}.

           @arg scalar: Factor (C{scalar}).

           @return: Product (L{Vector3d}).
        '''
        return self.times(scalar)

    def __ne__(self, other):
        '''Is this vector not equal to an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if so, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return not self.isequalTo(other, eps=EPS0)

    def __neg__(self):
        '''Return the opposite of this vector.

           @return: A negated copy (L{Vector3d})
        '''
        return self.negate()

    def __pos__(self):  # PYCHOK no cover
        '''Return this vector I{as-is} or a copy.

           @return: This instance (L{Vector3d})
        '''
        return self if _pos_self else self.copy()

    def __pow__(self, other, *mod):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other, *mod)

    __radd__ = __add__  # PYCHOK no cover

    def __rdivmod__ (self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

#   def __repr__(self):
#       '''Return the default C{repr(this)}.
#       '''
#       return self.toRepr()

    def __rfloordiv__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rmatmul__(self, other):  # PYCHOK Python 3.5+
        '''Compute the cross product of an other and this vector, C{B{other} @ this}.

           @arg other: The other vector (L{Vector3d}).

           @return: Cross product (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.others(other).cross(self)

    def __rmod__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    __rmul__ = __mul__

    def __round__(self, *ndigits):  # PYCHOK no cover
        '''Return a vector with these components C{rounded}.

           @arg ndigits: Optional number of digits (C{int}).

           @return: Rounded (L{Vector3d}).
        '''
        # <https://docs.Python.org/3.12/reference/datamodel.html?#object.__round__>
        return self.classof(*(round(_, *ndigits) for _ in self.xyz3))

    def __rpow__(self, other, *mod):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other, *mod)

    def __rsub__(self, other):  # PYCHOK no cover
        '''Subtract this vector from an other vector, C{B{other} - this}.

           @arg other: The other vector (L{Vector3d}).

           @return: Difference (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.others(other).minus(self)

    def __rtruediv__(self, scalar):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, scalar)

#   def __str__(self):
#       '''Return the default C{str(self)}.
#       '''
#       return self.toStr()

    def __sub__(self, other):
        '''Subtract an other vector from this vector, C{this - B{other}}.

           @arg other: The other vector (L{Vector3d}).

           @return: Difference (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.minus(other)

    def __truediv__(self, scalar):
        '''Divide this vector by a scalar, C{this / B{scalar}}.

           @arg scalar: The divisor (C{scalar}).

           @return: Quotient (L{Vector3d}).

           @raise TypeError: Non-scalar B{C{scalar}}.
        '''
        return self.dividedBy(scalar)

    def __trunc__(self):  # PYCHOK no cover
        '''Return a vector with the C{trunc} of these components.

           @return: Trunc-ed (L{Vector3d}).
        '''
        return self._mapped(_trunc)

    if _MODS.sys_version_info2 < (3, 0):  # PYCHOK no cover
        # <https://docs.Python.org/2/library/operator.html#mapping-operators-to-functions>
        __div__     = __truediv__
        __idiv__    = __itruediv__
        __long__    = __int__
        __nonzero__ = __bool__
        __rdiv__    = __rtruediv__

    def angleTo(self, other, vSign=None, wrap=False):
        '''Compute the angle between this and an other vector.

           @arg other: The other vector (L{Vector3d}).
           @kwarg vSign: Optional vector, if supplied (and out of the
                         plane of this and the other), angle is signed
                         positive if this->other is clockwise looking
                         along vSign or negative in opposite direction,
                         otherwise angle is unsigned.
           @kwarg wrap: If C{True}, wrap/unroll the angle to +/-PI (C{bool}).

           @return: Angle (C{radians}).

           @raise TypeError: If B{C{other}} or B{C{vSign}} not a L{Vector3d}.
        '''
        x = self.cross(other)
        s = x.length
        # use vSign as reference to set sign of s
        if s and vSign and x.dot(vSign) < 0:
            s = -s

        a = atan2(s, self.dot(other))
        if wrap and fabs(a) > PI:
            a -= _copysign(PI2, a)
        return a

    def apply(self, fun2, other_x, *y_z, **fun2_kwds):
        '''Apply a 2-argument function pairwise to the components
           of this and an other vector.

           @arg fun2: 2-Argument callable (C{any(scalar, scalar}),
                      return a C{scalar} or L{INT0} result.
           @arg other_x: Other X component (C{scalar}) or a vector
                         with X, Y and Z components (C{Cartesian},
                         L{Ecef9Tuple}, C{Nvector}, L{Vector3d},
                         L{Vector3Tuple} or L{Vector4Tuple}).
           @arg y_z: Other Y and Z components, positional (C{scalar}, C{scalar}).
           @kwarg fun2_kwds: Optional keyword arguments for B{C{fun2}}.

           @return: New, applied vector (L{Vector3d}).

           @raise ValueError: Invalid B{C{other_x}} or B{C{y_z}}.
        '''
        _xcallable(fun2=fun2)
        if fun2_kwds:
            def _f2(a, b):
                return fun2(a, b, **fun2_kwds)
        else:
            _f2 = fun2

        xyz = _xyz3(self.apply, other_x, *y_z)
        xyz = (_f2(a, b) for a, b in _zip(self.xyz3, xyz))  # strict=True
        return self.classof(*xyz)

    def bools(self):
        '''Return the vector with C{bool} components.
        '''
        return self._mapped(bool)

    def cross(self, other, raiser=None, eps0=EPS):  # raiser=NN
        '''Compute the cross product of this and an other vector.

           @arg other: The other vector (L{Vector3d}).
           @kwarg raiser: Optional, L{CrossError} label if raised (C{str}, non-L{NN}).
           @kwarg eps0: Near-zero tolerance (C{scalar}), same units as C{x}, C{y} and
                        C{z}.

           @return: Cross product (L{Vector3d}).

           @raise CrossError: Zero or near-zero cross product and if B{C{raiser}} and
                              L{crosserrors<pygeodesy.crosserrors>} are both C{True}.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        X, Y, Z = self.others(other).xyz3
        x, y, z = self.xyz3
        xyz = (fdot_(y, Z, -z, Y),
               fdot_(z, X, -x, Z),
               fdot_(x, Y, -y, X))

        if raiser and self.crosserrors and eps0 > 0 \
                  and max(map(fabs, xyz)) < eps0:
            r =  other._fromll or other
            s =   self._fromll or self
            t =  self.isequalTo(other, eps=eps0)
            t = _coincident_ if t else _colinear_
            raise CrossError(raiser, s, other=r, txt=t)

        return self.classof(*xyz)  # name__=self.cross

    @property_doc_('''raise or ignore L{CrossError} exceptions (C{bool}).''')
    def crosserrors(self):
        '''Get L{CrossError} exceptions (C{bool}).
        '''
        return self._crosserrors

    @crosserrors.setter  # PYCHOK setter!
    def crosserrors(self, raiser):
        '''Raise or ignore L{CrossError} exceptions (C{bool}).
        '''
        self._crosserrors = bool(raiser)

    def dividedBy(self, divisor):
        '''Divide this vector by a scalar.

           @arg divisor: The divisor (C{scalar}).

           @return: New, scaled vector (L{Vector3d}).

           @raise TypeError: Non-scalar B{C{divisor}}.

           @raise VectorError: Invalid or zero B{C{divisor}}.
        '''
        d = Scalar(divisor=divisor)
        try:
            return self._times(_1_0 / d)
        except (ValueError, ZeroDivisionError) as x:
            raise VectorError(divisor=divisor, cause=x)

    def dot(self, other):
        '''Compute the dot (scalar) product of this and an other vector.

           @arg other: The other vector (L{Vector3d}).

           @return: Dot product (C{float}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.length2 if other is self else fdot(
               self.xyz3, *self.others(other).xyz3)

    @deprecated_method
    def equals(self, other, units=False):  # PYCHOK no cover
        '''DEPRECATED, use method C{isequalTo}.
        '''
        return self.isequalTo(other, units=units)

    @Property_RO
    def euclid(self):
        '''I{Approximate} the length (norm, magnitude) of this vector (C{Float}).

           @see: Properties C{length} and C{length2} and function
                 L{pygeodesy.euclid_}.
        '''
        return Float(euclid=euclid_(*self.xyz3))

    def equirectangular(self, other):
        '''I{Approximate} the difference between this and an other vector.

           @arg other: Vector to subtract (C{Vector3dBase}).

           @return: The length I{squared} of the difference (C{Float}).

           @raise TypeError: Incompatible B{C{other}} C{type}.

           @see: Property C{length2}.
        '''
        d = self.minus(other)
        return Float(equirectangular=hypot2_(*d.xyz3))

    def fabs(self):
        '''Return the vector with C{fabs} components.
        '''
        return self._mapped(fabs)

    def floats(self):
        '''Return the vector with C{float} components.
        '''
        return self._mapped(_float0)

    @Property
    def _fromll(self):
        '''(INTERNAL) Get the latlon reference (C{LatLon}) or C{None}.
        '''
        return self._ll

    @_fromll.setter  # PYCHOK setter!
    def _fromll(self, ll):
        '''(INTERNAL) Set the latlon reference (C{LatLon}) or C{None}.
        '''
        self._ll = ll or None

    @property_RO
    def homogeneous(self):
        '''Get this vector's homogeneous representation (L{Vector3d}).
        '''
        x, y, z = self.xyz3
        if z:
            x =  x / z  # /= chokes PyChecker
            y =  y / z
#           z = _1_0
        else:
            if isneg0(z):
                x = -x
                y = -y
            x = _copysignINF(x)
            y = _copysignINF(y)
#           z =  NAN
        return self.classof(x, y, _1_0)

    def intermediateTo(self, other, fraction, **unused):  # height=None, wrap=False
        '''Locate the vector at a given fraction between (or along) this
           and an other vector.

           @arg other: The other vector (L{Vector3d}).
           @arg fraction: Fraction between both vectors (C{scalar},
                          0.0 for this and 1.0 for the other vector).

           @return: Intermediate vector (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        f = Scalar(fraction=fraction)
        if isnear0(f):  # PYCHOK no cover
            r = self
        else:
            r = self.others(other)
            if not isnear1(f):  # self * (1 - f) + r * f
                r = self.plus(r.minus(self)._times(f))
        return r

    def ints(self):
        '''Return the vector with C{int} components.
        '''
        return self._mapped(int)

    def isconjugateTo(self, other, minum=1, eps=EPS):
        '''Determine whether this and an other vector are conjugates.

           @arg other: The other vector (C{Cartesian}, L{Ecef9Tuple},
                       L{Vector3d}, C{Vector3Tuple} or C{Vector4Tuple}).
           @kwarg minum: Minimal number of conjugates required (C{int}, 0..3).
           @kwarg eps: Tolerance for equality and conjugation (C{scalar}),
                       same units as C{x}, C{y}, and C{z}.

           @return: C{True} if both vector's components either match
                    or at least C{B{minum}} have opposite signs.

           @raise TypeError: Incompatible B{C{other}} C{type}.

           @see: Method C{isequalTo}.
        '''
        self.others(other)
        n = 0
        for a, b in zip(self.xyz3, other.xyz3):
            if fabs(a + b) < eps and ((a < 0 and b > 0) or
                                      (a > 0 and b < 0)):
                n += 1  # conjugate
            elif fabs(a - b) > eps:
                return False  # unequal
        return bool(n >= minum)

    def isequalTo(self, other, units=False, eps=EPS):
        '''Check if this and an other vector are equal or equivalent.

           @arg other: The other vector (L{Vector3d}).
           @kwarg units: Optionally, compare the normalized, unit
                         version of both vectors.
           @kwarg eps: Tolerance for equality (C{scalar}), same units as
                       C{x}, C{y}, and C{z}.

           @return: C{True} if vectors are identical, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.

           @see: Method C{isconjugateTo}.
        '''
        if units:
            self.others(other)
            d = self.unit().minus(other.unit())
        else:
            d = self.minus(other)
        return max(map(fabs, d.xyz3)) < eps

    @Property_RO
    def length(self):  # __dict__ value overwritten by Property_RO C{_united}
        '''Get the length (norm, magnitude) of this vector (C{Float}).

           @see: Properties L{length2} and L{euclid}.
        '''
        return Float(length=hypot_(self.x, self.y, self.z))

    @Property_RO
    def length2(self):  # __dict__ value overwritten by Property_RO C{_united}
        '''Get the length I{squared} of this vector (C{Float}).

           @see: Property L{length} and method C{equirectangular}.
        '''
        return Float(length2=hypot2_(self.x, self.y, self.z))

    def _mapped(self, func):
        '''(INTERNAL) Apply C{func} to all components.
        '''
        return self.classof(*map2(func, self.xyz3))

    def minus(self, other):
        '''Subtract an other vector from this vector.

           @arg other: The other vector (L{Vector3d}).

           @return: New vector difference (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self._minus(*self.others(other).xyz3)

    def _minus(self, x, y, z):
        '''(INTERNAL) Helper for methods C{.minus} and C{.minus_}.
        '''
        return self.classof(self.x - x, self.y - y, self.z - z)

    def minus_(self, other_x, *y_z):
        '''Subtract separate X, Y and Z components from this vector.

           @arg other_x: X component (C{scalar}) or a vector's
                         X, Y, and Z components (C{Cartesian},
                         L{Ecef9Tuple}, C{Nvector}, L{Vector3d},
                         L{Vector3Tuple}, L{Vector4Tuple}).
           @arg y_z: Y and Z components (C{scalar}, C{scalar}),
                     ignored if B{C{other_x}} is not C{scalar}.

           @return: New, vectiorial vector (L{Vector3d}).

           @raise ValueError: Invalid B{C{other_x}} or B{C{y_z}}.
        '''
        return self._minus(*_xyz3(self.minus_, other_x, *y_z))

    def negate(self):
        '''Return the opposite of this vector.

           @return: A negated copy (L{Vector3d})
        '''
        return self.classof(-self.x, -self.y, -self.z)

    @Property_RO
    def _N_vector(self):
        '''(INTERNAL) Get the (C{nvectorBase._N_vector_})
        '''
        return _MODS.nvectorBase._N_vector_(*self.xyz3, name=self.name)

    def _other_cmp(self, other):
        '''(INTERNAL) Return the value for comparison.
        '''
        return other if isscalar(other) else self.others(other).length

    def others(self, *other, **name_other_up):
        '''Refined class comparison.

           @arg other: The other vector (L{Vector3d}).
           @kwarg name_other_up: Overriding C{name=other} and C{up=1}
                                 keyword arguments.

           @return: The B{C{other}} if compatible.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        other, name, up = _xother3(self, other, **name_other_up)
        if not isinstance(other, Vector3dBase):
            _NamedBase.others(self, other, name=name, up=up + 1)
        return other

    def plus(self, other):
        '''Add this vector and an other vector.

           @arg other: The other vector (L{Vector3d}).

           @return: Vectorial sum (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self._plus(*self.others(other).xyz3)

    sum = plus  # alternate name

    def _plus(self, x, y, z):
        '''(INTERNAL) Helper for methods C{.plus} and C{.plus_}.
        '''
        return self.classof(self.x + x, self.y + y, self.z + z)

    def plus_(self, other_x, *y_z):
        '''Sum of this vector and separate X, Y and Z components.

           @arg other_x: X component (C{scalar}) or a vector's
                         X, Y, and Z components (C{Cartesian},
                         L{Ecef9Tuple}, C{Nvector}, L{Vector3d},
                         L{Vector3Tuple}, L{Vector4Tuple}).
           @arg y_z: Y and Z components (C{scalar}, C{scalar}),
                     ignored if B{C{other_x}} is not C{scalar}.

           @return: New, vectiorial vector (L{Vector3d}).

           @raise ValueError: Invalid B{C{other_x}} or B{C{y_z}}.
        '''
        return self._plus(*_xyz3(self.plus_, other_x, *y_z))

    def rotate(self, axis, theta, fma=False):
        '''Rotate this vector around an axis by a specified angle.

           @arg axis: The axis being rotated around (L{Vector3d}).
           @arg theta: The angle of rotation (C{radians}).
           @kwarg fma: If C{True}, use fused-multiply-add (C{bool}).

           @return: New, rotated vector (L{Vector3d}).

           @see: U{Rotation matrix from axis and angle<https://WikiPedia.org/wiki/
                 Rotation_matrix#Rotation_matrix_from_axis_and_angle>} and
                 U{Quaternion-derived rotation matrix<https://WikiPedia.org/wiki/
                 Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix>}.
        '''
        s, c = sincos2(theta)  # rotation angle
        d = _1_0 - c
        if d or s:
            p = self.unit().xyz  # point being rotated
            r = self.others(axis=axis).unit()  # axis being rotated around

            ax, ay, az = r.xyz3  # quaternion-derived rotation matrix
            bx, by, bz = r.times(d).xyz3
            sx, sy, sz = r.times(s).xyz3

            if fma:
                _fma = _MODS.fmath.fma
            else:
                def _fma(a, b, c):
                    return a * b + c

            x = fdot(p, _fma(ax, bx,   c), _fma(ax, by, -sz), _fma(ax, bz,  sy))
            y = fdot(p, _fma(ay, bx,  sz), _fma(ay, by,   c), _fma(ay, bz, -sx))
            z = fdot(p, _fma(az, bx, -sy), _fma(az, by,  sx), _fma(az, bz,   c))
        else:  # unrotated
            x, y, z = self.xyz3
        return self.classof(x, y, z)

    @deprecated_method
    def rotateAround(self, axis, theta):
        '''DEPRECATED, use method C{rotate}.'''
        return self.rotate(axis, theta)  # PYCHOK no cover

    def times(self, factor):
        '''Multiply this vector by a scalar.

           @arg factor: Scale factor (C{scalar}).

           @return: New, scaled vector (L{Vector3d}).

           @raise TypeError: Non-scalar B{C{factor}}.
        '''
        return self._times(Scalar(factor=factor))

    def _times(self, s):
        '''(INTERNAL) Helper for C{.dividedBy} and C{.times}.
        '''
        return self.classof(self.x * s, self.y * s, self.z * s)

    def times_(self, other_x, *y_z):
        '''Multiply this vector's components by separate X, Y and Z factors.

           @arg other_x: X scale factor (C{scalar}) or a vector's
                         X, Y, and Z components as scale factors
                         (C{Cartesian}, L{Ecef9Tuple}, C{Nvector},
                         L{Vector3d}, L{Vector3Tuple}, L{Vector4Tuple}).
           @arg y_z: Y and Z scale factors (C{scalar}, C{scalar}),
                     ignored if B{C{other_x}} is not C{scalar}.

           @return: New, scaled vector (L{Vector3d}).

           @raise ValueError: Invalid B{C{other_x}} or B{C{y_z}}.
        '''
        x, y, z = _xyz3(self.times_, other_x, *y_z)
        return self.classof(self.x * x, self.y * y, self.z * z)

#   @deprecated_method
#   def to2ab(self):  # PYCHOK no cover
#       '''DEPRECATED, use property C{Nvector.philam}.
#
#          @return: A L{PhiLam2Tuple}C{(phi, lam)}.
#       '''
#       return _MODS.nvectorBase.n_xyz2philam(self.x, self.y, self.z)

#   @deprecated_method
#   def to2ll(self):  # PYCHOK no cover
#       '''DEPRECATED, use property C{Nvector.latlon}.
#
#          @return: A L{LatLon2Tuple}C{(lat, lon)}.
#       '''
#       return _MODS.nvectorBase.n_xyz2latlon(self.x, self.y, self.z)

    @deprecated_method
    def to3xyz(self):  # PYCHOK no cover
        '''DEPRECATED, use property L{xyz}.
        '''
        return self.xyz

    def toStr(self, prec=5, fmt=Fmt.PAREN, sep=_COMMASPACE_):  # PYCHOK expected
        '''Return a string representation of this vector.

           @kwarg prec: Number of decimal places (C{int}).
           @kwarg fmt: Enclosing format to use (C{str}).
           @kwarg sep: Separator between components (C{str}).

           @return: Vector as "(x, y, z)" (C{str}).
        '''
        t = sep.join(strs(self.xyz3, prec=prec))
        return (fmt % (t,)) if fmt else t

    def unit(self, ll=None):
        '''Normalize this vector to unit length.

           @kwarg ll: Optional, original location (C{LatLon}).

           @return: Normalized vector (L{Vector3d}).
        '''
        u = self._united
        if ll:
            u._fromll = ll
        return u

    @Property_RO
    def _united(self):  # __dict__ value overwritten below
        '''(INTERNAL) Get normalized vector (L{Vector3d}).
        '''
        n = self.length
        if n > EPS0 and fabs(n - _1_0) > EPS0:
            u = self._xnamed(self.dividedBy(n))
            u._update(False, length=_1_0, length2=_1_0, _united=u)
        else:
            u = self.copy()
            u._update(False, _united=u)
        if self._fromll:
            u._fromll = self._fromll
        return u

    @Property
    def x(self):
        '''Get the X component (C{float}).
        '''
        return self._x

    @x.setter  # PYCHOK setter!
    def x(self, x):
        '''Set the X component, if different (C{float}).
        '''
        x = Float(x=x)
        if self._x != x:
            _update_all(self, needed=3)
            self._x = x

    @Property
    def xyz(self):
        '''Get the X, Y and Z components (L{Vector3Tuple}C{(x, y, z)}).
        '''
        return _MODS.namedTuples.Vector3Tuple(*self.xyz3, name=self.name)

    @xyz.setter  # PYCHOK setter!
    def xyz(self, xyz):
        '''Set the X, Y and Z components (C{Cartesian}, L{Ecef9Tuple},
           C{Nvector}, L{Vector3d}, L{Vector3Tuple}, L{Vector4Tuple}
           or a C{tuple} or C{list} of 3+ C{scalar} items).
        '''
        self._xyz(xyz)

    def _xyz(self, x_xyz, *y_z):
        '''(INTERNAL) Set the C{_x}, C{_y} and C{_z} attributes.
        '''
        _update_all(self, needed=3)
        self._x, self._y, self._z = _xyz3(_xyz_, x_xyz, *y_z)
        return self

    @property_RO
    def xyz3(self):
        '''Get the X, Y and Z components as C{3-tuple}.
        '''
        return self.x, self.y, self.z

    @property_RO
    def x2y2z2(self):
        '''Get the X, Y and Z components I{squared} (3-tuple C{(x**2, y**2, z**2)}).
        '''
        return self.x**2, self.y**2, self.z**2

    @Property
    def y(self):
        '''Get the Y component (C{float}).
        '''
        return self._y

    @y.setter  # PYCHOK setter!
    def y(self, y):
        '''Set the Y component, if different (C{float}).
        '''
        y = Float(y=y)
        if self._y != y:
            _update_all(self, needed=3)
            self._y = y

    @Property
    def z(self):
        '''Get the Z component (C{float}).
        '''
        return self._z

    @z.setter  # PYCHOK setter!
    def z(self, z):
        '''Set the Z component, if different (C{float}).
        '''
        z = Float(z=z)
        if self._z != z:
            _update_all(self, needed=3)
            self._z = z


def _xyz3(where, x_xyz, *y_z):  # in .cartesianBase._rtp3
    '''(INTERNAL) Get , Y and Z as 3-tuple C{(x, y, z)}.
    '''
    try:
        xyz3 = map1(_float0, x_xyz, *y_z) if y_z else (  # islistuple for Vector*Tuple
               map2(_float0, x_xyz[:3]) if islistuple(x_xyz, minum=3) else
               x_xyz.xyz)  # .xyz3
    except (AttributeError, TypeError, ValueError) as x:
        raise _xError(x, unstr(where, x_xyz, *y_z))
    return xyz3


__all__ += _ALL_DOCS(Vector3dBase)

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
