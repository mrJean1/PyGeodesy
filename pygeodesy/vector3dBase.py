
# -*- coding: utf-8 -*-

u'''I{Veness}' 3-D vector base class C{Vector3dBase}.

Pure Python implementation of vector-based functions by I{(C) Chris Veness
2011-2015} published under the same MIT Licence**, see U{Vector-based geodesy
<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.
'''

from pygeodesy.basics import copysign0, isnear0, isnear1, isscalar, map1, \
                            _xinstanceof
from pygeodesy.errors import CrossError, _InvalidError, _IsnotError, VectorError
from pygeodesy.fmath import euclid_, fdot, hypot_, hypot2_, _sys_version_info2
from pygeodesy.formy import n_xyz2latlon, n_xyz2philam, sincos2
from pygeodesy.interns import EPS, EPS0, NN, PI, PI2, _coincident_, \
                             _colinear_, _COMMASPACE_, _1_0
# from pygeodesy.lazily import _sys_version_info2  # from .fmath
from pygeodesy.named import _NamedBase, _NotImplemented, _xother3
from pygeodesy.namedTuples import Vector3Tuple, Vector4Tuple
from pygeodesy.props import deprecated_method, Property, Property_RO, \
                            property_doc_
from pygeodesy.streprs import Fmt, strs
from pygeodesy.units import Float, Scalar
# from pygeodesy.utily import sincos2  # from .formy

from math import atan2

__all__ = ()
__version__ = '21.12.18'


class Vector3dBase(_NamedBase):  # XXX or _NamedTuple or Vector3Tuple?
    '''(INTERNAL) Generic 3-D vector base class.

       In a geodesy context, these may be used to represent:
        - n-vector representing a normal to point on earth's surface
        - earth-centered, earth-fixed cartesian (= spherical n-vector)
        - great circle normal to vector
        - motion vector on earth's surface
        - etc.
    '''
    _crosserrors = True  # un/set by .errors.crosserrors

    _fromll = None  # original latlon, '_fromll'
    _x      = 0     # X component
    _y      = 0     # Y component
    _z      = 0     # Z component

    def __init__(self, x_xyz, y=0, z=0, ll=None, name=NN):
        '''New L{Vector3d} or C{Vector3dBase} instance.

           The vector may be normalised or use x, y, z for position and
           distance from earth centre or height relative to the surface
           of the earth' sphere or ellipsoid.

           @arg x_xyz: X component of vector (C{scalar}) or (3-D) vector
                       (C{Cartesian}, C{Nvector}, L{Vector3d} or L{Vector3Tuple}).
           @kwarg y: Y component of vector (C{scalar}), ignored if B{C{x_xyz}}
                     is not C{scalar}, otherwise same units as B{C{x_xyz}}.
           @kwarg z: Z component of vector (C{scalar}), ignored if B{C{x_xyz}}
                     is not C{scalar}, otherwise same units as B{C{x_xyz}}.
           @kwarg ll: Optional latlon reference (C{LatLon}).
           @kwarg name: Optional name (C{str}).

           @raise VectorError: Invalid B{C{x_xyz}}.
        '''
        try:
            self._x, self._y, self._z = (x_xyz, y, z) if isscalar(x_xyz) else x_xyz.xyz
        except (AttributeError, TypeError, ValueError) as x:
            raise VectorError(x=x_xyz, y=y, z=z, txt=str(x))
        if ll:
            self._fromll = ll
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

    __radd__ = __add__  # PYCHOK no cover

    def __bool__(self):  # PYCHOK PyChecker
        '''Is this vector non-zero?
        '''
        return bool(self.x or self.y or self.z)

    def __cmp__(self, other):  # Python 2-
        '''Compare this and an other vector

           @arg other: The other vector (L{Vector3d}).

           @return: -1, 0 or +1 (C{int}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)
        return -1 if self.length < other.length else (
               +1 if self.length > other.length else 0)

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
#       self.others(other)
        return self.isequalTo(other, eps=EPS0)

    def __floordiv__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __format__(self, *other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, *other)

    def __ge__(self, other):
        '''Is this vector longer than or equal to an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if so, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)
        return self.length >= other.length

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
        self.others(other)
        return self.length > other.length

    def __iadd__(self, other):
        '''Add an other vector to this one I{in-place}, C{this += B{other}}.

           @arg other: The other vector (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.xyz = self.plus(other).xyz

    def __imatmul__(self, other):  # PYCHOK Python 3.5+
        '''Cross multiply an other vector and this one I{in-place}, C{this @= B{other}}.

           @arg other: The other vector (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.

           @see: Luciano Ramalho, "Fluent Python", page 397-398, O'Reilly 2016.
        '''
        self.xyz = self.cross(other).xyz

    def __imul__(self, scalar):
        '''Multiply this vector by a scalar I{in-place}, C{this *= B{scalar}}.

           @arg scalar: Factor (C{scalar}).

           @raise TypeError: Non-scalar B{C{scalar}}.
        '''
        self.xyz = self.times(scalar).xyz

    def __isub__(self, other):
        '''Subtract an other vector from this one I{in-place}, C{this -= B{other}}.

           @arg other: The other vector (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.xyz = self.minus(other).xyz

    def __itruediv__(self, scalar):
        '''Divide this vector by a scalar I{in-place}, C{this /= B{scalar}}.

           @arg scalar: The divisor (C{scalar}).

           @raise TypeError: Non-scalar B{C{scalar}}.
        '''
        self.xyz = self.dividedBy(scalar).xyz

    def __le__(self, other):  # Python 3+
        '''Is this vector shorter than or equal to an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if so, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)
        return self.length <= other.length

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
        self.others(other)
        return self.length < other.length

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
        '''Multiply this vector by a scalar.

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
        return not self.__eq__(other)

    def __pos__(self):  # PYCHOK no cover
        '''Return this vector I{as-is}.

           @return: This instance (L{Vector3d})
        '''
        return self

    def __pow__(self, other, *mod):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other, *mod)

    def __rdivmod__ (self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rfloordiv__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rmatmul__(self, other):  # PYCHOK Python 3.5+
        '''Compute the cross product of an other and this vector, C{B{other} @ this}.

           @arg other: The other vector (L{Vector3d}).

           @return: Cross product (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)
        return other.cross(self)

    def __rmod__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rmul__(self, scalar):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, scalar)

    def __round__(self, ndigits=None):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, ndigits=ndigits)

    def __rpow__(self, other, *mod):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other, *mod)

    def __rsub__(self, other):  # PYCHOK no cover
        '''Subtract this vector from an other vector.

           @arg other: The other vector (L{Vector3d}).

           @return: Difference (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)
        return other.minus(self)

    def __rtruediv__(self, scalar):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, scalar)

    def __sub__(self, other):
        '''Subtract an other vector from this vector.

           @arg other: The other vector (L{Vector3d}).

           @return: Difference (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.minus(other)

    def __truediv__(self, scalar):
        '''Divide this vector by a scalar.

           @arg scalar: The divisor (C{scalar}).

           @return: Quotient (L{Vector3d}).

           @raise TypeError: Non-scalar B{C{scalar}}.
        '''
        return self.dividedBy(scalar)

    if _sys_version_info2 < (3, 0):  # PYCHOK no cover
        # <https://docs.Python.org/2/library/operator.html#mapping-operators-to-functions>
        __div__     = __truediv__
        __idiv__    = __itruediv__
        __nonzero__ = __bool__

    def angleTo(self, other, vSign=None, wrap=False):
        '''Compute the angle between this and an other vector.

           @arg other: The other vector (L{Vector3d}).
           @kwarg vSign: Optional vector, if supplied (and out of the
                         plane of this and the other), angle is signed
                         positive if this->other is clockwise looking
                         along vSign or negative in opposite direction,
                         otherwise angle is unsigned.
           @kwarg wrap: Wrap/unroll the angle to +/-PI (C{bool}).

           @return: Angle (C{radians}).

           @raise TypeError: If B{C{other}} or B{C{vSign}} not a L{Vector3d}.
        '''
        x = self.cross(other)
        s = x.length
        # use vSign as reference to set sign of s
        if s and vSign and x.dot(vSign) < 0:
            s = -s

        a = atan2(s, self.dot(other))
        if wrap and abs(a) > PI:
            a -= copysign0(PI2, a)
        return a

    def apply(self, fun2, other_x, *y_z, **fun2_kwds):
        '''Apply a function component-wise to this and an other vector.

           @arg fun2: 2-Argument callable (C{any(scalar, scalar}).
           @arg other_x: An other vector factors (L{Vector3dBase},
                         L{Vector3Tuple}, L{Vector4Tuple},
                         L{Ecef9Tuple} cartesian) or X scale
                         factor (C{scalar}).
           @arg y_z: Y and Z scale factors (C{scalar}, C{scalar}).
           @kwarg fun2_kwds: Optional keyword arguments for B{C{fun2}}.

           @return: New, applied vector (L{Vector3d}).

           @raise ValueError: Invalid B{C{other_x}} or B{C{y_z}}.
        '''
        if not callable(fun2):
            raise _IsnotError(callable.__name__, fun2=fun2)

        if fun2_kwds:
            def _f2(a, b):
                return fun2(a, b, **fun2_kwds)
        else:
            _f2 = fun2

        xyz = _other_xyz3(other_x, y_z)
        xyz = (_f2(a, b) for a, b in zip(self.xyz, xyz))
        return self.classof(*xyz)

    def cross(self, other, raiser=None):  # raiser=NN
        '''Compute the cross product of this and an other vector.

           @arg other: The other vector (L{Vector3d}).
           @kwarg raiser: Optional, L{CrossError} label if raised (C{str}).

           @return: Cross product (L{Vector3d}).

           @raise CrossError: Zero or near-zero cross product and both
                              B{C{raiser}} and L{pygeodesy.crosserrors} set.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)

        x = self.y * other.z - self.z * other.y
        y = self.z * other.x - self.x * other.z
        z = self.x * other.y - self.y * other.x

        if raiser and self.crosserrors and max(map1(abs, x, y, z)) < EPS:
            t = _coincident_ if self.isequalTo(other) else _colinear_
            r = getattr(other, '_fromll', None) or other
            raise CrossError(raiser, r, txt=t)

        return self.classof(x, y, z)

    @property_doc_('''raise or ignore L{CrossError} exceptions (C{bool}).''')
    def crosserrors(self):
        '''Get L{CrossError} exceptions (C{bool}).
        '''
        return self._crosserrors

    @crosserrors.setter  # PYCHOK setter!
    def crosserrors(self, raiser):
        '''Raise L{CrossError} exceptions (C{bool}).
        '''
        self._crosserrors = bool(raiser)

    def dividedBy(self, factor):
        '''Divide this vector by a scalar.

           @arg factor: The divisor (C{scalar}).

           @return: New, scaled vector (L{Vector3d}).

           @raise TypeError: Non-scalar B{C{factor}}.

           @raise VectorError: Invalid or zero B{C{factor}}.
        '''
        f = Scalar(factor=factor)
        try:
            return self.times(_1_0 / f)
        except (ValueError, ZeroDivisionError) as x:
            raise VectorError(factor=factor, txt=str(x))

    def dot(self, other):
        '''Compute the dot (scalar) product of this and an other vector.

           @arg other: The other vector (L{Vector3d}).

           @return: Dot product (C{float}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        if other is self:
            d = self.length2
        else:
            self.others(other)
            d = fdot(self.xyz, *other.xyz)
        return d

    @deprecated_method
    def equals(self, other, units=False):  # PYCHOK no cover
        '''DEPRECATED, use method C{isequalTo}.
        '''
        return self.isequalTo(other, units=units)

    @Property_RO
    def euclid(self):
        '''Approximate the length (norm, magnitude) of this vector (C{Float}).

           @see: Properties C{length} and C{length2} and function
                 L{pygeodesy.euclid_}.
        '''
        return Float(euclid=euclid_(self.x, self.y, self.z))

    def equirectangular(self, other):
        '''Approximate the different between this and an other vector.

           @arg other: Vector to subtract (C{Vector3dBase}).

           @return: The lenght I{squared} of the difference (C{Float}).

           @raise TypeError: Incompatible B{C{other}} C{type}.

           @see: Property C{length2}.
        '''
        d = self.minus(other)
        return Float(equirectangular=hypot2_(d.x, d.y, d.z))

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
        if isnear0(f):
            r = self
        elif isnear1(f):
            r = self.others(other)
        else:
            d = self.others(other).minus(self)
            r = self.plus(d.times(f))
        return r

    def isconjugateTo(self, other, minum=1, eps=EPS):
        '''Determine whether this and an other vector are conjugates.

           @arg other: The other vector (L{Vector3d}, C{Vector3Tuple} or
                       C{Vector4Tuple}).
           @kwarg minum: Minimal number of conjugates (C{int}, 0..3).
           @kwarg eps: Tolerance for equality and conjugation (C{scalar}),
                       same units as C{x}, C{y}, and C{z}.

           @return: C{True} if C{x}, C{y} and C{z} of this match the other
                    vector's or at least C{B{minum}} have opposite signs.

           @raise TypeError: Incompatible B{C{other}} C{type}.

           @see: Method C{isequalTo}.
        '''
        self.others(other)

        n = 0
        for a, b in zip(self.xyz, other.xyz):
            if abs(a + b) < eps and (a * b) < 0:
                n += 1  # conjugate
            elif abs(a - b) > eps:
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
        return max(map(abs, d.xyz)) < eps

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

    def minus(self, other):
        '''Subtract an other vector from this vector.

           @arg other: The other vector (L{Vector3d}).

           @return: New vector difference (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)

        return self.classof(self.x - other.x,
                            self.y - other.y,
                            self.z - other.z)

    def negate(self):
        '''Return this vector in opposite direction.

           @return: New, opposite vector (L{Vector3d}).
        '''
        return self.classof(-self.x, -self.y, -self.z)

    __neg__ = negate  # PYCHOK no cover

    @Property_RO
    def _N_vector(self):
        '''(INTERNAL) Get the (C{nvectorBase._N_vector_})
        '''
        from pygeodesy.nvectorBase import _N_vector_
        return _N_vector_(*self.xyz, name=self.name)

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
        self.others(other)

        return self.classof(self.x + other.x,
                            self.y + other.y,
                            self.z + other.z)

    sum = plus  # alternate name

    def rotate(self, axis, theta):
        '''Rotate this vector around an axis by a specified angle.

           See U{Rotation matrix from axis and angle
           <https://WikiPedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle>}
           and U{Quaternion-derived rotation matrix
           <https://WikiPedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix>}.

           @arg axis: The axis being rotated around (L{Vector3d}).
           @arg theta: The angle of rotation (C{radians}).

           @return: New, rotated vector (L{Vector3d}).

           @JSname: I{rotateAround}.
        '''
        p = self.unit().xyz  # point being rotated
        r = self.others(axis=axis).unit()  # axis being rotated around

        s, c = sincos2(theta)  # rotation angle

        ax, ay, az = r.xyz  # quaternion-derived rotation matrix
        bx, by, bz = r.times(_1_0 - c).xyz
        sx, sy, sz = r.times(s).xyz
        return self.classof(fdot(p, ax * bx + c,  ax * by - sz, ax * bz + sy),
                            fdot(p, ay * bx + sz, ay * by + c,  ay * bz - sx),
                            fdot(p, az * bx - sy, az * by + sx, az * bz + c))

    @deprecated_method
    def rotateAround(self, axis, theta):  # PYCHOK no cover
        '''DEPRECATED, use method C{rotate}.'''
        return self.rotate(axis, theta)

    def times(self, factor):
        '''Multiply this vector by a scalar.

           @arg factor: Scale factor (C{scalar}).

           @return: New, scaled vector (L{Vector3d}).

           @raise TypeError: Non-scalar B{C{factor}}.
        '''
        s = Scalar(factor=factor)
        return self.classof(self.x * s, self.y * s, self.z * s)

    def times_(self, other_x, *y_z):
        '''Multiply this vector's components by separate scalars.

           @arg other_x: An other vector factors (L{Vector3dBase},
                         L{Vector3Tuple}, L{Vector4Tuple},
                         L{Ecef9Tuple} cartesian) or X scale
                         factor (C{scalar}).
           @arg y_z: Y and Z scale factors (C{scalar}, C{scalar}),
                     ignored if B{C{other_x}} is not C{scalar}.

           @return: New, scaled vector (L{Vector3d}).

           @raise ValueError: Invalid B{C{other_x}} or B{C{y_z}}.
        '''
        x, y, z = _other_xyz3(other_x, y_z)
        return self.classof(self.x * x, self.y * y, self.z * z)

    @deprecated_method
    def to2ab(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{Nvector.philam}.

           @return: A L{PhiLam2Tuple}C{(phi, lam)}.
        '''
        return n_xyz2philam(self.x, self.y, self.z)

    @deprecated_method
    def to2ll(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{Nvector.latlon}.

           @return: A L{LatLon2Tuple}C{(lat, lon)}.
        '''
        return n_xyz2latlon(self.x, self.y, self.z)

    @deprecated_method
    def to3xyz(self):  # PYCHOK no cover
        '''DEPRECATED, use property L{xyz}.
        '''
        return self.xyz

    def toStr(self, prec=5, fmt=Fmt.PAREN, sep=_COMMASPACE_):  # PYCHOK expected
        '''Return a string representation of this vector.

           @kwarg prec: Optional number of decimal places (C{int}).
           @kwarg fmt: Optional, enclosing format to use (C{str}).
           @kwarg sep: Optional separator between components (C{str}).

           @return: Vector as "(x, y, z)" (C{str}).
        '''
        t = sep.join(strs(self.xyz, prec=prec))
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
        if n > EPS0 and abs(n - _1_0) > EPS0:
            u = self._xnamed(self.dividedBy(n))
            u._overwrite(length=_1_0, length2=_1_0, _united=u)
        else:
            u = self.copy()
            u._overwrite(_united=u)
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
        if self.x != x:
            self._update(True)
            self._x = x

    @Property
    def xyz(self):
        '''Get the X, Y and Z components (L{Vector3Tuple}C{(x, y, z)}).
        '''
        return Vector3Tuple(self.x, self.y, self.z, name=self.name)

    @xyz.setter  # PYCHOK setter!
    def xyz(self, xyz):
        '''Set the X, Y and Z components, if different (L{Vector3dBase}, L{Vector3Tuple} or L{Vector4Tuple}).
        '''
        _xinstanceof(Vector3dBase, Vector3Tuple, Vector4Tuple, xyz=xyz)
        t3 = xyz.x,  xyz.y,  xyz.z
        if (self.x, self.y, self.z) != t3:
            self._update(True)
            self.x, self.y, self.z = t3

    @Property_RO
    def x2y2z2(self):
        '''Get the X, Y and Z components I{squared} (C{Vector3Tuple}C{(x2, y2, z2)}).
        '''
        return Vector3Tuple(self.x**2, self.y**2, self.z**2, name=self.name)

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
        if self.y != y:
            self._update(True)
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
        if self.z != z:
            self._update(True)
            self._z = z


def _other_xyz3(other_x, y_z):
    '''(INTERNAL) Helper.
    '''
    try:
        if y_z:
            y, z = y_z
            x = Scalar(x=other_x)
        else:
            x, y, z = other_x.x, other_x.y, other_x.z
    except (AttributeError, TypeError, ValueError) as x:
        raise _InvalidError(other_x=other_x, y_z=y_z, txt=str(x))
    return x, y, z

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
