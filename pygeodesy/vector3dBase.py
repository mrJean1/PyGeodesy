
# -*- coding: utf-8 -*-

u'''I{Veness}' 3-D vector base class C{Vector3dBase}.

Pure Python implementation of vector-based functions by I{(C) Chris Veness
2011-2015} published under the same MIT Licence**, see U{Vector-based geodesy
<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.
'''

from pygeodesy.basics import copysign0, map1
from pygeodesy.errors import CrossError, VectorError
from pygeodesy.fmath import euclid_, fdot, hypot_, hypot2_
from pygeodesy.formy import n_xyz2latlon, n_xyz2philam, sincos2
from pygeodesy.interns import EPS, EPS0, NN, PI, PI2, _coincident_, \
                             _colinear_, _COMMASPACE_, _1_0
from pygeodesy.named import _NamedBase, _xother3
from pygeodesy.namedTuples import Vector3Tuple  # Vector4Tuple
from pygeodesy.props import deprecated_method, Property_RO, property_doc_
from pygeodesy.streprs import Fmt, strs
from pygeodesy.units import Float, Scalar
# from pygeodesy.utily import sincos2  # from .formy

from math import atan2

__all__ = ()
__version__ = '21.07.02'


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

    def __init__(self, x, y, z, ll=None, name=NN):
        '''New L{Vector3d} or C{Vector3dBase} instance.

           The vector may be normalised or use x/y/z values
           for height relative to the surface of the sphere
           or ellipsoid, distance from earth centre, etc.

           @arg x: X component of vector (C{scalar}).
           @arg y: Y component of vector (C{scalar}).
           @arg z: Z component of vector (C{scalar}).
           @kwarg ll: Optional latlon reference (C{LatLon}).
           @kwarg name: Optional name (C{str}).
        '''
        self._x = x
        self._y = y
        self._z = z
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
#   __iadd__ = __add__
    __radd__ = __add__

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

    def __div__(self, scalar):
        '''Divide this vector by a scalar.

           @arg scalar: The divisor (C{scalar}).

           @return: Quotient (L{Vector3d}).

           @raise TypeError: Non-scalar B{C{scalar}}.
        '''
        return self.dividedBy(scalar)
#   __itruediv__ = __div__
    __truediv__ = __div__

    def __eq__(self, other):
        '''Is this vector equal to an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if equal, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)
        return self.isequalTo(other)

    def __ge__(self, other):
        '''Is this vector longer than or equal to an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if so, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)
        return self.length >= other.length

    def __gt__(self, other):
        '''Is this vector longer than an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if so, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)
        return self.length > other.length

    def __le__(self, other):  # Python 3+
        '''Is this vector shorter than or equal to an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if so, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)
        return self.length <= other.length

    def __lt__(self, other):  # Python 3+
        '''Is this vector shorter than an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if so, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)
        return self.length < other.length

    # Luciano Ramalho, "Fluent Python", page 397, O'Reilly 2016
    def __matmul__(self, other):  # PYCHOK Python 3.5+ ... c = a @ b
        '''Compute the cross product of this and an other vector.

           @arg other: The other vector (L{Vector3d}).

           @return: Cross product (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.cross(other)
#   __imatmul__ = __matmul__

    def __mul__(self, scalar):
        '''Multiply this vector by a scalar

           @arg scalar: Factor (C{scalar}).

           @return: Product (L{Vector3d}).
        '''
        return self.times(scalar)
#   __imul__ = __mul__
#   __rmul__ = __mul__

    def __ne__(self, other):
        '''Is this vector not equal to an other vector?

           @arg other: The other vector (L{Vector3d}).

           @return: C{True} if so, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)
        return not self.isequalTo(other)

    def __neg__(self):
        '''Negate this vector.

           @return: Negative (L{Vector3d})
        '''
        return self.negate()

    def __pos__(self):
        '''Copy this vector.

           @return: Positive (L{Vector3d})
        '''
        return self.copy()

    # Luciano Ramalho, "Fluent Python", page 397, O'Reilly 2016
    def __rmatmul__(self, other):  # PYCHOK Python 3.5+ ... c = a @ b
        '''Compute the cross product of an other and this vector.

           @arg other: The other vector (L{Vector3d}).

           @return: Cross product (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)
        return other.cross(self)

    def __rsub__(self, other):
        '''Subtract this vector from an other vector.

           @arg other: The other vector (L{Vector3d}).

           @return: Difference (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)
        return other.minus(self)

    def __sub__(self, other):
        '''Subtract an other vector from this vector.

           @arg other: The other vector (L{Vector3d}).

           @return: Difference (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.minus(other)
#   __isub__ = __sub__

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

    def cross(self, other, raiser=None):  # raiser=NN
        '''Compute the cross product of this and an other vector.

           @arg other: The other vector (L{Vector3d}).
           @kwarg raiser: Optional, L{CrossError} label if raised (C{str}).

           @return: Cross product (L{Vector3d}).

           @raise CrossError: Zero or near-zero cross product and both
                              B{C{raiser}} and L{crosserrors} set.

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
        '''Approximate the length (norm, magnitude) of this vector.

           @see: Function L{euclid_} and properties C{length} and C{length2}.
        '''
        return Float(euclid=euclid_(self.x, self.y, self.z))

    def isequalTo(self, other, units=False, eps=EPS):
        '''Check if this and an other vector are equal or equivalent.

           @arg other: The other vector (L{Vector3d}).
           @kwarg units: Optionally, compare the normalized, unit
                         version of both vectors.
           @kwarg eps: Tolerance (C{scalar}), same units as C{x},
                       C{y}, and C{z}.

           @return: C{True} if vectors are identical, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)

        if units:
            d = self.unit().minus(other.unit())
        else:
            d = self.minus(other)
        return max(map(abs, d.xyz)) < eps

    @Property_RO
    def length(self):  # __dict__ value overwritten by Property_RO C{_united}
        '''Get the length (norm, magnitude) of this vector (C{float}).

           @see: Properties L{length2} and L{euclid}.
        '''
        return Float(length=hypot_(self.x, self.y, self.z))

    @Property_RO
    def length2(self):  # __dict__ value overwritten by Property_RO C{_united}
        '''Get the length I{squared} of this vector (C{float}).

           @see: Property L{length}.
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
        s, c = sincos2(theta)

        a = self.others(axis=axis).unit()  # axis being rotated around
        b = a.times(_1_0 - c)
        s = a.times(s)

        p = self.unit().xyz  # point being rotated
        # multiply p by a quaternion-derived rotation matrix
        ax, ay, az = a.x, a.y, a.z
        bx, by, bz = b.x, b.y, b.z
        sx, sy, sz = s.x, s.y, s.z
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
        f = Scalar(factor=factor)
        return self.classof(self.x * f, self.y * f, self.z * f)

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

    @Property_RO
    def x(self):
        '''Get the X component (C{float}).
        '''
        return self._x

    @Property_RO
    def xyz(self):
        '''Get the X, Y and Z components (L{Vector3Tuple}C{(x, y, z)}).
        '''
        return Vector3Tuple(self.x, self.y, self.z, name=self.name)

    @Property_RO
    def y(self):
        '''Get the Y component (C{float}).
        '''
        return self._y

    @Property_RO
    def z(self):
        '''Get the Z component (C{float}).
        '''
        return self._z

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
