
# -*- coding: utf-8 -*-

u'''Generic 3-D vector base class L{Vector3d} and function L{sumOf}.

Pure Python implementation of vector-based functions by I{(C) Chris
Veness 2011-2015} published under the same MIT Licence**, see
U{Vector-based geodesy
<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.

@newfield example: Example, Examples
'''

from pygeodesy.basics import EPS, _isnotError, isscalar, len2, \
                             map1, property_doc_, property_RO
from pygeodesy.fmath import fdot, fsum, hypot_
from pygeodesy.formy import n_xyz2latlon, n_xyz2philam
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _NamedBase, Vector3Tuple
from pygeodesy.streprs import strs

from math import atan2, cos, sin

# all public constants, classes and functions
__all__ = _ALL_LAZY.vector3d + ('Vector3d', 'sumOf')
__version__ = '20.04.02'


def _xyzn4(xyz, y, z, Error=TypeError):  # imported by .ecef.py
    '''(INTERNAL) Get an C{(x, y, z, name)} 4-tuple.
    '''
    try:
        t = xyz.x, xyz.y, xyz.z
    except AttributeError:
        t = xyz, y, z
    try:
        x, y, z = map1(float, *t)
    except (TypeError, ValueError) as x:
        raise Error('%s invalid: %r, %s' % ('xyz, y or z', t, x))

    return x, y, z, getattr(xyz, 'name', '')


def _xyzhdn6(xyz, y, z, height, datum, ll, Error=TypeError):  # .cartesianBase.py, .nvectorBase.py
    '''(INTERNAL) Get an C{(x, y, z, h, d, name)} 6-tuple.
    '''
    x, y, z, n = _xyzn4(xyz, y, z, Error=Error)

    h = height or getattr(xyz, 'height', None) \
               or getattr(xyz, 'h', None) \
               or getattr(ll,  'height', None)

    d = datum or getattr(xyz, 'datum', None) \
              or getattr(ll,  'datum', None)

    return x, y, z, h, d, n


class CrossError(ValueError):
    '''Error raised for zero or near-zero vectorial cross products,
       occurring for coincident or colinear points, paths or bearings.
    '''
    pass


def crosserrors(raiser=None):
    '''Report or ignore vectorial cross product errors.

       @kwarg raiser: Use C{True} to throw or C{False} to ignore
                      L{CrossError} exceptions.  Use C{None} to
                      leave the setting unchanged.

       @return: Previous setting (C{bool}).
    '''
    t = Vector3d._crosserrors
    if raiser in (True, False):
        Vector3d._crosserrors = raiser
    return t


class VectorError(ValueError):
    '''Issue with class L{Vector3d} or C{*Nvector}.
    '''
    pass


class Vector3d(_NamedBase):  # XXX or _NamedTuple or Vector3Tuple?
    '''Generic 3-D vector manipulation.

       In a geodesy context, these may be used to represent:
        - n-vector representing a normal to point on earth's surface
        - earth-centered, earth-fixed vector (= n-vector for spherical model)
        - great circle normal to vector
        - motion vector on earth's surface
        - etc.
    '''
    _crosserrors = True  # set by function crosserrors above

    _fromll = None  #: (INTERNAL) original ll.
    _length = None  #: (INTERNAL) cached length.
    _united = None  #: (INTERNAL) cached norm, unit (L{Vector3d}).
    _xyz    = None  #: (INTERNAL) cached xyz (L{Vector3Tuple}).

    _x = 0  #: (INTERNAL) X component.
    _y = 0  #: (INTERNAL) Y component.
    _z = 0  #: (INTERNAL) Z component.

    def __init__(self, x, y, z, ll=None, name=''):
        '''New 3-D vector.

           The vector may be normalised, or use x/y/z values for
           height relative to the surface of the sphere or ellipsoid,
           distance from earth centre, etc.

           @arg x: X component of vector (C{scalar}).
           @arg y: Y component of vector (C{scalar}).
           @arg z: Z component of vector (C{scalar}).
           @kwarg ll: Optional, original latlon (C{LatLon}).
           @kwarg name: Optional name (C{str}).
        '''
        self._x = x
        self._y = y
        self._z = z
        if ll:
            self._fromll = ll
        if name:
            self.name = name

    def __add__(self, other):
        '''Add this to an other vector (L{Vector3d}).

           @return: Vectorial sum (L{Vector3d}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        return self.plus(other)
#   __iadd__ = __add__
    __radd__ = __add__

    def __abs__(self):
        '''Return the norm of this vector.

           @return: Norm, unit length (C{float});
        '''
        return self.length

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

    def _update(self, updated, *attrs):
        '''(INTERNAL) Zap cached attributes if updated.
        '''
        if updated:
            _NamedBase._update(self, updated, '_length',  # '_fromll'
                                   '_united', '_xyz', *attrs)

    def angleTo(self, other, vSign=None):
        '''Compute the angle between this and an other vector.

           @arg other: The other vector (L{Vector3d}).
           @kwarg vSign: Optional vector, if supplied (and out of the
                         plane of this and the other), angle is signed
                         positive if this->other is clockwise looking
                         along vSign or negative in opposite direction,
                         otherwise angle is unsigned.

           @return: Angle (C{radians}).

           @raise TypeError: If B{C{other}} or B{C{vSign}} not a L{Vector3d}.
        '''
        x = self.cross(other)
        s = x.length
        if s < EPS:
            return 0.0
        # use vSign as reference to get sign of s
        if vSign and x.dot(vSign) < 0:
            s = -s
        return atan2(s, self.dot(other))

    def cross(self, other, raiser=None):
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
            t = 'coincident' if self.isequalTo(other) else 'colinear'
            r = getattr(other, '_fromll', None) or other
            raise CrossError('%s %s: %r' % (t, raiser, r))

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
        if not isscalar(factor):
            raise _isnotError('scalar', factor=factor)
        try:
            return self.times(1.0 / factor)
        except (ValueError, ZeroDivisionError):
            raise VectorError('%s invalid: %r' % ('factor', factor))

    def dot(self, other):
        '''Compute the dot (scalar) product of this and an other vector.

           @arg other: The other vector (L{Vector3d}).

           @return: Dot product (C{float}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        self.others(other)

        return fdot(self.xyz, *other.xyz)

    def equals(self, other, units=False):  # PYCHOK no cover
        '''DEPRECATED, use method C{isequalTo}.
        '''
        return self.isequalTo(other, units=units)

    def isequalTo(self, other, units=False, eps=EPS):
        '''Check if this and an other vector are equal or equivalent.

           @arg other: The other vector (L{Vector3d}).
           @kwarg units: Optionally, compare the normalized, unit
                         version of both vectors.
           @kwarg eps: Tolerance for equality (C{scalar}).

           @return: C{True} if vectors are identical, C{False} otherwise.

           @raise TypeError: Incompatible B{C{other}} C{type}.

           @example:

           >>> v1 = Vector3d(52.205, 0.119)
           >>> v2 = Vector3d(52.205, 0.119)
           >>> e = v1.isequalTo(v2)  # True
        '''
        self.others(other)

        if units:
            d = self.unit().minus(other.unit())
        else:
            d = self.minus(other)
        return max(map(abs, d.xyz)) < eps

    @property_RO
    def length(self):
        '''Get the length (norm, magnitude) of this vector (C{float}).
        '''
        if self._length is None:
            self._length = hypot_(self.x, self.y, self.z)
        return self._length

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

    @property_RO
    def _N_vector(self):
        '''(INTERNAL) Get the (C{nvectorBase._N_vector_})
        '''
        from pygeodesy.nvectorBase import _N_vector_
        return _N_vector_(*(self._xyz or self.xyz))

    def others(self, other, name='other'):
        '''Refined class comparison.

           @arg other: The other vector (L{Vector3d}).
           @kwarg name: Optional, other's name (C{str}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        try:
            _NamedBase.others(self, other, name=name)
        except TypeError:
            if not isinstance(other, Vector3d):
                raise

    def parse(self, str3d, sep=','):
        '''Parse an C{"x, y, z"} string.

           @arg str3d: X, y and z values (C{str}).
           @kwarg sep: Optional separator (C{str}).

           @return: New vector (L{Vector3d}).

           @raise VectorError: Invalid B{C{str3d}}.
        '''
        try:
            v = [float(v.strip()) for v in str3d.split(sep)]
            if len(v) != 3:
                raise ValueError
        except (TypeError, ValueError):
            raise VectorError('%s invalid: %r' % ('str3d', str3d))

        return self.classof(*v)

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
        self.others(axis, name='axis')

        c = cos(theta)
        a = axis.unit()  # axis being rotated around
        b = a.times(1 - c)
        s = a.times(sin(theta))

        p = self.unit().xyz  # point being rotated

        # multiply p by a quaternion-derived rotation matrix
        return self.classof(fdot(p, a.x * b.x + c,   a.x * b.y - s.z, a.x * b.z + s.y),
                            fdot(p, a.y * b.x + s.z, a.y * b.y + c,   a.y * b.z - s.x),
                            fdot(p, a.z * b.x - s.y, a.z * b.y + s.x, a.z * b.z + c))

    def rotateAround(self, axis, theta):  # PYCHOK no cover
        '''DEPRECATED, use method C{rotate}.
        '''
        return self.rotate(axis, theta)

    def times(self, factor):
        '''Multiply this vector by a scalar.

           @arg factor: Scale factor (C{scalar}).

           @return: New, scaled vector (L{Vector3d}).

           @raise TypeError: Non-scalar B{C{factor}}.
        '''
        if not isscalar(factor):
            raise _isnotError('scalar', factor=factor)
        return self.classof(self.x * factor,
                            self.y * factor,
                            self.z * factor)

    def to2ab(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{Nvector.philam}.

           @return: A L{PhiLam2Tuple}C{(phi, lam)}.
        '''
        return n_xyz2philam(self.x, self.y, self.z)

    def to2ll(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{Nvector.latlon}.

           @return: A L{LatLon2Tuple}C{(lat, lon)}.
        '''
        return n_xyz2latlon(self.x, self.y, self.z)

    def to3xyz(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{xyz}.

           @return: A L{Vector3Tuple}C{(x, y, z)}.
        '''
        return self.xyz

    def toStr(self, prec=5, fmt='(%s)', sep=', '):  # PYCHOK expected
        '''Return a string representation of this vector.

           @kwarg prec: Optional number of decimal places (C{int}).
           @kwarg fmt: Optional, enclosing format to use (C{str}).
           @kwarg sep: Optional separator between components (C{str}).

           @return: Vector as "(x, y, z)" (C{str}).
        '''
        return fmt % (sep.join(strs(self.xyz, prec=prec)),)

    def unit(self, ll=None):
        '''Normalize this vector to unit length.

           @kwarg ll: Optional, original location (C{LatLon}).

           @return: Normalized vector (L{Vector3d}).
        '''
        if self._united is None:
            n = self.length
            if n > EPS and abs(n - 1) > EPS:
                u = self._xnamed(self.dividedBy(n))
                u._length = 1
            else:
                u = self.copy()
            u._fromll = ll or self._fromll
            self._united = u._united = u
        return self._united

    @property_RO
    def x(self):
        '''Get the X component (C{float}).
        '''
        return self._x

    @property_RO
    def xyz(self):
        '''Get the X, Y and Z components (L{Vector3Tuple}C{(x, y, z)}).
        '''
        if self._xyz is None:
            self._xyz = Vector3Tuple(self.x, self.y, self.z)
        return self._xnamed(self._xyz)

    @property_RO
    def y(self):
        '''Get the Y component (C{float}).
        '''
        return self._y

    @property_RO
    def z(self):
        '''Get the Z component (C{float}).
        '''
        return self._z


def sumOf(vectors, Vector=Vector3d, **Vector_kwds):
    '''Compute the vectorial sum of several vectors.

       @arg vectors: Vectors to be added (L{Vector3d}[]).
       @kwarg Vector: Optional class for the vectorial sum (L{Vector3d}).
       @kwarg Vector_kwds: Optional B{C{Vector}} keyword arguments,
                           ignored if B{C{Vector=None}}.

       @return: Vectorial sum (B{C{Vector}}).

       @raise VectorError: No B{C{vectors}}.
    '''
    n, vectors = len2(vectors)
    if n < 1:
        raise VectorError('no vectors: %r' & (n,))

    r = Vector3Tuple(fsum(v.x for v in vectors),
                     fsum(v.y for v in vectors),
                     fsum(v.z for v in vectors))
    if Vector is not None:
        r = Vector(r.x, r.y, r.z, **Vector_kwds)  # PYCHOK x, y, z
    return r

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
