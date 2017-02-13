
# -*- coding: utf-8 -*-

'''Generic 3-D vector base class L{Vector3d} and function L{sumOf}.

Python implementation of vector-based functions by I{(C) Chris Veness
2011-2015} published under the same MIT Licence**, see
U{http://www.movable-type.co.uk/scripts/latlong-vectors.html}.

@newfield example: Example, Examples
'''

from bases import VectorBase
from utils import EPS, degrees90, degrees180, fdot, fStr, fsum, \
                  hypot3, isscalar, len2

from math import atan2, cos, hypot, sin

# all public contants, classes and functions
__all__ = ('Vector3d',  # classes
           'sumOf')  # functions
__version__ = '17.02.11'

try:
    _cmp = cmp
except NameError:  # Python 3+
    def _cmp(a, b):
        if a < b:
            return -1
        elif a > b:
            return +1
        else:
            return 0


class Vector3d(VectorBase):
    '''Generic 3-D vector manipulation.

       In a geodesy context, these may be used to represent:
        - n-vector representing a normal to point on earth's surface
        - earth-centered, earth-fixed vector (= n-vector for spherical model)
        - great circle normal to vector
        - motion vector on earth's surface
        - etc.
    '''

    _length = None  #: (INTERNAL) cached length.
    _united = None  #: (INTERNAL) cached norm, unit.

    _x = 0  #: (INTERNAL) X component.
    _y = 0  #: (INTERNAL) Y component.
    _z = 0  #: (INTERNAL) Z component.

    def __init__(self, x, y, z):
        '''New 3-D vector.

           The vector may be normalised, or use x/y/z values for
           height relative to the surface of the sphere or ellipsoid,
           distance from earth centre, etc.

           @param x: X component of vector.
           @param y: Y component of vector.
           @param z: Z component of vector.
        '''
        self._x = x
        self._y = y
        self._z = z

    def __add__(self, other):
        '''This plus an other vector (L{Vector3d}).
        '''
        return self.plus(other)
    __radd__ = __add__

    def __abs__(self):
        '''Norm of this vector (scalar).
        '''
        return self.length()

    def __cmp__(self, other):  # Python 2-
        '''Compare this and an other vector?
           @return: -1, 0 or +1 (int).
        '''
        return _cmp(self.length(), other.length())

    def __div__(self, scalar):
        '''Divide this vector by a scalar.
           @return: Quotient (L{Vector3d}).
        '''
        return self.dividedBy(scalar)

    def __gt__(self, other):
        '''Is this vector longer than other vector?
           @return: True if so (bool).
        '''
        return self.length() > other.length()

    def __mul__(self, scalar):
        '''Multiply this vector by a scalar
           @return: Product (L{Vector3d}).
        '''
        return self.times(scalar)

    def __lt__(self, other):  # Python 3+
        '''Is this vector shorter than other vector?
           @return: True if so (bool).
        '''
        return self.length() < other.length()

    def __neg__(self):
        '''Negate this vector.
           @return: Opposite (L{Vector3d})
        '''
        return self.negate()

    def __sub__(self, other):
        '''This minus an other vector.
           @return: Difference (L{Vector3d}).
        '''
        return self.minus(other)

    def __rsub__(self, other):
        '''An other minus this vector.
           @return: Difference (L{Vector3d}).
        '''
        return other.minus(self)

    def _update(self, updated):
        '''(INTERNAL) Clear caches.
        '''
        if updated:  # reset caches
            self._length = self._united = None

    def angleTo(self, other, vSign=None):
        '''Angle between this and an other vector.

           @param other: The other vector (L{Vector3d}).
           @keyword vSign: Vector, if supplied (and out of the plane of
                           this and the other), angle is signed positive
                           if this->other is clockwise looking along vSign
                           or negative in opposite direction, otherwise
                           angle is unsigned.

           @return: Angle (radians).
        '''
        x = self.cross(other)
        s = x.length()
        # use vSign as reference to get sign of s
        if vSign is not None and x.dot(vSign) < 0:
            s = -s
        return atan2(s, self.dot(other))

    def copy(self):
        '''Copy this vector.

           @return: New, vector copy (Vector3d).
        '''
        v = self.topsub(self.x, self.y, self.z)
        v._length = self._length
        v._united = self._united
        return v

    def cross(self, other):
        '''Cross product of this and an other vector.

           @param other: The other vector (L{Vector3d}).

           @return: Cross product (L{Vector3d}).
        '''
        self.others(other)

        return self.topsub(self.y * other.z - self.z * other.y,
                           self.z * other.x - self.x * other.z,
                           self.x * other.y - self.y * other.x)

    def dividedBy(self, factor):
        '''Divide this vector by a scalar.

           @param factor: The divisor (scalar).

           @return: New, scaled vector (L{Vector3d}).

           @raise TypeError: If L{factor} not scalar'
        '''
        if not isscalar(factor):
            raise TypeError('%s not scalar: %r' % ('factor', factor))
        return self.times(1.0 / factor)

    def dot(self, other):
        '''Dot (scalar) product of this and an other vector.

           @param other: The other vector (L{Vector3d}).

           @return: Dot product (float).
        '''
        self.others(other)

        return fdot(self.to3xyz(), *other.to3xyz())

    def equals(self, other, units=False):
        '''Check if this and an other vector are equal or equivalent.

           @param other: The other vector (L{Vector3d}).
           @keyword units: Use units=True to compare the normalized,
                           unit version of both vectors.

           @return: True if vectors are identical.

           @example:

           >>> v1 = Vector3d(52.205, 0.119)
           >>> v2 = Vector3d(52.205, 0.119)
           >>> e = v1.equals(v2)  # True
        '''
        self.others(other)

        if units:
            d = self.unit().minus(other.unit())
        else:
            d = self.minus(other)
        return max(map(abs, d.to3xyz())) < EPS

    def length(self):
        '''Length (aka norm or magnitude) of this vector.

           @return: Norm (float).
        '''
        if self._length is None:
            self._length = hypot3(self.x, self.y, self.z)
        return self._length

    def minus(self, other):
        '''Subtract an other vector from this vector.

           @param other: The other vector (L{Vector3d}).

           @return: New vector difference (L{Vector3d}).
        '''
        self.others(other)

        return self.topsub(self.x - other.x,
                           self.y - other.y,
                           self.z - other.z)

    def negate(self):
        '''This vector in opposite direction.

           @return: New, opposite vector (L{Vector3d}).
        '''
        return self.topsub(-self.x, -self.y, -self.z)

    def others(self, other, name='other'):
        '''Refined class comparison.

           @param other: The other vector (L{Vector3d}).
           @keyword name: Other's name (string).

           @raise TypeError: Incompatible type(L{other}).
        '''
        try:
            VectorBase.others(self, other, name=name)
        except TypeError:
            if not isinstance(other, Vector3d):
                raise

    def parse(self, str3d):
        '''Parse an "x, y, z" string representing a L{Vector3d}.

           The x, y and z must be separated by a comma.

           @param str3d: X, y and z string.

           @return: New vector (L{Vector3d}).

           @raise ValueError: Invalid L{str3d}.
        '''
        try:
            v = [float(v.strip()) for v in str3d.split(',')]
            if len(v) != 3:
                raise ValueError
        except ValueError:
            raise ValueError('%s invalid: %r' % ('str3d', str3d))

        return self.topsub(*v)

    def plus(self, other):
        '''Add this and an other vector.

           @param other: The other vector (L{Vector3d}).

           @return: New vector sum (L{Vector3d}).
        '''
        self.others(other)

        return self.topsub(self.x + other.x,
                         self.y + other.y,
                         self.z + other.z)

    sum = plus  # alternate name

    def rotate(self, axis, theta):
        '''Rotate this vector by a specified angle around an axis.

           See U{Rotation matrix from axis and angle<http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix/Rodrigues'_rotation_formula...>}

           @param axis: The axis being rotated around (L{Vector3d}).
           @param theta: The angle of rotation (radians).

           @return: New, rotated vector (L{Vector3d}).
        '''
        self.others(axis, name='axis')

        c = cos(theta)
        a = axis.unit()  # axis being rotated around
        b = a.times(1 - c)
        s = a.times(sin(theta))

        p = self.unit().to3xyz()  # point being rotated
        # multiply p by a quaternion-derived rotation matrix
        return self.topsub(fdot(p, a.x * b.x + c,   a.x * b.y - s.z, a.x * b.z + s.y),
                           fdot(p, a.y * b.x + s.z, a.y * b.y + c,   a.y * b.z - s.x),
                           fdot(p, a.z * b.x - s.y, a.z * b.y + s.x, a.z * b.z + c))

    rotateAround = rotate  # alternate name

    def times(self, factor):
        '''Multiply this vector by a scalar.

           @param factor: Scale factor (scalar).

           @return: New, scaled vector (L{Vector3d}).
        '''
        if not isscalar(factor):
            raise TypeError('%s not scalar: %r' % ('factor', factor))
        return self.topsub(self.x * factor,
                           self.y * factor,
                           self.z * factor)

    def to2ll(self):
        '''Convert this vector to (geodetic) lat- and longitude.

           @return: 2-Tuple (lat, lon) in (degrees90, degrees180).

           @example:

           >>> v = Vector3d(0.500, 0.500, 0.707)
           >>> a, b = v.to2ll()  # 45.0, 45.0
        '''
        a = atan2(self.z, hypot(self.x, self.y))
        b = atan2(self.y, self.x)
        return degrees90(a), degrees180(b)

    def to3xyz(self):
        '''Return this vector as a 3-tuple.

           @return: 3-Tuple (x, y, z) as (scalar, ...).
        '''
        return self.x, self.y, self.z

    def toStr(self, prec=5, fmt='(%s)', sep=', '):  # PYCHOK expected
        '''String representation of this vector.

           @keyword prec: Number of decimal places (int).
           @keyword fmt: Enclosing format to use (string).
           @keyword sep: Separator between components (string).

           @return: Vector as "(x, y, z)" (string).
        '''
        return fmt % (fStr(self.to3xyz(), prec=prec, sep=sep),)

    def unit(self):
        '''Normalize this vector to unit length.

           @return: Normalised, unit vector (L{Vector3d}).
        '''
        if self._united is None:
            n = self.length()
            if n > EPS and abs(n - 1) > EPS:
                u = self.dividedBy(n)
                u._length = 1
            else:
                u = self.copy()
            self._united = u._united = u
        return self._united

    @property
    def x(self):
        '''Get the X component (scalar).
        '''
        return self._x

    @property
    def y(self):
        '''Get the Y component (scalar).
        '''
        return self._y

    @property
    def z(self):
        '''Get the Z component (scalar).
        '''
        return self._z


def sumOf(vectors, Vector=Vector3d, **kwds):
    '''Vectorially add a number of vectors.

       @param vectors: Vectors to be added (L{Vector3d}[]).
       @keyword Vector: Vector class to instantiate (L{Vector3d}).
       @keyword kwds: Optional, additional L{Vector} keyword argments.

       @return: Vectorial sum (L{Vector}).

       @raise ValueError: No L{vectors}.
    '''
    n, vectors = len2(vectors)
    if n < 1:
        raise ValueError('no vectors: %r' & (n,))
    return Vector(fsum(v.x for v in vectors),
                  fsum(v.y for v in vectors),
                  fsum(v.z for v in vectors), **kwds)

# **) MIT License
#
# Copyright (C) 2016-2017 -- mrJean1@Gmail.com
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
