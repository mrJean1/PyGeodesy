
# -*- coding: utf-8 -*-

# Python implementation of vector-based geodetic (lat-/longitude) functions
# by (C) Chris Veness 2011-2015 published under the same MIT Licence**,
# see <http://www.movable-type.co.uk/scripts/latlong-vectors.html>
#
# These functions work with
# a) geodesic (polar) lat-/longitude points on the earth's surface
# b) 3D vectors used as n-vectors representing points on the earth's
#    surface or vectors normal to the plane of a great circle

from bases import _VectorBase
from utils import EPS2, \
                  degrees90, degrees180, fdot, fsum, \
                  hypot3, isscalar, fStr, len2
from math import atan2, cos, hypot, sin

# all public contants, classes and functions
__all__ = ('Vector3d',  # classes
           'sumOf')  # functions
__version__ = '16.11.30'


class Vector3d(_VectorBase):
    '''Generic 3-d vector manipulation.

       In a geodesy context, these may be used to represent:
       - n-vector representing a normal to point on earth's surface
       - earth-centered, earth-fixed vector (= n-vector for spherical model)
       - great circle normal to vector
       - motion vector on earth's surface
       - etc.
    '''

    _length = None  # cache length
    _united = None  # cache normalised, unit

    _x = 0
    _y = 0
    _z = 0

    def __init__(self, x, y, z):
        '''Create a 3-d vector.

           The vector may be normalised, or use x/y/z values for
           height relative to the surface of the sphere or ellipsoid,
           distance from earth centre, etc.

           @param {number} x - X component of vector.
           @param {number} y - Y component of vector.
           @param {number} z - Z component of vector.
        '''
        self._x = x
        self._y = y
        self._z = z

    def __add__(self, other):
        return self.plus(other)
    __radd__ = __add__

    def __abs__(self):
        return self.length()

    def __cmp__(self, other):
        return cmp(self.length(), other.length())

    def __div__(self, scalar):
        return self.dividedBy(scalar)

    def __mul__(self, scalar):
        return self.times(scalar)

    def __neg__(self):
        return self.negate()

    def __sub__(self, other):
        return self.minus(other)

    def __rsub__(self, other):
        return other.minus(self)

    def _update(self, updated):
        if updated:  # reset caches
            self._length = self._united = None

    def angleTo(self, other, vSign=None):
        '''Calculates the angle between this and an other vector.

           @param {Vector3d} other - The other vector.
           @param {Vector3d} [vSign=None] - If supplied (and out of
                             plane of this and other), angle is signed
                             positive if this->other is clockwise
                             looking along vSign or negative in opposite
                             direction (otherwise unsigned angle).

           @returns {radians} Angle between this and other vector.
        '''
        self.others(other)

        x = self.cross(other)
        s = x.length()
        # use vSign as reference to get sign of s
        if vSign is not None and x.dot(vSign) < 0:
            s = -s
        return atan2(s, self.dot(other))

    def copy(self):
        '''Return a copy of this vector.

           @returns {Vector3d} Copy of this vector.
        '''
        v = self.Top(self.x, self.y, self.z)
        v._length = self._length
        v._united = self._united
        return v

    def cross(self, other):
        '''Return cross product of this and an other vector.

           @param {Vector3d} other - Vector to be crossed with this vector.

           @returns {Vector3d} Cross product of this and the other.
        '''
        self.others(other)

        return self.Top(self.y * other.z - self.z * other.y,
                        self.z * other.x - self.x * other.z,
                        self.x * other.y - self.y * other.x)

    def dividedBy(self, factor):
        '''Return this vector divided by a scalar.

           @param {number} factor - Scale factor.

           @returns {Vector3d} New vector scaled.
        '''
        if not isscalar(factor):
            raise TypeError('%s not scalar: %r' % ('factor', factor))
        return self.times(1.0 / factor)

    def dot(self, other):
        '''Return the dot (scalar) product of this and an other vector.

           @param {Vector3d} other - Vector to be dotted with this.
           @returns {number} Dot product of this and the other.
        '''
        self.others(other)

        return fdot(self.to3tuple(), *other.to3tuple())

    def equals(self, other, units=False):
        '''Check if this vector is equal or equivalent to an other.

           @param {Vector3d} other - Vector to be compared against this.
           @param {bool} [units=False] - Use units=True to compare the
                                         normalized version of both
                                         vectors.

           @returns {bool} True if vectors are identical.

           @example
           v1 = Vector3d(52.205, 0.119)
           v2 = Vector3d(52.205, 0.119)
           e = v1.equals(v2)  # True
        '''
        self.others(other)

        if units:
            d = self.unit().minus(other.unit())
            return max(map(abs, d.to3tuple())) < EPS2

        else:
            return self.x == other.x and \
                   self.y == other.y and \
                   self.z == other.z

    def length(self):
        '''Return the length (magnitude or norm) of this vector.

           @returns {number} Magnitude of this vector.
        '''
        if self._length is None:
            self._length = hypot3(self.x, self.y, self.z)
        return self._length

    def minus(self, other):
        '''Return the vectorial difference of this and an other vector.

           @param {Vector3d} other - Vector3d to be subtracted from this.

           @returns {Vector3d} New vector, difference of this and the other.
        '''
        self.others(other)

        return self.Top(self.x - other.x,
                        self.y - other.y,
                        self.z - other.z)

    def negate(self):
        '''Return the vector in opposite direction of this one.

           @returns {Vector3d} New vector.
        '''
        return self.Top(-self.x, -self.y, -self.z)

    def others(self, other, name='other'):
        '''Refine class comparison.
        '''
        try:
            _VectorBase.others(self, other, name=name)
        except TypeError:
            if not isinstance(other, Vector3d):
                raise

    def parse(self, str3d):
        '''Parse a string representing x, y and z and return a
           Vector3d.

           The x, y and z must be separated by a comma.

           @param {string} str3d - X, y and z string.

           @returns {Vector3d} Point for the location.

           @throws {ValueError} Invalid str3d.
        '''
        try:
            v = [float(v.strip()) for v in str3d.split(',')]
            if len(v) != 3:
                raise ValueError
        except ValueError:
            raise ValueError('parsing %r' % (str3d,))

        return self.Top(*v)

    def plus(self, other):
        '''Return the vectorial sum of this and an other vector.

           @param {Vector3d} other - Vector to be added to this.

           @returns {Vector3d} New vector, sum of this and the other.
        '''
        self.others(other)

        return self.Top(self.x + other.x,
                        self.y + other.y,
                        self.z + other.z)

    def rotate(self, axis, theta):
        '''Rotates this vector by a specified angle around an axis.

           @param {Vector3d} axis - The axis being rotated around.
           @param {number} theta - The angle of rotation (in radians).

           @returns {Vector3d} The rotated vector.

           http://en.wikipedia.org/wiki/
                  Rotation_matrix#Rotation_matrix_from_axis_and_angle
                  Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix
                  Rodrigues'_rotation_formula...
        '''
        self.others(axis, name='axis')

        c = cos(theta)
        a = axis.unit()  # axis being rotated around
        b = a.times(1 - c)
        s = a.times(sin(theta))

        p = self.unit().to3tuple()  # point being rotated
        # multiply p by a quaternion-derived rotation matrix
        return self.Top(fdot(p, a.x * b.x + c,   a.x * b.y - s.z, a.x * b.z + s.y),
                        fdot(p, a.y * b.x + s.z, a.y * b.y + c,   a.y * b.z - s.x),
                        fdot(p, a.z * b.x - s.y, a.z * b.y + s.x, a.z * b.z + c))

    sum = plus  # alternate name

    def times(self, factor):
        '''Return this vector multiplied by a scalar.

           @param {number} factor - Scale factor.

           @returns {Vector3d} New vector scaled.
        '''
        if not isscalar(factor):
            raise TypeError('%s not scalar: %r' % ('factor', factor))
        return self.Top(self.x * factor,
                        self.y * factor,
                        self.z * factor)

    def to2latlon(self):
        '''Convert this vector to (geodetic) lat- and longitude.

           @returns {(degrees90, degrees180)} 2-Tuple with (lat, lon).

           @example
           v = Vector3d(0.500, 0.500, 0.707)
           a, b = v.to2latlon()  # 45.0, 45.0
        '''
        a = atan2(self.z, hypot(self.x, self.y))
        b = atan2(self.y, self.x)
        return degrees90(a), degrees180(b)

    def to3tuple(self):
        '''Return this vector as a 3-tuple.

           @returns {(x, y, z)} 3-Tuple of x, y and z.
        '''
        return self.x, self.y, self.z

    def toStr(self, prec=5, fmt='(%s)', sep=', '):  # PYCHOK expected
        '''String representation of this vector.

           @param {number} [prec=6] - Number of decimal places.
           @param {string} [fmt='(%s)'] - Format to use.
           @param {string} [sep=', '] - Separator between NEDs.

           @returns {string} Vector represented as "(x, y, z)".
        '''
        return fmt % (fStr(self.to3tuple(), prec=prec, sep=sep),)

    def unit(self):
        '''Normalize this vector to unit length.

           @returns {Vector3d} Normalised vector.
        '''
        if self._united is None:
            n = self.length()
            if n > EPS2 and abs(n - 1) > EPS2:
                u = self.dividedBy(n)
                u._length = 1
            else:
                u = self.copy()
            self._united = u._united = u
        return self._united

    @property
    def x(self):
        '''The x component.
        '''
        return self._x

    @property
    def y(self):
        '''The y component.
        '''
        return self._y

    @property
    def z(self):
        '''The z component.
        '''
        return self._z


def sumOf(vectors):
    '''Return the vectorial sum of any number of vectors.

       @param {Vector3d[]} vectors - Array of Vector3d to be added.

       @returns {Vector3d} New vector, sum of the vectors.

       @throws {ValueError} No vectors.
    '''
    n, vectors = len2(vectors)
    if n < 1:
        raise ValueError('no vectors: %r' & (n,))
    return Vector3d(fsum(v.x for v in vectors),
                    fsum(v.y for v in vectors),
                    fsum(v.z for v in vectors))

# **) MIT License
#
# Copyright (c) 2016-2017 -- mrJean1@Gmail.com
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
