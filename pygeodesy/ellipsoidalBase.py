
# -*- coding: utf-8 -*-

u'''(INTERNAL) Ellipsoidal base classes.

Pure Python implementation of geodesy tools for ellipsoidal earth models,
transcribed in part from JavaScript originals by I{(C) Chris Veness 2005-2016}
and published under the same MIT Licence**, see for example
U{http://www.movable-type.co.uk/scripts/geodesy/docs/latlon-ellipsoidal.js.html}.

@newfield example: Example, Examples
'''

from bases import LatLonHeightBase
from datum import Datum, Datums
from dms import parse3llh
from utils import EPS, EPS1, degrees90, degrees180, hypot, hypot1
from vector3d import Vector3d

from math import atan2, copysign, cos, sin, sqrt

# XXX the following classes are listed only to get
# Epydoc to include class and method documentation
__all__ = ('CartesianBase', 'LatLonEllipsoidalBase')
__version__ = '17.09.22'


class CartesianBase(Vector3d):
    '''(INTERNAL) Base class for ellipsoidal Cartesians.
    '''

    def _applyHelmert(self, transform, inverse=False):
        '''(INTERNAL) Return a new (geocentric) Cartesian point
           by applying a Helmert transform to this point.

           @param transform: Transform to apply (L{Transform}).
           @keyword inverse: Optionally, apply the inverse of
                             Helmert transform (bool).

           @return: The transformed point (L{Cartesian}).
        '''
        x, y, z = self.to3xyz()
        return self.classof(*transform.transform(x, y, z, inverse))

    def to3llh(self, datum=Datums.WGS84):
        '''Convert this (geocentric) Cartesian (x/y/z) point to
           (ellipsoidal) geodetic lat-, longitude and height on
           the given datum.

           Uses Bowring’s (1985) formulation for μm precision in
           concise form: 'The accuracy of geodetic latitude and
           height equations', B. R. Bowring, Survey Review, Vol
           28, 218, Oct 1985.

           See also 'An Efficient Algorithm for Geocentric to
           Geodetic Coordinate Conversion', Ralph M. Toms,
           Lawrence Livermore National Laboratory, Sept 1995,
           U{http://www.osti.gov/scitech/biblio/110235} and 'An
           Improved Algorithm for Geocentric to Geodetic Coordinate
           Conversion', Ralph Toms, Apr 1996, U{http://www.osti.gov/
           scitech/servlets/purl/231228}.

           @keyword datum: Optional datum to use (L{Datum}).

           @return: 3-Tuple (lat, lon, heigth) in (degrees90,
                    degrees180, meter).
        '''
        E = datum.ellipsoid
        x, y, z = self.to3xyz()

        p = hypot(x, y)  # distance from minor axis
        r = hypot(p, z)  # polar radius

        if min(p, r) > EPS:
            # parametric latitude (Bowring eqn 17, replaced)
            t = (E.b * z) / (E.a * p) * (1 + E.e22 * E.b / r)
            c = 1 / hypot1(t)
            s = t * c

            # geodetic latitude (Bowring eqn 18)
            a = atan2(z + E.e22 * E.b * s**3,
                      p - E.e2  * E.a * c**3)
            b = atan2(y, x)  # ... and longitude

            # height above ellipsoid (Bowring eqn 7)
            ca, sa = cos(a), sin(a)
#           r = E.a / E.e2s2(sa)  # length of normal terminated by minor axis
#           h = p * ca + z * sa - (E.a * E.a / r)
            h = p * ca + z * sa - (E.a * E.e2s2(sa))

            a, b = degrees90(a), degrees180(b)

        # see <http://GIS.StackExchange.com/questions/28446/>
        elif p > EPS:  # latitude arbitrarily zero
            a, b, h = 0.0, degrees180(atan2(y, x)), p - E.a
        else:  # polar latitude, longitude arbitrarily zero
            a, b, h = copysign(90.0, z), 0.0, abs(z) - E.b

        return a, b, h

    def toStr(self, prec=3, fmt='[%s]', sep=', '):  # PYCHOK expected
        '''Return the string representation of this cartesian.

           @keyword prec: Optional number of decimals, unstripped (int).
           @keyword fmt: Optional enclosing backets format (string).
           @keyword sep: Optional separator to join (string).

           @return: Cartesian represented as "[x, y, z]" (string).
        '''
        return Vector3d.toStr(self, prec=prec, fmt=fmt, sep=sep)


class LatLonEllipsoidalBase(LatLonHeightBase):
    '''(INTERNAL) Base class for ellipsoidal LatLons.
    '''
    _datum = Datums.WGS84  #: (INTERNAL) Datum (L{Datum}).
    _osgr  = None  #: (INTERNAL) cache toOsgr ({Osgr}).
    _utm   = None  #: (INTERNAL) cache toUtm (L{Utm}).

    def __init__(self, lat, lon, height=0, datum=None):
        '''Create an (ellipsoidal) I{LatLon} point frome the given
           lat-, longitude and height (elevation, altitude) on the
           given datum.

           @param lat: Latitude (degrees or DMS string with N or S suffix).
           @param lon: Longitude (degrees or DMS string with E or W suffix).
           @keyword height: Optional elevation (meter or the same units
                            as the datum's half-axes).
           @keyword datum: Optional datum to use (L{Datum}).

           @example:

           >>> p = LatLon(51.4778, -0.0016)  # height=0, datum=Datums.WGS84
        '''
        LatLonHeightBase.__init__(self, lat, lon, height=height)
        if datum:  # check datum
            self.datum = datum

    def _update(self, updated):
        if updated:  # reset caches
            self._osgr = self._utm = None
            LatLonHeightBase._update(self, updated)

    def convertDatum(self, datum):
        '''Convert this I{LatLon} point to a new coordinate system.

           @param datum: Datum to convert to (L{Datum}).

           @return: The converted point (L{LatLonEllipsoidalBase}).

           @example:

           >>> pWGS84 = LatLon(51.4778, -0.0016)  # default Datums.WGS84
           >>> pOSGB  = pWGS84.convertDatum(Datums.OSGB36)  # 51.477284°N, 000.00002°E
        '''
        if self.datum == datum:
            return self.copy()

        elif self.datum == Datums.WGS84:
            # converting from WGS 84
            ll, t, i = self, datum.transform, False

        elif datum == Datums.WGS84:
            # converting to WGS84, use inverse transform
            ll, t, i = self, self.datum.transform, True

        else:  # neither self.datum nor datum is WGS84, convert to WGS84 first
            ll, t, i = self.convertDatum(Datums.WGS84), datum.transform, False

        return ll.toCartesian()._applyHelmert(t, i).toLatLon(datum=datum)

    def copy(self):
        '''Copy this I{LatLon} point.

           @return: Copy of this point (L{LatLonEllipsoidalBase}).
        '''
        p = LatLonHeightBase.copy(self)
        assert hasattr(p, 'datum')
        p.datum = self.datum
        return p

    @property
    def datum(self):
        '''Get this point's datum (L{Datum}).
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set this point's datum without conversion.

           @param datum: New datum (L{Datum}).

           @raise TypeError: The datum is not a L{Datum}.

           @raise ValueError: The datum is not ellipsoidal.
        '''
        if not isinstance(datum, Datum):
            raise TypeError('%r not a %s: %r' % ('datum', Datum.__name__, datum))
        if not datum.isEllipsoidal:
            raise ValueError('%r not %s: %r' % ('datum', 'ellipsoidal', datum))
        self._update(datum != self._datum)
        self._datum = datum

    def ellipsoid(self, datum=Datums.WGS84):
        '''Return the ellipsoid of this point's datum or the given datum.

           @keyword datum: Optional datum (L{Datum}).

           @return: The ellipsoid (L{Ellipsoid}).
        '''
        return getattr(self, 'datum', datum).ellipsoid

    def ellipsoids(self, other):
        '''Check the type and ellipsoid of this points' datum and an other datum.

           @param other: The other datum (L{Datum}).

           @return: This datum's ellipsoid (L{Ellipsoid}).

           @raise TypeError: The other point is not LatLon.

           @raise ValueError: If datum ellipsoids are incompatible.
        '''
        self.others(other)

        E = self.ellipsoid()
        try:  # other may be Sphere, etc.
            e = other.ellipsoid()
        except AttributeError:
            try:  # no ellipsoid method, try datum
                e = other.datum.ellipsoid
            except AttributeError:
                e = E  # no datum, XXX assume equivalent?
        if e != E:
            c = E.__class__.__name__
            raise ValueError('%s %s mistmatch: %ss.%s vs %ss.%s' %
                             ('other', c, c, e.name, c, E.name))
        return E

    @property
    def isEllipsoidal(self):
        '''Check whether this I{LatLon} point is ellipsoidal (bool).
        '''
        return self.datum.isEllipsoidal

    @property
    def isSpherical(self):
        '''Check whether this I{LatLon} point is spherical (bool).
        '''
        return self.datum.isSpherical

    def parse(self, strll, height=0, datum=None, sep=','):
        '''Parse a string representing this I{LatLon} point.

           The lat- and longitude must be separated by a sep[arator]
           character.  If height is present it must follow and be
           separated by another sep[arator].  Lat- and longitude
           may be swapped, provided at least one ends with the
           proper compass direction.

           For more details, see functions L{parse3llh} and L{parseDMS}
           in sub-module L{dms}.

           @param strll: Lat, lon [, height] (string).
           @keyword height: Optional, default height (meter or None).
           @keyword datum: Optional, default datum (L{Datum}).
           @keyword sep: Optional separator (string).

           @return: The point (L{LatLonEllipsoidalBase}).

           @raise ValueError: Invalid strll.
        '''
        a, b, h = parse3llh(strll, height=height, sep=sep)
        return self.classof(a, b, height=h, datum=datum or self.datum)

    def to3xyz(self):  # overloads _LatLonHeightBase.to3xyz
        '''Convert this (ellipsoidal) geodetic I{LatLon} point to
           (geocentric) cartesian x/y/z components.

           @return: 3-Tuple (x, y, z) in (meter).
        '''
        a, b = self.to2ab()
        sa = sin(a)

        E = self.ellipsoid()
        # radius of curvature in prime vertical
        t = E.e2 * sa * sa
        if t < EPS:
            r = E.a
        elif t < EPS1:
            r = E.a / sqrt(1 - t)
        else:
            r = 0

        h = self.height
        t = (h + r) * cos(a)
        return (t * cos(b),
                t * sin(b),
               (h + r * E.e12) * sa)

    def toOsgr(self):
        '''Convert this I{LatLon} point to an OSGR coordinate.

           See function L{toOsgr} in module L{osgr} for more details.

           @return: The OSGR coordinate (L{Osgr}).
        '''
        if self._osgr is None:
            from osgr import toOsgr  # PYCHOK recursive import
            self._osgr = toOsgr(self, datum=self.datum)
            self._osgr._latlon = self
        return self._osgr

    def toUtm(self):
        '''Convert this I{LatLon} point to a UTM coordinate.

           See function L{toUtm} in module L{utm} for more details.

           @return: The UTM coordinate (L{Utm}).
        '''
        if self._utm is None:
            from utm import toUtm  # PYCHOK recursive import
            self._utm = toUtm(self, datum=self.datum)
            self._utm._latlon = self
        return self._utm

# **) MIT License
#
# Copyright (C) 2016-2018 -- mrJean1 at Gmail dot com
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
