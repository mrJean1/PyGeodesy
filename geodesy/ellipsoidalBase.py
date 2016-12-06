
# -*- coding: utf-8 -*-

# Python implementation of geodesy tools for an ellipsoidal earth model.
# Transcribed from JavaScript originals by (C) Chris Veness 2005-2016
# and published under the same MIT Licence**, see for example
# <http://www.movable-type.co.uk/scripts/geodesy/docs/latlon-ellipsoidal.js.html>

from bases import _LatLonHeightBase
from datum import Datum, Datums
from dms import parse3llh
from utils import EPS2, degrees90, degrees180, hypot1
from vector3d import Vector3d
from math import atan2, copysign, cos, hypot, sin, sqrt

# all public constants, classes and functions
__all__ = ()  # none
__version__ = '16.12.02'


class _CartesianBase(Vector3d):
    '''Base class for ellipsoidal Cartesian.
    '''

    def _applyHelmert(self, transform, inverse=False):
        '''Return a new (geocentric) Cartesian point by
           applying a Helmert transform to this point.

           @private
           @param {Transform} transform - Transform to apply.
           @param (bool} inverse - Apply inverse Helmert transform.

           @returns {Cartesian} The transformed point.
        '''
        x, y, z = self.to3tuple()
        return self.Top(*transform.transform(x, y, z, inverse))

    def to3latlonheight(self, datum=Datums.WGS84):
        '''Converts this (geocentric) Cartesian (x/y/z) point to
           (ellipsoidal geodetic) lat-, longitude and height on
           the given datum.

           @param {Datum} [datum=Datums.WGS84] - Datum to use.

           @returns {(degrees90, degrees180, meter)} 3-Tuple
                                          (lat, lon, heigth).

           Uses Bowring’s (1985) formulation for μm precision in
           concise form: 'The accuracy of geodetic latitude and
           height equations', B. R. Bowring, Survey Review, Vol
           28, 218, Oct 1985.
        '''
        E = datum.ellipsoid
        x, y, z = self.to3tuple()

        p = hypot(x, y)  # distance from minor axis
        r = hypot(p, z)  # polar radius

        if min(p, r) > EPS2:
            # parametric latitude (Bowring eqn 17, replaced)
            t = (E.b * z) / (E.a * p) * (1 + E.e22 * E.b / r)
            s = t / hypot1(t)
            c = s / t

            # geodetic latitude (Bowring eqn 18)
            a = atan2(z + E.e22 * E.b * s * s * s,
                      p - E.e2  * E.a * c * c * c)
            b = atan2(y, x)  # ... and longitude

            # height above ellipsoid (Bowring eqn 7)
            ca, sa = cos(a), sin(a)
#           r = E.a / E.e2s2(sa)  # length of normal terminated by minor axis
#           h = p * ca + z * sa - (E.a * E.a / r)
            h = p * ca + z * sa - (E.a * E.e2s2(sa))

            a, b = degrees90(a), degrees180(b)

        # see <http://GIS.StackExchange.com/questions/28446/>
        elif p > EPS2:  # latitude arbitrarily zero
            a, b, h = 0.0, degrees180(atan2(y, x)), p - E.a
        else:  # polar latitude, longitude arbitrarily zero
            a, b, h = copysign(90.0, z), 0.0, abs(z) - E.b

        return a, b, h

    def toStr(self, prec=3, fmt='[%s]', sep=', '):  # PYCHOK expected
        '''String representation of this cartesion.

           @param {number} [prec=3] - Number of decimals, unstripped.
           @param {string} [fmt='[%s]'] - Enclosing backets format.
           @param {string} [sep=', '] - Separator to join.

           @returns {string} Cartesion represented as "[x, y, z]".
        '''
        return Vector3d.toStr(self, prec=prec, fmt=fmt, sep=sep)


class _LatLonHeightDatumBase(_LatLonHeightBase):
    '''Base class for ellipsoidal LatLon.
    '''
    _datum = Datums.WGS84
    _osgr  = None  # cache toOsgr
    _utm   = None  # cache toUtm

    def __init__(self, lat, lon, height=0, datum=None):
        '''Create an (ellipsoidal) LatLon point frome the given
           lat-, longitude and height (elevation, altitude) on
           a given datum.

           @constructor
           @param {degrees|strDMS} lat - Latitude as degrees or string.
           @param {degrees|strDMS} lon - Longitude as degrees or string.
           @param {number} [height=0] - Elevation in meter (must be
                          the same units as the datum's half-axes).
           @param {Datum} [datum=Datums.WGS84] - Datum to use.

           @example
           p = LatLon(51.4778, -0.0016)  # height=0, datum=Datums.WGS84
        '''
        _LatLonHeightBase.__init__(self, lat, lon, height=height)
        if datum:
            self.datum = datum

    def _update(self, updated):
        if updated:  # reset caches
            self._osgr = self._utm = None
#           _LatLonHeightBase._update(self, updated)

    def convertDatum(self, toDatum):
        '''Converts this LatLon instance to new coordinate system.

           @param {Datum} toDatum - Datum to convert to.

           @returns {LatLon} This point converted to toDatum.

           @example
           pWGS84 = LatLon(51.4778, -0.0016)  # default Datums.WGS84
           pOSGB  = pWGS84.convertDatum(Datums.OSGB36)  # 51.477284°N, 000.00002°E
        '''
        if self.datum == toDatum:
            return self.copy()

        elif self.datum == Datums.WGS84:
            # converting from WGS 84
            ll, t, i = self, toDatum.transform, False

        elif toDatum == Datums.WGS84:
            # converting to WGS84, use inverse transform
            ll, t, i = self, self.datum.transform, True

        else:  # neither self.datum nor toDatum is WGS84, convert to WGS84 first
            ll, t, i = self.convertDatum(Datums.WGS84), toDatum.transform, False

        return ll.toCartesian()._applyHelmert(t, i).toLatLon(datum=toDatum)

    toDatum = convertDatum  # alternate name

    def copy(self):
        '''Return a copy of this point.

           @returns {LatLon} Copy of this point.
        '''
        p = _LatLonHeightBase.copy(self)
        assert hasattr(p, 'datum')
        p.datum = self.datum
        return p

    @property
    def datum(self):
        '''Get this point's datum.
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set this point's datum without conversion.

           @param {Datum} datum - New datum.

           @throws {TypeError} If datum is not a Datum.
           @throws {ValueError} If datum.ellipsoid is not ellipsoidal.
        '''
        if not isinstance(datum, Datum):
            raise TypeError('%r not a %s: %r' % ('datum', Datum.__name__, datum))
        E = datum.ellipsoid
        if not E.a > E.R > E.b:
            raise ValueError('%r not %s: %r' % ('datum', 'ellipsoidal', E))
        self._update(datum != self._datum)
        self._datum = datum

    def ellipsoid(self, datum=Datums.WGS84):
        '''Return the ellipsoid of this or the given datum.
        '''
        return getattr(self, 'datum', datum).ellipsoid

    def ellipsoids(self, other):
        '''Check type and ellipsoid model match.
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
            raise ValueError('other %s mistmatch: %ss.%s vs %ss.%s' %
                             (c, c, e.name, c, E.name))
        return E

    def parse(self, strll, height=0, datum=None, sep=','):
        '''Parse a string representing lat-/longitude point and
           return a LatLon.

           The lat- and longitude must be separated by a sep[arator]
           character.  If height is present it must follow and be
           separated by another sep[arator].  Lat- and longitude
           may be swapped, provided at least one ends with the
           proper compass direction.

           For more details, see functions parse3llh and parseDMS
           in module dms.

           @param {string} strll - Latitude, longitude [, height].
           @param {meter|None} [height=0] - Default height in meter.
           @param {Datum} [datum=self.datum] - Default datum.
           @param {string} [sep=','] - Separator.

           @returns {LatLon} Point for the location.

           @throws {ValueError} Invalid strll.
        '''
        a, b, h = parse3llh(strll, height=height, sep=sep)
        return self.Top(a, b, height=h, datum=datum or self.datum)

    def to3xyz(self):  # overloads _LatLonHeightBase.to3xyz
        '''Converts this (ellipsoidal geodetic) LatLon point to
           (ellipsoidal geocentric) cartesian x/y/z components.

           @returns {(meter, meter, meter)} 3-Tuple (x, y, z).
        '''
        a, b = self.toradians()
        ca, sa = cos(a), sin(a)

        E = self.ellipsoid()
        # radius of curvature in prime vertical
        r = E.a / sqrt(1 - E.e2 * sa * sa)

        h = self.height
        return ((h + r) * ca * cos(b),
                (h + r) * ca * sin(b),
                (h + r * (1 - E.e2)) * sa)

    def toOsgr(self):
        '''Convert this lat-/longitude to an OSGR coordinate.

           See function toOsgr in module osgr for more details.

           @returns {Osgr} The OSGR coordinate.
        '''
        if self._osgr is None:
            from osgr import toOsgr  # PYCHOK recursive import
            self._osgr = toOsgr(self, datum=self.datum)
            self._osgr._latlon = self
        return self._osgr

    def toUtm(self):
        '''Convert this lat-/longitude to a UTM coordinate.

           See function toUtm in module utm for more details.

           @returns {Utm} The UTM coordinate.
        '''
        if self._utm is None:
            from utm import toUtm  # PYCHOK recursive import
            self._utm = toUtm(self, datum=self.datum)
            self._utm._latlon = self
        return self._utm

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
