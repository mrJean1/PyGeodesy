
# -*- coding: utf-8 -*-

# Python implementation of geodetic (lat-/longitude) functions.
# Transcribed from JavaScript originals by (C) Chris Veness 2011-2016
# and published under the same MIT Licence**, see
# <http://www.movable-type.co.uk/scripts/latlong.html>.

from bases import _LatLonHeightBase
from datum import R_M, Datum, Datums
from dms import parse3llh
from utils import EPS, EPS2, PI, PI2, \
                  degrees90, degrees180, degrees360, \
                  radians, radiansPI, tanPI_2_2
from math import acos, atan2, cos, hypot, log, sin

# all public contants, classes and functions
__all__ = ()  # classes
__version__ = '16.11.11'


class _LatLonSphericalBase(_LatLonHeightBase):
    '''Private class containing methods shared between the
       sphericalNvector and -Trig flavors.  Otherwise, such
       methods would have to be duplicated in both flavors.
    '''

    _datum = Datums.Sphere  # XXX TBD

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
           @throws {ValueError} If datum.ellipsoid is not spherical.
        '''
        if not isinstance(datum, Datum):
            raise TypeError('%r not a %s: %r' % ('datum', Datum.__name__, datum))
        E = datum.ellipsoid
        if not E.a == E.b == E.R:
            raise ValueError('%r not %s: %r' % ('datum', 'spherical', E))
        self._update(datum != self._datum)
        self._datum = datum

    def finalBearingTo(self, other):
        '''Return the final bearing (reverse azimuth) from this
           to an other point, in compass degrees from North.

           @param {LatLon} other - The other point.

           @returns {degrees360} Final bearing in degrees from North.

           @example
           p = LatLon(52.205, 0.119)
           q = LatLon(48.857, 2.351)
           b = p.finalBearingTo(q)  # 157.9
        '''
        self.others(other)

        # final bearing is reverse
        b = other.bearingTo(self) + 180
        if b > 360:
            b -= 360
        return b  # 0..360

    def maxLat(self, bearing):
        '''Returns maximum latitude reached when travelling
           on a great circle on given bearing from this
           point (based on 'Clairaut's formula').

           The maximum latitude is independent of longitude,
           it is the same for all points on a given latitude.

           Negate the result for the minimum latitude (in
           the Southern hemisphere).

           @param {degrees} bearing - Initial bearing in degrees.

           @returns {degrees90} Maximum latitude in degrees.
        '''
        m = acos(abs(sin(radians(bearing)) * cos(radians(self.lat))))
        return degrees90(m)

    maxLatitude = maxLat  # XXX original name

    def minLat(self, bearing):
        '''Returns minimum latitude reached when travelling
           on a great circle on given bearing from this
           point.  See method maxLat for more details.

           @param {degrees} bearing - Initial bearing in degrees.

           @returns {degrees90} Minimum latitude in degrees.
        '''
        return -self.maxLat(bearing)

    minLatitude = minLat  # XXX original name?

    def parse(self, strll, height=0, sep=','):
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
           @param {meter} [height=0] - Default height in meter.
           @param {string} [sep=','] - Separator.

           @returns {LatLon} Point for the location.

           @throws {ValueError} Invalid strll.
        '''
        return self.Top(*parse3llh(strll, height=height, sep=sep))

    def _rhumb3(self, other):
        # rhumb_ helper function
        self.others(other)

        a1 = radians(self.lat)
        a2 = radians(other.lat)
        # if |db| > 180 take shorter rhumb
        # line across the anti-meridian
        db = radiansPI(other.lon - self.lon)
        ds = log(tanPI_2_2(a2) / tanPI_2_2(a1))
        return (a2 - a1), db, ds

    def rhumbBearingTo(self, other):
        '''Return the initial bearing (forward azimuth) from this
           to an other point along a rhumb (loxodrome) line, in
           compass degrees from North.

           @param {LatLon} other - The other LatLon point.

           @returns {degrees360} Initial bearing in degrees from North.

           @example
           p = LatLon(51.127, 1.338)
           q = LatLon(50.964, 1.853)
           b = p.rhumbBearingTo(q)  # 116.7
        '''
        _, db, ds = self._rhumb3(other)

        return degrees360(atan2(db, ds))

    def rhumbDistanceTo(self, other, radius=R_M):
        '''Returns distance from this to an other point along
           a rhumb (loxodrome) line.

           @param {LatLon} other - the other LatLon point.
           @param {number} [radius=R_M] - Mean radius of earth,
                                          defaults to meter.

           @returns {number} Distance between this and the other
                             point (in the same units as radius).

           @example
           p = LatLon(51.127, 1.338)
           q = LatLon(50.964, 1.853)
           d = p.rhumbDistanceTo(q)  # 403100
        '''
        # see http://williams.best.vwh.net/avform.htm#Rhumb
        da, db, ds = self._rhumb3(other)

        # on Mercator projection, longitude distances shrink
        # by latitude; the 'stretch factor' q becomes ill-
        # conditioned along E-W line (0/0); use an empirical
        # tolerance to avoid it
        if abs(da) < EPS2:
            q = cos(radians(self.lat))
        else:
            q = da / ds

        return float(radius) * hypot(da, q * db)

    def rhumbMidpointTo(self, other):
        '''Return the loxodromic midpoint between this and
           an other point.

           @param {LatLon} other - the other LatLon point.

           @returns {LatLon} Midpoint between both points.

           @example
           p = LatLon(51.127, 1.338)
           q = LatLon(50.964, 1.853)
           m = p.rhumb_midpointTo(q)
           m.toStr()  # '51.0455°N, 001.5957°E'
        '''
        self.others(other)

        # see http://mathforum.org/library/drmath/view/51822.html
        a1, b1 = self.toradians()
        a2, b2 = other.toradians()
        if abs(b2 - b1) > PI:
            b1 += PI2  # crossing anti-meridian

        a3 = (a1 + a2) * 0.5
        b3 = (b1 + b2) * 0.5

        f1 = tanPI_2_2(a1)
        if abs(f1) > EPS:
            f2 = tanPI_2_2(a2)
            f = f2 / f1
            if abs(f) > EPS:
                f = log(f)
                if abs(f) > EPS:
                    f3 = tanPI_2_2(a3)
                    b3 = (b1 * log(f2) -
                          b2 * log(f1) + (b2 - b1) * log(f3)) / f

        return self.Top(degrees90(a3), degrees180(b3), height=self._alter(other))

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
