
# -*- coding: utf-8 -*-

u'''(INTERNAL) Spherical base classes.

Pure Python implementation of geodetic (lat-/longitude) functions,
transcribed in part from JavaScript originals by I{(C) Chris Veness 2011-2016}
and published under the same MIT Licence**, see
U{Latitude/Longitude<http://www.movable-type.co.uk/scripts/latlong.html>}.

@newfield example: Example, Examples
'''

from bases import LatLonHeightBase
from datum import R_M, R_MA, Datum, Datums
from dms   import parse3llh
from fmath import EPS, favg, fsum_, hypot
from utils import PI, PI2, PI_2, degrees90, degrees180, degrees360, \
                  radians, tanPI_2_2, wrapPI

from math import acos, atan2, cos, log, sin

# XXX the following classes are listed only to get
# Epydoc to include class and method documentation
__all__ = ('LatLonSphericalBase',)
__version__ = '18.02.05'


class LatLonSphericalBase(LatLonHeightBase):
    '''(INTERNAL) Base class for spherical Latlons.
    '''
    _datum = Datums.Sphere  #: (INTERNAL) XXX TBD

    @property
    def datum(self):
        '''Get this point's datum (L{Datum}).
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set this point's datum without conversion.

           @param datum: New datum (L{Datum}).

           @raise TypeError: If datum is not a L{Datum}.

           @raise ValueError: If datum is not spherical.
        '''
        if not isinstance(datum, Datum):
            raise TypeError('%r not a %s: %r' % ('datum', Datum.__name__, datum))
        if not datum.isSpherical:
            raise ValueError('%r not %s: %r' % ('datum', 'spherical', datum))
        self._update(datum != self._datum)
        self._datum = datum

    def finalBearingTo(self, other, wrap=False):
        '''Return the final bearing (reverse azimuth) from this
           to an other point.

           @param other: The other point (spherical LatLon).
           @keyword wrap: Optionally, unroll the longitudinal delta
                          within the M{-180..+180} range (bool).

           @return: Final bearing (compass degrees).

           @raise TypeError: The other point is not spherical.

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(48.857, 2.351)
           >>> b = p.finalBearingTo(q)  # 157.9
        '''
        self.others(other)

        # final bearing is the reverse of initial one
        b = other.initialBearingTo(self, wrap=wrap) + 180
        if b > 360:
            b -= 360
        return b  # 0..360

    @property
    def isEllipsoidal(self):
        '''Check whether this I{LatLon} is ellipsoidal (bool).
        '''
        return self.datum.isEllipsoidal

    @property
    def isSpherical(self):
        '''Check whether this I{LatLon} is spherical (bool).
        '''
        return self.datum.isSpherical

    def maxLat(self, bearing):
        '''Return the maximum latitude reached when travelling
           on a great circle on given bearing from this point
           (based on 'Clairaut's formula').

           The maximum latitude is independent of longitude,
           it is the same for all points on a given latitude.

           Negate the result for the minimum latitude (on
           the Southern hemisphere).

           @param bearing: Initial bearing (compass degrees).

           @return: Maximum latitude (degrees90).

           @JSname: I{maxLatitude}.
        '''
        m = acos(abs(sin(radians(bearing)) * cos(radians(self.lat))))
        return degrees90(m)

    def minLat(self, bearing):
        '''Return the minimum latitude reached when travelling
           on a great circle on given bearing from this point.
           See method L{maxLat} for more details.

           @param bearing: Initial bearing (compass degrees).

           @return: Minimum latitude (degrees90).

           @JSname: I{minLatitude}.
        '''
        return -self.maxLat(bearing)

    def parse(self, strll, height=0, sep=','):
        '''Parse a string representing lat-/longitude point and
           return a LatLon.

           The lat- and longitude must be separated by a sep[arator]
           character.  If height is present it must follow and be
           separated by another sep[arator].  Lat- and longitude
           may be swapped, provided at least one ends with the
           proper compass direction.

           For more details, see functions L{parse3llh} and L{parseDMS}
           in module L{dms}.

           @param strll: Lat, lon [, height] (string).
           @keyword height: Optional , default height (meter).
           @keyword sep: Optional separator (string).

           @return: The point (spherical LatLon).

           @raise ValueError: Invalid strll.
        '''
        return self.classof(*parse3llh(strll, height=height, sep=sep))

    def _rhumb3(self, other):
        '''(INTERNAL) Rhumb_ helper function.

           @param other: The other point (spherical LatLon).
        '''
        self.others(other)

        a1, b1 = self.to2ab()
        a2, b2 = other.to2ab()
        # if |db| > 180 take shorter rhumb
        # line across the anti-meridian
        db = wrapPI(b2 - b1)
        dp = log(tanPI_2_2(a2) / tanPI_2_2(a1))
        return (a2 - a1), db, dp

    def rhumbBearingTo(self, other):
        '''Return the initial bearing (forward azimuth) from this
           to an other point along a rhumb (loxodrome) line.

           @param other: The other point (spherical LatLon).

           @return: Initial bearing (compass degrees).

           @raise TypeError: The other point is not spherical.

           @example:

           >>> p = LatLon(51.127, 1.338)
           >>> q = LatLon(50.964, 1.853)
           >>> b = p.rhumbBearingTo(q)  # 116.7
        '''
        _, db, dp = self._rhumb3(other)

        return degrees360(atan2(db, dp))

    def rhumbDestination(self, distance, bearing, radius=R_M, height=None):
        '''Return the destination point having travelled along
           a rhumb (loxodrome) line from this point the given
           distance on the given bearing.

           @param distance: Distance travelled (same units as radius).
           @param bearing: Bearing from this point (compass degrees).
           @keyword radius: Optional, mean earth radius (meter).
           @keyword height: Optional height, overriding the default
                            height (meter or same unit as radius).

           @return: The destination point (spherical LatLon).

           @example:

           >>> p = LatLon(51.127, 1.338)
           >>> q = p.rhumbDestination(40300, 116.7)  # 50.9642°N, 001.8530°E

           @JSname: I{rhumbDestinationPoint}
        '''
        a1, b1 = self.to2ab()

        r = float(distance) / float(radius)  # angular distance in radians
        t = radians(bearing)

        da = r * cos(t)
        a2 = a1 + da
        # normalize latitude if past pole
        if a2 > PI_2:
            a2 =  PI - a2
        elif a2 < -PI_2:
            a2 = -PI - a2

        dp = log(tanPI_2_2(a2) / tanPI_2_2(a1))
        # E-W course becomes ill-conditioned with 0/0
        if abs(dp) > EPS:
            q = da / dp
        else:
            q = cos(a1)

        if abs(q) > EPS:
            b2 = b1 + r * sin(t) / q
        else:
            b2 = b1

        h = self.height if height is None else height
        return self.classof(degrees90(a2), degrees180(b2), height=h)

    def rhumbDistanceTo(self, other, radius=R_M):
        '''Return the distance from this to an other point along
           a rhumb (loxodrome) line.

           @param other: The other point (spherical LatLon).
           @keyword radius: Optional mean radius of earth (scalar, default meter).

           @return: Distance (in the same units as radius).

           @raise TypeError: The other point is not spherical.

           @example:

           >>> p = LatLon(51.127, 1.338)
           >>> q = LatLon(50.964, 1.853)
           >>> d = p.rhumbDistanceTo(q)  # 403100
        '''
        # see <http://www.edwilliams.org/avform.htm#Rhumb>
        da, db, dp = self._rhumb3(other)

        # on Mercator projection, longitude distances shrink
        # by latitude; the 'stretch factor' q becomes ill-
        # conditioned along E-W line (0/0); use an empirical
        # tolerance to avoid it
        if abs(dp) > EPS:
            q = da / dp
        else:
            a, _ = self.to2ab()
            q = cos(a)

        return float(radius) * hypot(da, q * db)

    def rhumbMidpointTo(self, other, height=None):
        '''Return the (loxodromic) midpoint between this and
           an other point.

           @param other: The other point (spherical LatLon).
           @keyword height: Optional height, overriding the mean height
                            (meter or same unit as radius).

           @return: The midpoint (spherical LatLon).

           @raise TypeError: The other point is not spherical.

           @example:

           >>> p = LatLon(51.127, 1.338)
           >>> q = LatLon(50.964, 1.853)
           >>> m = p.rhumb_midpointTo(q)
           >>> m.toStr()  # '51.0455°N, 001.5957°E'
        '''
        self.others(other)

        # see <http://mathforum.org/library/drmath/view/51822.html>
        a1, b1 = self.to2ab()
        a2, b2 = other.to2ab()
        if abs(b2 - b1) > PI:
            b1 += PI2  # crossing anti-meridian

        a3 = favg(a1, a2)
        b3 = favg(b1, b2)

        f1 = tanPI_2_2(a1)
        if abs(f1) > EPS:
            f2 = tanPI_2_2(a2)
            f = f2 / f1
            if abs(f) > EPS:
                f = log(f)
                if abs(f) > EPS:
                    f3 = tanPI_2_2(a3)
                    b3 = fsum_(b1 * log(f2),
                              -b2 * log(f1), (b2 - b1) * log(f3)) / f

        h = self._havg(other) if height is None else height
        return self.classof(degrees90(a3), degrees180(b3), height=h)

    def toWm(self, radius=R_MA):
        '''Convert this I{LatLon} point to a WM coordinate.

           See function L{toWm} in module L{webmercator} for details.

           @keyword radius: Optional earth radius (meter).

           @return: The WM coordinate (L{Wm}).
        '''
        from webmercator import toWm  # PYCHOK recursive import
        return toWm(self, radius=radius)

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
