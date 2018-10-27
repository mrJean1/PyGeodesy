
# -*- coding: utf-8 -*-

u'''(INTERNAL) Spherical base classes.

Pure Python implementation of geodetic (lat-/longitude) functions,
transcribed in part from JavaScript originals by I{(C) Chris Veness 2011-2016}
and published under the same MIT Licence**, see
U{Latitude/Longitude<http://www.Movable-Type.co.UK/scripts/latlong.html>}.

@newfield example: Example, Examples
'''

from bases import LatLonHeightBase
from datum import R_M, R_MA, Datum, Datums
from dms   import parse3llh
from fmath import EPS, acos1, favg, fsum_
from utily import PI, PI2, PI_2, degrees90, degrees180, degrees360, \
                  property_RO, tanPI_2_2, wrapPI

from math import atan2, cos, hypot, log, radians, sin

# XXX the following classes are listed only to get
# Epydoc to include class and method documentation
__all__ = ('LatLonSphericalBase',)  # for documentation
__version__ = '18.10.26'


class LatLonSphericalBase(LatLonHeightBase):
    '''(INTERNAL) Base class for spherical C{LatLon}.
    '''
    _datum = Datums.Sphere  #: (INTERNAL) XXX TBD

    def bearingTo2(self, other, wrap=False, raiser=False):
        '''Return the initial and final bearing (forward and reverse
           azimuth) from this to an other point.

           @param other: The other point (C{LatLon}).
           @keyword wrap: Wrap and unroll longitudes (C{bool}).
           @keyword raiser: Optionally, raise L{CrossError} (C{bool}).

           @return: 2-Tuple (initial, final) bearings (compass C{degrees360}).

           @raise TypeError: The I{other} point is not spherical.

           @see: Methods C{initialBearingTo} and C{finalBearingTo}.
        '''
        return (self.initialBearingTo(other, wrap=wrap, raiser=raiser),  # PYCHOK expected
                self.finalBearingTo(other, wrap=wrap, raiser=raiser))

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

    def finalBearingTo(self, other, wrap=False, raiser=False):
        '''Return the final bearing (reverse azimuth) from this to
           an other point.

           @param other: The other point (spherical C{LatLon}).
           @keyword wrap: Wrap and unroll longitudes (C{bool}).
           @keyword raiser: Optionally, raise L{CrossError} (C{bool}).

           @return: Final bearing (compass C{degrees360}).

           @raise TypeError: The I{other} point is not spherical.

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(48.857, 2.351)
           >>> b = p.finalBearingTo(q)  # 157.9
        '''
        self.others(other)

        # final bearing is the reverse of the other initial one
        b = other.initialBearingTo(self, wrap=wrap, raiser=raiser)
        return (b + 180) % 360

    @property_RO
    def isEllipsoidal(self):
        '''Check whether this C{LatLon} is ellipsoidal (C{bool}).
        '''
        return self.datum.isEllipsoidal

    @property_RO
    def isSpherical(self):
        '''Check whether this C{LatLon} is spherical (C{bool}).
        '''
        return self.datum.isSpherical

    def maxLat(self, bearing):
        '''Return the maximum latitude reached when travelling
           on a great circle on given bearing from this point
           based on Clairaut's formula.

           The maximum latitude is independent of longitude
           and the same for all points on a given latitude.

           Negate the result for the minimum latitude (on the
           Southern hemisphere).

           @param bearing: Initial bearing (compass C{degrees360}).

           @return: Maximum latitude (C{degrees90}).

           @JSname: I{maxLatitude}.
        '''
        a, _ = self.to2ab()
        m = acos1(abs(sin(radians(bearing)) * cos(a)))
        return degrees90(m)

    def minLat(self, bearing):
        '''Return the minimum latitude reached when travelling
           on a great circle on given bearing from this point.

           @param bearing: Initial bearing (compass C{degrees360}).

           @return: Minimum latitude (C{degrees90}).

           @see: Method L{maxLat} for more details.

           @JSname: I{minLatitude}.
        '''
        return -self.maxLat(bearing)

    def parse(self, strll, height=0, sep=','):
        '''Parse a string representing lat-/longitude point and
           return a C{LatLon}.

           The lat- and longitude must be separated by a sep[arator]
           character.  If height is present it must follow and be
           separated by another sep[arator].  Lat- and longitude
           may be swapped, provided at least one ends with the
           proper compass direction.

           For more details, see functions L{parse3llh} and L{parseDMS}
           in module L{dms}.

           @param strll: Lat, lon [, height] (C{str}).
           @keyword height: Optional , default height (C{meter}).
           @keyword sep: Optional separator (C{str}).

           @return: The point (spherical C{LatLon}).

           @raise ValueError: Invalid I{strll}.
        '''
        return self.classof(*parse3llh(strll, height=height, sep=sep))

    def _rhumb3(self, other):
        '''(INTERNAL) Rhumb_ helper function.

           @param other: The other point (spherical C{LatLon}).
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
        '''Return the initial bearing (forward azimuth) from this to
           an other point along a rhumb (loxodrome) line.

           @param other: The other point (spherical C{LatLon}).

           @return: Initial bearing (compass C{degrees360}).

           @raise TypeError: The I{other} point is not spherical.

           @example:

           >>> p = LatLon(51.127, 1.338)
           >>> q = LatLon(50.964, 1.853)
           >>> b = p.rhumbBearingTo(q)  # 116.7
        '''
        _, db, dp = self._rhumb3(other)

        return degrees360(atan2(db, dp))

    def rhumbDestination(self, distance, bearing, radius=R_M, height=None):
        '''Return the destination point having travelled along a rhumb
           (loxodrome) line from this point the given distance on the
           given bearing.

           @param distance: Distance travelled (C{meter}, same units as
                            I{radius}).
           @param bearing: Bearing from this point (compass C{degrees360}).
           @keyword radius: Optional, mean earth radius (C{meter}).
           @keyword height: Optional height, overriding the default
                            height (C{meter}, same unit as I{radius}).

           @return: The destination point (spherical C{LatLon}).

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
        q  = (da / dp) if abs(dp) > EPS else cos(a1)
        b2 = (b1 + r * sin(t) / q) if abs(q) > EPS else b1

        h = self.height if height is None else height
        return self.classof(degrees90(a2), degrees180(b2), height=h)

    def rhumbDistanceTo(self, other, radius=R_M):
        '''Return the distance from this to an other point along a rhumb
           (loxodrome) line.

           @param other: The other point (spherical C{LatLon}).
           @keyword radius: Optional mean radius of earth (scalar, default meter).

           @return: Distance (C{meter}, the same units as I{radius}).

           @raise TypeError: The I{other} point is not spherical.

           @example:

           >>> p = LatLon(51.127, 1.338)
           >>> q = LatLon(50.964, 1.853)
           >>> d = p.rhumbDistanceTo(q)  # 403100
        '''
        # see <http://www.EdWilliams.org/avform.htm#Rhumb>
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
                            (C{meter}).

           @return: The midpoint (spherical C{LatLon}).

           @raise TypeError: The I{other} point is not spherical.

           @example:

           >>> p = LatLon(51.127, 1.338)
           >>> q = LatLon(50.964, 1.853)
           >>> m = p.rhumb_midpointTo(q)
           >>> m.toStr()  # '51.0455°N, 001.5957°E'
        '''
        self.others(other)

        # see <http://MathForum.org/library/drmath/view/51822.html>
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
        '''Convert this C{LatLon} point to a I{WM} coordinate.

           @keyword radius: Optional earth radius (C{meter}).

           @return: The WM coordinate (L{Wm}).

           @see: Function L{toWm} in module L{webmercator} for details.
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
