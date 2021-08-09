
# -*- coding: utf-8 -*-

u'''(INTERNAL) Spherical geodesy bases.

Base classes C{CartesianSphericalBase} and C{LatLonSphericalBase}
imported by L{sphericalNvector} or L{sphericalTrigonometry}.

Pure Python implementation of geodetic (lat-/longitude) functions,
transcoded in part from JavaScript originals by I{(C) Chris Veness 2011-2016}
and published under the same MIT Licence**, see
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division

from pygeodesy.basics import isnear0, isnon0, map1
from pygeodesy.cartesianBase import CartesianBase
from pygeodesy.datums import Datums, _spherical_datum
from pygeodesy.ellipsoids import R_M, R_MA
from pygeodesy.errors import IntersectionError
from pygeodesy.fmath import favg, fdot, hypot
from pygeodesy.interns import EPS, NN, PI, PI2, PI_2, _COMMA_, \
                             _datum_, _distant_, _exceed_PI_radians_, \
                             _name_, _near_concentric_, _too_, \
                             _1_0, _180_0, _360_0
from pygeodesy.latlonBase import LatLonBase, _trilaterate5  # PYCHOK passed
from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.namedTuples import Bearing2Tuple
from pygeodesy.nvectorBase import NvectorBase, _xattrs  # streprs
from pygeodesy.props import property_doc_
from pygeodesy.units import Bearing_, Height, Radians_, Radius, Radius_
from pygeodesy.utily import acos1, atan2b, degrees90, degrees180, \
                            sincos2, tanPI_2_2, wrapPI

from math import cos, log, sin, sqrt

__all__ = ()
__version__ = '21.08.07'


def _angular(distance, radius):  # PYCHOK for export
    '''(INTERNAL) Return the angular distance in C{radians}.

       @raise UnitError: Invalid B{C{distance}} or B{C{radius}}.
    '''
    return Radians_(distance / Radius_(radius=radius), low=EPS)


def _rads3(rad1, rad2, radius):  # in .sphericalTrigonometry
    '''(INTERNAL) Convert radii to radians.
    '''
    r1 = Radius_(rad1=rad1)
    r2 = Radius_(rad2=rad2)
    if radius is not None:  # convert radii to radians
        r1 = _angular(r1, radius)
        r2 = _angular(r2, radius)

    x = r1 < r2
    if x:
        r1, r2 = r2, r1
    if r1 > PI:
        raise IntersectionError(rad1=rad1, rad2=rad2,
                                txt=_exceed_PI_radians_)
    return r1, r2, x


class CartesianSphericalBase(CartesianBase):
    '''(INTERNAL) Base class for spherical C{Cartesian}s.
    '''
    _datum = Datums.Sphere  # L{Datum}

    def intersections2(self, rad1, other, rad2, radius=R_M):
        '''Compute the intersection points of two circles each defined
           by a center point and a radius.

           @arg rad1: Radius of the this circle (C{meter} or C{radians},
                      see B{C{radius}}).
           @arg other: Center of the other circle (C{Cartesian}).
           @arg rad2: Radius of the other circle (C{meter} or C{radians},
                      see B{C{radius}}).
           @kwarg radius: Mean earth radius (C{meter} or C{None} if both
                          B{C{rad1}} and B{C{rad2}} are given in C{radians}).

           @return: 2-Tuple of the intersection points, each C{Cartesian}.
                    The intersection points are the same C{Cartesian}
                    instance for abutting circles, aka I{radical center}.

           @raise IntersectionError: Concentric, antipodal, invalid or
                                     non-intersecting circles.

           @raise TypeError: If B{C{other}} is not C{Cartesian}.

           @raise ValueError: Invalid B{C{rad1}}, B{C{rad2}} or B{C{radius}}.

           @see: U{Calculating intersection of two Circles
                 <https://GIS.StackExchange.com/questions/48937/
                 calculating-intersection-of-two-circles>} and method
                 or function C{trilaterate3d2}.
        '''
        x1, x2 = self, self.others(other)
        r1, r2, x = _rads3(rad1, rad2, radius)
        if x:
            x1, x2 = x2, x1
        try:
            n, q = x1.cross(x2), x1.dot(x2)
            n2, q1 = n.length2, (_1_0 - q**2)
            if n2 < EPS or isnear0(q1):
                raise ValueError(_near_concentric_)
            c1, c2 = cos(r1), cos(r2)
            x0 = x1.times((c1 - q * c2) / q1).plus(
                 x2.times((c2 - q * c1) / q1))
            n1 = _1_0 - x0.length2
            if n1 < EPS:
                raise ValueError(_too_(_distant_))
        except ValueError as x:
            raise IntersectionError(center=self, rad1=rad1,
                                    other=other, rad2=rad2, txt=str(x))
        n = n.times(sqrt(n1 / n2))
        if n.length > EPS:
            x1 = x0.plus(n)
            x2 = x0.minus(n)
        else:  # abutting circles
            x1 = x2 = x0

        return (_xattrs(x1, self, _datum_, _name_),
                _xattrs(x2, self, _datum_, _name_))


class LatLonSphericalBase(LatLonBase):
    '''(INTERNAL) Base class for spherical C{LatLon}s.
    '''
    _datum = Datums.Sphere  # spherical L{Datum}

    def __init__(self, lat, lon, height=0, datum=None, name=NN):
        '''Create a spherical C{LatLon} point frome the given
           lat-, longitude and height on the given datum.

           @arg lat: Latitude (C{degrees} or DMS C{[N|S]}).
           @arg lon: Longitude (C{degrees} or DMS C{str[E|W]}).
           @kwarg height: Optional elevation (C{meter}, the same units
                          as the datum's half-axes).
           @kwarg datum: Optional, spherical datum to use (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2}, L{a_f2Tuple})
                         or C{scalar} earth radius).
           @kwarg name: Optional name (string).

           @raise TypeError: If B{C{datum}} invalid or not
                             not spherical.

           @example:

            >>> p = LatLon(51.4778, -0.0016)  # height=0, datum=Datums.WGS84
        '''
        LatLonBase.__init__(self, lat, lon, height=height, name=name)
        if datum not in (None, self.datum):
            self._datum = _spherical_datum(datum, name=self.name, raiser=True)

    def bearingTo2(self, other, wrap=False, raiser=False):
        '''Return the initial and final bearing (forward and reverse
           azimuth) from this to an other point.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).
           @kwarg raiser: Optionally, raise L{CrossError} (C{bool}).

           @return: A L{Bearing2Tuple}C{(initial, final)}.

           @raise TypeError: The B{C{other}} point is not spherical.

           @see: Methods C{initialBearingTo} and C{finalBearingTo}.
        '''
        # .initialBearingTo is inside .-Nvector and .-Trigonometry
        i = self.initialBearingTo(other, wrap=wrap, raiser=raiser)  # PYCHOK .initialBearingTo
        f = self.finalBearingTo(  other, wrap=wrap, raiser=raiser)
        return Bearing2Tuple(i, f, name=self.name)

    @property_doc_(''' this point's datum (L{Datum}).''')
    def datum(self):
        '''Get this point's datum (L{Datum}).
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set this point's datum I{without conversion}.

           @arg datum: New spherical datum (L{Datum}, L{Ellipsoid},
                       L{Ellipsoid2}, L{a_f2Tuple}) or C{scalar}
                       earth radius).

           @raise TypeError: If B{C{datum}} invalid or not
                             not spherical.
        '''
        d = _spherical_datum(datum, name=self.name, raiser=True)
        self._update(d != self._datum)
        self._datum = d

    def finalBearingTo(self, other, wrap=False, raiser=False):
        '''Return the final bearing (reverse azimuth) from this to
           an other point.

           @arg other: The other point (spherical C{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).
           @kwarg raiser: Optionally, raise L{CrossError} (C{bool}).

           @return: Final bearing (compass C{degrees360}).

           @raise TypeError: The B{C{other}} point is not spherical.

           @example:

            >>> p = LatLon(52.205, 0.119)
            >>> q = LatLon(48.857, 2.351)
            >>> b = p.finalBearingTo(q)  # 157.9
        '''
        self.others(other)

        # final bearing is the reverse of the other, initial one;
        # .initialBearingTo is inside .-Nvector and .-Trigonometry
        b = other.initialBearingTo(self, wrap=wrap, raiser=raiser)
        return (b + _180_0) % _360_0  # == wrap360 since b >= 0

    def maxLat(self, bearing):
        '''Return the maximum latitude reached when travelling
           on a great circle on given bearing from this point
           based on Clairaut's formula.

           The maximum latitude is independent of longitude
           and the same for all points on a given latitude.

           Negate the result for the minimum latitude (on the
           Southern hemisphere).

           @arg bearing: Initial bearing (compass C{degrees360}).

           @return: Maximum latitude (C{degrees90}).

           @raise ValueError: Invalid B{C{bearing}}.

           @JSname: I{maxLatitude}.
        '''
        m = acos1(abs(sin(Bearing_(bearing)) * cos(self.phi)))
        return degrees90(m)

    def minLat(self, bearing):
        '''Return the minimum latitude reached when travelling
           on a great circle on given bearing from this point.

           @arg bearing: Initial bearing (compass C{degrees360}).

           @return: Minimum latitude (C{degrees90}).

           @see: Method L{maxLat} for more details.

           @raise ValueError: Invalid B{C{bearing}}.

           @JSname: I{minLatitude}.
        '''
        return -self.maxLat(bearing)

    def parse(self, strllh, height=0, sep=_COMMA_, name=NN):
        '''Parse a string representing a similar, spherical C{LatLon}
           point, consisting of C{"lat, lon[, height]"}.

           @arg strllh: Lat, lon and optional height (C{str}),
                        see function L{parse3llh}.
           @kwarg height: Optional, default height (C{meter}).
           @kwarg sep: Optional separator (C{str}).
           @kwarg name: Optional instance name (C{str}),
                        overriding this name.

           @return: The similar point (spherical C{LatLon}).

           @raise ParseError: Invalid B{C{strllh}}.
        '''
        from pygeodesy.dms import parse3llh
        r = self.classof(*parse3llh(strllh, height=height, sep=sep))
        if name:
            r.rename(name)
        return r

    def _rhumb3(self, other):
        '''(INTERNAL) Rhumb_ helper function.

           @arg other: The other point (spherical C{LatLon}).
        '''
        self.others(other)

        a1, b1 = self.philam
        a2, b2 = other.philam
        # if |db| > 180 take shorter rhumb
        # line across the anti-meridian
        db = wrapPI(b2 - b1)
        dp = log(tanPI_2_2(a2) / tanPI_2_2(a1))
        return (a2 - a1), db, dp

    def rhumbBearingTo(self, other):
        '''Return the initial bearing (forward azimuth) from this to
           an other point along a rhumb (loxodrome) line.

           @arg other: The other point (spherical C{LatLon}).

           @return: Initial bearing (compass C{degrees360}).

           @raise TypeError: The B{C{other}} point is not spherical.

           @example:

            >>> p = LatLon(51.127, 1.338)
            >>> q = LatLon(50.964, 1.853)
            >>> b = p.rhumbBearingTo(q)  # 116.7
        '''
        _, db, dp = self._rhumb3(other)
        return atan2b(db, dp)

    def rhumbDestination(self, distance, bearing, radius=R_M, height=None):
        '''Return the destination point having travelled along a rhumb
           (loxodrome) line from this point the given distance on the
           given bearing.

           @arg distance: Distance travelled (C{meter}, same units as
                          B{C{radius}}).
           @arg bearing: Bearing from this point (compass C{degrees360}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg height: Optional height, overriding the default height
                          (C{meter}, same unit as B{C{radius}}).

           @return: The destination point (spherical C{LatLon}).

           @raise ValueError: Invalid B{C{distance}}, B{C{bearing}},
                              B{C{radius}} or B{C{height}}.

           @example:

            >>> p = LatLon(51.127, 1.338)
            >>> q = p.rhumbDestination(40300, 116.7)  # 50.9642°N, 001.8530°E

           @JSname: I{rhumbDestinationPoint}
        '''
        r = _angular(distance, radius)

        a1, b1 = self.philam
        sb, cb = sincos2(Bearing_(bearing))

        da = r * cb
        a2 = a1 + da
        # normalize latitude if past pole
        if a2 > PI_2:
            a2 =  PI - a2
        elif a2 < -PI_2:
            a2 = -PI - a2

        dp = log(tanPI_2_2(a2) / tanPI_2_2(a1))
        # E-W course becomes ill-conditioned with 0/0
        q  = (da / dp) if abs(dp) > EPS else cos(a1)
        b2 = (b1 + r * sb / q) if abs(q) > EPS else b1

        h = self.height if height is None else Height(height)
        return self.classof(degrees90(a2), degrees180(b2), height=h)

    def rhumbDistanceTo(self, other, radius=R_M):
        '''Return the distance from this to an other point along a rhumb
           (loxodrome) line.

           @arg other: The other point (spherical C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None}.

           @return: Distance (C{meter}, the same units as B{C{radius}}
                    or C{radians} if B{C{radius}} is C{None}).

           @raise TypeError: The B{C{other}} point is not spherical.

           @raise ValueError: Invalid B{C{radius}}.

           @example:

            >>> p = LatLon(51.127, 1.338)
            >>> q = LatLon(50.964, 1.853)
            >>> d = p.rhumbDistanceTo(q)  # 403100
        '''
        # see <https://www.EdWilliams.org/avform.htm#Rhumb>
        da, db, dp = self._rhumb3(other)

        # on Mercator projection, longitude distances shrink
        # by latitude; the 'stretch factor' q becomes ill-
        # conditioned along E-W line (0/0); use an empirical
        # tolerance to avoid it
        q = (da / dp) if abs(dp) > EPS else cos(self.phi)
        r = hypot(da, q * db)
        return r if radius is None else (Radius(radius) * r)

    def rhumbMidpointTo(self, other, height=None):
        '''Return the (loxodromic) midpoint between this and
           an other point.

           @arg other: The other point (spherical LatLon).
           @kwarg height: Optional height, overriding the mean height
                          (C{meter}).

           @return: The midpoint (spherical C{LatLon}).

           @raise TypeError: The B{C{other}} point is not spherical.

           @raise ValueError: Invalid B{C{height}}.

           @example:

            >>> p = LatLon(51.127, 1.338)
            >>> q = LatLon(50.964, 1.853)
            >>> m = p.rhumb_midpointTo(q)
            >>> m.toStr()  # '51.0455°N, 001.5957°E'
        '''
        self.others(other)

        # see <https://MathForum.org/library/drmath/view/51822.html>
        a1, b1 = self.philam
        a2, b2 = other.philam
        if abs(b2 - b1) > PI:
            b1 += PI2  # crossing anti-meridian

        a3 = favg(a1, a2)
        b3 = favg(b1, b2)

        f1 = tanPI_2_2(a1)
        if isnon0(f1):
            f2 = tanPI_2_2(a2)
            f = f2 / f1
            if isnon0(f):
                f = log(f)
                if isnon0(f):
                    f3 = tanPI_2_2(a3)
                    b3 = fdot(map1(log, f1, f2, f3),
                                       -b2, b1, b2 - b1) / f

        h = self._havg(other) if height is None else Height(height)
        return self.classof(degrees90(a3), degrees180(b3), height=h)

    def toNvector(self, Nvector=NvectorBase, **Nvector_kwds):  # PYCHOK signature
        '''Convert this point to C{Nvector} components, I{including
           height}.

           @kwarg Nvector_kwds: Optional, additional B{C{Nvector}}
                                keyword arguments, ignored if
                                C{B{Nvector} is None}.

           @return: An B{C{Nvector}} or a L{Vector4Tuple}C{(x, y, z, h)}
                    if B{C{Nvector}} is C{None}.

           @raise TypeError: Invalid B{C{Nvector}} or B{C{Nvector_kwds}}.
        '''
        return LatLonBase.toNvector(self, Nvector=Nvector, **Nvector_kwds)

    def toWm(self, radius=R_MA):
        '''Convert this point to a I{WM} coordinate.

           @kwarg radius: Optional earth radius (C{meter}).

           @return: The WM coordinate (L{Wm}).

           @see: Function L{toWm} in module L{webmercator} for details.
        '''
        from pygeodesy.webmercator import toWm
        return toWm(self, radius=radius)


__all__ += _ALL_DOCS(CartesianSphericalBase, LatLonSphericalBase)

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
