# -*- coding: utf-8 -*-

u'''Slightly enhanced versions of classes U{PolygonArea
<https://GeographicLib.SourceForge.io/1.52/python/code.html#
module-geographiclib.polygonarea>} and C{Accumulator} from
I{Karney}'s Python U{geographiclib
<https://GeographicLib.SourceForge.io/1.52/python/index.html>}.

Class L{GeodesicAreaExact} is intended to work with instances
of class L{GeodesicExact} and of I{wrapped} class C{Geodesic},
see module L{pygeodesy.karney}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2008-2023)
and licensed under the MIT/X11 License.  For more information, see the
U{GeographicLib<https://GeographicLib.SourceForge.io>} documentation.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # PYCHOK semicolon

# from pygeodesy.basics import isodd, unsigned0  # from .karney
from pygeodesy.constants import NAN, _0_0, _0_5, _720_0
from pygeodesy.interns import NN, _COMMASPACE_
from pygeodesy.karney import Area3Tuple, _diff182, GeodesicError, isodd, \
                            _norm180, _remainder, _sum2_, unsigned0
from pygeodesy.lazily import _ALL_DOCS, printf
from pygeodesy.named import callername, _Dict, _NamedBase, pairs
from pygeodesy.props import Property, Property_RO, property_RO
# from pygeodesy.streprs import pairs  # from .named

from math import fmod

__all__ = ()
__version__ = '23.08.19'


class GeodesicAreaExact(_NamedBase):
    '''Area and perimeter of a geodesic polygon, an enhanced
       version of I{Karney}'s Python class U{PolygonArea
       <https://GeographicLib.SourceForge.io/html/python/
       code.html#module-geographiclib.polygonarea>} using
       the more accurate surface area.

       @note: The name of this class C{*Exact} is a misnomer, see
              I{Karney}'s comments at C++ attribute U{GeodesicExact._c2
              <https://GeographicLib.SourceForge.io/C++/doc/
              GeodesicExact_8cpp_source.html>}.
    '''
    _Area    =  None
    _g_gX    =  None  # Exact or not
    _lat0    = _lon0 = \
    _lat1    = _lon1 = NAN
    _mask    =  0
    _num     =  0
    _Peri    =  None
    _verbose =  False
    _xings   =  0

    def __init__(self, geodesic, polyline=False, name=NN):
        '''New L{GeodesicAreaExact} instance.

           @arg geodesic: A geodesic (L{GeodesicExact}, I{wrapped}
                          C{Geodesic} or L{GeodesicSolve}).
           @kwarg polyline: If C{True}, compute the perimeter only,
                            otherwise area and perimeter (C{bool}).
           @kwarg name: Optional name (C{str}).

           @raise GeodesicError: Invalid B{C{geodesic}}.
        '''
        try:  # results returned as L{GDict}
            if not (callable(geodesic._GDictDirect) and
                    callable(geodesic._GDictInverse)):
                raise TypeError
        except (AttributeError, TypeError):
            raise GeodesicError(geodesic=geodesic)

        self._g_gX = g = geodesic
        # use the class-level Caps since the values
        # differ between GeodesicExact and Geodesic
        self._mask = g.DISTANCE | g.LATITUDE | g.LONGITUDE
        self._Peri = _Accumulator(name='_Peri')
        if not polyline:  # perimeter and area
            self._mask |=  g.AREA | g.LONG_UNROLL
            self._Area  = _Accumulator(name='_Area')
        if g.debug:  # PYCHOK no cover
            self.verbose = True  # debug as verbosity
        if name:
            self.name = name

    def AddEdge(self, azi, s):
        '''Add another polygon edge.

           @arg azi: Azimuth at the current point (compass
                     C{degrees360}).
           @arg s: Length of the edge (C{meter}).
        '''
        if self.num < 1:
            raise GeodesicError(num=self.num)
        r = self._Direct(azi, s)
        p = self._Peri.Add(s)
        if self._Area:
            a = self._Area.Add(r.S12)
            self._xings += r.xing
        else:
            a = NAN
        self._lat1 = r.lat2
        self._lon1 = r.lon2
        self._num += 1
        if self.verbose:  # PYCHOK no cover
            self._print(self.num, p, a, r, lat1=r.lat2, lon1=r.lon2,
                                           azi=azi, s=s)
        return self.num

    def AddPoint(self, lat, lon):
        '''Add another polygon point.

           @arg lat: Latitude of the point (C{degrees}).
           @arg lon: Longitude of the point (C{degrees}).
        '''
        if self.num > 0:
            r = self._Inverse(self.lat1, self.lon1, lat, lon)
            s = r.s12
            p = self._Peri.Add(s)
            if self._Area:
                a = self._Area.Add(r.S12)
                self._xings += r.xing
            else:
                a = NAN
        else:
            self._lat0 = lat
            self._lon0 = lon
            a = p = s = _0_0
            r = None
        self._lat1 = lat
        self._lon1 = lon
        self._num += 1
        if self.verbose:  # PYCHOK no cover
            self._print(self.num, p, a, r, lat1=lat, lon1=lon, s=s)
        return self.num

    @Property_RO
    def area0x(self):
        '''Get the ellipsoid's surface area (C{meter} I{squared}),
           more accurate for very I{oblate} ellipsoids.
        '''
        return self.ellipsoid.areax  # not .area!

    area0 = area0x  # for C{geographiclib} compatibility

    def Compute(self, reverse=False, sign=True):
        '''Compute the accumulated perimeter and area.

           @kwarg reverse: If C{True}, clockwise traversal counts as a
                           positive area instead of counter-clockwise
                           (C{bool}).
           @kwarg sign: If C{True}, return a signed result for the area if
                        the polygon is traversed in the "wrong" direction
                        instead of returning the area for the rest of the
                        earth.

           @return: L{Area3Tuple}C{(number, perimeter, area)} with the
                    number of points, the perimeter in C{meter} and the
                    area in C{meter**2}.  The perimeter includes the
                    length of a final edge, connecting the current to
                    the initial point, if this polygon was initialized
                    with C{polyline=False}.  For perimeter only, i.e.
                    C{polyline=True}, area is C{NAN}.

           @note: Arbitrarily complex polygons are allowed.  In the case
                  of self-intersecting polygons, the area is accumulated
                  "algebraically".  E.g., the areas of the 2 loops in a
                  I{figure-8} polygon will partially cancel.

           @note: More points and edges can be added after this call.
        '''
        r, n = None, self.num
        if n < 2:
            p = _0_0
            a = NAN if self.polyline else p
        elif self._Area:
            r = self._Inverse(self.lat1, self.lon1, self.lat0, self.lon0)
            a = self._reduced(r.S12, reverse, sign, r.xing)
            p = self._Peri.Sum(r.s12)
        else:
            p = self._Peri.Sum()
            a = NAN
        if self.verbose:  # PYCHOK no cover
            self._print(n, p, a, r, lat0=self.lat0, lon0=self.lon0)
        return Area3Tuple(n, p, a)

    def _Direct(self, azi, s):
        '''(INTERNAL) Edge helper.
        '''
        lon1 = self.lon1
        r = self._g_gX._GDictDirect(self.lat1, lon1, azi, False, s, self._mask)
        if self._Area:  # aka transitDirect
            # Count crossings of prime meridian exactly as
            # int(ceil(lon2 / 360)) - int(ceil(lon1 / 360))
            # Since we only need the parity of the result we
            # can use std::remquo but this is buggy with g++
            # 4.8.3 and requires C++11.  So instead we do:
            lon1 = fmod(  lon1, _720_0)  # r.lon1
            lon2 = fmod(r.lon2, _720_0)
            # int(True) == 1, int(False) == 0
            r.set_(xing=int(lon2 > 360 or -360 < lon2 <= 0) -
                        int(lon1 > 360 or -360 < lon1 <= 0))
        return r

    @Property_RO
    def ellipsoid(self):
        '''Get this area's ellipsoid (C{Ellipsoid[2]}).
        '''
        return self._g_gX.ellipsoid

    @Property_RO
    def geodesic(self):
        '''Get this area's geodesic object (C{Geodesic[Exact]}).
        '''
        return self._g_gX

    earth = geodesic  # for C{geographiclib} compatibility

    def _Inverse(self, lat1, lon1, lat2, lon2):
        '''(INTERNAL) Point helper.
        '''
        r = self._g_gX._GDictInverse(lat1, lon1, lat2, lon2, self._mask)
        if self._Area:  # aka transit
            # count crossings of prime meridian as +1 or -1
            # if in east or west direction, otherwise 0
            lon1 = _norm180(lon1)
            lon2 = _norm180(lon2)
            lon12, _ = _diff182(lon1, lon2)
            r.set_(xing=int(lon12 > 0 and lon1 <= 0 and lon2 > 0) or
                       -int(lon12 < 0 and lon2 <= 0 and lon1 > 0))
        return r

    @property_RO
    def lat0(self):
        '''Get the first point's latitude (C{degrees}).
        '''
        return self._lat0

    @property_RO
    def lat1(self):
        '''Get the most recent point's latitude (C{degrees}).
        '''
        return self._lat1

    @property_RO
    def lon0(self):
        '''Get the first point's longitude (C{degrees}).
        '''
        return self._lon0

    @property_RO
    def lon1(self):
        '''Get the most recent point's longitude (C{degrees}).
        '''
        return self._lon1

    @property_RO
    def num(self):
        '''Get the current number of points (C{int}).
        '''
        return self._num

    @Property_RO
    def polyline(self):
        '''Is this perimeter only (C{bool}), area NAN?
        '''
        return self._Area is None

    def _print(self, n, p, a, r, **kwds):  # PYCHOK no cover
        '''(INTERNAL) Print a verbose line.
        '''
        d = _Dict(p=p, s12=r.s12 if r else NAN, **kwds)
        if self._Area:
            d.set_(a=a, S12=r.S12 if r else NAN)
        t = _COMMASPACE_.join(pairs(d, prec=10))
        printf('%s %s: %s (%s)', self.named2, n, t, callername(up=2))

    def _reduced(self, S12, reverse, sign, xing):
        '''(INTERNAL) Accumulate and reduce area to allowed range.
        '''
        a0 =  self.area0x
        A  = _Accumulator(self._Area)
        _  =  A.Add(S12)
        a  =  A.Remainder(a0)  # clockwise
        if isodd(self._xings + xing):
            a = A.Add((a0 if a < 0 else -a0) * _0_5)
        if not reverse:
            a = A.Negate()  # counter-clockwise
        # (-area0x/2, area0x/2] if sign else [0, area0x)
        a0_ = a0 if sign else (a0 * _0_5)
        if a > a0_:
            a = A.Add(-a0)
        elif a <= -a0_:
            a = A.Add( a0)
        return unsigned0(a)

    def Reset(self):
        '''Reset this polygon to empty.
        '''
        if self._Area:
            self._Area.Reset()
        self._Peri.Reset()
        self._lat0 = self._lon0 = \
        self._lat1 = self._lon1 = NAN
        self._num  = self._xings = n = 0
        if self.verbose:  # PYCHOK no cover
            printf('%s %s: (%s)', self.named2, n, self.Reset.__name__)
        return n

    Clear = Reset

    def TestEdge(self, azi, s, reverse=False, sign=True):
        '''Compute the properties for a tentative, additional edge

           @arg azi: Azimuth at the current the point (compass C{degrees}).
           @arg s: Length of the edge (C{meter}).
           @kwarg reverse: If C{True}, clockwise traversal counts as a
                           positive area instead of counter-clockwise
                           (C{bool}).
           @kwarg sign: If C{True}, return a signed result for the area if
                        the polygon is traversed in the "wrong" direction
                        instead of returning the area for the rest of the
                        earth.

           @return: L{Area3Tuple}C{(number, perimeter, area)}.

           @raise GeodesicError: No points.
        '''
        n = self.num + 1
        p = self._Peri.Sum(s)
        if self.polyline:
            a, r = NAN, None
        elif n < 2:
            raise GeodesicError(num=self.num)
        else:
            d  = self._Direct(azi, s)
            r  = self._Inverse(d.lat2, d.lon2, self.lat0, self.lon0)
            a  = self._reduced(d.S12 + r.S12, reverse, sign, d.xing + r.xing)
            p += r.s12
        if self.verbose:  # PYCHOK no cover
            self._print(n, p, a, r, azi=azi, s=s)
        return Area3Tuple(n, p, a)

    def TestPoint(self, lat, lon, reverse=False, sign=True):
        '''Compute the properties for a tentative, additional vertex

           @arg lat: Latitude of the point (C{degrees}).
           @arg lon: Longitude of the point (C{degrees}).
           @kwarg reverse: If C{True}, clockwise traversal counts as a
                           positive area instead of counter-clockwise
                           (C{bool}).
           @kwarg sign: If C{True}, return a signed result for the area if
                        the polygon is traversed in the "wrong" direction
                        instead of returning the area for the rest of the
                        earth.

           @return: L{Area3Tuple}C{(number, perimeter, area)}.
        '''
        r, n = None, self.num + 1
        if n < 2:
            p = _0_0
            a = NAN if self.polyline else p
        else:
            i = self._Inverse(self.lat1, self.lon1, lat, lon)
            p = self._Peri.Sum(i.s12)
            if self._Area:
                r  = self._Inverse(lat, lon, self.lat0, self.lon0)
                a  = self._reduced(i.S12 + r.S12, reverse, sign, i.xing + r.xing)
                p += r.s12
            else:
                a = NAN
        if self.verbose:  # PYCHOK no cover
            self._print(n, p, a, r, lat=lat, lon=lon)
        return Area3Tuple(n, p, a)

    def toStr(self, prec=6, sep=_COMMASPACE_, **unused):  # PYCHOK signature
        '''Return this C{GeodesicExactArea} as string.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Separator to join (C{str}).

           @return: Area items (C{str}).
        '''
        n, p, a = self.Compute()
        d = dict(geodesic=self.geodesic, num=n, area=a,
                 perimeter=p, polyline=self.polyline)
        return sep.join(pairs(d, prec=prec))

    @Property
    def verbose(self):
        '''Get the C{verbose} option (C{bool}).
        '''
        return self._verbose

    @verbose.setter  # PYCHOK setter!
    def verbose(self, verbose):  # PYCHOK no cover
        '''Set the C{verbose} option (C{bool}) to print
           a message after each method invokation.
        '''
        self._verbose = bool(verbose)


class PolygonArea(GeodesicAreaExact):
    '''For C{geographiclib} compatibility, sub-class of L{GeodesicAreaExact}.
    '''
    def __init__(self, earth, polyline=False, name=NN):
        '''New L{PolygonArea} instance.

           @arg earth: A geodesic (L{GeodesicExact}, I{wrapped}
                       C{Geodesic} or L{GeodesicSolve}).
           @kwarg polyline: If C{True} perimeter only, otherwise
                            area and perimeter (C{bool}).
           @kwarg name: Optional name (C{str}).

           @raise GeodesicError: Invalid B{C{earth}}.
        '''
        GeodesicAreaExact.__init__(self, earth, polyline=polyline, name=name)


class _Accumulator(_NamedBase):
    '''Like C{math.fsum}, but allowing a running sum.

       Original from I{Karney}'s U{geographiclib
       <https://PyPI.org/project/geographiclib>}C{.accumulator},
       enhanced to return the current sum by most methods.
    '''
    _n =  0  # len()
    _s = _t = _0_0

    def __init__(self, y=0, name=NN):
        '''New L{_Accumulator}.
        '''
        if isinstance(y, _Accumulator):
            self._s, self._t, self._n = y._s, y._t, 1
        elif y:
            self._s, self._n = float(y), 1
        if name:
            self.name = name

    def Add(self, y):
        '''Add a value.

           @return: Current C{sum}.
        '''
        self._n += 1
        self._s, self._t = _sum2_(self._s, self._t, y)
        return self._s  # current .Sum()

    def Negate(self):
        '''Negate sum.

           @return: Current C{sum}.
        '''
        self._s = s = -self._s
        self._t =     -self._t
        return s  # current .Sum()

    @property_RO
    def num(self):
        '''Get the current number of C{Add}itions (C{int}).
        '''
        return self._n

    def Remainder(self, y):
        '''Remainder on division by B{C{y}}.

           @return: Remainder of C{sum} / B{C{y}}.
        '''
        self._s = _remainder(self._s, y)
#       self._t = _remainder(self._t, y)
        self._n = -1
        return self.Add(_0_0)

    def Reset(self, y=0):
        '''Set value from argument.
        '''
        self._n, self._s, self._t = 0, float(y), _0_0

    Set = Reset

    def Sum(self, y=0):
        '''Return C{sum + B{y}}.

           @note: B{C{y}} is included in the returned
                  result, but I{not} accumulated.
        '''
        if y:
            s = _Accumulator(self, name='_Sum')
            s.Add(y)
        else:
            s = self
        return s._s

    def toStr(self, prec=6, sep=_COMMASPACE_, **unused):  # PYCHOK signature
        '''Return this C{_Accumulator} as string.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Separator to join (C{str}).

           @return: Accumulator (C{str}).
        '''
        d = dict(n=self.num, s=self._s, t=self._t)
        return sep.join(pairs(d, prec=prec))


__all__ += _ALL_DOCS(GeodesicAreaExact, PolygonArea)

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
