
# -*- coding: utf-8 -*-

u'''(INTERNAL) Common base classes.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**,
see U{http://www.Movable-Type.co.UK/scripts/latlong.html}
and U{http://www.Movable-Type.co.UK/scripts/latlong-vectors.html}.

@newfield example: Example, Examples
'''
from dms import F_D, F_DMS, latDMS, lonDMS, parseDMS, parseDMS2
from fmath import EPS, favg, map1, scalar
from formy import antipode, compassAngle, equirectangular, \
                  haversine, isantipode
from utily import R_M, points2, property_RO, unStr

from math import asin, cos, degrees, radians, sin

# XXX the following classes are listed only to get
# Epydoc to include class and method documentation
__all__ = (  # 'Based', 'Named', 'VectorBased',
           'LatLonHeightBase',  # for documentation
           'classname', 'classnaming',
           'inStr')
__version__ = '18.10.26'

__X = object()  # unique instance


class Named(object):
    '''(INTERNAL) Base class for object with a name.
    '''
    _name        = ''     #: (INTERNAL) name (C{str})
    _classnaming = False  #: (INTERNAL) prefixed (C{bool})

    def _xcopy(self, *attrs):
        '''(INTERNAL) Must be overloaded.
        '''
        raise AssertionError(unStr(self.classname + '._xcopy', *attrs))

    @property_RO
    def classname(self):
        '''Get this object's C{[module.]class} name (C{str}), see L{classnaming}.
        '''
        return classname(self, prefixed=self._classnaming)

    @property
    def classnaming(self):
        '''Get the class naming (C{bool}).
        '''
        return self._classnaming

    @classnaming.setter  # PYCHOK setter!
    def classnaming(self, prefixed):
        '''Set the class naming for C{[module.].class} names.

           @param prefixed: Include the module name (C{bool}).
        '''
        self._classnaming = bool(prefixed)

    def copy(self):
        '''Make a copy of this instance.

           @return: The copy (C{This class} or subclass thereof).
        '''
        return self._xcopy()

    @property
    def name(self):
        '''Get the name (C{str}).
        '''
        return self._name

    @name.setter  # PYCHOK setter!
    def name(self, name):
        '''Set the name.

           @param name: New name (C{str}).
        '''
        self._name = str(name)


class Based(Named):
    '''(INTERNAL) Base class with name.
    '''

    def __repr__(self):
        return self.toStr2()

    def __str__(self):
        return self.toStr()

    def _update(self, unused):
        '''(INTERNAL) To be overloaded.
        '''
        pass

    def classof(self, *args, **kwds):
        '''Instantiate this very class.

           @param args: Optional, positional arguments.
           @keyword kwds: Optional, keyword arguments.

           @return: New instance (I{self.__class__}).
        '''
        return _xnamed(self.__class__(*args, **kwds), self.name)

#   def notImplemented(self, attr):
#       '''Raise error for a missing method, function or attribute.
#
#          @param attr: Attribute name (C{str}).
#
#          @raise NotImplementedError: No such attribute.
#       '''
#       c = self.__class__.__name__
#       return NotImplementedError('%s.%s' % (c, attr))

    def others(self, other, name='other'):
        '''Check this and an other instance for type compatiblility.

           @param other: The other instance (any C{type}).
           @keyword name: Optional, name for other (C{str}).

           @return: C{None}.

           @raise TypeError: Mismatch of this and I{other} C{type}.
        '''
        if not (isinstance(self, other.__class__) or
                isinstance(other, self.__class__)):
            raise TypeError('type(%s) mismatch: %s vs %s' % (name,
                             classname(other), self.classname))

    def toStr(self, **kwds):
        '''(INTERNAL) Must be overloaded.

           @param kwds: Optional, keyword arguments.
        '''
        raise AssertionError(unStr(self.classname + '.toStr', **kwds))

    def toStr2(self, **kwds):
        '''(INTERNAL) To be overloaded.

           @keyword kwds: Optional, keyword arguments.

           @return: L{toStr}() plus keyword arguments (as C{str}).
        '''
        t = self.toStr(**kwds).lstrip('([{').rstrip('}])')
        return '%s(%s)' % (self.classname, t)


class LatLonHeightBase(Based):
    '''(INTERNAL) Base class for C{LatLon} points on
       spherical or ellipsiodal earth models.
    '''
    _ab     = ()    #: (INTERNAL) Cache (lat, lon) radians (2-tuple)
    _datum  = None  #: (INTERNAL) Datum, overriden
    _height = 0     #: (INTERNAL) Height (C{meter})
    _lat    = 0     #: (INTERNAL) Latitude (C{degrees})
    _lon    = 0     #: (INTERNAL) Longitude (C{degrees})
    _name   = ''    #: (INTERNAL) name (C{str})

    def __init__(self, lat, lon, height=0, name=''):
        '''New C{LatLon}.

           @param lat: Latitude (C{degrees} or DMS C{str} with N or S suffix).
           @param lon: Longitude (C{degrees} or DMS C{str} with E or W suffix).
           @keyword height: Optional height (C{meter} above or below the earth surface).
           @keyword name: Optional name (C{str}).

           @return: New instance (C{LatLon}).

           @raise RangeError: Value of I{lat} or I{lon} outside the valid
                              range and I{rangerrrors} set to C{True}.

           @raise ValueError: Invalid I{lat} or I{lon}.

           @example:

           >>> p = LatLon(50.06632, -5.71475)
           >>> q = LatLon('50°03′59″N', """005°42'53"W""")
        '''
        self._lat, self._lon = parseDMS2(lat, lon)
        if height:  # elevation
            self._height = scalar(height, None, name='height')
        if name:
            self.name = name

    def __eq__(self, other):
        return self.isequalTo(other)

    def __ne__(self, other):
        return not self.isequalTo(other)

    def __str__(self):
        return self.toStr(form=F_D, prec=6)

    def _havg(self, other, f=0.5):
        '''(INTERNAL) Weighted, average height.

           @param other: An other point (C{LatLon}).
           @keyword f: Optional fraction (C{float}).

           @return: Average, fractional height (C{float}).
        '''
        return favg(self.height, other.height, f=f)

    def _update(self, updated):
        '''(INTERNAL) Reset caches if updated.
        '''
        if updated:  # reset caches
            self._ab = self._wm = None

    def _xcopy(self, *attrs):
        '''(INTERNAL) Make copy with add'l, subclass attributes.
        '''
        return _xattrs(self.classof(self.lat, self.lon),
                       self, '_height', *attrs)

    def antipode(self, height=None):
        '''Return the antipode, the point diametrically opposite
           to this point.

           @keyword height: Optional height of the antipode, height
                            of this point otherwise (C{meter}).

           @return: The antipodal point (C{LatLon}).
        '''
        a, b = antipode(self.lat, self.lon)
        h = self.height if height is None else height
        return self.classof(a, b, height=h)

    def bounds(self, wide, high, radius=R_M):
        '''Return the SE and NW lat-/longitude of a great circle
           bounding box centered at this location.

           @param wide: Longitudinal box width (C{meter}, same unts as I{radius}).
           @param high: Latitudinal box height (C{meter}, same unts as I{radius}).
           @keyword radius: Optional, mean earth radius (C{meter}).

           @return: 2-Tuple (LatLonSW, LatLonNE) of (C{LatLon}[]).

           @see: U{http://www.Movable-Type.co.UK/scripts/latlong-db.html}
        '''
        a, _ = self.to2ab()
        ca = cos(a)

        if ca > EPS:
            w = abs(degrees(asin(wide * 0.5 / radius) / ca))
        else:
            w = 0  # XXX
        h = abs(degrees(high * 0.5 / radius))

        return self.classof(self.lat - h, self.lon - w, height=self.height), \
               self.classof(self.lat + h, self.lon + w, height=self.height)

    def compassAngle(self, other):
        '''DEPRECATED, use method C{compassAngleTo}.
        '''
        return self.compassAngleTo(other)

    def compassAngleTo(self, other, adjust=True, wrap=False):
        '''Return the angle from North for the direction vector between
           this and an other point.

           Suitable only for short, non-near-polar vectors up to a few
           hundred Km or Miles.  Use method I{initialBearingTo} for
           larger distances.

           @param other: The other point (C{LatLon}).
           @keyword adjust: Adjust the longitudinal delta by the
                            cosine of the mean latitude (C{bool}).
           @keyword wrap: Wrap and L{unroll180} longitudes and longitudinal
                          delta (C{bool}).

           @return: Compass angle from North (C{degrees360}).

           @raise TypeError: The I{other} point is not C{LatLon}.

           @note: Courtesy Martin Schultz.

           @see: U{Local, flat earth approximation
                 <http://www.EdWilliams.org/avform.htm#flat>}.
        '''
        self.others(other)
        return compassAngle(self.lat, self.lon, other.lat, other.lon,
                            adjust=adjust, wrap=wrap)

    def copy(self):
        '''Copy this point.

           @return: The copy (C{LatLon} or subclass thereof).
        '''
        return self._xcopy()

    def equals(self, other, eps=None):
        '''DEPRECATED, use method C{isequalTo}.
        '''
        return self.isequalTo(other, eps=eps)

    def equals3(self, other, eps=None):
        '''DEPRECATED, use method C{isequalTo3}.
        '''
        return self.isequalTo3(other, eps=eps)

    def equirectangularTo(self, other, radius=None, **options):
        '''Compute the distance between this and an other point
           using the U{Equirectangular Approximation / Projection
           <http://www.Movable-Type.co.UK/scripts/latlong.html>}.

           Suitable only for short, non-near-polar distances up to a
           few hundred Km or Miles.  Use method C{haversineTo} or
           C{distanceTo*} for more accurate and/or larger distances.

           See function L{equirectangular_} for more details, the
           available I{options} and errors raised.

           @param other: The other point (C{LatLon}).
           @keyword radius: Optional, mean earth radius (C{meter}) or
                            C{None} for the mean radius of this
                            point's datum ellipsoid.
           @keyword options: Optional keyword arguments for function
                             L{equirectangular}.

           @return: Distance (C{meter}, same units as I{radius}).

           @raise TypeError: The I{other} point is not C{LatLon}.
        '''
        self.others(other)

        r = radius or (self._datum.ellipsoid.R1 if self._datum else R_M)
        return equirectangular(self.lat, self.lon, other.lat, other.lon,
                               radius=r, **options)

    def haversineTo(self, other, radius=None, wrap=False):
        '''Compute the distance between this and an other point using the
           U{Haversine<http://www.Movable-Type.co.UK/scripts/latlong.html>}
           formula.

           @param other: The other point (C{LatLon}).
           @keyword radius: Optional, mean earth radius (C{meter}) or
                            C{None} for the mean radius of this
                            point's datum ellipsoid.
           @keyword wrap: Wrap and L{unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as I{radius}).

           @raise TypeError: The I{other} point is not C{LatLon}.

           @see: Function L{haversine}, methods C{equirectangularTo}
                 and C{distanceTo*}.
        '''
        self.others(other)

        r = radius or (self._datum.ellipsoid.R1 if self._datum else R_M)
        return haversine(self.lat, self.lon, other.lat, other.lon,
                                             radius=r, wrap=wrap)

    @property
    def height(self):
        '''Get the height (C{meter}).
        '''
        return self._height

    @height.setter  # PYCHOK setter!
    def height(self, height):
        '''Set the height.

           @param height: New height (C{meter}).

           @raise TypeError: Invalid I{height} C{type}.

           @raise ValueError: Invalid I{height}.
        '''
        h = scalar(height, None, name='height')
        self._update(h != self._height)
        self._height = h

    def isantipodeTo(self, other, eps=EPS):
        '''Check whether this and an other point are antipodal,
           on diametrically opposite sides of the earth.

           @param other: The other point (C{LatLon}).
           @keyword eps: Tolerance for near-equality (C{degrees}).

           @return: C{True} if points are antipodal within the given
                    tolerance, C{False} otherwise.
        '''
        return isantipode(self.lat,  self.lon,
                         other.lat, other.lon, eps=eps)

    def isantipode(self, other, eps=EPS):
        '''DEPRECATED, use method C{isantipodeTo}.
        '''
        return self.isantipodeTo(other, eps=eps)

    def isequalTo(self, other, eps=None):
        '''Compare this point with an other point.

           @param other: The other point (C{LatLon}).
           @keyword eps: Tolerance for equality (C{degrees}).

           @return: C{True} if both points are identical,
                    I{ignoring} height, C{False} otherwise.

           @raise TypeError: The I{other} point is not C{LatLon}.

           @see: Method L{isequalTo3}.

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(52.205, 0.119)
           >>> e = p.isequalTo(q)  # True
        '''
        self.others(other)

        if eps and eps > 0:
            return abs(self.lat - other.lat) < eps and \
                   abs(self.lon - other.lon) < eps
        else:
            return self.lat == other.lat and \
                   self.lon == other.lon

    def isequalTo3(self, other, eps=None):
        '''Compare this point with an other point.

           @param other: The other point (C{LatLon}).
           @keyword eps: Tolerance for equality (C{degrees}).

           @return: C{True} if both points are identical
                    I{including} height, C{False} otherwise.

           @raise TypeError: The I{other} point is not C{LatLon}.

           @see: Method L{isequalTo}.

           @example:

           >>> p = LatLon(52.205, 0.119, 42)
           >>> q = LatLon(52.205, 0.119)
           >>> e = p.isequalTo3(q)  # False
        '''
        return self.isequalTo(other, eps=eps) and self.height == other.height

    @property
    def lat(self):
        '''Get the latitude (C{degrees90}).
        '''
        return self._lat

    @lat.setter  # PYCHOK setter!
    def lat(self, lat):
        '''Set the latitude.

           @param lat: New latitude (C{str[N|S]} or C{degrees}).

           @raise ValueError: Invalid I{lat}.
        '''
        lat = parseDMS(lat, suffix='NS', clip=90)
        self._update(lat != self._lat)
        self._lat = lat

    @property
    def latlon(self):
        '''Get the lat- and longitude (2-tuple of C{degrees90}, C{degrees180}).
        '''
        return self._lat, self._lon

    @latlon.setter  # PYCHOK setter!
    def latlon(self, latlonh):
        '''Set the lat- and longitude and optionally the height.

           @param latlonh: New lat-, longitude and height (2- or
                           3-tuple of C{degrees} and C{meter}).

           @raise TypeError: Height of I{latlonh} not C{scalar} or
                             I{latlonh} not C{list} or C{tuple}.

           @raise ValueError: Invalid I{latlonh} or M{len(latlonh)}.

           @see: Function L{parse3llh} to parse a I{latlonh} string
                 into a 3-tuple (lat, lon, h).
        '''
        if not isinstance(latlonh, (list, tuple)):
            raise TypeError('%s invalid: %r' % ('latlonh', latlonh))

        if len(latlonh) == 3:
            h = scalar(latlonh[2], None, name='latlonh')
        elif len(latlonh) != 2:
            raise ValueError('%s invalid: %r' % ('latlonh', latlonh))
        else:
            h = self._height

        lat, lon = parseDMS2(latlonh[0], latlonh[1])
        self._update(lat != self._lat or
                     lon != self._lon or h != self._height)
        self._lat, self._lon, self._height = lat, lon, h

    def latlon_(self, ndigits=0):
        '''DEPRECATED, use method C{latlon2}.
        '''
        return self.latlon2(ndigits)

    def latlon2(self, ndigits=0):
        '''Return this point's lat- and longitude, rounded.

           @keyword ndigits: Number of decimal digits (C{int}).

           @return: 2-Tuple (lat, lon) in (C{degrees90}, C{degrees180}),
                    rounded away from zero.

           @see: Built-in function C{round}.
        '''
        return round(self.lat, ndigits), round(self.lon, ndigits)

    def latlon2round(self, ndigits=0):
        '''DEPRECATED, use method C{latlon2}.
        '''
        return self.latlon2(ndigits)

    @property
    def lon(self):
        '''Get the longitude (C{degrees180}).
        '''
        return self._lon

    @lon.setter  # PYCHOK setter!
    def lon(self, lon):
        '''Set the longitude.

           @param lon: New longitude (C{str[E|W]} or C{degrees}).

           @raise ValueError: Invalid I{lon}.
        '''
        lon = parseDMS(lon, suffix='EW', clip=180)
        self._update(lon != self._lon)
        self._lon = lon

    def points(self, points, closed=True):
        '''DEPRECATED, use method C{points2}.
        '''
        return self.points2(points, closed=closed)

    def points2(self, points, closed=True):
        '''Check a polygon represented by points.

           @param points: The polygon points (C{LatLon}[])
           @keyword closed: Optionally, consider the polygon closed,
                            ignoring any duplicate or closing final
                            I{points} (C{bool}).

           @return: 2-Tuple (number, ...) of points (C{int}, C{list} or
                    C{tuple}).

           @raise TypeError: Some I{points} are not C{LatLon}.

           @raise ValueError: Insufficient number of I{points}.
        '''
        return points2(points, closed=closed, base=self)

    def to2ab(self):
        '''Return this point's lat- and longitude in radians.

           @return: 2-Tuple (lat, lon) in (C{radiansPI_2}, C{radiansPI}).
        '''
        if not self._ab:
            self._ab = map1(radians, self.lat, self.lon)
        return self._ab

    def to3llh(self):
        '''Return this point's lat-, longitude and height.

           @return: 3-Tuple (lat, lon, h) in (C{degrees90}, C{degrees180}, C{meter}).
        '''
        return self.lat, self.lon, self.height

    def to3xyz(self):
        '''Convert this (geodetic) point to (n-)vector (normal
           to the earth's surface) x/y/z components, ignoring
           the height.

           @return: 3-Tuple (x, y, z) in (units, NOT meter).
        '''
        # Kenneth Gade eqn 3, but using right-handed
        # vector x -> 0°E,0°N, y -> 90°E,0°N, z -> 90°N
        a, b = self.to2ab()
        ca = cos(a)
        return ca * cos(b), ca * sin(b), sin(a)

    def toStr(self, form=F_DMS, prec=None, m='m', sep=', '):  # PYCHOK expected
        '''Convert this point to a "lat, lon [+/-height]" string,
           formatted in the given form.

           @keyword form: Optional format, F_D, F_DM, F_DMS for
                          deg°, deg°min′, deg°min′sec″ (C{str}).
           @keyword prec: Optional number of decimal digits (0..8 or C{None}).
           @keyword m: Optional unit of the height (C{str}).
           @keyword sep: Optional separator to join (C{str}).

           @return: Point in the specified form (C{str}).

           @example:

           >>> LatLon(51.4778, -0.0016).toStr()  # 51°28′40″N, 000°00′06″W
           >>> LatLon(51.4778, -0.0016).toStr(F_D)  # 51.4778°N, 000.0016°W
           >>> LatLon(51.4778, -0.0016, 42).toStr()  # 51°28′40″N, 000°00′06″W, +42.00m

        '''
        t = [latDMS(self.lat, form=form, prec=prec),
             lonDMS(self.lon, form=form, prec=prec)]
        if self.height:
            t += ['%+.2f%s' % (self.height, m)]
        return sep.join(t)


class VectorBased(Based):
    '''(INTERNAL) Base class for I{Vector3d}.
    '''
    def __init__(self, name='', **unused):
        if name:
            self.name = name


def _nameof(inst):
    '''(INTERNAL) Get the instance' name or C{''}.
    '''
    try:
        return inst.name
    except AttributeError:
        return ''


def _xattrs(inst, other, *attrs):
    '''(INTERNAL) Copy attribute values from I{other} to I{inst}.

       @param inst: Instance to copy atrribute values to.
       @param other: Instance to copy atrribute values from.
       @param attrs: Attribute names (C{str})s.

       @return: The I{inst}, updated.

       @raise AttributeError: If attribute doesn't exist or is not settable.
    '''
    for a in attrs:
        g = getattr(inst, a, __X)
        if g is __X:
            raise AttributeError('invalid %r.%s' % (inst, a))
        s = getattr(other, a, __X)
        if s is __X:
            raise AttributeError('invalid %r.%s' % (other, a))
        elif s != g:
            setattr(inst, a, s)  # not settable?
    return inst


def _xnamed(inst, name):
    '''(INTERNAL) Set the instance' C{.name = }C{name}.

       @param name: The name (C{str}).

       @return: The I{inst}, named if unnamed before.
    '''
    if name and isinstance(inst, Named):
        try:
            if not inst.name:
                inst.name = name
        except AttributeError:
            pass
    return inst


def classname(inst, prefixed=None):
    '''Return an instance' module and class name.

       @param inst: The object (any type).
       @keyword prefixed: Prefix the module name (C{bool}), see
                          function L{classnaming}.

       @return: The I{inst}'s C{[module.]class} name (C{str}).
    '''
    try:
        n = inst.__class__.__name__
    except AttributeError:
        n = 'Nn'
    if prefixed or (getattr(inst, 'classnaming', Named._classnaming)
                    if prefixed is None else False):
        try:
            m = inst.__module__
            n = '.'.join(m.split('.')[-1:] + [n])
        except AttributeError:
            pass
    return n


def classnaming(prefixed=None):
    '''Set the default class naming for C{[module.]class} names.

       @keyword prefixed: Include the module name (C{bool}).

       @return: Previous class naming setting (C{bool}).
    '''
    t = Named._classnaming
    if prefixed in (True, False):
        Named._classnaming = bool(prefixed)
    return t


def inStr(inst, *args, **kwds):
    '''Return the string representation of an instance.

       @param inst: The instance (any type).
       @param args: Optional positional arguments.
       @keyword kwds: Optional keyword arguments.

       @return: The I{inst}'s representation (C{str}).
    '''
    return unStr(classname(inst), *args, **kwds)

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
