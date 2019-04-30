
# -*- coding: utf-8 -*-

u'''(INTERNAL) Ellipsoidal base classes.

Pure Python implementation of geodesy tools for ellipsoidal earth models,
transcribed in part from JavaScript originals by I{(C) Chris Veness 2005-2016}
and published under the same MIT Licence**, see for example U{latlon-ellipsoidal
<http://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}.

@newfield example: Example, Examples
'''

from bases import LatLonHeightBase, _nameof, _xnamed
from datum import Datum, Datums
from dms import parseDMS2
from fmath import EPS, EPS1, fsum_, hypot, hypot1, isscalar
from utily import degrees90, degrees180, _for_docs, property_RO, \
                  sincos2, _Strs, wrap90, wrap360
from vector3d import Vector3d

from math import atan2, copysign, sqrt

__all__ = _for_docs('CartesianBase', 'LatLonEllipsoidalBase')
__version__ = '19.04.24'

# some INTERNAL constants, to avoid circular ups, utm, utmups imports
_UTM_LAT_MAX      =  84  # PYCHOK for export (C{degrees})
_UTM_LAT_MIN      = -80  # PYCHOK for export (C{degrees})
_UTM_ZONE_MAX     =  60  # PYCHOK for export
_UTM_ZONE_MIN     =   1  # PYCHOK for export
_UTM_ZONE_OFF_MAX =  60  # PYCHOK max Central meridian offset (C{degrees})

_UPS_LAT_MAX  = _UTM_LAT_MAX - 0.5     # PYCHOK includes 30' UTM overlap
_UPS_LAT_MIN  = _UTM_LAT_MIN + 0.5     # PYCHOK includes 30' UTM overlap
_UPS_ZONE     = _UTM_ZONE_MIN - 1      # PYCHOK for export
_UPS_ZONE_STR = '%02d' % (_UPS_ZONE,)  # PYCHOK for export

_UTMUPS_ZONE_INVALID = -4             # PYCHOK for export too
_UTMUPS_ZONE_MIN     = _UPS_ZONE      # PYCHOK for export too
_UTMUPS_ZONE_MAX     = _UTM_ZONE_MAX  # PYCHOK for export too

# _MAX_PSEUDO_ZONE      = -1
# _MIN_PSEUDO_ZONE      = -4
# _UTMUPS_ZONE_MATCH    = -3
# _UTMUPS_ZONE_STANDARD = -1
# _UTM                  = -2


def _hemi(lat):  # imported by .ups, .utm
    '''Return the hemisphere letter.

       @param lat: Latitude (C{degrees} or C{radians}).

       @return: C{'N'|'S'} for north-/southern hemisphere.
    '''
    return 'S' if lat < 0 else 'N'


def _to4lldn(latlon, lon, datum, name):
    '''(INTERNAL) Return 4-tuple (C{lat, lon, datum, name}).
    '''
    try:
        # if lon is not None:
        #     raise AttributeError
        lat, lon = latlon.lat, latlon.lon
        if not isinstance(latlon, LatLonEllipsoidalBase):
            raise TypeError('%s not %s: %r' % ('latlon', 'ellipsoidal', latlon))
        d = datum or latlon.datum
    except AttributeError:
        lat, lon = parseDMS2(latlon, lon)
        d = datum or Datums.WGS84
    return lat, lon, d, (name or _nameof(latlon))


def _to3zBhp(zone, band, hemipole=''):  # imported by .epsg, .ups, .utm, .utmups
    '''Parse UTM/UPS zone, Band letter and hemisphere/pole letter.

       @param zone: Zone with/-out Band (C{scalar} or C{str}).
       @keyword band: Optional (longitudinal/polar) Band letter (C{str}).
       @keyword hemipole: Optional hemisphere/pole letter (C{str}).

       @return: 3-Tuple (C{zone, Band, hemisphere/pole}) as (C{int,
                str, 'N'|'S'}) where C{zone} is C{0} for UPS or
                C{1..60} for UTM and C{Band} is C{'A'..'Z'} I{NOT}
                checked for valid UTM/UPS bands.

       @raise ValueError: Invalid I{zone}, I{band} or I{hemipole}.
    '''
    B = band
    try:
        z = _UTMUPS_ZONE_INVALID
        if isscalar(zone) or zone.isdigit():
            z = int(zone)
        elif zone and isinstance(zone, _Strs):
            if len(zone) > 1:
                B = zone[-1:]
                z = int(zone[:-1])
            elif zone in 'AaBbYyZz':  # single letter
                B = zone
                z = _UPS_ZONE

        if _UTMUPS_ZONE_MIN <= z <= _UTMUPS_ZONE_MAX:
            hp = hemipole[:1].upper()
            if hp in ('N', 'S') or not hp:
                B = B.upper()
                if B.isalpha():
                    return z, B, (hp or ('S' if B < 'N' else 'N'))
                elif not B:
                    return z, B, hp

    except (AttributeError, TypeError, ValueError):
        pass
    raise ValueError('%s, %s or %s invalid: %r' %
                     ('zone', 'band', 'hemipole', (zone, B, hemipole)))


def _to3zll(lat, lon):  # imported by .ups, .utm
    '''Wrap lat- and longitude and determine UTM zone.

       @param lat: Latitude (C{degrees}).
       @param lon: Longitude (C{degrees}).

       @return: 3-Tuple (C{zone, lat, lon}) as (C{int}, C{degrees90},
                C{degrees180}) where C{zone} is C{1..60} for UTM.
    '''
    x = wrap360(lon + 180)  # use wrap360 to get ...
    z = int(x) // 6 + 1  # ... longitudinal UTM zone [1, 60] and ...
    lon = x - 180  # ... lon [-180, 180) i.e. -180 <= lon < 180
    return z, wrap90(lat), lon


class CartesianBase(Vector3d):
    '''(INTERNAL) Base class for ellipsoidal I{Cartesian}.
    '''

    def _applyHelmert(self, transform, inverse=False):
        '''(INTERNAL) Return a new (geocentric) Cartesian point
           by applying a Helmert transform to this point.

           @param transform: Transform to apply (L{Transform}).
           @keyword inverse: Optionally, apply the inverse of
                             Helmert transform (C{bool}).

           @return: The transformed point (L{Cartesian}).
        '''
        x, y, z = self.to3xyz()
        return self.classof(*transform.transform(x, y, z, inverse))

    def to3llh(self, datum=Datums.WGS84):
        '''Convert this (geocentric) Cartesian (x/y/z) point to
           (ellipsoidal, geodetic) lat-, longitude and height on
           the given datum.

           Uses Bowring’s (1985) formulation for μm precision in concise
           form: U{'The accuracy of geodetic latitude and height equations'
           <http://www.ResearchGate.net/publication/
           233668213_The_Accuracy_of_Geodetic_Latitude_and_Height_Equations>},
           B. R. Bowring, Survey Review, Vol 28, 218, Oct 1985.

           See also Ralph M. Toms U{'An Efficient Algorithm for Geocentric to
           Geodetic Coordinate Conversion'<http://www.OSTI.gov/scitech/biblio/110235>},
           Sept 1995 and U{'An Improved Algorithm for Geocentric to Geodetic Coordinate
           Conversion'<http://www.OSTI.gov/scitech/servlets/purl/231228>},
           Apr 1996, from Lawrence Livermore National Laboratory.

           @keyword datum: Optional datum to use (L{Datum}).

           @return: 3-Tuple (lat, lon, heigth) in (C{degrees90},
                    C{degrees180}, C{meter}).
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
            sa, ca = sincos2(a)
#           r = E.a / E.e2s(sa)  # length of normal terminated by minor axis
#           h = p * ca + z * sa - (E.a * E.a / r)
            h = fsum_(p * ca, z * sa, -E.a * E.e2s(sa))

            a, b = degrees90(a), degrees180(b)

        # see <http://GIS.StackExchange.com/questions/28446>
        elif p > EPS:  # latitude arbitrarily zero
            a, b, h = 0.0, degrees180(atan2(y, x)), p - E.a
        else:  # polar latitude, longitude arbitrarily zero
            a, b, h = copysign(90.0, z), 0.0, abs(z) - E.b

        return a, b, h

    def _toLLhd(self, LL, datum):
        '''(INTERNAL) Helper for I{subclass.toLatLon}.
        '''
        a, b, h = self.to3llh(datum)
        return (a, b, h) if LL is None else _xnamed(LL(
                a, b, height=h, datum=datum), self.name)

    def toStr(self, prec=3, fmt='[%s]', sep=', '):  # PYCHOK expected
        '''Return the string representation of this cartesian.

           @keyword prec: Optional number of decimals, unstripped (C{int}).
           @keyword fmt: Optional enclosing backets format (string).
           @keyword sep: Optional separator to join (string).

           @return: Cartesian represented as "[x, y, z]" (string).
        '''
        return Vector3d.toStr(self, prec=prec, fmt=fmt, sep=sep)


class LatLonEllipsoidalBase(LatLonHeightBase):
    '''(INTERNAL) Base class for ellipsoidal C{LatLon}.
    '''
    _convergence  = None  #: (INTERNAL) UTM/UPS meridian convergence (C{degrees}).
    _datum        = Datums.WGS84  #: (INTERNAL) Datum (L{Datum}).
    _elevation2   = ()    #: (INTERNAL) cached C{elevation2} result.
    _geoidHeight2 = ()    #: (INTERNAL) cached C{geoidHeight2} result.
    _osgr         = None  #: (INTERNAL) cache toOsgr (C{Osgr}).
    _scale        = None  #: (INTERNAL) UTM/UPS scale factor (C{float}).
    _ups          = None  #: (INTERNAL) cache toUps (L{Ups}).
    _utm          = None  #: (INTERNAL) cache toUtm (L{Utm}).
    _wm           = None  #: (INTERNAL) cache toWm (webmercator.Wm instance).

    def __init__(self, lat, lon, height=0, datum=None, name=''):
        '''Create an (ellipsoidal) C{LatLon} point frome the given
           lat-, longitude and height on the given datum.

           @param lat: Latitude (C{degrees} or DMS C{[N|S]}).
           @param lon: Longitude (C{degrees} or DMS C{str[E|W]}).
           @keyword height: Optional elevation (C{meter}, the same units
                            as the datum's half-axes).
           @keyword datum: Optional datum to use (L{Datum}).
           @keyword name: Optional name (string).

           @example:

           >>> p = LatLon(51.4778, -0.0016)  # height=0, datum=Datums.WGS84
        '''
        LatLonHeightBase.__init__(self, lat, lon, height=height, name=name)
        if datum:  # check datum
            self.datum = datum

    def _Radjust2(self, adjust, datum, meter, model):
        '''(INTERNAL) Adjust elevation or geoidHeight with diference
           in Gaussian radii of curvature of given datum and NAD83.
           This is an arbitrary, possibly incorrect adjustment.
        '''
        if adjust and isinstance(meter, float):
            n = Datums.NAD83.ellipsoid.rocGauss(self.lat)
            if min(abs(meter), n) > EPS:
                # use ratio, datum and NAD83 units may differ
                e = self.ellipsoid(datum).rocGauss(self.lat)
                if min(abs(e - n), e) > EPS:
                    meter *= e / n
        return meter, model

    def _update(self, updated):
        if updated:  # reset caches
            self._osgr = self._utm = self._wm = None
            self._elevation2 = self._geoidHeight2 = ()
            LatLonHeightBase._update(self, updated)

    def _xcopy(self, *attrs):
        '''(INTERNAL) Make copy with add'l, subclass attributes.
        '''
        return LatLonHeightBase._xcopy(self, '_datum', *attrs)

    def antipode(self, height=None):
        '''Return the antipode, the point diametrically opposite
           to this point.

           @keyword height: Optional height of the antipode, height
                            of this point otherwise (C{meter}).

           @return: The antipodal point (C{LatLon}).
        '''
        lla = LatLonHeightBase.antipode(self, height=height)
        if lla.datum != self.datum:
            lla.datum = self.datum
        return lla

    @property_RO
    def convergence(self):
        '''Get this point's UTM or UPS meridian convergence (C{degrees})
           or C{None} if not converted from L{Utm} ot L{Ups}.
        '''
        return self._convergence

    def convertDatum(self, datum):
        '''Convert this C{LatLon} point to a new datum.

           @param datum: Datum to convert to (L{Datum}).

           @return: The converted point (ellipsoidal C{LatLon}).

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

    @property
    def datum(self):
        '''Get this point's datum (L{Datum}).
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set this point's datum without conversion.

           @param datum: New datum (L{Datum}).

           @raise TypeError: The I{datum} is not a L{Datum}.

           @raise ValueError: The I{datum} is not ellipsoidal.
        '''
        if not isinstance(datum, Datum):
            raise TypeError('%r not a %s: %r' % ('datum', Datum.__name__, datum))
        if not datum.isEllipsoidal:
            raise ValueError('%r not %s: %r' % ('datum', 'ellipsoidal', datum))
        self._update(datum != self._datum)
        self._datum = datum

    def distanceTo2(self, other):
        '''Approximate the distance and bearing between this and an
           other (ellipsoidal) point based on the radii of curvature.

           Suitable only for short distances up to a few hundred Km
           or Miles and only between non-near-polar points.

           @param other: The other point (C{LatLon}).

           @return: 2-Tuple (distance, bearing) in (C{meter}, C{degrees360}).

           @raise TypeError: The I{other} point is not C{LatLon}.

           @raise ValueError: Incompatible datum ellipsoids.

           @see: Method L{Ellipsoid.distance2} and U{Local, flat earth
                 approximation<http://www.EdWilliams.org/avform.htm#flat>}.
        '''
        return self.ellipsoids(other).distance2(self.lat,  self.lon,
                                               other.lat, other.lon)

    def elevation2(self, adjust=True, datum=Datums.WGS84, timeout=2):
        '''Return elevation of this point for its or the given datum.

           @keyword adjust: Adjust the elevation for a I{datum} other
                            than C{NAD83}.
           @keyword datum: Optional datum (L{Datum}).
           @keyword timeout: Optional query timeout (seconds).

           @return: 2-Tuple (elevation, data_source) in (C{meter}, C{str})
                    or in case of errors (C{None}, I{<error>}).

           @note: The adjustment applied is the difference in geocentric
                  earth radius for the I{datum} used and the C{NAV83}
                  datum upon which L{elevations.elevation2} is based.

           @note: NED elevation is only available for locations in the
                  U{Conterminous US (CONUS)
                  <http://WikiPedia.org/wiki/Contiguous_United_States>}.

           @see: Function L{elevations.elevation2} and method
                 L{Ellipsoid.Rgeocentric} for further details.
        '''
        if not self._elevation2:  # get elevation and data source
            from elevations import elevation2
            self._elevation2 = elevation2(self.lat, self.lon,
                                          timeout=timeout)
        return self._Radjust2(adjust, datum, *self._elevation2)

    def ellipsoid(self, datum=Datums.WGS84):
        '''Return the ellipsoid of this point's datum or the given datum.

           @keyword datum: Optional datum (L{Datum}).

           @return: The ellipsoid (L{Ellipsoid}).
        '''
        return getattr(self, 'datum', datum).ellipsoid

    def ellipsoids(self, other):
        '''Check the type and ellipsoid of this and an other point's datum.

           @param other: The other point (C{LatLon}).

           @return: This point's datum ellipsoid (L{Ellipsoid}).

           @raise TypeError: The I{other} point is not C{LatLon}.

           @raise ValueError: Incompatible datum ellipsoids.
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

    def geoidHeight2(self, adjust=False, datum=Datums.WGS84, timeout=2):
        '''Return geoid height of this point for its or the given datum.

           @keyword adjust: Adjust the geoid height for a I{datum} other
                            than C{NAD83/NADV88}.
           @keyword datum: Optional datum (L{Datum}).
           @keyword timeout: Optional query timeout (seconds).

           @return: 2-Tuple (height, model_name) in (C{meter}, C{str})
                    or in case of errors (C{None}, I{<error>}).

           @note: The adjustment applied is the difference in geocentric
                  earth radius for the given I{datum} and the C{NAV83/NADV88}
                  datum upon which L{elevations.geoidHeight2} is based.

           @note: Ths geoid height is only available for locations in
                  the U{Conterminous US (CONUS)
                  <http://WikiPedia.org/wiki/Contiguous_United_States>}.

           @see: Function L{elevations.geoidHeight2} and method
                 L{Ellipsoid.Rgeocentric} for further details.
        '''
        if not self._geoidHeight2:  # get elevation and data source
            from elevations import geoidHeight2
            self._geoidHeight2 = geoidHeight2(self.lat, self.lon,
                                              model=0, timeout=timeout)
        return self._Radjust2(adjust, datum, *self._geoidHeight2)

    @property_RO
    def isEllipsoidal(self):
        '''Check whether this C{LatLon} point is ellipsoidal (C{bool}).
        '''
        return self.datum.isEllipsoidal

    @property_RO
    def isSpherical(self):
        '''Check whether this C{LatLon} point is spherical (C{bool}).
        '''
        return self.datum.isSpherical

    def parse(self, strll, height=0, datum=None, sep=','):
        '''Parse a string representing this C{LatLon} point.

           The lat- and longitude must be separated by a sep[arator]
           character.  If height is present it must follow and be
           separated by another sep[arator].  Lat- and longitude
           may be swapped, provided at least one ends with the
           proper compass direction.

           For more details, see functions L{parse3llh} and L{parseDMS}
           in sub-module L{dms}.

           @param strll: Lat, lon [, height] (string).
           @keyword height: Optional, default height (C{meter} or C{None}).
           @keyword datum: Optional, default datum (L{Datum}).
           @keyword sep: Optional separator (string).

           @return: The point (L{LatLonEllipsoidalBase}).

           @raise ValueError: Invalid I{strll}.
        '''
        from dms import parse3llh
        a, b, h = parse3llh(strll, height=height, sep=sep)
        return self.classof(a, b, height=h, datum=datum or self.datum)

    @property_RO
    def scale(self):
        '''Get this point's UTM grid or UPS point scale factor (C{float})
           or C{None} if not converted from L{Utm} or L{Ups}.
        '''
        return self._scale

    def to3xyz(self):  # overloads _LatLonHeightBase.to3xyz
        '''Convert this (ellipsoidal) geodetic C{LatLon} point to
           (geocentric) cartesian x/y/z components.

           @return: 3-Tuple (x, y, z) in (C{meter}).
        '''
        a, b = self.to2ab()
        sa, ca, sb, cb = sincos2(a, b)

        E = self.ellipsoid()
        # radius of curvature in prime vertical
        t = E.e2s2(sa)
        if t > EPS1:
            r = E.a
        elif t > EPS:
            r = E.a / sqrt(t)
        else:
            r = 0

        h = self.height
        t = (h + r) * ca
        return (t * cb,
                t * sb,
               (h + r * E.e12) * sa)

    def toOsgr(self):
        '''Convert this C{LatLon} point to an OSGR coordinate.

           @see: Function L{toOsgr} in module L{osgr}.

           @return: The OSGR coordinate (L{Osgr}).
        '''
        if self._osgr is None:
            from osgr import toOsgr  # PYCHOK recursive import
            self._osgr = toOsgr(self, datum=self.datum)
            self._osgr._latlon = self
        return self._osgr

    def toUps(self, pole='N', falsed=True):
        '''Convert this C{LatLon} point to a UPS coordinate.

           @keyword pole: Optional top/center of (stereographic)
                          projection (C{str}, 'N[orth]' or 'S[outh]').
           @keyword falsed: False easting and northing (C{bool}).

           @return: The UPS coordinate (L{Ups}).

           @see: Function L{toUps8}.
        '''
        if self._ups is None:
            from ups import toUps8, Ups  # PYCHOK recursive import
            self._ups = toUps8(self, datum=self.datum, Ups=Ups,
                                     pole=pole, falsed=falsed)
            self._ups._latlon = self
        return self._ups

    def toUtm(self):
        '''Convert this C{LatLon} point to a UTM coordinate.

           @return: The UTM coordinate (L{Utm}).

           @see: Function L{toUtm8}.
        '''
        if self._utm is None:
            from utm import toUtm8, Utm  # PYCHOK recursive import
            self._utm = toUtm8(self, datum=self.datum, Utm=Utm)
            self._utm._latlon = self
        return self._utm

    def toUtmUps(self, pole=''):
        '''Convert this C{LatLon} point to a UTM or UPS coordinate.

           @keyword pole: Optional top/center of UPS (stereographic)
                          projection (C{str}, 'N[orth]' or 'S[outh]').

           @return: The UTM or UPS coordinate (L{Utm} or L{Ups}).

           @raise TypeError: Result in L{Utm} or L{Ups}.

           @see: Function L{toUtmUps}.
        '''
        if self._utm:
            u = self._utm
        elif self._ups and (self._utm.pole == pole or not pole):
            u = self._ups
        else:
            from utmups import toUtmUps8, Utm, Ups  # PYCHOK recursive import
            u = toUtmUps8(self, datum=self.datum, Utm=Utm, Ups=Ups, pole=pole)
            u._latlon = self
            if isinstance(u, Utm):
                self._utm = u
            elif isinstance(u, Ups):
                self._ups = u
            else:
                raise TypeError('%s: %r' % ('toUtmUps8', u))
        return u

    def toWm(self):
        '''Convert this C{LatLon} point to a WM coordinate.

           @see: Function L{toWm} in module L{webmercator}.

           @return: The WM coordinate (L{Wm}).
        '''
        if self._wm is None:
            from webmercator import toWm  # PYCHOK recursive import
            self._wm = toWm(self)
        return self._wm

# **) MIT License
#
# Copyright (C) 2016-2019 -- mrJean1 at Gmail dot com
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
