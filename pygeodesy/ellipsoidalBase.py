
# -*- coding: utf-8 -*-

u'''(INTERNAL) Ellipsoidal base classes C{CartesianEllipsoidalBase} and
C{LatLonEllipsoidalBase}.

Pure Python implementation of geodesy tools for ellipsoidal earth models,
transcribed in part from JavaScript originals by I{(C) Chris Veness 2005-2016}
and published under the same MIT Licence**, see for example U{latlon-ellipsoidal
<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}.

@newfield example: Example, Examples
'''

from pygeodesy.basics import EPS, _isnotError, property_doc_, \
                             property_RO, _TypeError
from pygeodesy.cartesianBase import CartesianBase
from pygeodesy.datum import Datum, Datums
from pygeodesy.ecef import EcefVeness
from pygeodesy.latlonBase import LatLonBase
from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.named import Vector3Tuple
from pygeodesy.trf import _2epoch, RefFrame, TRFError, _reframeTransforms

__all__ = _ALL_DOCS('CartesianEllipsoidalBase', 'LatLonEllipsoidalBase')
__version__ = '20.03.23'


class CartesianEllipsoidalBase(CartesianBase):
    '''(INTERNAL) Base class for ellipsoidal C{Cartesian}s.
    '''
    _datum = Datums.WGS84  #: (INTERNAL) L{Datum}.
    _Ecef  = EcefVeness    #: (INTERNAL) Preferred C{Ecef...} class, backward compatible.

    def convertRefFrame(self, reframe2, reframe, epoch=None):
        '''Convert this cartesian point from one to an other reference frame.

           @arg reframe2: Reference frame to convert I{to} (L{RefFrame}).
           @arg reframe: Reference frame to convert I{from} (L{RefFrame}).
           @kwarg epoch: Optional epoch to observe for B{C{reframe}}, a
                         fractional calendar year (C{scalar}).

           @return: The converted point (C{Cartesian}) or this point if
                    conversion is C{nil}.

           @raise TRFError: No conversion available from B{C{reframe}}
                            to B{C{reframe2}}.

           @raise TypeError: B{C{reframe2}} or B{C{reframe}} not a
                             L{RefFrame} or B{C{epoch}} not C{scalar}.
        '''
        _TypeError(RefFrame, reframe2=reframe2, reframe=reframe)

        c, d = self, self.datum
        for t in _reframeTransforms(reframe2, reframe, reframe.epoch if
                                    epoch is None else _2epoch(epoch)):
            c = c._applyHelmert(t, False, datum=d)
        return c


class LatLonEllipsoidalBase(LatLonBase):
    '''(INTERNAL) Base class for ellipsoidal C{LatLon}s.
    '''
    _convergence  = None  #: (INTERNAL) UTM/UPS meridian convergence (C{degrees}).
    _datum        = Datums.WGS84  #: (INTERNAL) Datum (L{Datum}).
    _elevation2   = ()    #: (INTERNAL) Cached C{elevation2} result.
    _epoch        = None  #: (INTERNAL) overriding .reframe.epoch (C{float}).
    _etm          = None  #: (INTERNAL) Cached toEtm (L{Etm}).
    _geoidHeight2 = ()    #: (INTERNAL) Cached C{geoidHeight2} result.
    _lcc          = None  #: (INTERNAL) Cached toLcc (C{Lcc}).
    _osgr         = None  #: (INTERNAL) Cached toOsgr (C{Osgr}).
    _reframe      = None  #: (INTERNAL) reference frame (L{RefFrame}).
    _scale        = None  #: (INTERNAL) UTM/UPS scale factor (C{float}).
    _ups          = None  #: (INTERNAL) Cached toUps (L{Ups}).
    _utm          = None  #: (INTERNAL) Cached toUtm (L{Utm}).
    _wm           = None  #: (INTERNAL) Cached toWm (webmercator.Wm instance).

    def __init__(self, lat, lon, height=0, datum=None, reframe=None,
                                           epoch=None, name=''):
        '''Create an ellipsoidal C{LatLon} point frome the given
           lat-, longitude and height on the given datum and with
           the given reference frame and epoch.

           @arg lat: Latitude (C{degrees} or DMS C{[N|S]}).
           @arg lon: Longitude (C{degrees} or DMS C{str[E|W]}).
           @kwarg height: Optional elevation (C{meter}, the same units
                          as the datum's half-axes).
           @kwarg datum: Optional, ellipsoidal datum to use (L{Datum}).
           @kwarg reframe: Optional reference frame (L{RefFrame}).
           @kwarg epoch: Optional epoch to observe for B{C{reframe}}
                         (C{scalar}), a non-zero, fractional calendar year.
           @kwarg name: Optional name (string).

           @raise TypeError: B{C{datum}} is not a L{datum}, B{C{reframe}}
                             is not a L{RefFrame} or B{C{epoch}} is not
                             C{scalar} non-zero.

           @example:

           >>> p = LatLon(51.4778, -0.0016)  # height=0, datum=Datums.WGS84
        '''
        LatLonBase.__init__(self, lat, lon, height=height, name=name)
        if datum:
            self.datum = datum
        if reframe:
            self.reframe = reframe
            self.epoch = epoch

    def _Radjust2(self, adjust, datum, meter_text2):
        '''(INTERNAL) Adjust elevation or geoidHeight with difference
           in Gaussian radii of curvature of given datum and NAD83.

           @note: This is an arbitrary, possibly incorrect adjustment.
        '''
        if adjust:  # Elevation2Tuple or GeoidHeight2Tuple
            m, t = meter_text2
            if isinstance(m, float):
                n = Datums.NAD83.ellipsoid.rocGauss(self.lat)
                if min(abs(m), n) > EPS:
                    # use ratio, datum and NAD83 units may differ
                    e = self.ellipsoid(datum).rocGauss(self.lat)
                    if min(abs(e - n), e) > EPS:
                        m *= e / n
                        meter_text2 = meter_text2.classof(m, t)
        return self._xnamed(meter_text2)

    def _update(self, updated, *attrs):
        '''(INTERNAL) Zap cached attributes if updated.
        '''
        if updated:
            LatLonBase._update(self, updated, '_etm', '_lcc', '_osgr',
                                      '_ups', '_utm', '_wm', *attrs)
            if self._elevation2:
                self._elevation2 = ()
            if self._geoidHeight2:
                self._geoidHeight2 = ()

    def antipode(self, height=None):
        '''Return the antipode, the point diametrically opposite
           to this point.

           @kwarg height: Optional height of the antipode, height
                          of this point otherwise (C{meter}).

           @return: The antipodal point (C{LatLon}).
        '''
        lla = LatLonBase.antipode(self, height=height)
        if lla.datum != self.datum:
            lla.datum = self.datum
        return lla

    @property_RO
    def convergence(self):
        '''Get this point's UTM or UPS meridian convergence (C{degrees})
           or C{None} if not converted from L{Utm} or L{Ups}.
        '''
        return self._convergence

    def convertDatum(self, datum2):
        '''Convert this point to an other datum.

           @arg datum2: Datum to convert I{to} (L{Datum}).

           @return: The converted point (ellipsoidal C{LatLon}).

           @raise TypeError: The B{C{datum2}} is not a L{Datum}.

           @example:

           >>> p = LatLon(51.4778, -0.0016)  # default Datums.WGS84
           >>> p.convertDatum(Datums.OSGB36)  # 51.477284°N, 000.00002°E
        '''
        if self.datum == datum2:
            return self.copy()

        c = self.toCartesian().convertDatum(datum2)
        return c.toLatLon(datum=datum2, LatLon=self.classof)

    def convertRefFrame(self, reframe2):
        '''Convert this point to an other reference frame.

           @arg reframe2: Reference frame to convert I{to} (L{RefFrame}).

           @return: The converted point (ellipsoidal C{LatLon}) or
                    this point if conversion is C{nil}.

           @raise TRFError: No B{C{.reframe}} or no conversion
                            available from B{C{.reframe}} to
                            B{C{reframe2}}.

           @raise TypeError: The B{C{reframe2}} is not a L{RefFrame}.

           @example:

           >>> p = LatLon(51.4778, -0.0016, reframe=RefFrames.ETRF2000)  # default Datums.WGS84
           >>> p.convertRefFrame(RefFrames.ITRF2014)  # 51.477803°N, 000.001597°W, +0.01m
        '''
        _TypeError(RefFrame, reframe2=reframe2)

        if not self.reframe:
            raise TRFError('no %r.%s' % (self, 'reframe'))

        ts = _reframeTransforms(reframe2, self.reframe, self.epoch)
        if ts:
            c = self.toCartesian()
            for t in ts:
                c = c._applyHelmert(t, False)
            ll = c.toLatLon(datum=self.datum, LatLon=self.classof,
                            epoch=self.epoch, reframe=reframe2)
            # ll.reframe, ll.epoch = reframe2, self.epoch
        else:
            ll = self
        return ll

    @property_doc_(''' this points's datum (L{Datum}).''')
    def datum(self):
        '''Get this point's datum (L{Datum}).
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set this point's datum I{without conversion}.

           @arg datum: New datum (L{Datum}).

           @raise TypeError: The B{C{datum}} is not a L{Datum}
                             or not ellipsoidal.
        '''
        _TypeError(Datum, datum=datum)
        if not datum.isEllipsoidal:
            raise _isnotError('ellipsoidal', datum=datum)
        self._update(datum != self._datum)
        self._datum = datum

    def distanceTo2(self, other):
        '''Approximate the distance and (initial) bearing between this
           and an other (ellipsoidal) point based on the radii of curvature.

           Suitable only for short distances up to a few hundred Km
           or Miles and only between non-near-polar points.

           @arg other: The other point (C{LatLon}).

           @return: An L{Distance2Tuple}C{(distance, initial)}.

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: Incompatible datum ellipsoids.

           @see: Method L{Ellipsoid.distance2} and U{Local, flat earth
                 approximation<https://www.EdWilliams.org/avform.htm#flat>}.
        '''
        return self.ellipsoids(other).distance2(self.lat,  self.lon,
                                               other.lat, other.lon)

    def elevation2(self, adjust=True, datum=Datums.WGS84, timeout=2):
        '''Return elevation of this point for its or the given datum.

           @kwarg adjust: Adjust the elevation for a B{C{datum}} other
                          than C{NAD83}.
           @kwarg datum: Optional datum (L{Datum}).
           @kwarg timeout: Optional query timeout (seconds).

           @return: An L{Elevation2Tuple}C{(elevation, data_source)}
                    or C{(None, error)} in case of errors.

           @note: The adjustment applied is the difference in geocentric
                  earth radius for the B{C{datum}} used and the C{NAV83}
                  datum upon which L{elevations.elevation2} is based.

           @note: NED elevation is only available for locations within
                  the U{Conterminous US (CONUS)
                  <https://WikiPedia.org/wiki/Contiguous_United_States>}.

           @see: Function L{elevations.elevation2} and method
                 L{Ellipsoid.Rgeocentric} for further details and
                 possible C{error}s.
        '''
        if not self._elevation2:  # get elevation and data source
            from pygeodesy.elevations import elevation2
            self._elevation2 = elevation2(self.lat, self.lon,
                                          timeout=timeout)
        return self._Radjust2(adjust, datum, self._elevation2)

    def ellipsoid(self, datum=Datums.WGS84):
        '''Return the ellipsoid of this point's datum or the given datum.

           @kwarg datum: Default datum (L{Datum}).

           @return: The ellipsoid (L{Ellipsoid}).
        '''
        return getattr(self, 'datum', datum).ellipsoid

    def ellipsoids(self, other):
        '''Check the type and ellipsoid of this and an other point's datum.

           @arg other: The other point (C{LatLon}).

           @return: This point's datum ellipsoid (L{Ellipsoid}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

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
            c = E.classname
            raise ValueError('%s %s mistmatch: %ss.%s vs %ss.%s' %
                             ('other', c, c, e.name, c, E.name))
        return E

    @property_doc_(''' this point's observed epoch (C{float}).''')
    def epoch(self):
        '''Get this point's observed epoch (C{float}) or C{None}.
        '''
        return self._epoch or (self.reframe.epoch if self.reframe else None)

    @epoch.setter  # PYCHOK setter!
    def epoch(self, epoch):
        '''Set or clear this point's observed epoch.

           @arg epoch: Observed epoch, a fractional calendar year
                       (C{scalar}) or C{None}.

           @raise TypeError: The B{C{epoch}} is not C{scalar}.
        '''
        self._epoch = None if epoch is None else _2epoch(epoch)

    def geoidHeight2(self, adjust=False, datum=Datums.WGS84, timeout=2):
        '''Return geoid height of this point for its or the given datum.

           @kwarg adjust: Adjust the geoid height for a B{C{datum}} other
                          than C{NAD83/NADV88}.
           @kwarg datum: Optional datum (L{Datum}).
           @kwarg timeout: Optional query timeout (seconds).

           @return: An L{GeoidHeight2Tuple}C{(height, model_name)} or
                    C{(None, error)} in case of errors.

           @note: The adjustment applied is the difference in geocentric
                  earth radius for the given B{C{datum}} and the C{NAV83/NADV88}
                  datum of the L{elevations.geoidHeight2}.

           @note: The geoid height is only available for locations within
                  the U{Conterminous US (CONUS)
                  <https://WikiPedia.org/wiki/Contiguous_United_States>}.

           @see: Function L{elevations.geoidHeight2} and method
                 L{Ellipsoid.Rgeocentric} for further details and
                 possible C{error}s.
        '''
        if not self._geoidHeight2:  # get elevation and data source
            from pygeodesy.elevations import geoidHeight2
            self._geoidHeight2 = geoidHeight2(self.lat, self.lon,
                                              model=0, timeout=timeout)
        return self._Radjust2(adjust, datum, self._geoidHeight2)

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

    @property_RO
    def iteration(self):
        '''Get the iteration number (C{int}) or C{None} if not available/applicable.
        '''
        return None

    def parse(self, strll, height=0, datum=None, sep=','):
        '''Parse a string representing this C{LatLon} point.

           The lat- and longitude must be separated by a sep[arator]
           character.  If height is present it must follow and be
           separated by another sep[arator].  Lat- and longitude
           may be swapped, provided at least one ends with the
           proper compass direction.

           For more details, see functions L{parse3llh} and L{parseDMS}
           in sub-module L{dms}.

           @arg strll: Lat, lon [, height] (string).
           @kwarg height: Optional, default height (C{meter} or C{None}).
           @kwarg datum: Optional, default datum (L{Datum}).
           @kwarg sep: Optional separator (string).

           @return: The point (L{LatLonEllipsoidalBase}).

           @raise ValueError: Invalid B{C{strll}}.
        '''
        from pygeodesy.dms import parse3llh
        a, b, h = parse3llh(strll, height=height, sep=sep)
        return self.classof(a, b, height=h, datum=datum or self.datum)

    @property_doc_(''' this point's reference frame (L{RefFrame}).''')
    def reframe(self):
        '''Get this point's reference frame (L{RefFrame}) or C{None}.
        '''
        return self._reframe

    @reframe.setter  # PYCHOK setter!
    def reframe(self, reframe):
        '''Set or clear this point's reference frame.

           @arg reframe: Reference frame (L{RefFrame}) or C{None}.

           @raise TypeError: The B{C{reframe}} is not a L{RefFrame}.
        '''
        if isinstance(reframe, RefFrame):
            self._reframe = reframe
        elif reframe is not None:
            _TypeError(RefFrame, reframe=reframe)
        elif self.reframe is not None:
            self._reframe = None

    @property_RO
    def scale(self):
        '''Get this point's UTM grid or UPS point scale factor (C{float})
           or C{None} if not converted from L{Utm} or L{Ups}.
        '''
        return self._scale

    def to3xyz(self):  # PYCHOK no cover
        '''DEPRECATED, use method C{toEcef}.

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @note: Overloads C{LatLonBase.to3xyz}
        '''
        r = self.toEcef()
        return self._xnamed(Vector3Tuple(r.x, r.y, r.z))

    def toEtm(self):
        '''Convert this C{LatLon} point to an ETM coordinate.

           @return: The ETM coordinate (L{Etm}).

           @see: Function L{toEtm8}.
        '''
        if self._etm is None:
            from pygeodesy.etm import toEtm8, Etm  # PYCHOK recursive import
            self._etm = toEtm8(self, datum=self.datum, Etm=Etm)
        return self._etm

    def toLcc(self):
        '''Convert this C{LatLon} point to a Lambert location.

           @see: Function L{toLcc} in module L{lcc}.

           @return: The Lambert location (L{Lcc}).
        '''
        if self._lcc is None:
            from pygeodesy.lcc import Lcc, toLcc  # PYCHOK recursive import
            self._lcc = toLcc(self, height=self.height, Lcc=Lcc,
                                                        name=self.name)
        return self._lcc

    def toOsgr(self):
        '''Convert this C{LatLon} point to an OSGR coordinate.

           @see: Function L{toOsgr} in module L{osgr}.

           @return: The OSGR coordinate (L{Osgr}).
        '''
        if self._osgr is None:
            from pygeodesy.osgr import Osgr, toOsgr  # PYCHOK recursive import
            self._osgr = toOsgr(self, datum=self.datum, Osgr=Osgr,
                                                        name=self.name)
        return self._osgr

    def toUps(self, pole='N', falsed=True):
        '''Convert this C{LatLon} point to a UPS coordinate.

           @kwarg pole: Optional top/center of (stereographic)
                        projection (C{str}, 'N[orth]' or 'S[outh]').
           @kwarg falsed: False easting and northing (C{bool}).

           @return: The UPS coordinate (L{Ups}).

           @see: Function L{toUps8}.
        '''
        if self._ups is None:
            from pygeodesy.ups import toUps8, Ups  # PYCHOK recursive import
            self._ups = toUps8(self, datum=self.datum, Ups=Ups,
                                     pole=pole, falsed=falsed)
        return self._ups

    def toUtm(self):
        '''Convert this C{LatLon} point to a UTM coordinate.

           @return: The UTM coordinate (L{Utm}).

           @see: Function L{toUtm8}.
        '''
        if self._utm is None:
            from pygeodesy.utm import toUtm8, Utm  # PYCHOK recursive import
            self._utm = toUtm8(self, datum=self.datum, Utm=Utm)
        return self._utm

    def toUtmUps(self, pole=''):
        '''Convert this C{LatLon} point to a UTM or UPS coordinate.

           @kwarg pole: Optional top/center of UPS (stereographic)
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
            from pygeodesy.utmups import toUtmUps8, Utm, Ups  # PYCHOK recursive import
            u = toUtmUps8(self, datum=self.datum, Utm=Utm, Ups=Ups, pole=pole)
            if isinstance(u, Utm):
                self._utm = u
            elif isinstance(u, Ups):
                self._ups = u
            else:
                _TypeError(Utm, Ups, toUtmUps8=u)
        return u

    def toWm(self):
        '''Convert this C{LatLon} point to a WM coordinate.

           @see: Function L{toWm} in module L{webmercator}.

           @return: The WM coordinate (L{Wm}).
        '''
        if self._wm is None:
            from pygeodesy.webmercator import toWm  # PYCHOK recursive import
            self._wm = toWm(self)
        return self._wm

# **) MIT License
#
# Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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
