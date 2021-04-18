
# -*- coding: utf-8 -*-

u'''(INTERNAL) Ellipsoidal base classes C{CartesianEllipsoidalBase} and
C{LatLonEllipsoidalBase} and several ellipsoidal functions, all used by
C{.ellipsoidalKarney}, C{.ellipsoidalNvector} or C{.ellipsoidalVincenty}.

Pure Python implementation of geodesy tools for ellipsoidal earth models,
transcribed in part from JavaScript originals by I{(C) Chris Veness 2005-2016}
and published under the same MIT Licence**, see for example U{latlon-ellipsoidal
<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}.

@newfield example: Example, Examples
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division

from pygeodesy.basics import issubclassof, _xinstanceof
from pygeodesy.cartesianBase import CartesianBase
from pygeodesy.datums import Datum, Datums, _ellipsoidal_datum, _WGS84
from pygeodesy.errors import _AssertionError, _incompatible, IntersectionError, \
                             _IsnotError, RangeError, TRFError, _ValueError, \
                             _xellipsoidal
from pygeodesy.fmath import euclid, favg, fmean_, fsum_
from pygeodesy.formy import _radical2
from pygeodesy.interns import _ellipsoidal_  # PYCHOK used!
from pygeodesy.interns import EPS, EPS0, EPS1, MISSING, NN, PI, _COMMA_, \
                             _conversion_, _datum_, _DOT_, _exceed_PI_radians_, \
                             _N_, _no_, _near_concentric_, _SPACE_, _too_
from pygeodesy.latlonBase import LatLonBase, _trilaterate5
from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.named import _xnamed
from pygeodesy.namedTuples import _LL4Tuple, Vector3Tuple
from pygeodesy.props import deprecated_method, Property_RO, \
                            property_doc_, property_RO
from pygeodesy.streprs import Fmt
from pygeodesy.units import Epoch, Height, Radius_, Scalar, _1mm as _TOL_M
from pygeodesy.utily import m2degrees, unroll180

__all__ = ()
__version__ = '21.04.15'

_reframe_ = 'reframe'
_TRIPS    =  17  # _intersects2, _nearestOn interations, 6 is sufficient


class CartesianEllipsoidalBase(CartesianBase):
    '''(INTERNAL) Base class for ellipsoidal C{Cartesian}s.
    '''
    _datum = _WGS84  # L{Datum}

    def _applyHelmerts(self, *transforms):
        '''(INTERNAL) Apply one I{or more} Helmert transforms.
        '''
        xyz = self.xyz
        for t in transforms:
            xyz = t.transform(*xyz)
        return self.classof(xyz, datum=self.datum)

    @deprecated_method
    def convertRefFrame(self, reframe2, reframe, epoch=None):
        '''DEPRECATED, use method L{toRefFrame}.'''
        return self.toRefFrame(reframe2, reframe, epoch=epoch)

    def toRefFrame(self, reframe2, reframe, epoch=None):
        '''Convert this cartesian point from one to an other reference frame.

           @arg reframe2: Reference frame to convert I{to} (L{RefFrame}).
           @arg reframe: Reference frame to convert I{from} (L{RefFrame}).
           @kwarg epoch: Optional epoch to observe (C{scalar}, fractional
                         calendar year), overriding B{C{reframe}}'s epoch.

           @return: The converted point (C{Cartesian}) or this point if
                    conversion is C{nil}.

           @raise TRFError: No conversion available from B{C{reframe}}
                            to B{C{reframe2}} or invalid B{C{epoch}}.

           @raise TypeError: B{C{reframe2}} or B{C{reframe}} not a
                             L{RefFrame}.
        '''
        from pygeodesy.trf import RefFrame, _reframeTransforms2
        _xinstanceof(RefFrame, reframe2=reframe2, reframe=reframe)

        _, xs = _reframeTransforms2(reframe2, reframe, epoch)
        return self._applyHelmerts(*xs) if xs else self


class LatLonEllipsoidalBase(LatLonBase):
    '''(INTERNAL) Base class for ellipsoidal C{LatLon}s.
    '''
    _convergence    =  None   # UTM/UPS meridian convergence (C{degrees})
    _datum          = _WGS84  # L{Datum}
    _elevation2to   =  None   # _elevation2 timeout (C{secs})
    _epoch          =  None   # overriding .reframe.epoch (C{float})
    _geoidHeight2to =  None   # _geoidHeight2 timeout (C{secs})
    _iteration      =  None   # iteration number (C{int} or C{None})
    _reframe        =  None   # reference frame (L{RefFrame})
    _scale          =  None   # UTM/UPS scale factor (C{float})

    def __init__(self, lat, lon, height=0, datum=None, reframe=None,
                                           epoch=None, name=NN):
        '''Create an ellipsoidal C{LatLon} point frome the given
           lat-, longitude and height on the given datum and with
           the given reference frame and epoch.

           @arg lat: Latitude (C{degrees} or DMS C{[N|S]}).
           @arg lon: Longitude (C{degrees} or DMS C{str[E|W]}).
           @kwarg height: Optional elevation (C{meter}, the same units
                          as the datum's half-axes).
           @kwarg datum: Optional, ellipsoidal datum to use (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg reframe: Optional reference frame (L{RefFrame}).
           @kwarg epoch: Optional epoch to observe for B{C{reframe}}
                         (C{scalar}), a non-zero, fractional calendar
                         year; silently ignored if C{B{reframe}=None}.
           @kwarg name: Optional name (string).

           @raise RangeError: Value of B{C{lat}} or B{C{lon}} outside the valid
                              range and C{rangerrors} set to C{True}.

           @raise TypeError: B{C{datum}} is not a L{datum}, B{C{reframe}}
                             is not a L{RefFrame} or B{C{epoch}} is not
                             C{scalar} non-zero.

           @raise UnitError: Invalid B{C{lat}}, B{C{lon}} or B{C{height}}.

           @example:

            >>> p = LatLon(51.4778, -0.0016)  # height=0, datum=Datums.WGS84
        '''
        LatLonBase.__init__(self, lat, lon, height=height, name=name)
        if datum not in (None, self._datum):
            self.datum = _ellipsoidal_datum(datum, name=name)
        if reframe:
            self.reframe = reframe
            self.epoch = epoch

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
        '''Get this point's UTM or UPS meridian convergence (C{degrees}) or
           C{None} if not available or not converted from L{Utm} or L{Ups}.
        '''
        return self._convergence

    @deprecated_method
    def convertDatum(self, datum2):
        '''DEPRECATED, use method L{toDatum}.'''
        return self.toDatum(datum2)

    @deprecated_method
    def convertRefFrame(self, reframe2):
        '''DEPRECATED, use method L{toRefFrame}.'''
        return self.toRefFrame(reframe2)

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
        _xinstanceof(Datum, datum=datum)
        if not datum.isEllipsoidal:
            raise _IsnotError(_ellipsoidal_, datum=datum)
        self._update(datum != self._datum)
        self._datum = datum

    def distanceTo2(self, other):
        '''I{Approximate} the distance and (initial) bearing between this
           and an other (ellipsoidal) point based on the radii of curvature.

           I{Suitable only for short distances up to a few hundred Km
           or Miles and only between points not near-polar}.

           @arg other: The other point (C{LatLon}).

           @return: An L{Distance2Tuple}C{(distance, initial)}.

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: Incompatible datum ellipsoids.

           @see: Method L{Ellipsoid.distance2} and U{Local, flat earth
                 approximation<https://www.EdWilliams.org/avform.htm#flat>}
                 aka U{Hubeny<https://www.OVG.AT/de/vgi/files/pdf/3781/>}
                 formula.
        '''
        return self.ellipsoids(other).distance2(self.lat,  self.lon,
                                               other.lat, other.lon)

    @Property_RO
    def _elevation2(self):
        '''(INTERNAL) Get elevation and data source.
        '''
        from pygeodesy.elevations import elevation2
        return elevation2(self.lat, self.lon, timeout=self._elevation2to)

    def elevation2(self, adjust=True, datum=_WGS84, timeout=2):
        '''Return elevation of this point for its or the given datum.

           @kwarg adjust: Adjust the elevation for a B{C{datum}} other
                          than C{NAD83} (C{bool}).
           @kwarg datum: Optional datum (L{Datum}).
           @kwarg timeout: Optional query timeout (C{seconds}).

           @return: An L{Elevation2Tuple}C{(elevation, data_source)}
                    or C{(None, error)} in case of errors.

           @note: The adjustment applied is the difference in geocentric
                  earth radius between the B{C{datum}} and C{NAV83}
                  upon which the L{elevations.elevation2} is based.

           @note: NED elevation is only available for locations within
                  the U{Conterminous US (CONUS)
                  <https://WikiPedia.org/wiki/Contiguous_United_States>}.

           @see: Function L{elevations.elevation2} and method
                 L{Ellipsoid.Rgeocentric} for further details and
                 possible C{error}s.
        '''
        if self._elevation2to != timeout:
            self._elevation2to = timeout
            LatLonEllipsoidalBase._elevation2._update(self)
        return self._Radjust2(adjust, datum, self._elevation2)

    def ellipsoid(self, datum=_WGS84):
        '''Return the ellipsoid of this point's datum or the given datum.

           @kwarg datum: Default datum (L{Datum}).

           @return: The ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).
        '''
        return getattr(self, _datum_, datum).ellipsoid

    def ellipsoids(self, other):
        '''Check the type and ellipsoid of this and an other point's datum.

           @arg other: The other point (C{LatLon}).

           @return: This point's datum ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: Incompatible datum ellipsoids.
        '''
        self.others(other, up=2)  # ellipsoids' caller

        E = self.ellipsoid()
        try:  # other may be Sphere, etc.
            e = other.ellipsoid()
        except AttributeError:  # PYCHOK no cover
            try:  # no ellipsoid method, try datum
                e = other.datum.ellipsoid
            except AttributeError:
                e = E  # no datum, XXX assume equivalent?
        if e != E:
            raise _ValueError(e.named2, txt=_incompatible(E.named2))
        return E

    @property_doc_(''' this point's observed or C{reframe} epoch (C{float}).''')
    def epoch(self):
        '''Get this point's observed or C{reframe} epoch (C{float}) or C{None}.
        '''
        return self._epoch or (self.reframe.epoch if self.reframe else None)

    @epoch.setter  # PYCHOK setter!
    def epoch(self, epoch):
        '''Set or clear this point's observed epoch.

           @arg epoch: Observed epoch, a fractional calendar year
                       (L{Epoch}, C{scalar}) or C{None}.

           @raise TRFError: Invalid B{C{epoch}}.
        '''
        self._epoch = None if epoch is None else Epoch(epoch)

    @Property_RO
    def _etm(self):
        '''(INTERNAL) Get this C{LatLon} point as an ETM coordinate (L{toEtm8}).
        '''
        from pygeodesy.etm import toEtm8, Etm
        return toEtm8(self, datum=self.datum, Etm=Etm)

    @Property_RO
    def _geoidHeight2(self):
        '''(INTERNAL) Get geoid height and model.
        '''
        from pygeodesy.elevations import geoidHeight2
        return geoidHeight2(self.lat, self.lon, model=0, timeout=self._geoidHeight2to)

    def geoidHeight2(self, adjust=False, datum=_WGS84, timeout=2):
        '''Return geoid height of this point for its or the given datum.

           @kwarg adjust: Adjust the geoid height for a B{C{datum}}
                          other than C{NAD83/NADV88} (C{bool}).
           @kwarg datum: Optional datum (L{Datum}).
           @kwarg timeout: Optional query timeout (C{seconds}).

           @return: A L{GeoidHeight2Tuple}C{(height, model_name)} or
                    C{(None, error)} in case of errors.

           @note: The adjustment applied is the difference in geocentric
                  earth radius between the B{C{datum}} and C{NAV83/NADV88}
                  upon which the L{elevations.geoidHeight2} is based.

           @note: The geoid height is only available for locations within
                  the U{Conterminous US (CONUS)
                  <https://WikiPedia.org/wiki/Contiguous_United_States>}.

           @see: Function L{elevations.geoidHeight2} and method
                 L{Ellipsoid.Rgeocentric} for further details and
                 possible C{error}s.
        '''
        if self._geoidHeight2to != timeout:
            self._geoidHeight2to = timeout
            LatLonEllipsoidalBase._geoidHeight2._update(self)
        return self._Radjust2(adjust, datum, self._geoidHeight2)

    def intersections2(self, radius1, other, radius2, height=None, wrap=True,
                                                 equidistant=None, tol=_TOL_M):
        '''Compute the intersection points of two circles each defined
           by a center point and a radius.

           @arg radius1: Radius of the this circle (C{meter}, conventionally).
           @arg other: Center of the other circle (C{LatLon}).
           @arg radius2: Radius of the other circle (C{meter}, same units as
                         B{C{radius1}}).
           @kwarg height: Optional height for the intersection points,
                          overriding the "radical height" at the "radical
                          line" between both centers (C{meter}) or C{None}.
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).
           @kwarg equidistant: An azimuthal equidistant projection class
                               (L{Equidistant} or L{EquidistantKarney}),
                               function L{equidistant} will be invoked
                               if left unspecified.
           @kwarg tol: Convergence tolerance (C{meter}, same units as B{C{radius1}}
                       and B{C{radius2}}).

           @return: 2-Tuple of the intersection points, each a C{LatLon}
                    instance.  For abutting circles, both intersection
                    points are the same instance.

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found.

           @raise IntersectionError: Concentric, antipodal, invalid or
                                     non-intersecting circles or no
                                     convergence for B{C{tol}}.

           @raise TypeError: Invalid B{C{other}} or B{C{equidistant}}.

           @raise UnitError: Invalid B{C{radius1}}, B{C{radius2}} or B{C{height}}.

           @see: U{The B{ellipsoidal} case<https://GIS.StackExchange.com/questions/48937/
                 calculating-intersection-of-two-circles>}, U{Karney's paper
                 <https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section B{14. MARITIME BOUNDARIES},
                 U{circle-circle<https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>} and
                 U{sphere-sphere<https://MathWorld.Wolfram.com/Sphere-SphereIntersection.html>}
                 intersections.
        '''
        self.others(other)
        return _intersections2(self, radius1, other, radius2, height=height, wrap=wrap,
                                     equidistant=equidistant, tol=tol,
                                     LatLon=self.classof, datum=self.datum)

    @property_RO
    def iteration(self):
        '''Get the most recent C{intersections2} or C{nearestOn} iteration
           number (C{int}) or C{None} if not available/applicable.
        '''
        return self._iteration

    @Property_RO
    def _lcc(self):
        '''(INTERNAL) Get this C{LatLon} point to a Lambert location (L{Lcc}).
        '''
        from pygeodesy.lcc import Lcc, toLcc
        return toLcc(self, height=self.height, Lcc=Lcc, name=self.name)

    def nearestOn(self, point1, point2, within=True, height=None, wrap=True,
                                        equidistant=None, tol=_TOL_M):
        '''Locate the closest point between two other points.

           @arg point1: Start point (C{LatLon}).
           @arg point2: End point (C{LatLon}).
           @kwarg within: If C{True} return the closest point I{between}
                          B{C{point1}} and B{C{point2}}, otherwise the
                          closest point elsewhere on the arc (C{bool}).
           @kwarg height: Optional height for the closest point (C{meter})
                          or C{None} to interpolate the height.
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).
           @kwarg equidistant: An azimuthal equidistant projection class
                               (L{Equidistant} or L{EquidistantKarney}),
                               function L{equidistant} will be invoked
                               if left unspecified.
           @kwarg tol: Convergence tolerance (C{meter}).

           @return: Closest point (C{LatLon}).

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found.

           @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{equidistant}}.
        '''
        p1 = self.others(point1=point1)
        p2 = self.others(point2=point2)
        return _nearestOn(self, p1, p2, within=within, height=height, wrap=wrap,
                                equidistant=equidistant, tol=tol,
                                LatLon=self.classof, datum=self.datum)

    @Property_RO
    def _osgr(self):
        '''(INTERNAL) Get this C{LatLon} point to an OSGR coordinate (L{Osgr}).
        '''
        from pygeodesy.osgr import Osgr, toOsgr
        return toOsgr(self, datum=self.datum, Osgr=Osgr, name=self.name)

    def parse(self, strllh, height=0, datum=None, sep=_COMMA_, name=NN):
        '''Parse a string representing a similar, ellipsoidal C{LatLon}
           point, consisting of C{"lat, lon[, height]"}.

           @arg strllh: Lat, lon and optional height (C{str}),
                        see function L{parse3llh}.
           @kwarg height: Optional, default height (C{meter} or
                          C{None}).
           @kwarg datum: Optional datum (L{Datum}), overriding this
                         datum I{without conversion}.
           @kwarg sep: Optional separator (C{str}).
           @kwarg name: Optional instance name (C{str}), overriding
                        this name.

           @return: The similar point (ellipsoidal C{LatLon}).

           @raise ParseError: Invalid B{C{strllh}}.
        '''
        from pygeodesy.dms import parse3llh
        a, b, h = parse3llh(strllh, height=height, sep=sep)
        r = self.classof(a, b, height=h, datum=self.datum)
        if datum not in (None, self.datum):
            r.datum = datum
        if name:
            r = _xnamed(r, name, force=True)
        return r

    def _Radjust2(self, adjust, datum, meter_text2):
        '''(INTERNAL) Adjust elevation or geoidHeight with difference
           in Gaussian radii of curvature of given datum and NAD83.

           @note: This is an arbitrary, possibly incorrect adjustment.
        '''
        if adjust:  # Elevation2Tuple or GeoidHeight2Tuple
            m, t = meter_text2
            if isinstance(m, float) and abs(m) > EPS:
                n = Datums.NAD83.ellipsoid.rocGauss(self.lat)
                if n > EPS0:
                    # use ratio, datum and NAD83 units may differ
                    e = self.ellipsoid(datum).rocGauss(self.lat)
                    if e > EPS0 and abs(e - n) > EPS:  # EPS1
                        m *= e / n
                        meter_text2 = meter_text2.classof(m, t)
        return self._xnamed(meter_text2)

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
        if reframe is not None:
            from pygeodesy.trf import RefFrame
            _xinstanceof(RefFrame, reframe=reframe)
            self._reframe = reframe
        elif self.reframe is not None:
            self._reframe = None

    @Property_RO
    def scale(self):
        '''Get this point's UTM grid or UPS point scale factor (C{float})
           or C{None} if not converted from L{Utm} or L{Ups}.
        '''
        return self._scale

    def toDatum(self, datum2):
        '''Convert this point to an other datum.

           @arg datum2: Datum to convert I{to} (L{Datum}).

           @return: The converted point (ellipsoidal C{LatLon}).

           @raise TypeError: The B{C{datum2}} invalid.

           @example:

            >>> p = LatLon(51.4778, -0.0016)  # default Datums.WGS84
            >>> p.toDatum(Datums.OSGB36)  # 51.477284°N, 000.00002°E
        '''
        d2 = _ellipsoidal_datum(datum2, name=self.name)
        if self.datum == d2:
            return self.copy()

        c = self.toCartesian().toDatum(d2)
        return c.toLatLon(datum=d2, LatLon=self.classof)

    def toEtm(self):
        '''Convert this C{LatLon} point to an ETM coordinate.

           @return: The ETM coordinate (L{Etm}).

           @see: Function L{toEtm8}.
        '''
        return self._etm

    def toLcc(self):
        '''Convert this C{LatLon} point to a Lambert location.

           @see: Function L{toLcc} in module L{lcc}.

           @return: The Lambert location (L{Lcc}).
        '''
        return self._lcc

    def toOsgr(self):
        '''Convert this C{LatLon} point to an OSGR coordinate.

           @see: Function L{toOsgr} in module L{osgr}.

           @return: The OSGR coordinate (L{Osgr}).
        '''
        return self._osgr

    def toRefFrame(self, reframe2):
        '''Convert this point to an other reference frame.

           @arg reframe2: Reference frame to convert I{to} (L{RefFrame}).

           @return: The converted point (ellipsoidal C{LatLon}) or this
                    point if conversion is C{nil}.

           @raise TRFError: This B{C{reframe}} not defined or no
                            conversion available from this B{C{reframe}}
                            to B{C{reframe2}}.

           @raise TypeError: The B{C{reframe2}} is not a L{RefFrame}.

           @example:

            >>> p = LatLon(51.4778, -0.0016, reframe=RefFrames.ETRF2000)  # default Datums.WGS84
            >>> p.toRefFrame(RefFrames.ITRF2014)  # 51.477803°N, 000.001597°W, +0.01m
        '''
        from pygeodesy.trf import RefFrame, _reframeTransforms2
        _xinstanceof(RefFrame, reframe2=reframe2)

        if not self.reframe:
            t = _SPACE_(_DOT_(repr(self), _reframe_), MISSING)
            raise TRFError(_no_(_conversion_), txt=t)

        e, xs = _reframeTransforms2(reframe2, self.reframe, self.epoch)
        if xs:
            c = self.toCartesian()._applyHelmerts(*xs)
            ll = c.toLatLon(datum=self.datum, epoch=e, reframe=reframe2,
                                              LatLon=self.classof)
        else:
            ll = self
        return ll

    def toUps(self, pole=_N_, falsed=True):
        '''Convert this C{LatLon} point to a UPS coordinate.

           @kwarg pole: Optional top/center of (stereographic)
                        projection (C{str}, 'N[orth]' or 'S[outh]').
           @kwarg falsed: False easting and northing (C{bool}).

           @return: The UPS coordinate (L{Ups}).

           @see: Function L{toUps8}.
        '''
        if self._upsOK(pole, falsed):
            u = self._ups
        else:
            from pygeodesy.ups import toUps8, Ups
            u = toUps8(self, datum=self.datum, Ups=Ups,
                              pole=pole, falsed=falsed)
        return u

    def toUtm(self):
        '''Convert this C{LatLon} point to a UTM coordinate.

           @return: The UTM coordinate (L{Utm}).

           @see: Function L{toUtm8}.
        '''
        return self._utm

    def toUtmUps(self, pole=NN):
        '''Convert this C{LatLon} point to a UTM or UPS coordinate.

           @kwarg pole: Optional top/center of UPS (stereographic)
                        projection (C{str}, 'N[orth]' or 'S[outh]').

           @return: The UTM or UPS coordinate (L{Utm} or L{Ups}).

           @raise TypeError: Result in L{Utm} or L{Ups}.

           @see: Function L{toUtmUps}.
        '''
        if self._utmOK():
            u = self._utm
        elif self._upsOK(pole):
            u = self._ups
        else:  # no cover
            from pygeodesy.utmups import toUtmUps8, Utm, Ups
            u = toUtmUps8(self, datum=self.datum, Utm=Utm, Ups=Ups,
                                 pole=pole, name=self.name)
            if isinstance(u, Utm):
                self._overwrite(_utm=u)
            elif isinstance(u, Ups):
                self._overwrite(_ups=u)
            else:
                _xinstanceof(Utm, Ups, toUtmUps8=u)
        return u

    def toWm(self):
        '''Convert this C{LatLon} point to a WM coordinate.

           @see: Function L{toWm} in module L{webmercator}.

           @return: The WM coordinate (L{Wm}).
        '''
        return self._wm

    @deprecated_method
    def to3xyz(self):  # PYCHOK no cover
        '''DEPRECATED, use method L{toEcef}.

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @note: Overloads C{LatLonBase.to3xyz}
        '''
        r = self.toEcef()
        return Vector3Tuple(r.x, r.y, r.z, name=self.name)

    def trilaterate5(self, distance1, point2, distance2, point3, distance3,
                           area=True, eps=EPS1, wrap=False):
        '''Trilaterate three points by area overlap or perimeter intersection
           of three intersecting circles.

           @arg distance1: Distance to this point (C{meter}), same units
                           as B{C{eps}}).
           @arg point2: Second center point (C{LatLon}).
           @arg distance2: Distance to point2 (C{meter}, same units as
                           B{C{eps}}).
           @arg point3: Third center point (C{LatLon}).
           @arg distance3: Distance to point3 (C{meter}, same units as
                           B{C{eps}}).
           @kwarg area: If C{True} compute the area overlap, otherwise the
                        perimeter intersection of the circles (C{bool}).
           @kwarg eps: The required I{minimal overlap} for C{B{area}=True}
                       or the I{intersection margin} for C{B{area}=False}
                       (C{meter}, conventionally).
           @kwarg wrap: Wrap/unroll angular distances (C{bool}).

           @return: A L{Trilaterate5Tuple}C{(min, minPoint, max, maxPoint, n)}
                    with C{min} and C{max} in C{meter}, same units as B{C{eps}},
                    the corresponding trilaterated points C{minPoint} and
                    C{maxPoint} as I{ellipsoidal} C{LatLon} and C{n}, the number
                    of trilatered points found for the given B{C{eps}}.

                    If only a single trilaterated point is found, C{min I{is}
                    max}, C{minPoint I{is} maxPoint} and C{n = 1}.

                    For C{B{area}=True}, C{min} and C{max} are the smallest
                    respectively largest I{radial} overlap found.

                    For C{B{area}=False}, C{min} and C{max} represent the
                    nearest respectively farthest intersection margin.

                    If C{B{area}=True} and all 3 circles are concentric, C{n=0}
                    and C{minPoint} and C{maxPoint} are the B{C{point#}} with
                    the smallest B{C{distance#}} C{min} respectively C{max} the
                    largest B{C{distance#}}.

           @raise IntersectionError: Trilateration failed for the given B{C{eps}},
                                     insufficient overlap for C{B{area}=True} or
                                     no intersection or all (near-)concentric for
                                     C{B{area}=False}.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @raise ValueError: Coincident B{C{points}} or invalid B{C{distance1}},
                              B{C{distance2}} or B{C{distance3}}.

           @note: Ellipsoidal trilateration invokes methods C{LatLon.intersections2}
                  and C{LatLon.nearestOn}.  Install Karney's Python package
                  U{geographiclib<https://PyPI.org/project/geographiclib>} to obtain
                  the most accurate results for both C{ellipsoidalVincenty.LatLon}
                  and C{ellipsoidalKarney.LatLon} points.
        '''
        return _trilaterate5(self, distance1,
                             self.others(point2=point2), distance2,
                             self.others(point3=point3), distance3,
                             area=area, eps=eps, wrap=wrap)

    @Property_RO
    def _ups(self):  # __dict__ value overwritten by method C{toUtmUps}
        '''(INTERNAL) Get this C{LatLon} point as UPS coordinate (L{Ups}), see L{toUps8}.
        '''
        from pygeodesy.ups import toUps8, Ups
        return toUps8(self, datum=self.datum, Ups=Ups,
                             pole=NN, falsed=True, name=self.name)

    def _upsOK(self, pole=NN, falsed=True):
        '''(INTERNAL) Check matching C{Ups}.
        '''
        try:
            u = self._ups
        except RangeError:
            return False
        return falsed and (u.pole == pole[:1].upper() or not pole)

    @Property_RO
    def _utm(self):  # __dict__ value overwritten by method C{toUtmUps}
        '''(INTERNAL) Get this C{LatLon} point as UTM coordinate (L{Utm}), see L{toUtm8}.
        '''
        from pygeodesy.utm import toUtm8, Utm
        return toUtm8(self, datum=self.datum, Utm=Utm, name=self.name)

    def _utmOK(self):
        '''(INTERNAL) Check C{Utm}.
        '''
        try:
            _ = self._utm
        except RangeError:
            return False
        return True

    @Property_RO
    def _wm(self):
        '''(INTERNAL) Get this C{LatLon} point as webmercator (L{WM}).
        '''
        from pygeodesy.webmercator import toWm
        return toWm(self)


def _Equidistant2(equidistant, datum):
    # (INTERNAL) Get an C{Equidistant} or C{EquidistantKarney} instance
    import pygeodesy.azimuthal as _az

    if equidistant is None or not callable(equidistant):
        equidistant = _az.equidistant
    elif not (issubclassof(equidistant, _az.Equidistant) or
              issubclassof(equidistant, _az.EquidistantKarney)):
        raise _IsnotError(_az.Equidistant.__name__,
                          _az.EquidistantKarney.__name__,
                           equidistant=equidistant)
    return equidistant(0, 0, datum)


def _intermediateTo(p1, p2, fraction, height, wrap):
    # (INTERNAL) Helper for C{ellipsoidalKarney.LatLon.intermediateTo}
    # and C{ellipsoidalVincenty.LatLon.intermediateTo}.
    t = p1.distanceTo3(p2, wrap=wrap)
    f = Scalar(fraction=fraction)
    h = p1._havg(p2, f=f) if height is None else Height(height)
    return p1.destination(t.distance * f, t.initial, height=h)


def _intersections2(center1, radius1, center2, radius2, height=None, wrap=True,
                    equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    # (INTERNAL) Iteratively compute the intersection points of two circles
    # each defined by an (ellipsoidal) center point and a radius, imported
    # by .ellipsoidalKarney and -Vincenty

    c1 = _xellipsoidal(center1=center1)
    c2 = c1.others(center2=center2)

    r1 = Radius_(radius1=radius1)
    r2 = Radius_(radius2=radius2)

    try:
        return _intersects2(c1, r1, c2, r2, height=height, wrap=wrap,
                                       equidistant=equidistant, tol=tol,
                                            LatLon=LatLon, **LatLon_kwds)
    except (TypeError, ValueError) as x:
        raise IntersectionError(center1=center1, radius1=radius1,
                                center2=center2, radius2=radius2, txt=str(x))


def _intersects2(c1, r1, c2, r2, height=None, wrap=True,  # MCCABE 17
                 equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    # (INTERNAL) Intersect two (ellipsoidal) circles, see L{_intersections2}
    # above, separated to allow callers to embellish any exceptions

    from pygeodesy.sphericalTrigonometry import _intersects2 as _si2, LatLon as _LLS
    from pygeodesy.vector3d import _intersects2 as _vi2

    def _latlon4(t, h, n):
        r = _LL4Tuple(t.lat, t.lon, h, t.datum, LatLon, LatLon_kwds, name=n)
        r._iteration = t.iteration  # ._iteration for tests
        return r

    if r1 < r2:
        c1, c2 = c2, c1
        r1, r2 = r2, r1

    E = c1.ellipsoids(c2)
    if r1 > (min(E.b, E.a) * PI):
        raise ValueError(_exceed_PI_radians_)

    if wrap:  # unroll180 == .karney._unroll2
        c2 = _unrollon(c1, c2)

    # distance between centers and radii are
    # measured along the ellipsoid's surface
    m = c1.distanceTo(c2, wrap=False)  # meter
    if m < max(r1 - r2, EPS):
        raise ValueError(_near_concentric_)
    if fsum_(r1, r2, -m) < 0:
        raise ValueError(_too_(Fmt.distant(m)))

    f = _radical2(m, r1, r2).ratio  # "radical fraction"
    r = E.rocMean(favg(c1.lat, c2.lat, f=f))
    e = max(m2degrees(tol, radius=r), EPS)

    # get the azimuthal equidistant projection
    A = _Equidistant2(equidistant, c1.datum)

    # gu-/estimate initial intersections, spherically ...
    t1, t2 = _si2(_LLS(c1.lat, c1.lon, height=c1.height), r1,
                  _LLS(c2.lat, c2.lon, height=c2.height), r2,
                   radius=r, height=height, wrap=False, too_d=m)
    h, n = t1.height, t1.name

    # ... and then iterate like Karney suggests to find
    # tri-points of median lines, @see: references under
    # method LatLonEllipsoidalBase.intersections2 above
    ts, ta = [], None
    for t in ((t1,) if t1 is t2 else (t1, t2)):
        p = None  # force first d == p to False
        for i in range(1, _TRIPS):
            A.reset(t.lat, t.lon)  # gu-/estimate as origin
            # convert centers to projection space
            t1 = A.forward(c1.lat, c1.lon)
            t2 = A.forward(c2.lat, c2.lon)
            # compute intersections in projection space
            v1, v2 = _vi2(t1, r1,  # XXX * t1.scale?,
                          t2, r2,  # XXX * t2.scale?,
                          sphere=False, too_d=m)
            # convert intersections back to geodetic
            t1 = A.reverse(v1.x, v1.y)
            d1 = euclid(t1.lat - t.lat, t1.lon - t.lon)
            if v1 is v2:  # abutting
                t, d = t1, d1
            else:
                t2 = A.reverse(v2.x, v2.y)
                d2 = euclid(t2.lat - t.lat, t2.lon - t.lon)
                # consider only the closer intersection
                t, d = (t1, d1) if d1 < d2 else (t2, d2)
            # break if below tolerance or if unchanged
            if d < e or d == p:
                t._iteration = i  # _NamedTuple._iteration
                ts.append(t)
                if v1 is v2:  # abutting
                    ta = t
                break
            p = d
        else:
            raise ValueError(_no_(Fmt.convergence(tol)))

    if ta:  # abutting circles
        pass
    elif len(ts) == 2:
        return (_latlon4(ts[0], h, n),
                _latlon4(ts[1], h, n))
    elif len(ts) == 1:  # PYCHOK no cover
        ta = ts[0]  # assume abutting
    else:
        raise _AssertionError(ts=ts)
    r = _latlon4(ta, h, n)
    return r, r


def _nearestOn(p, p1, p2, within=True, height=None, wrap=True,
               equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    # (INTERNAL) Get closet point, like L{_intersects2} above,
    # separated to allow callers to embellish any exceptions

    from pygeodesy.sphericalNvector import LatLon as _LLS
    from pygeodesy.vector3d import _nearestOn as _vnOn, Vector3d

    def _v(t, h):
        return Vector3d(t.x, t.y, h)

    _ = p.ellipsoids(p1)
    E = p.ellipsoids(p2)

    if wrap:
        p1 = _unrollon(p,  p1)
        p2 = _unrollon(p,  p2)
        p2 = _unrollon(p1, p2)

    r = E.rocMean(fmean_(p.lat, p1.lat, p2.lat))
    e = max(m2degrees(tol, radius=r), EPS)

    # get the azimuthal equidistant projection
    A = _Equidistant2(equidistant, p.datum)

    # gu-/estimate initial nearestOn, spherically ... wrap=False
    t = _LLS(p.lat,  p.lon,  height=p.height).nearestOn(
        _LLS(p1.lat, p1.lon, height=p1.height),
        _LLS(p2.lat, p2.lon, height=p2.height), within=within, height=height)
    n = t.name

    if height is False:  # use height as Z component
        h  = t.height
        h1 = p1.height
        h2 = p2.height
    else:
        h = h1 = h2 = 0

    # ... and then iterate like Karney suggests to find
    # tri-points of median lines, @see: references under
    # method LatLonEllipsoidalBase.intersections2 above
    c = None  # force first d == c to False
    # closest to origin, .z to interpolate height
    p = Vector3d(0, 0, h)
    for i in range(1, _TRIPS):
        A.reset(t.lat, t.lon)  # gu-/estimate as origin
        # convert points to projection space
        t1 = A.forward(p1.lat, p1.lon)
        t2 = A.forward(p2.lat, p2.lon)
        # compute nearestOn in projection space
        v = _vnOn(p, _v(t1, h1), _v(t2, h2), within=within)
        # convert nearestOn back to geodetic
        r = A.reverse(v.x, v.y)
        d = euclid(r.lat - t.lat, r.lon - t.lon)
        # break if below tolerance or if unchanged
        t = r
        if d < e or d == c:
            t._iteration = i  # _NamedTuple._iteration
            if height is False:
                h = v.z  # nearest interpolated
            break
        c = d
    else:
        raise ValueError(_no_(Fmt.convergence(tol)))

    r = _LL4Tuple(t.lat, t.lon, h, t.datum, LatLon, LatLon_kwds, name=n)
    r._iteration = t.iteration  # ._iteration for tests
    return r


def _unrollon(p1, p2):  # unroll180 == .karney._unroll2
    # wrap, unroll and replace longitude if different
    _, lon = unroll180(p1.lon, p2.lon, wrap=True)
    if abs(lon - p2.lon) > EPS:
        p2 = p2.classof(p2.lat, lon, p2.height, datum=p2.datum)
    return p2


__all__ += _ALL_DOCS(CartesianEllipsoidalBase, LatLonEllipsoidalBase)

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
