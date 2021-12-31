
# -*- coding: utf-8 -*-

u'''(INTERNAL) Ellipsoidal geodesy base classes C{CartesianEllipsoidalBase}
and C{LatLonEllipsoidalBase}.

Pure Python implementation of geodesy tools for ellipsoidal earth models,
transcoded in part from JavaScript originals by I{(C) Chris Veness 2005-2016}
and published under the same MIT Licence**, see for example U{latlon-ellipsoidal
<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import _xinstanceof
from pygeodesy.cartesianBase import CartesianBase, Vector3Tuple
from pygeodesy.datums import Datum, Datums, _ellipsoidal_datum, \
                            _spherical_datum, _WGS84
from pygeodesy.errors import _incompatible, _IsnotError, RangeError, TRFError, \
                             _ValueError, _xellipsoidal, _xError, _xkwds, \
                             _xkwds_get, _xkwds_not
from pygeodesy.interns import _ellipsoidal_  # PYCHOK used!
from pygeodesy.interns import EPS, EPS0, EPS1, MISSING, NN, _COMMA_, \
                             _conversion_, _datum_, _DOT_, _N_, _no_, \
                             _reframe_, _SPACE_, _0_0
from pygeodesy.latlonBase import LatLonBase, _trilaterate5
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
# from pygeodesy.namedTuples import Vector3Tuple  # from .cartesianBase
from pygeodesy.props import deprecated_method, Property_RO, \
                            property_doc_, property_RO
from pygeodesy.units import Epoch, _1mm as _TOL_M, Radius_

__all__ = _ALL_LAZY.ellipsoidalBase
__version__ = '21.12.29'


class CartesianEllipsoidalBase(CartesianBase):
    '''(INTERNAL) Base class for ellipsoidal C{Cartesian}s.
    '''
    _datum   = _WGS84  # L{Datum}
    _reframe =  None

#   def __matmul__(self, other):  # PYCHOK Python 3.5+
#       '''Return C{NotImplemented} for C{c_ = c @ datum}, C{c_ = c @ reframe} and C{c_ = c @ Transform}.
#       '''
#       from pygeodesy.trf import RefFrame
#       return NotImplemented if isinstance(other, (Datum, RefFrame, Transform)) else \
#             _NotImplemented(self, other)

    @deprecated_method
    def convertRefFrame(self, reframe2, reframe, epoch=None):
        '''DEPRECATED, use method L{toRefFrame}.'''
        return self.toRefFrame(reframe2, reframe, epoch=epoch)

    def intersections2(self, radius, center2, radius2, sphere=True,
                                                       Vector=None, **Vector_kwds):
        '''Compute the intersection of two spheres or circles, each defined by a
           cartesian center point and a radius.

           @arg radius: Radius of this sphere or circle (same units as this point's
                        coordinates).
           @arg center2: Center of the second sphere or circle (C{Cartesian}, L{Vector3d},
                         C{Vector3Tuple} or C{Vector4Tuple}).
           @arg radius2: Radius of the second sphere or circle (same units as this and
                         the B{C{other}} point's coordinates).
           @kwarg sphere: If C{True} compute the center and radius of the intersection
                          of two I{spheres}.  If C{False}, ignore the C{z}-component and
                          compute the intersection of two I{circles} (C{bool}).
           @kwarg Vector: Class to return intersections (C{Cartesian}, L{Vector3d} or
                          C{Vector3Tuple}) or C{None} for an instance of this (sub-)class.
           @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword arguments,
                               ignored if C{B{Vector} is None}.

           @return: If B{C{sphere}} is C{True}, a 2-tuple of the C{center} and C{radius}
                    of the intersection of the I{spheres}.  The C{radius} is C{0.0} for
                    abutting spheres (and the C{center} is aka I{radical center}).

                    If B{C{sphere}} is C{False}, a 2-tuple with the two intersection
                    points of the I{circles}.  For abutting circles, both points are
                    the same instance, aka I{radical center}.

           @raise IntersectionError: Concentric, invalid or non-intersecting spheres or circles.

           @raise TypeError: Invalid B{C{center2}}.

           @raise UnitError: Invalid B{C{radius}} or B{C{radius2}}.

           @see: U{Sphere-Sphere<https://MathWorld.Wolfram.com/Sphere-SphereIntersection.html>},
                 U{Circle-Circle<https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>}
                 Intersection and function L{pygeodesy.radical2}.
        '''
        try:
            return _MODS.vector3d._intersects2(self,    Radius_(radius=radius),
                                               center2, Radius_(radius2=radius2),
                                               sphere=sphere, clas=self.classof,
                                               Vector=Vector, **Vector_kwds)
        except (TypeError, ValueError) as x:
            raise _xError(x, center=self, radius=radius, center2=center2, radius2=radius2)

    @property_doc_(''' this cartesian's reference frame (L{RefFrame}).''')
    def reframe(self):
        '''Get this cartesian's reference frame (L{RefFrame}) or C{None}.
        '''
        return self._reframe

    @reframe.setter  # PYCHOK setter!
    def reframe(self, reframe):
        '''Set or clear this cartesian's reference frame.

           @arg reframe: Reference frame (L{RefFrame}) or C{None}.

           @raise TypeError: The B{C{reframe}} is not a L{RefFrame}.
        '''
        _set_reframe(self, reframe)

    def toRefFrame(self, reframe2, reframe=None, epoch=None):
        '''Convert this cartesian point from one to an other reference frame.

           @arg reframe2: Reference frame to convert I{to} (L{RefFrame}).
           @arg reframe: Reference frame to convert I{from} (L{RefFrame}),
                         overriding this cartesian's C{reframe}.
           @kwarg epoch: Optional epoch to observe (C{scalar}, fractional
                         calendar year), overriding B{C{reframe}}'s epoch.

           @return: The converted point (C{Cartesian}) or this point if
                    conversion is C{nil}.

           @raise TRFError: No conversion available from B{C{reframe}}
                            to B{C{reframe2}} or invalid B{C{epoch}}.

           @raise TypeError: B{C{reframe2}} or B{C{reframe}} not a
                             L{RefFrame}.
        '''
        r = self.reframe if reframe is None else reframe
        if r in (None, reframe2):
            xs = None  # XXX _set_reframe(self, reframe2)?
        else:
            trf = _MODS.trf
            _xinstanceof(trf.RefFrame, reframe2=reframe2, reframe=r)
            _, xs = trf._reframeTransforms2(reframe2, r, epoch)
        return self.toTransforms_(*xs) if xs else self

    def toTransforms_(self, *transforms, **datum):
        '''Apply none, one or several Helmert transforms.

           @arg transforms: Transforms to apply, in order (L{Transform}s).
           @kwarg datum: Datum for the transformed point (L{Datum}),
                         overriding this point's datum.

           @return: The transformed point (C{Cartesian}) or this point
                    if the B{C{transforms}} produce the same point.
        '''
        r = self
        if transforms:
            xyz = r.xyz
            for t in transforms:
                xyz = t.transform(*xyz)
            d = _xkwds_get(datum, datum=r.datum)
            if d != r.datum or xyz != r.xyz:
                r = r.classof(xyz, datum=d)
        return r


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
    _toLLEB_args    = ()      # Etm/Utm/Ups._toLLEB arguments

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
                         year; silently ignored if C{B{reframe} is None}.
           @kwarg name: Optional name (string).

           @raise RangeError: Value of B{C{lat}} or B{C{lon}} outside the valid
                              range and C{rangerrors} set to C{True}.

           @raise TypeError: B{C{datum}} is not a L{Datum}, B{C{reframe}}
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

#   def __matmul__(self, other):  # PYCHOK Python 3.5+
#       '''Return C{NotImplemented} for C{ll_ = ll @ datum} and C{ll_ = ll @ reframe}.
#       '''
#       from pygeodesy.trf import RefFrame
#       return NotImplemented if isinstance(other, (Datum, RefFrame)) else \
#             _NotImplemented(self, other)

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

    @Property_RO
    def _css(self):
        '''(INTERNAL) Get this C{LatLon} point as a Cassini-Soldner location (L{Css}).
        '''
        css = _MODS.css
        return css.toCss(self, height=self.height, Css=css.Css, name=self.name)

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
        return _MODS.elevations.elevation2(self.lat, self.lon,
                                           timeout=self._elevation2to)

    def elevation2(self, adjust=True, datum=None, timeout=2):
        '''Return elevation of this point for its or the given datum, ellipsoid
           or sphere.

           @kwarg adjust: Adjust the elevation for a B{C{datum}} other than
                          C{NAD83} (C{bool}).
           @kwarg datum: Optional datum overriding this point's datum (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2}, L{a_f2Tuple} or C{scalar}
                         radius).
           @kwarg timeout: Optional query timeout (C{seconds}).

           @return: An L{Elevation2Tuple}C{(elevation, data_source)} or
                    C{(None, error)} in case of errors.

           @note: The adjustment applied is the difference in geocentric earth
                  radius between the B{C{datum}} and C{NAV83} upon which the
                  L{elevations.elevation2} is based.

           @note: NED elevation is only available for locations within the
                  U{Conterminous US (CONUS)
                  <https://WikiPedia.org/wiki/Contiguous_United_States>}.

           @see: Function L{elevations.elevation2} and method C{Ellipsoid.Rgeocentric}
                 for further details and possible C{error}s.
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
        except AttributeError:
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
    def Equidistant(self):
        '''Get the prefered azimuthal equidistant projection I{class} (L{EquidistantKarney} or L{EquidistantExact}).
        '''
        try:
            _ = self.datum.ellipsoid.geodesic
            return _MODS.azimuthal.EquidistantKarney
        except ImportError:  # no geographiclib
            return _MODS.azimuthal.EquidistantExact  # XXX no longer L{azimuthal.Equidistant}

    @Property_RO
    def _etm(self):
        '''(INTERNAL) Get this C{LatLon} point as an ETM coordinate (L{pygeodesy.toEtm8}).
        '''
        etm = _MODS.etm
        return etm.toEtm8(self, datum=self.datum, Etm=etm.Etm)

    @Property_RO
    def _geoidHeight2(self):
        '''(INTERNAL) Get geoid height and model.
        '''
        return _MODS.elevations.geoidHeight2(self.lat, self.lon, model=0,
                                             timeout=self._geoidHeight2to)

    def geoidHeight2(self, adjust=False, datum=None, timeout=2):
        '''Return geoid height of this point for its or the given datum, ellipsoid
           or sphere.

           @kwarg adjust: Adjust the geoid height for a B{C{datum}} other than
                          C{NAD83/NADV88} (C{bool}).
           @kwarg datum: Optional datum overriding this point's datum (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2}, L{a_f2Tuple} or C{scalar}
                         radius).
           @kwarg timeout: Optional query timeout (C{seconds}).

           @return: A L{GeoidHeight2Tuple}C{(height, model_name)} or
                    C{(None, error)} in case of errors.

           @note: The adjustment applied is the difference in geocentric earth
                  radius between the B{C{datum}} and C{NAV83/NADV88} upon which
                  the L{elevations.geoidHeight2} is based.

           @note: The geoid height is only available for locations within the
                  U{Conterminous US (CONUS)
                  <https://WikiPedia.org/wiki/Contiguous_United_States>}.

           @see: Function L{elevations.geoidHeight2} and method C{Ellipsoid.Rgeocentric}
                 for further details and possible C{error}s.
        '''
        if self._geoidHeight2to != timeout:
            self._geoidHeight2to = timeout
            LatLonEllipsoidalBase._geoidHeight2._update(self)
        return self._Radjust2(adjust, datum, self._geoidHeight2)

    def intersection3(self, end1, other, end2, height=None, wrap=True,
                                          equidistant=None, tol=_TOL_M):
        '''Interatively compute the intersection point of two paths, each
           defined by two points or a start point and bearing from North.

           @arg end1: End point of this path (C{LatLon}) or the initial
                      bearing at this point (compass C{degrees360}).
           @arg other: Start point of the other path (C{LatLon}).
           @arg end2: End point of the other path (C{LatLon}) or the
                      initial bearing at the other point (compass
                      C{degrees360}).
           @kwarg height: Optional height at the intersection (C{meter},
                          conventionally) or C{None} for the mean height.
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).
           @kwarg equidistant: An azimuthal equidistant projection (I{class} or
                               function L{pygeodesy.equidistant}), or C{None}
                               for this point's preferred C{.Equidistant}.
           @kwarg tol: Tolerance for skew line distance and length and for
                       convergence (C{meter}, conventionally).

           @return: An L{Intersection3Tuple}C{(point, outside1, outside2)}
                    with C{point} a C{LatLon} instance.

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found, but only if
                               C{B{equidistant}=}L{EquidistantKarney}.

           @raise IntersectionError: Skew, colinear, parallel or otherwise
                                     non-intersecting paths or no convergence
                                     for the given B{C{tol}}.

           @raise TypeError: If B{C{end1}}, B{C{other}} or B{C{end2}} point
                             is not C{LatLon}.

           @note: For each path specified with an initial bearing, a pseudo-end
                  point is computed as the C{destination} along that bearing at
                  about 1.5 times the distance from the start point to an initial
                  gu-/estimate of the intersection point (and between 1/8 and 3/8
                  of the authalic earth perimeter).

           @see: U{The B{ellipsoidal} case<https://GIS.StackExchange.com/questions/48937/
                 calculating-intersection-of-two-circles>} and U{Karney's paper
                 <https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section B{14. MARITIME
                 BOUNDARIES} for more details about the iteration algorithm.
        '''
        try:
            s2 = self.others(other)
            return _MODS.ellipsoidalBaseDI._intersect3(self, end1,
                                                       s2,   end2,
                                                       height=height, wrap=wrap,
                                                       equidistant=equidistant, tol=tol,
                                                       LatLon=self.classof, datum=self.datum)
        except (TypeError, ValueError) as x:
            raise _xError(x, start1=self, end1=end1, other=other, end2=end2,
                                          height=height, wrap=wrap, tol=tol)

    def intersections2(self, radius1, other, radius2, height=None, wrap=True,
                                                 equidistant=None, tol=_TOL_M):
        '''Interatively compute the intersection points of two circles,
           each defined by a center point and a radius.

           @arg radius1: Radius of this circle (C{meter}, conventionally).
           @arg other: Center of the other circle (C{LatLon}).
           @arg radius2: Radius of the other circle (C{meter}, same units as
                         B{C{radius1}}).
           @kwarg height: Optional height for the intersection points (C{meter},
                          conventionally) or C{None} for the I{"radical height"}
                          at the I{radical line} between both centers.
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).
           @kwarg equidistant: An azimuthal equidistant projection (I{class} or
                               function L{pygeodesy.equidistant}), or C{None}
                               for this point's preferred C{.Equidistant}.
           @kwarg tol: Convergence tolerance (C{meter}, same units as
                       B{C{radius1}} and B{C{radius2}}).

           @return: 2-Tuple of the intersection points, each a C{LatLon}
                    instance.  For abutting circles, both intersection
                    points are the same instance, aka I{radical center}.

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found, but only if
                               C{B{equidistant}=}L{EquidistantKarney}.

           @raise IntersectionError: Concentric, antipodal, invalid or
                                     non-intersecting circles or no
                                     convergence for the given B{C{tol}}.

           @raise TypeError: Invalid B{C{other}} or B{C{equidistant}}.

           @raise UnitError: Invalid B{C{radius1}}, B{C{radius2}} or B{C{height}}.

           @see: U{The B{ellipsoidal} case<https://GIS.StackExchange.com/questions/48937/
                 calculating-intersection-of-two-circles>}, U{Karney's paper
                 <https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section B{14. MARITIME BOUNDARIES},
                 U{circle-circle<https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>} and
                 U{sphere-sphere<https://MathWorld.Wolfram.com/Sphere-SphereIntersection.html>}
                 intersections.
        '''
        try:
            c2 = self.others(other)
            return _MODS.ellipsoidalBaseDI._intersections2(self, radius1,
                                                           c2,   radius2,
                                                           height=height, wrap=wrap,
                                                           equidistant=equidistant, tol=tol,
                                                           LatLon=self.classof, datum=self.datum)
        except (AssertionError, TypeError, ValueError) as x:
            raise _xError(x, center=self, radius1=radius1, other=other, radius2=radius2,
                                          height=height, wrap=wrap, tol=tol)

    @property_RO
    def iteration(self):
        '''Get the most recent C{intersections2} or C{nearestOn} iteration
           number (C{int}) or C{None} if not available/applicable.
        '''
        return self._iteration

    @Property_RO
    def _lcc(self):
        '''(INTERNAL) Get this C{LatLon} point as a Lambert location (L{Lcc}).
        '''
        lcc = _MODS.lcc
        return lcc.toLcc(self, height=self.height, Lcc=lcc.Lcc, name=self.name)

    def nearestOn(self, point1, point2, within=True, height=None, wrap=True,
                                        equidistant=None, tol=_TOL_M):
        '''Interatively locate the closest point on the geodesic between
           two other (ellipsoidal) points.

           @arg point1: Start point (C{LatLon}).
           @arg point2: End point (C{LatLon}).
           @kwarg within: If C{True} return the closest point I{between}
                          B{C{point1}} and B{C{point2}}, otherwise the
                          closest point elsewhere on the geodesic (C{bool}).
           @kwarg height: Optional height for the closest point (C{meter},
                          conventionally) or C{None} or C{False} for the
                          interpolated height.  If C{False}, the closest
                          takes the heights of the points into account.
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).
           @kwarg equidistant: An azimuthal equidistant projection (I{class} or
                               function L{pygeodesy.equidistant}), or C{None}
                               for this point's preferred C{.Equidistant}.
           @kwarg tol: Convergence tolerance (C{meter}, conventionally).

           @return: Closest point (C{LatLon}).

           @raise ImportError: Package U{geographiclib
                               <https://PyPI.org/project/geographiclib>}
                               not installed or not found, but only if
                               C{B{equidistant}=}L{EquidistantKarney}.

           @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or
                             B{C{equidistant}}.

           @raise ValueError: Datum or ellipsoid of B{C{point1}} or B{C{point2}} is
                              incompatible or no convergence for the given B{C{tol}}.

           @see: U{The B{ellipsoidal} case<https://GIS.StackExchange.com/questions/48937/
                 calculating-intersection-of-two-circles>} and U{Karney's paper
                 <https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section B{14. MARITIME
                 BOUNDARIES} for details about the iteration algorithm.
        '''
        try:
            t = _MODS.ellipsoidalBaseDI._nearestOn2(self, point1, point2, within=within,
                                                          height=height, wrap=wrap,
                                                          equidistant=equidistant,
                                                          tol=tol, LatLon=self.classof)
        except (TypeError, ValueError) as x:
            raise _xError(x, point=self, point1=point1, point2=point2, within=within,
                                         height=height, wrap=wrap, tol=tol)
        return t.closest

    @Property_RO
    def _osgr(self):
        '''(INTERNAL) Get this C{LatLon} point to an OSGR coordinate (L{Osgr}).
        '''
        osgr = _MODS.osgr
        return osgr.toOsgr(self, datum=self.datum, Osgr=osgr.Osgr, name=self.name)

    def parse(self, strllh, height=0, datum=None, epoch=None, reframe=None,
                                                  sep=_COMMA_, name=NN):
        '''Parse a string representing a similar, ellipsoidal C{LatLon}
           point, consisting of C{"lat, lon[, height]"}.

           @arg strllh: Lat, lon and optional height (C{str}),
                        see function L{pygeodesy.parse3llh}.
           @kwarg height: Optional, default height (C{meter} or
                          C{None}).
           @kwarg datum: Optional datum (L{Datum}), overriding this
                         datum I{without conversion}.
           @kwarg epoch: Optional datum (L{Epoch}), overriding this
                         epoch I{without conversion}.
           @kwarg reframe: Optional datum (L{RefFrame}), overriding
                           this reframe I{without conversion}.
           @kwarg sep: Optional separator (C{str}).
           @kwarg name: Optional instance name (C{str}), overriding
                        this name.

           @return: The similar point (ellipsoidal C{LatLon}).

           @raise ParseError: Invalid B{C{strllh}}.
        '''
        a, b, h = _MODS.dms.parse3llh(strllh, height=height, sep=sep)
        r = self.classof(a, b, height=h, datum=self.datum)
        if datum not in (None, self.datum):
            r.datum = datum
        if epoch not in (None, self.epoch):
            r.epoch = epoch
        if reframe not in (None, self.reframe):
            r.reframe = reframe
        return self._xnamed(r, name=name, force=True) if name else r

    def _Radjust2(self, adjust, datum, meter_text2):
        '''(INTERNAL) Adjust an C{elevation} or C{geoidHeight} with
           difference in Gaussian radii of curvature of the given
           datum and NAD83 ellipsoids at this point's latitude.

           @note: This is an arbitrary, possibly incorrect adjustment.
        '''
        if adjust:  # Elevation2Tuple or GeoidHeight2Tuple
            m, t = meter_text2
            if isinstance(m, float) and abs(m) > EPS:
                n = Datums.NAD83.ellipsoid.rocGauss(self.lat)
                if n > EPS0:
                    # use ratio, datum and NAD83 units may differ
                    E = self.ellipsoid() if datum in (None, self.datum) else \
                       _spherical_datum(datum).ellipsoid
                    r = E.rocGauss(self.lat)
                    if r > EPS0 and abs(r - n) > EPS:  # EPS1
                        m *= r / n
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
        _set_reframe(self, reframe)

    @Property_RO
    def scale(self):
        '''Get this point's UTM grid or UPS point scale factor (C{float})
           or C{None} if not converted from L{Utm} or L{Ups}.
        '''
        return self._scale

    def toCss(self, **toCss_kwds):
        '''Convert this C{LatLon} point to a Cassini-Soldner location.

           @kwarg toCss_kwds: Optional L{pygeodesy.toCss} keyword arguments.

           @return: The Cassini-Soldner location (L{Css}).

           @see: Function L{pygeodesy.toCss}.
        '''
        return self._css if not toCss_kwds else _MODS.css.toCss(
               self, **_xkwds(toCss_kwds, name=self.name))

    def toDatum(self, datum2, height=None, name=NN):
        '''Convert this point to an other datum.

           @arg datum2: Datum to convert I{to} (L{Datum}).
           @kwarg height: Optional height, overriding the
                          converted height (C{meter}).
           @kwarg name: Optional name (C{str}), iff converted.

           @return: The converted point (ellipsoidal C{LatLon})
                    or a copy of this point if B{C{datum2}}
                    matches this point's C{datum}.

           @raise TypeError: Invalid B{C{datum2}}.

           @example:

            >>> p = LatLon(51.4778, -0.0016)  # default Datums.WGS84
            >>> p.toDatum(Datums.OSGB36)  # 51.477284°N, 000.00002°E
        '''
        n  =  name or self.name
        d2 = _ellipsoidal_datum(datum2, name=n)
        if self.datum == d2:
            r = self.copy(name=name)
        else:
            kwds = _xkwds_not(None, LatLon=self.classof, name=n,
                                    epoch=self.epoch, reframe=self.reframe)
            c = self.toCartesian().toDatum(d2)
            r = c.toLatLon(datum=d2, height=height, **kwds)
        return r

    def toEtm(self, **toEtm8_kwds):
        '''Convert this C{LatLon} point to an ETM coordinate.

           @kwarg toEtm8_kwds: Optional L{pygeodesy.toEtm8} keyword arguments.

           @return: The ETM coordinate (L{Etm}).

           @see: Function L{pygeodesy.toEtm8}.
        '''
        return self._etm if not toEtm8_kwds else _MODS.etm.toEtm8(
               self, **_xkwds(toEtm8_kwds, name=self.name))

    def toLcc(self, **toLcc_kwds):
        '''Convert this C{LatLon} point to a Lambert location.

           @kwarg toLcc_kwds: Optional L{pygeodesy.toLcc} keyword arguments.

           @return: The Lambert location (L{Lcc}).

           @see: Function L{pygeodesy.toLcc}.
        '''
        return self._lcc if not toLcc_kwds else _MODS.lcc.toLcc(
               self, **_xkwds(toLcc_kwds, name=self.name))

    def toMgrs(self, center=False):
        '''Convert this C{LatLon} point to an MGRS coordinate.

           @kwarg center: If C{True}, try to I{un}-center MGRS
                          to its C{lowerleft} (C{bool}) or by
                          C{B{center} meter} (C{scalar}).

           @return: The MGRS coordinate (L{Mgrs}).

           @see: Method L{toUtm} and L{Mgrs.toLatLon}.
        '''
        return self.toUtm(center=center).toMgrs(center=False)

    def toOsgr(self, **toOsgr_kwds):
        '''Convert this C{LatLon} point to an OSGR coordinate.

           @kwarg toOsgr_kwds: Optional L{pygeodesy.toOsgr} keyword arguments.

           @return: The OSGR coordinate (L{Osgr}).

           @see: Function L{pygeodesy.toOsgr}.
        '''
        return self._osgr if not toOsgr_kwds else _MODS.osgr.toOsgr(
               self, **_xkwds(toOsgr_kwds, name=self.name))

    def toRefFrame(self, reframe2, height=None, name=NN):
        '''Convert this point to an other reference frame.

           @arg reframe2: Reference frame to convert I{to} (L{RefFrame}).
           @kwarg height: Optional height, overriding the converted
                          height (C{meter}).
           @kwarg name: Optional name (C{str}), iff converted.

           @return: The converted point (ellipsoidal C{LatLon}) or this
                    point if conversion is C{nil}, or a copy of this
                    point if the B{C{name}} is non-empty.

           @raise TRFError: This point's C{reframe} is not defined or
                            conversion from this point's C{reframe} to
                            B{C{reframe2}} is not available.

           @raise TypeError: Invalid B{C{reframe2}}, not a L{RefFrame}.

           @example:

            >>> p = LatLon(51.4778, -0.0016, reframe=RefFrames.ETRF2000)  # default Datums.WGS84
            >>> p.toRefFrame(RefFrames.ITRF2014)  # 51.477803°N, 000.001597°W, +0.01m
            >>> p.toRefFrame(RefFrames.ITRF2014, height=0)  # 51.477803°N, 000.001597°W
        '''
        if not self.reframe:
            t = _SPACE_(_DOT_(repr(self), _reframe_), MISSING)
            raise TRFError(_no_(_conversion_), txt=t)

        trf = _MODS.trf
        trf._xinstanceof(trf.RefFrame, reframe2=reframe2)

        e, xs = trf._reframeTransforms2(reframe2, self.reframe, self.epoch)
        if xs:
            c = self.toCartesian().toTransforms_(*xs)
            n = name or self.name
            ll = c.toLatLon(datum=self.datum, epoch=e, height=height,
                            LatLon=self.classof, name=n, reframe=reframe2)
        else:
            ll = self.copy(name=name) if name else self
        return ll

    def toUps(self, pole=_N_, falsed=True):
        '''Convert this C{LatLon} point to a UPS coordinate.

           @kwarg pole: Optional top/center of (stereographic)
                        projection (C{str}, 'N[orth]' or 'S[outh]').
           @kwarg falsed: False easting and northing (C{bool}).

           @return: The UPS coordinate (L{Ups}).

           @see: Function L{pygeodesy.toUps8}.
        '''
        if self._upsOK(pole, falsed):
            u = self._ups
        else:
            ups = _MODS.ups
            u = ups.toUps8(self, datum=self.datum, Ups=ups.Ups,
                                 pole=pole, falsed=falsed)
        return u

    def toUtm(self, center=False):
        '''Convert this C{LatLon} point to a UTM coordinate.

           @kwarg center: If C{True}, I{un}-center the UTM
                          to its C{lowerleft} (C{bool}) or
                          by C{B{center} meter} (C{scalar}).

           @return: The UTM coordinate (L{Utm}).

           @see: Method L{Mgrs.toUtm} and function L{pygeodesy.toUtm8}.
        '''
        if center in (False, 0, _0_0):
            u = self._utm
        elif center in (True,):
            u = self._utm._lowerleft
        else:
            u = _MODS.utm._lowerleft(self._utm, center)
        return u

    def toUtmUps(self, pole=NN):
        '''Convert this C{LatLon} point to a UTM or UPS coordinate.

           @kwarg pole: Optional top/center of UPS (stereographic)
                        projection (C{str}, 'N[orth]' or 'S[outh]').

           @return: The UTM or UPS coordinate (L{Utm} or L{Ups}).

           @raise TypeError: Result in L{Utm} or L{Ups}.

           @see: Function L{pygeodesy.toUtmUps8}.
        '''
        if self._utmOK():
            u = self._utm
        elif self._upsOK(pole):
            u = self._ups
        else:  # no cover
            from pygeodesy.utmups import toUtmUps8, Utm, Ups  # utmups = _MODS.utmups
            u = toUtmUps8(self, datum=self.datum, Utm=Utm, Ups=Ups,
                                pole=pole, name=self.name)
            if isinstance(u, Utm):
                self._overwrite(_utm=u)
            elif isinstance(u, Ups):
                self._overwrite(_ups=u)
            else:
                _xinstanceof(Utm, Ups, toUtmUps8=u)
        return u

    def toWm(self, **toWm_kwds):
        '''Convert this C{LatLon} point to a WM coordinate.

           @kwarg toWm_kwds: Optional L{pygeodesy.toWm} keyword arguments.

           @return: The WM coordinate (L{Wm}).

           @see: Function L{pygeodesy.toWm}.
        '''
        return self._wm if not toWm_kwds else _MODS.webmercator.toWm(
               self, **_xkwds(toWm_kwds, name=self.name))

    @deprecated_method
    def to3xyz(self):  # PYCHOK no cover
        '''DEPRECATED, use method C{toEcef}.

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
                  and C{LatLon.nearestOn} based on I{Karney}'s Python U{geographiclib
                  <https://PyPI.org/project/geographiclib>} if installed, otherwise
                  uses the accurate (but slower) C{ellipsoidalExact.LatLon} methods.
        '''
        return _trilaterate5(self, distance1,
                             self.others(point2=point2), distance2,
                             self.others(point3=point3), distance3,
                             area=area, eps=eps, wrap=wrap)

    @Property_RO
    def _ups(self):  # __dict__ value overwritten by method C{toUtmUps}
        '''(INTERNAL) Get this C{LatLon} point as UPS coordinate (L{Ups}),
           see L{pygeodesy.toUps8}.
        '''
        ups = _MODS.ups
        return ups.toUps8(self, datum=self.datum, Ups=ups.Ups,
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
        '''(INTERNAL) Get this C{LatLon} point as UTM coordinate (L{Utm}),
           see L{pygeodesy.toUtm8}.
        '''
        utm = _MODS.utm
        return utm.toUtm8(self, datum=self.datum, Utm=utm.Utm, name=self.name)

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
        '''(INTERNAL) Get this C{LatLon} point as webmercator (L{Wm}).
        '''
        return _MODS.webmercator.toWm(self)


def _nearestOn(point, point1, point2, within=True, height=None, wrap=True,
               equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    '''(INTERNAL) Get closest point, imported by .ellipsoidalExact,
       -GeodSolve, -Karney and -Vincenty to embellish exceptions.
    '''
    try:
        p = _xellipsoidal(point=point)
        t = _MODS.ellipsoidalBaseDI._nearestOn2(p, point1, point2, within=within,
                                                   height=height, wrap=wrap,
                                                   equidistant=equidistant, tol=tol,
                                                   LatLon=LatLon, **LatLon_kwds)
    except (TypeError, ValueError) as x:
        raise _xError(x, point=point, point1=point1, point2=point2)
    return t.closest


def _set_reframe(inst, reframe):
    '''(INTERNAL) Set or clear an instance's reference frame.
    '''
    if reframe is not None:
        _xinstanceof(_MODS.trf.RefFrame, reframe=reframe)
        inst._reframe = reframe
    elif inst.reframe is not None:
        inst._reframe = None


__all__ += _ALL_DOCS(CartesianEllipsoidalBase, LatLonEllipsoidalBase)

# **) MIT License
#
# Copyright (C) 2016-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
