
# -*- coding: utf-8 -*-

u'''(INTERNAL) Base class L{LatLonBase} for all elliposiodal, spherical and N-vectorial C{LatLon} classes.

@see: I{(C) Chris Veness}' U{latlong<https://www.Movable-Type.co.UK/scripts/latlong.html>},
      U{-ellipsoidal<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>} and
      U{-vectors<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>} and I{Charles Karney}'s
      U{Rhumb<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Rhumb.html>} and
      U{RhumbLine<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1RhumbLine.html>} classes.
'''

from pygeodesy.basics import isscalar, isstr, map1, _xinstanceof
from pygeodesy.constants import EPS, EPS0, EPS1, EPS4, INT0, R_M, \
                               _0_0, _0_5, _1_0
# from pygeodesy.datums import _spherical_datum  # from .formy
from pygeodesy.dms import F_D, F_DMS, latDMS, lonDMS, parse3llh
# from pygeodesy.ecef import EcefKarney  # _MODS
from pygeodesy.errors import _incompatible, IntersectionError, _IsnotError, \
                             _TypeError, _ValueError, _xdatum, _xError, \
                             _xkwds, _xkwds_not
# from pygeodesy.fmath import favg  # _MODS
from pygeodesy.formy import antipode, compassAngle, cosineAndoyerLambert_, \
                            cosineForsytheAndoyerLambert_, cosineLaw, \
                            equirectangular, euclidean, flatLocal_, \
                            flatPolar, hartzell, haversine, isantipode, \
                            _isequalTo, isnormal, normal, philam2n_xyz, \
                            thomas_, vincentys,  _spherical_datum
from pygeodesy.interns import NN, _COMMASPACE_, _concentric_, _height_, \
                             _intersection_, _m_, _LatLon_, _no_, \
                             _overlap_,  _point_  # PYCHOK used!
# from pygeodesy.iters import PointsIter, points2  # from .vector3d, _MODS
# from pygeodesy.karney import Caps  # _MODS
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
# from pygeodesy.ltp import Ltp, _xLtp  # _MODS
from pygeodesy.named import _NamedBase, notOverloaded,  Fmt
from pygeodesy.namedTuples import Bounds2Tuple, LatLon2Tuple, PhiLam2Tuple, \
                                  Trilaterate5Tuple, Vector3Tuple
# from pygeodesy.nvectorBase import _N_vector_  # _MODS
from pygeodesy.props import deprecated_method, Property, Property_RO, \
                            property_RO, _update_all
# from pygeodesy.streprs import Fmt, hstr  # from .named, _MODS
from pygeodesy.units import Distance_, Lat, Lon, Height, Radius, Radius_, \
                            Scalar, Scalar_
from pygeodesy.utily import _unrollon, _unrollon3, _Wrap
from pygeodesy.vector2d import _circin6,  Circin6Tuple, _circum3, circum4_, \
                                Circum3Tuple, _radii11ABC
from pygeodesy.vector3d import nearestOn6, Vector3d,  PointsIter

from contextlib import contextmanager
from math import asin, cos, degrees, fabs, radians

__all__ = _ALL_LAZY.latlonBase
__version__ = '23.09.09'


class LatLonBase(_NamedBase):
    '''(INTERNAL) Base class for C{LatLon} points on spherical or
       ellipsoidal earth models.
    '''
    _clipid = INT0  # polygonal clip, see .booleans
    _datum  = None  # L{Datum}, to be overriden
    _height = INT0  # height (C{meter}), default
    _lat    = 0     # latitude (C{degrees})
    _lon    = 0     # longitude (C{degrees})

    def __init__(self, latlonh, lon=None, height=0, wrap=False, name=NN):
        '''New C{LatLon}.

           @arg latlonh: Latitude (C{degrees} or DMS C{str} with N or S suffix) or
                         a previous C{LatLon} instance provided C{B{lon}=None}.
           @kwarg lon: Longitude (C{degrees} or DMS C{str} with E or W suffix) or
                       C(None), indicating B{C{latlonh}} is a C{LatLon}.
           @kwarg height: Optional height above (or below) the earth surface
                          (C{meter}, conventionally).
           @kwarg wrap: If C{True}, wrap or I{normalize} B{C{lat}} and B{C{lon}}
                        (C{bool}).
           @kwarg name: Optional name (C{str}).

           @return: New instance (C{LatLon}).

           @raise RangeError: A B{C{lon}} or C{lat} value outside the valid
                              range and L{rangerrors} set to C{True}.

           @raise TypeError: If B{C{latlonh}} is not a C{LatLon}.

           @raise UnitError: Invalid B{C{lat}}, B{C{lon}} or B{C{height}}.

           @example:

            >>> p = LatLon(50.06632, -5.71475)
            >>> q = LatLon('50°03′59″N', """005°42'53"W""")
            >>> r = LatLon(p)
        '''
        if name:
            self.name = name

        if lon is None:
            try:
                lat, lon = latlonh.lat, latlonh.lon
                height = latlonh.get(_height_, height)
            except AttributeError:
                raise _IsnotError(_LatLon_, latlonh=latlonh)
            if wrap:
                lat, lon = _Wrap.latlon(lat, lon)
        elif wrap:
            lat, lon = _Wrap.latlonDMS2(latlonh, lon)
        else:
            lat = latlonh

        self._lat = Lat(lat)  # parseDMS2(lat, lon)
        self._lon = Lon(lon)  # PYCHOK LatLon2Tuple
        if height:  # elevation
            self._height = Height(height)

    def __eq__(self, other):
        return self.isequalTo(other)

    def __ne__(self, other):
        return not self.isequalTo(other)

    def __str__(self):
        return self.toStr(form=F_D, prec=6)

    def antipode(self, height=None):
        '''Return the antipode, the point diametrically opposite
           to this point.

           @kwarg height: Optional height of the antipode (C{meter}),
                          this point's height otherwise.

           @return: The antipodal point (C{LatLon}).
        '''
        h = self._heigHt(height)
        return self.classof(*antipode(*self.latlon), height=h)

    @deprecated_method
    def bounds(self, wide, tall, radius=R_M):  # PYCHOK no cover
        '''DEPRECATED, use method C{boundsOf}.'''
        return self.boundsOf(wide, tall, radius=radius)

    def boundsOf(self, wide, tall, radius=R_M, height=None):
        '''Return the SW and NE lat-/longitude of a great circle
           bounding box centered at this location.

           @arg wide: Longitudinal box width (C{meter}, same units as
                      B{C{radius}} or C{degrees} if B{C{radius}} is C{None}).
           @arg tall: Latitudinal box size (C{meter}, same units as
                      B{C{radius}} or C{degrees} if B{C{radius}} is C{None}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None} if I{both}
                          B{C{wide}} and B{C{tall}} are in C{degrees}.
           @kwarg height: Height for C{latlonSW} and C{latlonNE} (C{meter}),
                          overriding the point's height.

           @return: A L{Bounds2Tuple}C{(latlonSW, latlonNE)}, the
                    lower-left and upper-right corner (C{LatLon}).

           @see: U{https://www.Movable-Type.co.UK/scripts/latlong-db.html}
        '''
        w = Scalar_(wide=wide) * _0_5
        t = Scalar_(tall=tall) * _0_5
        if radius is not None:
            r = Radius_(radius)
            c = cos(self.phi)
            w = degrees(asin(w / r) / c) if fabs(c) > EPS0 else _0_0  # XXX
            t = degrees(t / r)
        y, t = self.lat, fabs(t)
        x, w = self.lon, fabs(w)

        h  = self._heigHt(height)
        sw = self.classof(y - t, x - w, height=h)
        ne = self.classof(y + t, x + w, height=h)
        return Bounds2Tuple(sw, ne, name=self.name)

    def chordTo(self, other, height=None, wrap=False):
        '''Compute the length of the chord through the earth between
           this and an other point.

           @arg other: The other point (C{LatLon}).
           @kwarg height: Overriding height for both points (C{meter})
                          or C{None} for each point's height.
           @kwarg wrap: If C{True}, wrap or I{normalize} the B{C{other}}
                        point (C{bool}).

           @return: The chord length (conventionally C{meter}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.
        '''
        def _v3d(ll):
            t = ll.toEcef(height=height)  # .toVector(Vector=Vector3d)
            return Vector3d(t.x, t.y, t.z)

        p = self.others(other)
        if wrap:
            p = _Wrap.point(p)
        return _v3d(self).minus(_v3d(p)).length

    def circin6(self, point2, point3, eps=EPS4, wrap=False):
        '''Return the radius and center of the I{inscribed} aka I{In-}circle
           of the (planar) triangle formed by this and two other points.

           @arg point2: Second point (C{LatLon}).
           @arg point3: Third point (C{LatLon}).
           @kwarg eps: Tolerance for function L{pygeodesy.trilaterate3d2}.
           @kwarg wrap: If C{True}, wrap or I{normalize} B{C{point2}} and
                        B{C{point3}} (C{bool}).

           @return: L{Circin6Tuple}C{(radius, center, deltas, cA, cB, cC)}.  The
                    C{center} and contact points C{cA}, C{cB} and C{cC}, each an
                    instance of this (sub-)class, are co-planar with this and the
                    two given points, see the B{Note} below.

           @raise ImportError: Package C{numpy} not found, not installed or older
                               than version 1.10.

           @raise IntersectionError: Near-coincident or -colinear points or
                                     a trilateration or C{numpy} issue.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @note: The C{center} is trilaterated in cartesian (ECEF) space and converted
                  back to geodetic lat-, longitude and height.  The latter, conventionally
                  in C{meter} indicates whether the C{center} is above, below or on the
                  surface of the earth model.  If C{deltas} is C{None}, the C{center} is
                  I{un}ambigous.  Otherwise C{deltas} is a L{LatLon3Tuple}C{(lat, lon,
                  height)} representing the differences between both results from
                  L{pygeodesy.trilaterate3d2} and C{center} is the mean thereof.

           @see: Function L{pygeodesy.circin6}, method L{circum3}, U{Incircle
                 <https://MathWorld.Wolfram.com/Incircle.html>} and U{Contact Triangle
                 <https://MathWorld.Wolfram.com/ContactTriangle.html>}.
        '''
        with _toCartesian3(self, point2, point3, wrap) as cs:
            r, c, d, cA, cB, cC = _circin6(*cs, eps=eps, useZ=True, dLL3=True,
                                                datum=self.datum)  # PYCHOK unpack
            return Circin6Tuple(r, c.toLatLon(), d, cA.toLatLon(), cB.toLatLon(), cC.toLatLon())

    def circum3(self, point2, point3, circum=True, eps=EPS4, wrap=False):
        '''Return the radius and center of the smallest circle I{through} or I{containing}
           this and two other points.

           @arg point2: Second point (C{LatLon}).
           @arg point3: Third point (C{LatLon}).
           @kwarg circum: If C{True} return the C{circumradius} and C{circumcenter},
                          always, ignoring the I{Meeus}' Type I case (C{bool}).
           @kwarg eps: Tolerance for function L{pygeodesy.trilaterate3d2}.
           @kwarg wrap: If C{True}, wrap or I{normalize} B{C{point2}} and
                        B{C{point3}} (C{bool}).

           @return: A L{Circum3Tuple}C{(radius, center, deltas)}.  The C{center}, an
                    instance of this (sub-)class, is co-planar with this and the two
                    given points.  If C{deltas} is C{None}, the C{center} is
                    I{un}ambigous.  Otherwise C{deltas} is a L{LatLon3Tuple}C{(lat,
                    lon, height)} representing the difference between both results
                    from L{pygeodesy.trilaterate3d2} and C{center} is the mean thereof.

           @raise ImportError: Package C{numpy} not found, not installed or older than
                               version 1.10.

           @raise IntersectionError: Near-concentric, -coincident or -colinear points,
                                     incompatible C{Ecef} classes or a trilateration
                                     or C{numpy} issue.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @note: The C{center} is trilaterated in cartesian (ECEF) space and converted
                  back to geodetic lat-, longitude and height.  The latter, conventionally
                  in C{meter} indicates whether the C{center} is above, below or on the
                  surface of the earth model.  If C{deltas} is C{None}, the C{center} is
                  I{un}ambigous.  Otherwise C{deltas} is a L{LatLon3Tuple}C{(lat, lon,
                  height)} representing the difference between both results from
                  L{pygeodesy.trilaterate3d2} and C{center} is the mean thereof.

           @see: Function L{pygeodesy.circum3} and methods L{circin6} and L{circum4_}.
        '''
        with _toCartesian3(self, point2, point3, wrap, circum=circum) as cs:
            r, c, d = _circum3(*cs, circum=circum, eps=eps, useZ=True, dLL3=True,  # XXX -3d2
                                    clas=cs[0].classof, datum=self.datum)  # PYCHOK unpack
            return Circum3Tuple(r, c.toLatLon(), d)

    def circum4_(self, *points, **wrap):
        '''Best-fit a sphere through this and two or more other points.

           @arg points: The other points (each a C{LatLon}).
           @kwarg wrap: If C{True}, wrap or I{normalize} the B{C{points}}
                        (C{bool}), default C{False}.

           @return: L{Circum4Tuple}C{(radius, center, rank, residuals)} with C{center}
                    an instance of this (sub-)class.

           @raise ImportError: Package C{numpy} not found, not installed or older than
                               version 1.10.

           @raise NumPyError: Some C{numpy} issue.

           @raise TypeError: One of the B{C{points}} invalid.

           @raise ValueError: Too few B{C{points}}.

           @see: Function L{pygeodesy.circum4_} and L{circum3}.
        '''
        def _cs(ps, C, wrap=False):
            _wp = _Wrap.point if wrap else (lambda p: p)
            for i, p in enumerate(ps):
                yield C(i=i, points=_wp(p))

        C = self._toCartesianEcef
        c = C(point=self)
        t = circum4_(c, Vector=c.classof, *_cs(points, C, **wrap))
        c = t.center.toLatLon(LatLon=self.classof)
        return t.dup(center=c)

    @property
    def clipid(self):
        '''Get the (polygonal) clip (C{int}).
        '''
        return self._clipid

    @clipid.setter  # PYCHOK setter!
    def clipid(self, clipid):
        '''Get the (polygonal) clip (C{int}).
        '''
        self._clipid = int(clipid)

    @deprecated_method
    def compassAngle(self, other, **adjust_wrap):  # PYCHOK no cover
        '''DEPRECATED, use method L{compassAngleTo}.'''
        return self.compassAngleTo(other, **adjust_wrap)

    def compassAngleTo(self, other, **adjust_wrap):
        '''Return the angle from North for the direction vector between
           this and an other point.

           Suitable only for short, non-near-polar vectors up to a few
           hundred Km or Miles.  Use method C{initialBearingTo} for
           larger distances.

           @arg other: The other point (C{LatLon}).
           @kwarg adjust_wrap: Optional keyword arguments for function
                               L{pygeodesy.compassAngle}.

           @return: Compass angle from North (C{degrees360}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @note: Courtesy of Martin Schultz.

           @see: U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>}.
        '''
        p = self.others(other)
        return compassAngle(self.lat, self.lon, p.lat, p.lon, **adjust_wrap)

    def cosineAndoyerLambertTo(self, other, wrap=False):
        '''Compute the distance between this and an other point using the U{Andoyer-Lambert correction<https://
           navlib.net/wp-content/uploads/2013/10/admiralty-manual-of-navigation-vol-1-1964-english501c.pdf>}
           of the U{Law of Cosines<https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>} formula.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{other}} point (C{bool}).

           @return: Distance (C{meter}, same units as the axes of this
                    point's datum ellipsoid).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.cosineAndoyerLambert} and methods
                 L{cosineForsytheAndoyerLambertTo}, L{cosineLawTo},
                 C{distanceTo*}, L{equirectangularTo}, L{euclideanTo},
                 L{flatLocalTo}/L{hubenyTo}, L{flatPolarTo}, L{haversineTo},
                 L{thomasTo} and L{vincentysTo}.
        '''
        return self._distanceTo_(cosineAndoyerLambert_, other, wrap=wrap)

    def cosineForsytheAndoyerLambertTo(self, other, wrap=False):
        '''Compute the distance between this and an other point using
           the U{Forsythe-Andoyer-Lambert correction
           <https://www2.UNB.Ca/gge/Pubs/TR77.pdf>} of the U{Law of Cosines
           <https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
           formula.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{other}} point (C{bool}).

           @return: Distance (C{meter}, same units as the axes of
                    this point's datum ellipsoid).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.cosineForsytheAndoyerLambert} and methods
                 L{cosineAndoyerLambertTo}, L{cosineLawTo}, C{distanceTo*},
                 L{equirectangularTo}, L{euclideanTo}, L{flatLocalTo}/L{hubenyTo},
                 L{flatPolarTo}, L{haversineTo}, L{thomasTo} and L{vincentysTo}.
        '''
        return self._distanceTo_(cosineForsytheAndoyerLambert_, other, wrap=wrap)

    def cosineLawTo(self, other, radius=None, wrap=False):
        '''Compute the distance between this and an other point using the
           U{spherical Law of Cosines
           <https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
           formula.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None}
                          for the mean radius of this point's datum
                          ellipsoid.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{other}} point (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.cosineLaw} and methods L{cosineAndoyerLambertTo},
                 L{cosineForsytheAndoyerLambertTo}, C{distanceTo*}, L{equirectangularTo},
                 L{euclideanTo}, L{flatLocalTo}/L{hubenyTo}, L{flatPolarTo},
                 L{haversineTo}, L{thomasTo} and L{vincentysTo}.
        '''
        return self._distanceTo(cosineLaw, other, radius, wrap=wrap)

    @property_RO
    def datum(self):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self)

    def destinationXyz(self, delta, LatLon=None, **LatLon_kwds):
        '''Calculate the destination using a I{local} delta from this point.

           @arg delta: Local delta to the destination (L{XyzLocal}, L{Enu},
                       L{Ned} or L{Local9Tuple}).
           @kwarg LatLon: Optional (geodetic) class to return the destination
                          or C{None}.
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: Destination as a C{B{LatLon}(lat, lon, **B{LatLon_kwds})}
                    instance or if C{B{LatLon} is None}, a L{LatLon3Tuple}C{(lat,
                    lon, height)} respectively L{LatLon4Tuple}C{(lat, lon,
                    height, datum)} depending on whether a C{datum} keyword
                    is un-/specified.

           @raise TypeError: Invalid B{C{delta}}, B{C{LatLon}} or B{C{LatLon_kwds}}.
        '''
        t = self._ltp._local2ecef(delta, nine=True)
        return t.toLatLon(LatLon=LatLon, **_xkwds(LatLon_kwds, name=self.name))

    def _distanceTo(self, func, other, radius=None, **kwds):
        '''(INTERNAL) Helper for distance methods C{<func>To}.
        '''
        p, r = self.others(other, up=2), radius
        if r is None:
            r = self._datum.ellipsoid.R1 if self._datum else R_M
        return func(self.lat, self.lon, p.lat, p.lon, radius=r, **kwds)

    def _distanceTo_(self, func_, other, wrap=False, radius=None):
        '''(INTERNAL) Helper for (ellipsoidal) methods C{<func>To}.
        '''
        p = self.others(other, up=2)
        D = self.datum
        lam21, phi2, _ = _Wrap.philam3(self.lam, p.phi, p.lam, wrap)
        r = func_(phi2, self.phi, lam21, datum=D)
        return r * (D.ellipsoid.a if radius is None else radius)

    @Property_RO
    def Ecef(self):
        '''Get the ECEF I{class} (L{EcefKarney}), I{lazily}.
        '''
        return _MODS.ecef.EcefKarney  # default

    @Property_RO
    def _Ecef_forward(self):
        '''(INTERNAL) Helper for L{_ecef9} and L{toEcef} (C{callable}).
        '''
        return self.Ecef(self.datum, name=self.name).forward

    @Property_RO
    def _ecef9(self):
        '''(INTERNAL) Helper for L{toCartesian}, L{toEcef} and L{toCartesian} (L{Ecef9Tuple}).
        '''
        return self._Ecef_forward(self, M=True)

    @deprecated_method
    def equals(self, other, eps=None):  # PYCHOK no cover
        '''DEPRECATED, use method L{isequalTo}.'''
        return self.isequalTo(other, eps=eps)

    @deprecated_method
    def equals3(self, other, eps=None):  # PYCHOK no cover
        '''DEPRECATED, use method L{isequalTo3}.'''
        return self.isequalTo3(other, eps=eps)

    def equirectangularTo(self, other, **radius_adjust_limit_wrap):
        '''Compute the distance between this and an other point
           using the U{Equirectangular Approximation / Projection
           <https://www.Movable-Type.co.UK/scripts/latlong.html#equirectangular>}.

           Suitable only for short, non-near-polar distances up to a
           few hundred Km or Miles.  Use method L{haversineTo} or
           C{distanceTo*} for more accurate and/or larger distances.

           @arg other: The other point (C{LatLon}).
           @kwarg radius_adjust_limit_wrap: Optional keyword arguments
                         for function L{pygeodesy.equirectangular},
                         overriding the default mean C{radius} of this
                         point's datum ellipsoid.

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.equirectangular} and methods L{cosineAndoyerLambertTo},
                 L{cosineForsytheAndoyerLambertTo}, L{cosineLawTo}, C{distanceTo*},
                 C{euclideanTo}, L{flatLocalTo}/L{hubenyTo}, L{flatPolarTo},
                 L{haversineTo}, L{thomasTo} and L{vincentysTo}.
        '''
        return self._distanceTo(equirectangular, other, **radius_adjust_limit_wrap)

    def euclideanTo(self, other, **radius_adjust_wrap):
        '''Approximate the C{Euclidian} distance between this and
           an other point.

           See function L{pygeodesy.euclidean} for the available B{C{options}}.

           @arg other: The other point (C{LatLon}).
           @kwarg radius_adjust_wrap: Optional keyword arguments for function
                         L{pygeodesy.euclidean}, overriding the default mean
                         C{radius} of this point's datum ellipsoid.

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.euclidean} and methods L{cosineAndoyerLambertTo},
                 L{cosineForsytheAndoyerLambertTo}, L{cosineLawTo}, C{distanceTo*},
                 L{equirectangularTo}, L{flatLocalTo}/L{hubenyTo}, L{flatPolarTo},
                 L{haversineTo}, L{thomasTo} and L{vincentysTo}.
        '''
        return self._distanceTo(euclidean, other, **radius_adjust_wrap)

    def flatLocalTo(self, other, radius=None, wrap=False):
        '''Compute the distance between this and an other point using the
           U{ellipsoidal Earth to plane projection
           <https://WikiPedia.org/wiki/Geographical_distance#Ellipsoidal_Earth_projected_to_a_plane>}
           aka U{Hubeny<https://www.OVG.AT/de/vgi/files/pdf/3781/>} formula.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None} for
                          the I{equatorial radius} of this point's
                          datum ellipsoid.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{other}} point (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: Invalid B{C{radius}}.

           @see: Function L{pygeodesy.flatLocal}/L{pygeodesy.hubeny}, methods
                 L{cosineAndoyerLambertTo}, L{cosineForsytheAndoyerLambertTo},
                 L{cosineLawTo}, C{distanceTo*}, L{equirectangularTo}, L{euclideanTo},
                 L{flatPolarTo}, L{haversineTo}, L{thomasTo} and L{vincentysTo} and
                 U{local, flat Earth approximation<https://www.edwilliams.org/avform.htm#flat>}.
        '''
        return self._distanceTo_(flatLocal_, other, wrap=wrap, radius=
                     radius if radius in (None, R_M, _1_0, 1) else Radius(radius))  # PYCHOK kwargs

    hubenyTo = flatLocalTo  # for Karl Hubeny

    def flatPolarTo(self, other, **radius_wrap):
        '''Compute the distance between this and an other point using
           the U{polar coordinate flat-Earth<https://WikiPedia.org/wiki/
           Geographical_distance#Polar_coordinate_flat-Earth_formula>} formula.

           @arg other: The other point (C{LatLon}).
           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.flatPolar}, overriding the
                               default mean C{radius} of this point's
                               datum ellipsoid.

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.flatPolar} and methods L{cosineAndoyerLambertTo},
                 L{cosineForsytheAndoyerLambertTo}, L{cosineLawTo}, C{distanceTo*},
                 L{equirectangularTo}, L{euclideanTo}, L{flatLocalTo}/L{hubenyTo},
                 L{haversineTo}, L{thomasTo} and L{vincentysTo}.
        '''
        return self._distanceTo(flatPolar, other, **radius_wrap)

    def hartzell(self, los=None, earth=None):
        '''Compute the intersection of a Line-Of-Sight (los) from this Point-Of-View
           (pov) with this point's ellipsoid surface.

           @kwarg los: Line-Of-Sight, I{direction} to earth (L{Vector3d}) or
                       C{None} to point to the ellipsoid's center.
           @kwarg earth: The earth model (L{Datum}, L{Ellipsoid}, L{Ellipsoid2},
                         L{a_f2Tuple} or C{scalar} radius in C{meter}) overriding
                         this point's C{datum} ellipsoid.

           @return: The ellipsoid intersection (C{LatLon}) or this very instance
                    if this C{pov's height} is C{0}.

           @raise IntersectionError: Null C{pov} or B{C{los}} vector, this
                                     C{pov's height} is negative or B{C{los}}
                                     points outside the ellipsoid or in an
                                     opposite direction.

           @raise TypeError: Invalid B{C{los}}.

           @see: Function C{hartzell} for further details.
        '''
        h = self.height
        if not h:
            r = self
        elif h < 0:
            raise IntersectionError(pov=self, los=los, height=h, txt=_no_(_height_))
        elif los is None:
            d = self.datum if earth is None else _spherical_datum(earth)
            r = self.dup(datum=d, height=0, name=self.hartzell.__name__)
        else:
            c = self.toCartesian()
            r = hartzell(c, los=los, earth=earth or self.datum, LatLon=self.classof)
        return r

    def haversineTo(self, other, **radius_wrap):
        '''Compute the distance between this and an other point using the
           U{Haversine<https://www.Movable-Type.co.UK/scripts/latlong.html>}
           formula.

           @arg other: The other point (C{LatLon}).
           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.haversine}, overriding the
                               default mean C{radius} of this point's
                               datum ellipsoid.

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.haversine} and methods L{cosineAndoyerLambertTo},
                 L{cosineForsytheAndoyerLambertTo}, L{cosineLawTo}, C{distanceTo*},
                 L{equirectangularTo}, L{euclideanTo}, L{flatLocalTo}/L{hubenyTo},
                 L{flatPolarTo}, L{thomasTo} and L{vincentysTo}.
        '''
        return self._distanceTo(haversine, other, **radius_wrap)

    def _havg(self, other, f=_0_5, h=None):
        '''(INTERNAL) Weighted, average height.

           @arg other: An other point (C{LatLon}).
           @kwarg f: Optional fraction (C{float}).
           @kwarg h: Overriding height (C{meter}).

           @return: Average, fractional height (C{float}) or
                    the overriding B{C{height}} (C{Height}).
        '''
        return Height(h) if h is not None else \
              _MODS.fmath.favg(self.height, other.height, f=f)

    @Property
    def height(self):
        '''Get the height (C{meter}).
        '''
        return self._height

    @height.setter  # PYCHOK setter!
    def height(self, height):
        '''Set the height (C{meter}).

           @raise TypeError: Invalid B{C{height}} C{type}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        h = Height(height)
        if self._height != h:
            _update_all(self)
            self._height = h

    def _heigHt(self, height):
        '''(INTERNAL) Overriding this C{height}.
        '''
        return self.height if height is None else Height(height)

    def height4(self, earth=None, normal=True, LatLon=None, **LatLon_kwds):
        '''Compute the height above or below and the projection of this point
           on this datum's or on an other earth's ellipsoid surface.

           @kwarg earth: A datum, ellipsoid, triaxial ellipsoid or earth radius
                         I{overriding} this datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2}, L{a_f2Tuple}, L{Triaxial}, L{Triaxial_},
                         L{JacobiConformal} or C{meter}, conventionally).
           @kwarg normal: If C{True} the projection is the nearest point on the
                          ellipsoid's surface, otherwise the intersection of the
                          radial line to the center and the ellipsoid's surface.
           @kwarg LatLon: Optional class to return the  height and projection
                          (C{LatLon}) or C{None}.
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword arguments,
                               ignored if C{B{LatLon} is None}.

           @note: Use keyword argument C{height=0} to override C{B{LatLon}.height}
                  to {0} or any other C{scalar}, conventionally in C{meter}.

           @return: An instance of B{C{LatLon}} or if C{B{LatLon} is None}, a
                    L{Vector4Tuple}C{(x, y, z, h)} with the I{projection} C{x}, C{y}
                    and C{z} coordinates and height C{h} in C{meter}, conventionally.

           @raise TriaxialError: No convergence in triaxial root finding.

           @raise TypeError: Invalid B{C{earth}}.

           @see: L{Ellipsoid.height4} and L{Triaxial_.height4} for more information.
        '''
        c = self.toCartesian()
        if LatLon is None:
            r = c.height4(earth=earth, normal=normal)
        else:
            r = c.height4(earth=earth, normal=normal, Cartesian=c.classof, height=0)
            r = r.toLatLon(LatLon=LatLon, **_xkwds(LatLon_kwds, height=r.height))
        return r

    def heightStr(self, prec=-2, m=_m_):
        '''Return this point's B{C{height}} as C{str}ing.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).
           @kwarg m: Optional unit of the height (C{str}).

           @see: Function L{pygeodesy.hstr}.
        '''
        return _MODS.streprs.hstr(self.height, prec=prec, m=m)

    @deprecated_method
    def isantipode(self, other, eps=EPS):  # PYCHOK no cover
        '''DEPRECATED, use method L{isantipodeTo}.'''
        return self.isantipodeTo(other, eps=eps)

    def isantipodeTo(self, other, eps=EPS):
        '''Check whether this and an other point are antipodal,
           on diametrically opposite sides of the earth.

           @arg other: The other point (C{LatLon}).
           @kwarg eps: Tolerance for near-equality (C{degrees}).

           @return: C{True} if points are antipodal within the given
                    tolerance, C{False} otherwise.
        '''
        p = self.others(other)
        return isantipode(*(self.latlon + p.latlon), eps=eps)

    @Property_RO
    def isEllipsoidal(self):
        '''Check whether this point is ellipsoidal (C{bool} or C{None} if unknown).
        '''
        return self.datum.isEllipsoidal if self._datum else None

    @Property_RO
    def isEllipsoidalLatLon(self):
        '''Get C{LatLon} base.
        '''
        return False

    def isequalTo(self, other, eps=None):
        '''Compare this point with an other point, I{ignoring} height.

           @arg other: The other point (C{LatLon}).
           @kwarg eps: Tolerance for equality (C{degrees}).

           @return: C{True} if both points are identical,
                    I{ignoring} height, C{False} otherwise.

           @raise TypeError: The B{C{other}} point is not C{LatLon}
                             or mismatch of the B{C{other}} and
                             this C{class} or C{type}.

           @raise UnitError: Invalid B{C{eps}}.

           @see: Method L{isequalTo3}.
        '''
        return _isequalTo(self, self.others(other), eps=eps)

    def isequalTo3(self, other, eps=None):
        '''Compare this point with an other point, I{including} height.

           @arg other: The other point (C{LatLon}).
           @kwarg eps: Tolerance for equality (C{degrees}).

           @return: C{True} if both points are identical
                    I{including} height, C{False} otherwise.

           @raise TypeError: The B{C{other}} point is not C{LatLon}
                             or mismatch of the B{C{other}} and
                             this C{class} or C{type}.

           @see: Method L{isequalTo}.
        '''
        return self.height == self.others(other).height and \
              _isequalTo(self, other, eps=eps)

    @Property_RO
    def isnormal(self):
        '''Return C{True} if this point is normal (C{bool}),
           meaning C{abs(lat) <= 90} and C{abs(lon) <= 180}.

           @see: Methods L{normal}, L{toNormal} and functions
                 L{pygeodesy.isnormal} and L{pygeodesy.normal}.
        '''
        return isnormal(self.lat, self.lon, eps=0)

    @Property_RO
    def isSpherical(self):
        '''Check whether this point is spherical (C{bool} or C{None} if unknown).
        '''
        return self.datum.isSpherical if self._datum else None

    @Property_RO
    def lam(self):
        '''Get the longitude (B{C{radians}}).
        '''
        return radians(self.lon)

    @Property
    def lat(self):
        '''Get the latitude (C{degrees90}).
        '''
        return self._lat

    @lat.setter  # PYCHOK setter!
    def lat(self, lat):
        '''Set the latitude (C{str[N|S]} or C{degrees}).

           @raise ValueError: Invalid B{C{lat}}.
        '''
        lat = Lat(lat)  # parseDMS(lat, suffix=_NS_, clip=90)
        if self._lat != lat:
            _update_all(self)
            self._lat = lat

    @Property
    def latlon(self):
        '''Get the lat- and longitude (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return LatLon2Tuple(self._lat, self._lon, name=self.name)

    @latlon.setter  # PYCHOK setter!
    def latlon(self, latlonh):
        '''Set the lat- and longitude and optionally the height
           (2- or 3-tuple or comma- or space-separated C{str}
           of C{degrees90}, C{degrees180} and C{meter}).

           @raise TypeError: Height of B{C{latlonh}} not C{scalar} or
                             B{C{latlonh}} not C{list} or C{tuple}.

           @raise ValueError: Invalid B{C{latlonh}} or M{len(latlonh)}.

           @see: Function L{pygeodesy.parse3llh} to parse a B{C{latlonh}}
                 string into a 3-tuple C{(lat, lon, h)}.
        '''
        if isstr(latlonh):
            latlonh = parse3llh(latlonh, height=self.height)
        else:
            _xinstanceof(list, tuple, latlonh=latlonh)
        if len(latlonh) == 3:
            h = Height(latlonh[2], name=Fmt.SQUARE(latlonh=2))
        elif len(latlonh) != 2:
            raise _ValueError(latlonh=latlonh)
        else:
            h = self.height

        llh = Lat(latlonh[0]), Lon(latlonh[1]), h  # parseDMS2(latlonh[0], latlonh[1])
        if (self._lat, self._lon, self._height) != llh:
            _update_all(self)
            self._lat, self._lon, self._height   = llh

    def latlon2(self, ndigits=0):
        '''Return this point's lat- and longitude in C{degrees}, rounded.

           @kwarg ndigits: Number of (decimal) digits (C{int}).

           @return: A L{LatLon2Tuple}C{(lat, lon)}, both C{float}
                    and rounded away from zero.

           @note: The C{round}ed values are always C{float}, also
                  if B{C{ndigits}} is omitted.
        '''
        return LatLon2Tuple(round(self.lat, ndigits),
                            round(self.lon, ndigits), name=self.name)

    @deprecated_method
    def latlon_(self, ndigits=0):  # PYCHOK no cover
        '''DEPRECATED, use method L{latlon2}.'''
        return self.latlon2(ndigits=ndigits)

    latlon2round = latlon_  # PYCHOK no cover

    @Property
    def latlonheight(self):
        '''Get the lat-, longitude and height (L{LatLon3Tuple}C{(lat, lon, height)}).
        '''
        return self.latlon.to3Tuple(self.height)

    @latlonheight.setter  # PYCHOK setter!
    def latlonheight(self, latlonh):
        '''Set the lat- and longitude and optionally the height
           (2- or 3-tuple or comma- or space-separated C{str}
           of C{degrees90}, C{degrees180} and C{meter}).

           @see: Property L{latlon} for more details.
        '''
        self.latlon = latlonh

    @Property
    def lon(self):
        '''Get the longitude (C{degrees180}).
        '''
        return self._lon

    @lon.setter  # PYCHOK setter!
    def lon(self, lon):
        '''Set the longitude (C{str[E|W]} or C{degrees}).

           @raise ValueError: Invalid B{C{lon}}.
        '''
        lon = Lon(lon)  # parseDMS(lon, suffix=_EW_, clip=180)
        if self._lon != lon:
            _update_all(self)
            self._lon = lon

    @Property_RO
    def _ltp(self):
        '''(INTERNAL) Cache for L{toLtp}.
        '''
        return _MODS.ltp.Ltp(self, ecef=self.Ecef(self.datum), name=self.name)

    def nearestOn6(self, points, closed=False, height=None, wrap=False):
        '''Locate the point on a path or polygon closest to this point.

           Points are converted to and distances are computed in
           I{geocentric}, cartesian space.

           @arg points: The path or polygon points (C{LatLon}[]).
           @kwarg closed: Optionally, close the polygon (C{bool}).
           @kwarg height: Optional height, overriding the height of
                          this and all other points (C{meter}).  If
                          C{None}, take the height of points into
                          account for distances.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{points}} (C{bool}).

           @return: A L{NearestOn6Tuple}C{(closest, distance, fi, j,
                    start, end)} with the C{closest}, the C{start}
                    and the C{end} point each an instance of this
                    C{LatLon} and C{distance} in C{meter}, same
                    units as the cartesian axes.

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} or some B{C{points}}'
                             C{Ecef} invalid.

           @raise ValueError: Some B{C{points}}' C{Ecef} is incompatible.

           @see: Function L{pygeodesy.nearestOn6}.
        '''
        def _cs(Ps, h, w, C):
            p = None  # not used
            for i, q in Ps.enumerate():
                if w and i:
                    q = _unrollon(p, q)
                yield C(height=h, i=i, up=3, points=q)
                p = q

        C  = self._toCartesianEcef  # to verify datum and Ecef
        Ps = self.PointsIter(points, wrap=wrap)

        c = C(height=height, this=self)  # this Cartesian
        t = nearestOn6(c, _cs(Ps, height, wrap, C), closed=closed)
        c, s, e = t.closest, t.start, t.end

        kwds = _xkwds_not(None, LatLon=self.classof,  # this LatLon
                                height=height)
        _r =  self.Ecef(self.datum).reverse
        p  = _r(c).toLatLon(**kwds)
        s  = _r(s).toLatLon(**kwds) if s is not c else p
        e  = _r(e).toLatLon(**kwds) if e is not c else p
        return t.dup(closest=p, start=s, end=e)

    def normal(self):
        '''Normalize this point I{in-place} to C{abs(lat) <= 90} and
           C{abs(lon) <= 180}.

           @return: C{True} if this point was I{normal}, C{False} if it
                    wasn't (but is now).

           @see: Property L{isnormal} and method L{toNormal}.
        '''
        n = self.isnormal
        if not n:
            self.latlon = normal(*self.latlon)
        return n

    @Property_RO
    def _N_vector(self):
        '''(INTERNAL) Get the (C{nvectorBase._N_vector_})
        '''
        return _MODS.nvectorBase._N_vector_(*self.xyzh)

    @Property_RO
    def phi(self):
        '''Get the latitude (B{C{radians}}).
        '''
        return radians(self.lat)

    @Property_RO
    def philam(self):
        '''Get the lat- and longitude (L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        return PhiLam2Tuple(self.phi, self.lam, name=self.name)

    def philam2(self, ndigits=0):
        '''Return this point's lat- and longitude in C{radians}, rounded.

           @kwarg ndigits: Number of (decimal) digits (C{int}).

           @return: A L{PhiLam2Tuple}C{(phi, lam)}, both C{float}
                    and rounded away from zero.

           @note: The C{round}ed values are always C{float}, also
                  if B{C{ndigits}} is omitted.
        '''
        return PhiLam2Tuple(round(self.phi, ndigits),
                            round(self.lam, ndigits), name=self.name)

    @Property_RO
    def philamheight(self):
        '''Get the lat-, longitude in C{radians} and height (L{PhiLam3Tuple}C{(phi, lam, height)}).
        '''
        return self.philam.to3Tuple(self.height)

    @deprecated_method
    def points(self, points, closed=True):  # PYCHOK no cover
        '''DEPRECATED, use method L{points2}.'''
        return self.points2(points, closed=closed)

    def points2(self, points, closed=True):
        '''Check a path or polygon represented by points.

           @arg points: The path or polygon points (C{LatLon}[])
           @kwarg closed: Optionally, consider the polygon closed,
                          ignoring any duplicate or closing final
                          B{C{points}} (C{bool}).

           @return: A L{Points2Tuple}C{(number, points)}, an C{int}
                    and C{list} or C{tuple}.

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} are not C{LatLon}.
        '''
        return _MODS.iters.points2(points, closed=closed, base=self)

    def PointsIter(self, points, loop=0, dedup=False, wrap=False):
        '''Return a C{PointsIter} iterator.

           @arg points: The path or polygon points (C{LatLon}[])
           @kwarg loop: Number of loop-back points (non-negative C{int}).
           @kwarg dedup: Skip duplicate points (C{bool}).
           @kwarg wrap: If C{True}, wrap or I{normalize} the
                        enum-/iterated B{C{points}} (C{bool}).

           @return: A new C{PointsIter} iterator.

           @raise PointsError: Insufficient number of B{C{points}}.
        '''
        return PointsIter(points, base=self, loop=loop, dedup=dedup, wrap=wrap)

    def radii11(self, point2, point3, wrap=False):
        '''Return the radii of the C{Circum-}, C{In-}, I{Soddy} and C{Tangent}
           circles of a (planar) triangle formed by this and two other points.

           @arg point2: Second point (C{LatLon}).
           @arg point3: Third point (C{LatLon}).
           @kwarg wrap: If C{True}, wrap or I{normalize} B{C{point2}} and
                        B{C{point3}} (C{bool}).

           @return: L{Radii11Tuple}C{(rA, rB, rC, cR, rIn, riS, roS, a, b, c, s)}.

           @raise IntersectionError: Near-coincident or -colinear points.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @see: Function L{pygeodesy.radii11}, U{Incircle
                 <https://MathWorld.Wolfram.com/Incircle.html>}, U{Soddy Circles
                 <https://MathWorld.Wolfram.com/SoddyCircles.html>} and U{Tangent
                 Circles<https://MathWorld.Wolfram.com/TangentCircles.html>}.
        '''
        with _toCartesian3(self, point2, point3, wrap) as cs:
            return _radii11ABC(*cs, useZ=True)[0]

    def _rhumb3(self, exact, radius):  # != .sphericalBase._rhumbs3
        '''(INTERNAL) Get the C{rhumb} for this point's datum or for
           the B{C{radius}}' earth model iff non-C{None}.
        '''
        try:
            d = self._rhumb3dict
            t = d[(exact, radius)]
        except KeyError:
            D = self.datum if radius is None else _spherical_datum(radius)  # ellipsoidal OK
            r = D.ellipsoid.rhumb_(exact=exact)  # or D.isSpherical)
            t = r, D, _MODS.karney.Caps
            while d:
                d.popitem()
            d[(exact, radius)] = t  # cache 3-tuple
        return t

    @Property_RO
    def _rhumb3dict(self):
        return {}  # single-item cache

    def rhumbAzimuthTo(self, other, exact=False, radius=None, wrap=False):
        '''Return the azimuth (bearing) of a rhumb line (loxodrome)
           between this and an other (ellipsoidal) point.

           @arg other: The other point (C{LatLon}).
           @kwarg exact: Exact C{Rhumb...} to use (C{bool} or C{Rhumb...}),
                         see method L{Ellipsoid.rhumb_}.
           @kwarg radius: Optional earth radius (C{meter}) or earth model
                          (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                          L{a_f2Tuple}), overriding this point's datum.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: Rhumb azimuth (compass C{degrees360}).

           @raise TypeError: The B{C{other}} point is incompatible or
                             B{C{radius}} is invalid.
        '''
        r, _, Cs = self._rhumb3(exact, radius)
        return r._Inverse(self, other, wrap, outmask=Cs.AZIMUTH).azi12

    def rhumbDestination(self, distance, azimuth, exact=False, radius=None, height=None):
        '''Return the destination point having travelled the given distance
           from this point along a rhumb line (loxodrome) at the given azimuth.

           @arg distance: Distance travelled (C{meter}, same units as this
                          point's datum (ellipsoid) axes or B{C{radius}},
                          may be negative.
           @arg azimuth: Azimuth (bearing) at this point (compass C{degrees}).
           @kwarg exact: Exact C{Rhumb...} to use (C{bool} or C{Rhumb...}),
                         see method L{Ellipsoid.rhumb_}.
           @kwarg radius: Optional earth radius (C{meter}) or earth model
                          (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                          L{a_f2Tuple}), overriding this point's datum.
           @kwarg height: Optional height, overriding the default height
                          (C{meter}).

           @return: The destination point (ellipsoidal C{LatLon}).

           @raise TypeError: Invalid B{C{radius}}.

           @raise ValueError: Invalid B{C{distance}}, B{C{azimuth}},
                              B{C{radius}} or B{C{height}}.
        '''
        r, D, _ = self._rhumb3(exact, radius)
        d = r._Direct(self, azimuth, distance)
        h = self._heigHt(height)
        return self.classof(d.lat2, d.lon2, datum=D, height=h)

    def rhumbDistanceTo(self, other, exact=False, radius=None, wrap=False):
        '''Return the distance from this to an other point along
           a rhumb line (loxodrome).

           @arg other: The other point (C{LatLon}).
           @kwarg exact: Exact C{Rhumb...} to use (C{bool} or C{Rhumb...}),
                         see method L{Ellipsoid.rhumb_}.
           @kwarg radius: Optional earth radius (C{meter}) or earth model
                          (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                          L{a_f2Tuple}), overriding this point's datum.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: Distance (C{meter}, the same units as this point's
                    datum (ellipsoid) axes or B{C{radius}}.

           @raise TypeError: The B{C{other}} point is incompatible or
                             B{C{radius}} is invalid.

           @raise ValueError: Invalid B{C{radius}}.
        '''
        r, _, Cs = self._rhumb3(exact, radius)
        return r._Inverse(self, other, wrap, outmask=Cs.DISTANCE).s12

    def rhumbLine(self, azimuth_other, exact=False, radius=None, wrap=False,
                                                               **name_caps):
        '''Get a rhumb line through this point at a given azimuth or
           through this and an other point.

           @arg azimuth_other: The azimuth of the rhumb line (compass
                               C{degrees}) or the other point (C{LatLon}).
           @kwarg exact: Exact C{Rhumb...} to use (C{bool} or C{Rhumb...}),
                         see method L{Ellipsoid.rhumb_}.
           @kwarg radius: Optional earth radius (C{meter}) or earth model
                          (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                          L{a_f2Tuple}), overriding this point's datum.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        C{azimuth_B{other}} point (C{bool}).
           @kwarg name_caps: Optional C{B{name}=str} and C{caps}, see
                             L{RhumbLine} C{B{caps}}.

           @return: A C{RhumbLine} instance.

           @raise TypeError: Invalid B{C{radius}} or BC{C{azimuth_other}}
                             not a C{scalar} nor a C{LatLon}.

           @see: Modules L{rhumbaux} and L{rhumbx}.
        '''
        r, _, _ = self._rhumb3(exact, radius)
        a, kwds = azimuth_other, _xkwds(name_caps, name=self.name)
        if isscalar(a):
            r = r._DirectLine(self, a, **kwds)
        elif isinstance(a, LatLonBase):
            r = r._InverseLine(self, a, wrap, **kwds)
        else:
            raise _TypeError(azimuth_other=a)
        return r

    def rhumbMidpointTo(self, other, exact=False, radius=None,
                                     height=None, fraction=_0_5, wrap=False):
        '''Return the (loxodromic) midpoint on the rhumb line between
           this and an other point.

           @arg other: The other point (C{LatLon}).
           @kwarg exact: Exact C{Rhumb...} to use (C{bool} or C{Rhumb...}),
                         see method L{Ellipsoid.rhumb_}.
           @kwarg radius: Optional earth radius (C{meter}) or earth model
                          (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                          L{a_f2Tuple}), overriding this point's datum.
           @kwarg height: Optional height, overriding the mean height
                          (C{meter}).
           @kwarg fraction: Midpoint location from this point (C{scalar}), 0
                            for this, 1 for the B{C{other}}, 0.5 for halfway
                            between this and the B{C{other}} point, may be
                            negative or greater than 1.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: The midpoint at the given B{C{fraction}} along the
                    rhumb line (C{LatLon}).

           @raise TypeError: The B{C{other}} point is incompatible or
                             B{C{radius}} is invalid.

           @raise ValueError: Invalid B{C{height}} or B{C{fraction}}.
        '''
        r, D, _ = self._rhumb3(exact, radius)
        f = Scalar(fraction=fraction)
        d = r._Inverse(self, other, wrap)  # C.AZIMUTH_DISTANCE
        d = r._Direct( self, d.azi12, d.s12 * f)
        h = self._havg(other, f=f, h=height)
        return self.classof(d.lat2, d.lon2, datum=D, height=h)

    def thomasTo(self, other, wrap=False):
        '''Compute the distance between this and an other point using
           U{Thomas'<https://apps.DTIC.mil/dtic/tr/fulltext/u2/703541.pdf>}
           formula.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{other}} point (C{bool}).

           @return: Distance (C{meter}, same units as the axes of
                    this point's datum ellipsoid).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.thomas} and methods L{cosineAndoyerLambertTo},
                 L{cosineForsytheAndoyerLambertTo}, L{cosineLawTo}, C{distanceTo*},
                 L{equirectangularTo}, L{euclideanTo}, L{flatLocalTo}/L{hubenyTo},
                 L{flatPolarTo}, L{haversineTo} and L{vincentysTo}.
        '''
        return self._distanceTo_(thomas_, other, wrap=wrap)

    @deprecated_method
    def to2ab(self):  # PYCHOK no cover
        '''DEPRECATED, use property L{philam}.'''
        return self.philam

    def toCartesian(self, height=None, Cartesian=None, **Cartesian_kwds):
        '''Convert this point to cartesian, I{geocentric} coordinates,
           also known as I{Earth-Centered, Earth-Fixed} (ECEF).

           @kwarg height: Optional height, overriding this point's height
                          (C{meter}, conventionally).
           @kwarg Cartesian: Optional class to return the geocentric
                             coordinates (C{Cartesian}) or C{None}.
           @kwarg Cartesian_kwds: Optional, additional B{C{Cartesian}}
                                  keyword arguments, ignored if
                                  C{B{Cartesian} is None}.

           @return: A B{C{Cartesian}} or if B{C{Cartesian}} is C{None},
                    an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M,
                    datum)} with C{C=0} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian}} or B{C{Cartesian_kwds}}.
        '''
        r = self._ecef9 if height is None else self.toEcef(height=height)
        if Cartesian is not None:  # class or .classof
            r = self._xnamed(Cartesian(r, **Cartesian_kwds))
        _xdatum(r.datum, self.datum)
        return r

    def _toCartesianEcef(self, height=None, i=None, up=2, **name_point):
        '''(INTERNAL) Convert to cartesian and check Ecef's before and after.
        '''
        p = self.others(up=up, **name_point)
        c = p.toCartesian(height=height)
        E = self.Ecef
        if E:
            for p in (p, c):
                e = getattr(p, LatLonBase.Ecef.name, None)
                if e not in (None, E):  # PYCHOK no cover
                    n, _ = name_point.popitem()
                    if i is not None:
                        Fmt.SQUARE(n, i)
                    raise _ValueError(n, e, txt=_incompatible(E.__name__))
        return c

    def toEcef(self, height=None, M=False):
        '''Convert this point to I{geocentric} coordinates, also known as
           I{Earth-Centered, Earth-Fixed} (U{ECEF<https://WikiPedia.org/wiki/ECEF>}).

           @kwarg height: Optional height, overriding this point's height
                          (C{meter}, conventionally).
           @kwarg M: Optionally, include the rotation L{EcefMatrix} (C{bool}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M, datum)}
                    with C{C=0} and C{M} if available.

           @raise EcefError: A C{.datum} or an ECEF issue.
        '''
        return self._ecef9 if height in (None, self.height) else \
               self._Ecef_forward(self.lat, self.lon, height=height, M=M)

    @deprecated_method
    def to3llh(self, height=None):  # PYCHOK no cover
        '''DEPRECATED, use property L{latlonheight} or C{latlon.to3Tuple(B{height})}.'''
        return self.latlonheight if height in (None, self.height) else \
               self.latlon.to3Tuple(height)

    def toLocal(self, Xyz=None, ltp=None, **Xyz_kwds):
        '''Convert this I{geodetic} point to I{local} C{X}, C{Y} and C{Z}.

           @kwarg Xyz: Optional class to return C{X}, C{Y} and C{Z}
                       (L{XyzLocal}, L{Enu}, L{Ned}) or C{None}.
           @kwarg ltp: The I{local tangent plane} (LTP) to use,
                       overriding this point's LTP (L{Ltp}).
           @kwarg Xyz_kwds: Optional, additional B{C{Xyz}} keyword
                            arguments, ignored if C{B{Xyz} is None}.

           @return: An B{C{Xyz}} instance or if C{B{Xyz} is None},
                    a L{Local9Tuple}C{(x, y, z, lat, lon, height,
                    ltp, ecef, M)} with C{M=None}, always.

           @raise TypeError: Invalid B{C{ltp}}.
        '''
        p = _MODS.ltp._xLtp(ltp, self._ltp)
        return p._ecef2local(self._ecef9, Xyz, Xyz_kwds)

    def toLtp(self, Ecef=None):
        '''Return the I{local tangent plane} (LTP) for this point.

           @kwarg Ecef: Optional ECEF I{class} (L{EcefKarney}, ...
                        L{EcefYou}), overriding this point's C{Ecef}.
        '''
        return self._ltp if Ecef in (None, self.Ecef) else _MODS.ltp.Ltp(
               self, ecef=Ecef(self.datum), name=self.name)

    def toNormal(self, deep=False, name=NN):
        '''Get this point I{normalized} to C{abs(lat) <= 90}
           and C{abs(lon) <= 180}.

           @kwarg deep: If C{True} make a deep, otherwise a
                        shallow copy (C{bool}).
           @kwarg name: Optional name of the copy (C{str}).

           @return: A copy of this point, I{normalized} and
                    optionally renamed (C{LatLon}).

           @see: Property L{isnormal}, method L{normal} and function
                 L{pygeodesy.normal}.
        '''
        ll = self.copy(deep=deep)
        _  = ll.normal()
        if name:
            ll.rename(name)
        return ll

    def toNvector(self, h=None, Nvector=None, **Nvector_kwds):
        '''Convert this point to C{n-vector} (normal to the earth's surface)
           components, I{including height}.

           @kwarg h: Optional height, overriding this point's
                     height (C{meter}).
           @kwarg Nvector: Optional class to return the C{n-vector}
                           components (C{Nvector}) or C{None}.
           @kwarg Nvector_kwds_wrap: Optional, additional B{C{Nvector}}
                               keyword arguments, ignored if C{B{Nvector}
                               is None}.

           @return: A B{C{Nvector}} or a L{Vector4Tuple}C{(x, y, z, h)}
                    if B{C{Nvector}} is C{None}.

           @raise TypeError: Invalid B{C{Nvector}} or B{C{Nvector_kwds}}.
        '''
        return self.toVector(Vector=Nvector, h=self.height if h is None else h,
                                            ll=self, **Nvector_kwds)

    def toStr(self, form=F_DMS, joined=_COMMASPACE_, m=_m_, **prec_sep_s_D_M_S):  # PYCHOK expected
        '''Convert this point to a "lat, lon[, +/-height]" string, formatted
           in the given C{B{form}at}.

           @kwarg form: The lat-/longitude C{B{form}at} to use (C{str}), see
                        functions L{pygeodesy.latDMS} or L{pygeodesy.lonDMS}.
           @kwarg joined: Separator to join the lat-, longitude and heigth
                          strings (C{str} or C{None} or C{NN} for non-joined).
           @kwarg m: Optional unit of the height (C{str}), use C{None} to
                     exclude height from the returned string.
           @kwarg prec_sep_s_D_M_S: Optional C{B{prec}ision}, C{B{sep}arator},
                      B{C{s_D}}, B{C{s_M}}, B{C{s_S}} and B{C{s_DMS}} keyword
                      arguments, see function L{pygeodesy.latDMS} or
                      L{pygeodesy.lonDMS}.

           @return: This point in the specified C{B{form}at}, etc. (C{str} or
                    a 2- or 3-tuple C{(lat_str, lon_str[, height_str])} if
                    C{B{joined}=NN} or C{B{joined}=None}).

           @see: Function L{pygeodesy.latDMS} or L{pygeodesy.lonDMS} for more
                 details about keyword arguments C{B{form}at}, C{B{prec}ision},
                 C{B{sep}arator}, B{C{s_D}}, B{C{s_M}}, B{C{s_S}} and B{C{s_DMS}}.

           @example:

            >>> LatLon(51.4778, -0.0016).toStr()  # 51°28′40″N, 000°00′06″W
            >>> LatLon(51.4778, -0.0016).toStr(F_D)  # 51.4778°N, 000.0016°W
            >>> LatLon(51.4778, -0.0016, 42).toStr()  # 51°28′40″N, 000°00′06″W, +42.00m
        '''
        t = (latDMS(self.lat, form=form, **prec_sep_s_D_M_S),
             lonDMS(self.lon, form=form, **prec_sep_s_D_M_S))
        if self.height and m is not None:
            t += (self.heightStr(m=m),)
        return joined.join(t) if joined else t

    def toVector(self, Vector=None, **Vector_kwds):
        '''Convert this point to C{n-vector} (normal to the earth's
           surface) components, I{ignoring height}.

           @kwarg Vector: Optional class to return the C{n-vector}
                          components (L{Vector3d}) or C{None}.
           @kwarg Vector_kwds: Optional, additional B{C{Vector}}
                               keyword arguments, ignored if
                               C{B{Vector} is None}.

           @return: A B{C{Vector}} or a L{Vector3Tuple}C{(x, y, z)}
                    if B{C{Vector}} is C{None}.

           @raise TypeError: Invalid B{C{Vector}} or B{C{kwds}}.

           @note: These are C{n-vector} x, y and z components,
                  I{NOT} geocentric (ECEF) x, y and z coordinates!
        '''
        r = self._vector3tuple
        if Vector is not None:
            r = Vector(*r, **_xkwds(Vector_kwds, name=self.name))
        return r

    def toVector3d(self):
        '''Convert this point to C{n-vector} (normal to the earth's
           surface) components, I{ignoring height}.

           @return: Unit vector (L{Vector3d}).

           @note: These are C{n-vector} x, y and z components,
                  I{NOT} geocentric (ECEF) x, y and z coordinates!
        '''
        return self._vector3d  # XXX .unit()

    def toWm(self, **toWm_kwds):
        '''Convert this point to a WM coordinate.

           @kwarg toWm_kwds: Optional L{pygeodesy.toWm} keyword arguments.

           @return: The WM coordinate (L{Wm}).

           @see: Function L{pygeodesy.toWm}.
        '''
        return self._wm if not toWm_kwds else _MODS.webmercator.toWm(
               self, **_xkwds(toWm_kwds, name=self.name))

    @deprecated_method
    def to3xyz(self):  # PYCHOK no cover
        '''DEPRECATED, use property L{xyz} or method L{toNvector}, L{toVector},
           L{toVector3d} or perhaps (geocentric) L{toEcef}.'''
        return self.xyz  # self.toVector()

    @Property_RO
    def _vector3d(self):
        '''(INTERNAL) Cache for L{toVector3d}.
        '''
        return self.toVector(Vector=Vector3d)  # XXX .unit()

    @Property_RO
    def _vector3tuple(self):
        '''(INTERNAL) Cache for L{toVector}.
        '''
        return philam2n_xyz(self.phi, self.lam, name=self.name)

    def vincentysTo(self, other, **radius_wrap):
        '''Compute the distance between this and an other point using
           U{Vincenty's<https://WikiPedia.org/wiki/Great-circle_distance>}
           spherical formula.

           @arg other: The other point (C{LatLon}).
           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.vincentys}, overriding the
                               default mean C{radius} of this point's
                               datum ellipsoid.

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.vincentys} and methods L{cosineAndoyerLambertTo},
                 L{cosineForsytheAndoyerLambertTo}, L{cosineLawTo}, C{distanceTo*},
                 L{equirectangularTo}, L{euclideanTo}, L{flatLocalTo}/L{hubenyTo},
                 L{flatPolarTo}, L{haversineTo} and L{thomasTo}.
        '''
        return self._distanceTo(vincentys, other, **_xkwds(radius_wrap, radius=None))

    @Property_RO
    def _wm(self):
        '''(INTERNAL) Get this point as webmercator (L{Wm}).
        '''
        return _MODS.webmercator.toWm(self)

    @Property_RO
    def xyz(self):
        '''Get the C{n-vector} X, Y and Z components (L{Vector3Tuple}C{(x, y, z)})

           @note: These are C{n-vector} x, y and z components, I{NOT}
                  geocentric (ECEF) x, y and z coordinates!
        '''
        return self.toVector(Vector=Vector3Tuple)

    @Property_RO
    def xyzh(self):
        '''Get the C{n-vector} X, Y, Z and H components (L{Vector4Tuple}C{(x, y, z, h)})

           @note: These are C{n-vector} x, y and z components, I{NOT}
                  geocentric (ECEF) x, y and z coordinates!
        '''
        return self.xyz.to4Tuple(self.height)


class _toCartesian3(object):  # see also .geodesicw._wargs, .vector2d._numpy
    '''(INTERNAL) Wrapper to convert 2 other points.
    '''
    @contextmanager  # <https://www.python.org/dev/peps/pep-0343/> Examples
    def __call__(self, p, p2, p3, wrap, **kwds):
        try:
            if wrap:
                p2, p3 =  map1(_Wrap.point, p2, p3)
                kwds   = _xkwds(kwds, wrap=wrap)
            yield (p. toCartesian().copy(name=_point_),  # copy to rename
                   p._toCartesianEcef(up=4, point2=p2),
                   p._toCartesianEcef(up=4, point3=p3))
        except (AssertionError, TypeError, ValueError) as x:
            raise _xError(x, point=p, point2=p2, point3=p3, **kwds)

_toCartesian3 = _toCartesian3()  # PYCHOK singleton


def _trilaterate5(p1, d1, p2, d2, p3, d3, area=True, eps=EPS1,  # MCCABE 13
                                          radius=R_M, wrap=False):
    '''(INTERNAL) Trilaterate three points by I{area overlap} or by
       I{perimeter intersection} of three circles.

       @note: The B{C{radius}} is only needed for both the n-vectorial
              and C{sphericalTrigonometry.LatLon.distanceTo} methods and
              silently ignored by the C{ellipsoidalExact}, C{-GeodSolve},
              C{-Karney} and C{-Vincenty.LatLon.distanceTo} methods.
    '''
    p2, p3, w = _unrollon3(p1, p2, p3, wrap)

    r1 = Distance_(distance1=d1)
    r2 = Distance_(distance2=d2)
    r3 = Distance_(distance3=d3)
    m  = 0 if area else (r1 + r2 + r3)
    pc = 0
    t  = []
    for _ in range(3):
        try:  # intersection of circle (p1, r1) and (p2, r2)
            c1, c2 = p1.intersections2(r1, p2, r2, wrap=w)

            if area:  # check overlap
                if c1 is c2:  # abutting
                    c = c1
                else:  # nearest point on radical
                    c = p3.nearestOn(c1, c2, within=True, wrap=w)
                d = r3 - p3.distanceTo(c, radius=radius, wrap=w)
                if d > eps:  # sufficient overlap
                    t.append((d, c))
                m = max(m, d)

            else:  # check intersection
                for c in ((c1,) if c1 is c2 else (c1, c2)):
                    d = fabs(r3 - p3.distanceTo(c, radius=radius, wrap=w))
                    if d < eps:  # below margin
                        t.append((d, c))
                    m = min(m, d)

        except IntersectionError as x:
            if _concentric_ in str(x):  # XXX ConcentricError?
                pc += 1

        p1, r1, p2, r2, p3, r3 = p2, r2, p3, r3, p1, r1  # rotate

    if t:  # get min, max, points and count ...
        t = tuple(sorted(t))
        n = len(t),  # as 1-tuple
        # ... or for a single trilaterated result,
        # min *is* max, min- *is* maxPoint and n=1, 2 or 3
        return Trilaterate5Tuple(t[0] + t[-1] + n)  # *(t[0] + ...)

    elif area and pc == 3:  # all pairwise concentric ...
        r, p = min((r1, p1), (r2, p2), (r3, p3))
        m = max(r1, r2, r3)
        # ... return "smallest" point twice, the smallest
        # and largest distance and n=0 for concentric
        return Trilaterate5Tuple(float(r), p, float(m), p, 0)

    n, f = (_overlap_, max) if area else (_intersection_, min)
    t = _COMMASPACE_(_no_(n), '%s %.3g' % (f.__name__, m))
    raise IntersectionError(area=area, eps=eps, wrap=wrap, txt=t)


__all__ += _ALL_DOCS(LatLonBase)

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
