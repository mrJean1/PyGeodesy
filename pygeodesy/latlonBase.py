
# -*- coding: utf-8 -*-

u'''(INTERNAL) Private base class L{LatLonBase} for elliposiodal, spherical
and N-vectorial C{LatLon}s.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**,
see U{https://www.Movable-Type.co.UK/scripts/latlong.html},
U{<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}
and U{https://www.Movable-Type.co.UK/scripts/latlong-vectors.html}.
'''

from pygeodesy.basics import isscalar, isstr, _xinstanceof
from pygeodesy.constants import EPS, EPS0, EPS1, EPS4, R_M, _0_0, _0_5, _1_0
# from pygeodesy.datums import _spherical_datum  # from .formy
from pygeodesy.dms import F_D, F_DMS, latDMS, lonDMS, parse3llh
from pygeodesy.errors import _incompatible, IntersectionError, _TypeError, \
                             _ValueError, _xdatum, _xError, _xkwds, _xkwds_not
from pygeodesy.formy import antipode, compassAngle, cosineAndoyerLambert_, \
                            cosineForsytheAndoyerLambert_, cosineLaw, \
                            equirectangular, euclidean, flatLocal_, \
                            flatPolar, hartzell, haversine, isantipode, \
                            isnormal, normal, philam2n_xyz, \
                            _spherical_datum, thomas_, vincentys
from pygeodesy.interns import NN, _COMMASPACE_, _concentric_, _height_, \
                             _intersection_, _m_, _no_, _overlap_, _point_
from pygeodesy.iters import PointsIter, points2
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _NamedBase, notOverloaded
from pygeodesy.namedTuples import Bounds2Tuple, LatLon2Tuple, PhiLam2Tuple, \
                                  Trilaterate5Tuple, Vector3Tuple
from pygeodesy.props import deprecated_method, Property, Property_RO, \
                            property_RO, _update_all
from pygeodesy.streprs import Fmt, hstr
from pygeodesy.units import Distance_, Lat, Lon, Height, Radius, Radius_, \
                            Scalar, Scalar_
from pygeodesy.utily import _unrollon, unrollPI
from pygeodesy.vector2d import _circin6,  Circin6Tuple, _circum3, Circum3Tuple, \
                                circum4_, Circum4Tuple, _radii11ABC
from pygeodesy.vector3d import nearestOn6, Vector3d

from math import asin, cos, degrees, radians

__all__ = _ALL_LAZY.latlonBase
__version__ = '22.10.12'


class LatLonBase(_NamedBase):
    '''(INTERNAL) Base class for C{LatLon} points on spherical or
       ellipsoidal earth models.
    '''
    _datum  = None  # L{Datum}, to be overriden
    _height = 0  # height (C{meter}), default
    _lat    = 0  # latitude (C{degrees})
    _lon    = 0  # longitude (C{degrees})

    def __init__(self, lat, lon, height=0, name=NN):
        '''New C{LatLon}.

           @arg lat: Latitude (C{degrees} or DMS C{str} with N or S suffix).
           @arg lon: Longitude (C{degrees} or DMS C{str} with E or W suffix).
           @kwarg height: Optional height (C{meter} above or below the earth surface).
           @kwarg name: Optional name (C{str}).

           @return: New instance (C{LatLon}).

           @raise RangeError: Value of B{C{lat}} or B{C{lon}} outside the valid
                              range and C{rangerrors} set to C{True}.

           @raise UnitError: Invalid B{C{lat}}, B{C{lon}} or B{C{height}}.

           @example:

            >>> p = LatLon(50.06632, -5.71475)
            >>> q = LatLon('50°03′59″N', """005°42'53"W""")
        '''
        if name:
            self.name = name

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
        a, b = antipode(self.lat, self.lon)  # PYCHOK LatLon2Tuple
        h = self.height if height is None else Height(height)
        return self.classof(a, b, height=h)

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
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg height: Height for C{latlonSW} and C{latlonNE} (C{meter}),
                          overriding the point's height.

           @return: A L{Bounds2Tuple}C{(latlonSW, latlonNE)}, the
                    lower-left and upper-right corner (C{LatLon}).

           @see: U{https://www.Movable-Type.co.UK/scripts/latlong-db.html}
        '''
        x = Scalar_(wide=wide) * _0_5
        y = Scalar_(tall=tall) * _0_5
        if radius is not None:
            r = Radius_(radius)
            c = cos(self.phi)
            x = degrees(asin(x / r) / c) if abs(c) > EPS0 else _0_0  # XXX
            y = degrees(y / r)
        x, y = abs(x), abs(y)

        h  = self.height if height is None else Height(height)
        sw = self.classof(self.lat - y, self.lon - x, height=h)
        ne = self.classof(self.lat + y, self.lon + x, height=h)
        return Bounds2Tuple(sw, ne, name=self.name)

    def chordTo(self, other, height=None):
        '''Compute the length of the chord through the earth between
           this and an other point.

           @arg other: The other point (C{LatLon}).
           @kwarg height: Overriding height for both points (C{meter})
                          or C{None} for each point's height.

           @return: The chord length (conventionally C{meter}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.
        '''
        def _v3d(ll):
            t = ll.toEcef(height=height)  # .toVector(Vector=Vector3d)
            return Vector3d(t.x, t.y, t.z)

        self.others(other)
        return _v3d(self).minus(_v3d(other)).length

    def circin6(self, point2, point3, eps=EPS4):
        '''Return the radius and center of the I{inscribed} aka I{In-}circle
           of the (planar) triangle formed by this and two other points.

           @arg point2: Second point (C{LatLon}).
           @arg point3: Third point (C{LatLon}).
           @kwarg eps: Tolerance for function L{pygeodesy.trilaterate3d2}.

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
        try:
            cs = self._toCartesian3(point2, point3)
            r, c, d, cA, cB, cC = _circin6(*cs, eps=eps, useZ=True, dLL3=True,
                                                datum=self.datum)  # PYCHOK unpack
            return Circin6Tuple(r, c.toLatLon(), d, cA.toLatLon(), cB.toLatLon(), cC.toLatLon())
        except (AssertionError, TypeError, ValueError) as x:
            raise _xError(x, point=self, point2=point2, point3=point3)

    def circum3(self, point2, point3, circum=True, eps=EPS4):
        '''Return the radius and center of the smallest circle I{through} or I{containing}
           this and two other points.

           @arg point2: Second point (C{LatLon}).
           @arg point3: Third point (C{LatLon}).
           @kwarg circum: If C{True} return the C{circumradius} and C{circumcenter},
                          always, ignoring the I{Meeus}' Type I case (C{bool}).
           @kwarg eps: Tolerance for function L{pygeodesy.trilaterate3d2}.

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
        try:
            cs = self._toCartesian3(point2, point3)
            r, c, d = _circum3(*cs, circum=circum, eps=eps, useZ=True, dLL3=True,  # XXX -3d2
                                    clas=cs[0].classof, datum=self.datum)  # PYCHOK unpack
            return Circum3Tuple(r, c.toLatLon(), d)
        except (AssertionError, TypeError, ValueError) as x:
            raise _xError(x, point=self, point2=point2, point3=point3, circum=circum)

    def circum4_(self, *points):
        '''Best-fit a sphere through this and two or more other points.

           @arg points: The other points (each a C{LatLon}).

           @return: L{Circum4Tuple}C{(radius, center, rank, residuals)} with C{center}
                    an instance of this (sub-)class.

           @raise ImportError: Package C{numpy} not found, not installed or older than
                               version 1.10.

           @raise NumPyError: Some C{numpy} issue.

           @raise TypeError: One of the B{C{points}} invalid.

           @raise ValueError: Too few B{C{points}}.

           @see: Function L{pygeodesy.circum4_} and L{circum3}.
        '''
        C = self._toCartesianEcef
        c = C(point=self)
        t = circum4_(c, Vector=c.classof, *(C(i=i, points=p) for i, p in enumerate(points)))
        c = t.center.toLatLon(LatLon=self.classof, name=t.name)
        return Circum4Tuple(t.radius, c, t.rank, t.residuals, name=c.name)

    @deprecated_method
    def compassAngle(self, other, adjust=True, wrap=False):  # PYCHOK no cover
        '''DEPRECATED, use method L{compassAngleTo}.'''
        return self.compassAngleTo(other, adjust=adjust, wrap=wrap)

    def compassAngleTo(self, other, adjust=True, wrap=False):
        '''Return the angle from North for the direction vector between
           this and an other point.

           Suitable only for short, non-near-polar vectors up to a few
           hundred Km or Miles.  Use method C{initialBearingTo} for
           larger distances.

           @arg other: The other point (C{LatLon}).
           @kwarg adjust: Adjust the longitudinal delta by the
                          cosine of the mean latitude (C{bool}).
           @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes and
                        longitudinal delta (C{bool}).

           @return: Compass angle from North (C{degrees360}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @note: Courtesy of Martin Schultz.

           @see: U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>}.
        '''
        self.others(other)
        return compassAngle(self.lat, self.lon, other.lat, other.lon,
                            adjust=adjust, wrap=wrap)

    def cosineAndoyerLambertTo(self, other, wrap=False):
        '''Compute the distance between this and an other point using
           the U{Andoyer-Lambert correction<https://navlib.net/wp-content/uploads/
           2013/10/admiralty-manual-of-navigation-vol-1-1964-english501c.pdf>} of
           the U{Law of Cosines<https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
           formula.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as the axes of
                    this point's datum ellipsoid).

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
           @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

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
           @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

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

    def _distanceTo(self, func, other, radius, **options):
        '''(INTERNAL) Helper for methods C{<func>To}.
        '''
        self.others(other)  # up=2
        if radius is None:
            radius = self._datum.ellipsoid.R1 if self._datum else R_M
        return func(self.lat, self.lon, other.lat, other.lon,
                                        radius=radius, **options)

    def _distanceTo_(self, func_, other, wrap=False):
        '''(INTERNAL) Helper for (ellipsoidal) methods C{<func>To}.
        '''
        self.others(other)  # up=2
        r, _ = unrollPI(self.lam, other.lam, wrap=wrap)
        r = func_(other.phi, self.phi, r, datum=self.datum)
        return r * self.datum.ellipsoid.a

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

    def equirectangularTo(self, other, radius=None, **options):
        '''Compute the distance between this and an other point
           using the U{Equirectangular Approximation / Projection
           <https://www.Movable-Type.co.UK/scripts/latlong.html#equirectangular>}.

           Suitable only for short, non-near-polar distances up to a
           few hundred Km or Miles.  Use method L{haversineTo} or
           C{distanceTo*} for more accurate and/or larger distances.

           See function L{pygeodesy.equirectangular_} for more details,
           the available B{C{options}} and errors raised.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None} for
                          the mean radius of this point's datum ellipsoid.
           @kwarg options: Optional keyword arguments for function
                           L{pygeodesy.equirectangular}.

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.equirectangular} and methods L{cosineAndoyerLambertTo},
                 L{cosineForsytheAndoyerLambertTo}, L{cosineLawTo}, C{distanceTo*},
                 C{euclideanTo}, L{flatLocalTo}/L{hubenyTo}, L{flatPolarTo},
                 L{haversineTo}, L{thomasTo} and L{vincentysTo}.
        '''
        return self._distanceTo(equirectangular, other, radius, **options)

    def euclideanTo(self, other, radius=None, **options):
        '''Approximate the C{Euclidian} distance between this and
           an other point.

           See function L{pygeodesy.euclidean} for the available B{C{options}}.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None} for
                          the mean radius of this point's datum ellipsoid.
           @kwarg options: Optional keyword arguments for function
                           L{pygeodesy.euclidean}.

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.euclidean} and methods L{cosineAndoyerLambertTo},
                 L{cosineForsytheAndoyerLambertTo}, L{cosineLawTo}, C{distanceTo*},
                 L{equirectangularTo}, L{flatLocalTo}/L{hubenyTo}, L{flatPolarTo},
                 L{haversineTo}, L{thomasTo} and L{vincentysTo}.
        '''
        return self._distanceTo(euclidean, other, radius, **options)

    def flatLocalTo(self, other, radius=None, wrap=False):
        '''Compute the distance between this and an other point using the
           U{ellipsoidal Earth to plane projection
           <https://WikiPedia.org/wiki/Geographical_distance#Ellipsoidal_Earth_projected_to_a_plane>}
           aka U{Hubeny<https://www.OVG.AT/de/vgi/files/pdf/3781/>} formula.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None} for the
                          equatorial radius of this point's datum/ellipsoid.
           @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: Invalid B{C{radius}}.

           @see: Function L{pygeodesy.flatLocal}/L{pygeodesy.hubeny}, methods
                 L{cosineAndoyerLambertTo}, L{cosineForsytheAndoyerLambertTo},
                 L{cosineLawTo}, C{distanceTo*}, L{equirectangularTo}, L{euclideanTo},
                 L{flatPolarTo}, L{haversineTo}, L{thomasTo} and L{vincentysTo} and
                 U{local, flat Earth approximation<https://www.edwilliams.org/avform.htm#flat>}.
        '''
        E = self.datum.ellipsoid
        r = self._distanceTo_(flatLocal_, other, wrap=wrap) * E.a2_
        a = E.a if radius in (None, 1, _1_0) else Radius(radius)
        return r * a

    hubenyTo = flatLocalTo  # for Karl Hubeny

    def flatPolarTo(self, other, radius=None, wrap=False):
        '''Compute the distance between this and an other point using
           the U{polar coordinate flat-Earth<https://WikiPedia.org/wiki/
           Geographical_distance#Polar_coordinate_flat-Earth_formula>}formula.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None} for the
                          mean radius of this point's datum ellipsoid.
           @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.flatPolar} and methods L{cosineAndoyerLambertTo},
                 L{cosineForsytheAndoyerLambertTo}, L{cosineLawTo}, C{distanceTo*},
                 L{equirectangularTo}, L{euclideanTo}, L{flatLocalTo}/L{hubenyTo},
                 L{haversineTo}, L{thomasTo} and L{vincentysTo}.
        '''
        return self._distanceTo(flatPolar, other, radius, wrap=wrap)

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

    def haversineTo(self, other, radius=None, wrap=False):
        '''Compute the distance between this and an other point using the
           U{Haversine<https://www.Movable-Type.co.UK/scripts/latlong.html>}
           formula.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None} for
                          the mean radius of this point's datum ellipsoid.
           @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.haversine} and methods L{cosineAndoyerLambertTo},
                 L{cosineForsytheAndoyerLambertTo}, L{cosineLawTo}, C{distanceTo*},
                 L{equirectangularTo}, L{euclideanTo}, L{flatLocalTo}/L{hubenyTo},
                 L{flatPolarTo}, L{thomasTo} and L{vincentysTo}.
        '''
        return self._distanceTo(haversine, other, radius, wrap=wrap)

    def _havg(self, other, f=_0_5):
        '''(INTERNAL) Weighted, average height.

           @arg other: An other point (C{LatLon}).
           @kwarg f: Optional fraction (C{float}).

           @return: Average, fractional height (C{float}).
        '''
        return _MODS.fmath.favg(self.height, other.height, f=f)

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

    def height4(self, earth=None, normal=True, LatLon=None, **LatLon_kwds):
        '''Compute the height above or below and the projection on this datum's
           ellipsoid surface.

           @kwarg earth: A datum, ellipsoid or earth radius I{overriding} this
                         datum (L{Datum}, L{Ellipsoid}, L{Ellipsoid2}, L{a_f2Tuple}
                         or C{meter}, conventionally).
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

           @raise TypeError: Invalid B{C{earth}}.

           @see: L{Ellipsoid.height4} for more information.
        '''
        if LatLon is None:
            r = self.toCartesian().height4(earth=earth, normal=normal)
        else:
            c = self.toCartesian()
            r = c.height4(earth=earth, normal=normal, Cartesian=c.classof, height=0)
            r = r.toLatLon(LatLon=LatLon, **_xkwds(LatLon_kwds, height=r.height))
        return r

    def heightStr(self, prec=-2, m=_m_):
        '''Return this B{C{height}} as C{str}ing.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).
           @kwarg m: Optional unit of the height (C{str}).

           @see: Function L{pygeodesy.hstr}.
        '''
        return hstr(self.height, prec=prec, m=m)

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
        return isantipode(self.lat,  self.lon,
                         other.lat, other.lon, eps=eps)

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
        self.others(other)

        return _isequalTo(self, other, eps=Scalar_(eps=eps)) if eps else \
                         (self.lat == other.lat and self.lon == other.lon)

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
        return self.height == other.height and self.isequalTo(other, eps=eps)

    @Property_RO
    def isnormal(self):
        '''Return C{True} if this point is normal (C{bool}),
           meaning C{abs(lat) <= 90} and C{abs(lon) <= 180}.

           @see: Functions L{pygeodesy.isnormal} and
                 L{pygeodesy.normal}.
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
                 string into a 3-tuple (lat, lon, h).
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
            self._lat, self._lon, self._height = llh

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

    @Property_RO
    def latlonheight(self):
        '''Get the lat-, longitude and height (L{LatLon3Tuple}C{(lat, lon, height)}).
        '''
        return self.latlon.to3Tuple(self.height)

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
           @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes
                        (C{bool}).

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
                if w and i != 0:
                    q = _unrollon(p, q)
                yield C(height=h, i=i, up=3, points=q)
                p = q

        C  = self._toCartesianEcef  # to verify datum and Ecef
        Ps = self.PointsIter(points)

        c = C(height=height, this=self)  # this Cartesian
        t = nearestOn6(c, _cs(Ps, height, wrap, C), closed=closed)
        c, s, e = t.closest, t.start, t.end

        kwds = _xkwds_not(None, LatLon=self.classof,  # this LatLon
                                height=height)
        r = self.Ecef(self.datum).reverse
        p = r(c).toLatLon(**kwds)
        s = r(s).toLatLon(**kwds) if s is not c else p
        e = r(e).toLatLon(**kwds) if e is not c else p
        return t.dup(closest=p, start=s, end=e)

    def normal(self):
        '''Normalize this point to C{abs(lat) <= 90} and C{abs(lon) <= 180}.

           @return: C{True} if this point was I{normal}, C{False}
                    if it wasn't (but is now).

           @see: Property L{isnormal} and function L{pygeodesy.normal}.
        '''
        n = self.isnormal
        if not n:
            self.latlon = normal(self.lat, self.lon)
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

           @return: A L{Points2Tuple}C{(number, points)}, C{int}
                    and C{list} or C{tuple}.

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} are not C{LatLon}.
        '''
        return points2(points, closed=closed, base=self)

    def PointsIter(self, points, loop=0, dedup=False):
        '''Return a C{PointsIter} iterator.

           @arg points: The path or polygon points (C{LatLon}[])
           @kwarg loop: Number of loop-back points (non-negative C{int}).
           @kwarg dedup: Skip duplicate points (C{bool}).

           @return: A new C{PointsIter} iterator.

           @raise PointsError: Insufficient number of B{C{points}}.
        '''
        return PointsIter(points, base=self, loop=loop, dedup=dedup)

    def radii11(self, point2, point3):
        '''Return the radii of the C{Circum-}, C{In-}, I{Soddy} and C{Tangent}
           circles of a (planar) triangle formed by this and two other points.

           @arg point2: Second point (C{LatLon}).
           @arg point3: Third point (C{LatLon}).

           @return: L{Radii11Tuple}C{(rA, rB, rC, cR, rIn, riS, roS, a, b, c, s)}.

           @raise IntersectionError: Near-coincident or -colinear points.

           @raise TypeError: Invalid B{C{point2}} or B{C{point3}}.

           @see: Function L{pygeodesy.radii11}, U{Incircle
                 <https://MathWorld.Wolfram.com/Incircle.html>}, U{Soddy Circles
                 <https://MathWorld.Wolfram.com/SoddyCircles.html>} and U{Tangent
                 Circles<https://MathWorld.Wolfram.com/TangentCircles.html>}.
        '''
        try:
            cs = self._toCartesian3(point2, point3)
            return _radii11ABC(*cs, useZ=True)[0]
        except (TypeError, ValueError) as x:
            raise _xError(x, point=self, point2=point2, point3=point3)

    def _rhumb2(self, exact, radius):
        '''(INTERNAL) Get the C{rhumb} for this point's datum  or for
           the earth model or earth B{C{radius}} if not C{None}.
        '''
        D = self.datum if radius is None else _spherical_datum(radius)  # ellipsoidal OK
        r = D.ellipsoid.rhumbx if exact else \
           _MODS.rhumbx.Rhumb(D, exact=False, name=D.name)
        return r, D

    def rhumbAzimuthTo(self, other, exact=False, radius=None):
        '''Return the azimuth (bearing) of a rhumb line (loxodrome)
           between this and an other (ellipsoidal) point.

           @arg other: The other point (C{LatLon}).
           @kwarg exact: If C{True}, use the I{exact} L{Rhumb} (C{bool}),
                         default C{False}.
           @kwarg radius: Optional earth radius (C{meter}) or earth model
                          (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                          L{a_f2Tuple}), overriding this point's datum.

           @return: Rhumb azimuth (compass C{degrees360}).

           @raise TypeError: The B{C{other}} point is incompatible or
                             B{C{radius}} is invalid.
        '''
        self.others(other)
        r, _ =  self._rhumb2(exact, radius)
        C = _MODS.rhumbx.Caps
        return r.Inverse(self.lat, self.lon, other.lat, other.lon,
                                             outmask=C.AZIMUTH).azi12

    def rhumbDestination(self, distance, azimuth, exact=False, radius=None, height=None):
        '''Return the destination point having travelled the given distance
           from this point along a rhumb line (loxodrome) at the given azimuth.

           @arg distance: Distance travelled (C{meter}, same units as this
                          point's datum (ellipsoid) axes or B{C{radius}},
                          may be negative.
           @arg azimuth: Azimuth (bearing) at this point (compass C{degrees}).
           @kwarg exact: If C{True}, use the I{exact} L{Rhumb} (C{bool}),
                         default C{False}.
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
        r, D = self._rhumb2(exact, radius)
        d = r.Direct(self.lat, self.lon, azimuth, distance)
        h = self.height if height is None else Height(height)
        return self.classof(d.lat2, d.lon2, datum=D, height=h)

    def rhumbDistanceTo(self, other, exact=False, radius=None):
        '''Return the distance from this to an other point along
           a rhumb line (loxodrome).

           @arg other: The other point (C{LatLon}).
           @kwarg exact: If C{True}, use the I{exact} L{Rhumb} (C{bool}),
                         default C{False}.
           @kwarg radius: Optional earth radius (C{meter}) or earth model
                          (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                          L{a_f2Tuple}), overriding this point's datum.

           @return: Distance (C{meter}, the same units as this point's
                    datum (ellipsoid) axes or B{C{radius}}.

           @raise TypeError: The B{C{other}} point is incompatible or
                             B{C{radius}} is invalid.

           @raise ValueError: Invalid B{C{radius}}.
        '''
        self.others(other)
        r, _ = self._rhumb2(exact, radius)
        C = _MODS.rhumbx.Caps
        return r.Inverse(self.lat, self.lon, other.lat, other.lon,
                                             outmask=C.DISTANCE).s12

    def rhumbLine(self, azimuth_other, exact=False, radius=None, name=NN, **caps):
        '''Get a rhumb line through this point at a given azimuth or
           through this and an other point.

           @arg azimuth_other: The azimuth of the rhumb line (compass)
                               C{degrees} or the other point (C{LatLon}).
           @kwarg exact: If C{True}, use the I{exact} L{Rhumb} (C{bool}),
                         default C{False}.
           @kwarg radius: Optional earth radius (C{meter}) or earth model
                          (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                          L{a_f2Tuple}), overriding this point's datum.
           @kwarg name: Optional name (C{str}).
           @kwarg caps: Optional C{caps}, see L{RhumbLine} C{B{caps}}.

           @return: A L{RhumbLine} instance.

           @raise TypeError: Invalid B{C{radius}} or BC{C{azimuth_other}}
                             not a C{scalar} nor a C{LatLon}.

           @see: Classes L{RhumbLine} and L{Rhumb}, property L{Rhumb.exact}
                 and methods L{Rhumb.DirectLine} and L{Rhumb.InverseLine}.
        '''
        r, _ = self._rhumb2(exact, radius)
        a = azimuth_other
        if isscalar(a):
            r = r.DirectLine(self.lat, self.lon, a,
                             name=name or self.name, **caps)
        elif isinstance(a, LatLonBase):
            self.others(a)
            r = r.InverseLine(self.lat, self.lon, a.lat, a.lon,
                              name=name or self.name, **caps)
        else:
            raise _TypeError(azimuth_other=a)
        return r

    def rhumbMidpointTo(self, other, exact=False, radius=None,
                                     height=None, fraction=_0_5):
        '''Return the (loxodromic) midpoint on the rhumb line between
           this and an other point.

           @arg other: The other point (C{LatLon}).
           @kwarg exact: If C{True}, use the I{exact} L{Rhumb} (C{bool}),
                         default C{False}.
           @kwarg radius: Optional earth radius (C{meter}) or earth model
                          (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                          L{a_f2Tuple}), overriding this point's datum.
           @kwarg height: Optional height, overriding the mean height
                          (C{meter}).
           @kwarg fraction: Midpoint location from this point (C{scalar}),
                            may be negative or greater than 1.0.

           @return: The midpoint at the given B{C{fraction}} along the
                    rhumb line (C{LatLon}).

           @raise TypeError: The B{C{other}} point is incompatible or
                             B{C{radius}} is invalid.

           @raise ValueError: Invalid B{C{height}} or B{C{fraction}}.
        '''
        self.others(other)
        r, D = self._rhumb2(exact, radius)
        f = Scalar(fraction=fraction)
        d = r.Inverse(self.lat, self.lon, other.lat, other.lon)
        d = r.Direct( self.lat, self.lon, d.azi12, d.s12 * f)
        h = self._havg(other, f=f) if height is None else Height(height)
        return self.classof(d.lat2, d.lon2, datum=D, height=h)

    def thomasTo(self, other, wrap=False):
        '''Compute the distance between this and an other point using
           U{Thomas'<https://apps.DTIC.mil/dtic/tr/fulltext/u2/703541.pdf>}
           formula.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: Wrap and L{pygeodesy.unrollPI} longitudes (C{bool}).

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

    def _toCartesian3(self, point2, point3):
        '''(INTERNAL) Convert this and 2 other points.
        '''
        return (self. toCartesian().copy(name=_point_),  # copy to rename
                self._toCartesianEcef(up=3, point2=point2),
                self._toCartesianEcef(up=3, point3=point3))

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

    def toNvector(self, h=None, Nvector=None, **Nvector_kwds):
        '''Convert this point to C{n-vector} (normal to the earth's
           surface) components, I{including height}.

           @kwarg h: Optional height, overriding this point's
                     height (C{meter}).
           @kwarg Nvector: Optional class to return the C{n-vector}
                           components (C{Nvector}) or C{None}.
           @kwarg Nvector_kwds: Optional, additional B{C{Nvector}}
                                keyword arguments, ignored if
                                C{B{Nvector} is None}.

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

    def vincentysTo(self, other, radius=None, wrap=False):
        '''Compute the distance between this and an other point using
           U{Vincenty's<https://WikiPedia.org/wiki/Great-circle_distance>}
           spherical formula.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None}
                          for the mean radius of this point's datum
                          ellipsoid.
           @kwarg wrap: Wrap and L{pygeodesy.unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{pygeodesy.vincentys} and methods L{cosineAndoyerLambertTo},
                 L{cosineForsytheAndoyerLambertTo}, L{cosineLawTo}, C{distanceTo*},
                 L{equirectangularTo}, L{euclideanTo}, L{flatLocalTo}/L{hubenyTo},
                 L{flatPolarTo}, L{haversineTo} and L{thomasTo}.
        '''
        return self._distanceTo(vincentys, other, radius, wrap=wrap)

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


def _isequalTo(point1, point2, eps=EPS):  # in .ellipsoidalBaseDI._intersect3._on, .formy
    '''(INTERNAL) Compare point lat-/lon without type.
    '''
    return abs(point1.lat - point2.lat) <= eps and \
           abs(point1.lon - point2.lon) <= eps


def _isequalTo_(point1, point2, eps=EPS):  # PYCHOK in .formy
    '''(INTERNAL) Compare point phi-/lam without type.
    '''
    return abs(point1.phi - point2.phi) <= eps and \
           abs(point1.lam - point2.lam) <= eps


def _trilaterate5(p1, d1, p2, d2, p3, d3, area=True, eps=EPS1,
                                          radius=R_M, wrap=False):
    '''(INTERNAL) Trilaterate three points by area overlap or by
       perimeter intersection of three circles.

       @note: The B{C{radius}} is only needed for both the n-vectorial
              and C{sphericalTrigonometry.LatLon.distanceTo} methods and
              silently ignored by the C{ellipsoidalExact}, C{-GeodSolve},
              C{-Karney} and C{-Vincenty.LatLon.distanceTo} methods.
    '''
    r1 = Distance_(distance1=d1)
    r2 = Distance_(distance2=d2)
    r3 = Distance_(distance3=d3)

    m  = 0 if area else (r1 + r2 + r3)
    pc = 0
    t  = []
    for _ in range(3):
        try:  # intersection of circle (p1, r1) and (p2, r2)
            c1, c2 = p1.intersections2(r1, p2, r2, wrap=wrap)

            if area:  # check overlap
                if c1 is c2:  # abutting
                    c = c1
                else:  # nearest point on radical
                    c = p3.nearestOn(c1, c2, within=True, wrap=wrap)
                d = r3 - p3.distanceTo(c, radius=radius, wrap=wrap)
                if d > eps:  # sufficient overlap
                    t.append((d, c))
                m = max(m, d)

            else:  # check intersection
                for c in ((c1,) if c1 is c2 else (c1, c2)):
                    d = abs(r3 - p3.distanceTo(c, radius=radius, wrap=wrap))
                    if d < eps:  # below margin
                        t.append((d, c))
                    m = min(m, d)

        except IntersectionError as x:
            if _concentric_ in str(x):  # XXX ConcentricError?
                pc += 1

        p1, r1, p2, r2, p3, r3 = p2, r2, p3, r3, p1, r1  # rotate

    if t:  # get min, max, points and count ...
        t =  tuple(sorted(t))
        n = (len(t),)  # as 1-tuple
        # ... or for a single trilaterated result,
        # min *is* max, min- *is* maxPoint and n=1
        return Trilaterate5Tuple(t[0] + t[-1] + n)  # *(t[0] + ...)

    elif area and pc == 3:  # all pairwise concentric ...
        r, p = min((r1, p1), (r2, p2), (r3, p3))
        m = max(r1, r2, r3)
        # ... return "smallest" point twice, the smallest
        # and largest distance and n=0 for concentric
        return Trilaterate5Tuple(float(r), p, float(m), p, 0)

    n, f = (_overlap_, max) if area else (_intersection_, min)
    t = '%s (%s %.3f)' % (_no_(n), f.__name__, m)
    raise IntersectionError(area=area, eps=eps, wrap=wrap, txt=t)


__all__ += _ALL_DOCS(LatLonBase)

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
