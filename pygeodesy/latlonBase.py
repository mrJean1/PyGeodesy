
# -*- coding: utf-8 -*-

u'''(INTERNAL) Base class L{LatLonBase} for elliposiodal, spherical and
N-vectorial C{LatLon}s.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**,
see U{https://www.Movable-Type.co.UK/scripts/latlong.html},
U{<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}
and U{https://www.Movable-Type.co.UK/scripts/latlong-vectors.html}.
'''

from pygeodesy.basics import isstr, _xinstanceof
from pygeodesy.dms import F_D, F_DMS, latDMS, lonDMS, parse3llh
from pygeodesy.errors import _datum_datum, IntersectionError, \
                             _ValueError, _xkwds
from pygeodesy.fmath import favg
from pygeodesy.formy import antipode, compassAngle, cosineAndoyerLambert_, \
                            cosineForsytheAndoyerLambert_, cosineLaw, \
                            equirectangular, euclidean, flatLocal_, \
                            flatPolar, haversine, isantipode, \
                            latlon2n_xyz, thomas_, vincentys
from pygeodesy.interns import EPS, EPS0, EPS1, NN, R_M, _COMMASPACE_, \
                             _intersection_, _m_, _near_concentric_, \
                             _no_, _overlap_, _0_0, _0_5, _1_0
from pygeodesy.iters import PointsIter, points2
from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.named import _NamedBase, notOverloaded
from pygeodesy.namedTuples import Bounds2Tuple, LatLon2Tuple, PhiLam2Tuple, \
                                  Trilaterate5Tuple, Vector3Tuple
from pygeodesy.props import deprecated_method, Property, Property_RO, \
                            property_doc_, property_RO
from pygeodesy.streprs import Fmt, hstr
from pygeodesy.units import Distance_, Lat, Lon, Height, Radius, Radius_, Scalar_
from pygeodesy.utily import unrollPI
from pygeodesy.vector3d import Vector3d

from math import asin, cos, degrees, radians

__all__ = ()
__version__ = '21.06.24'


class LatLonBase(_NamedBase):
    '''(INTERNAL) Base class for C{LatLon} points on spherical or
       ellipsoidal earth models.
    '''
    _datum  = None  # L{Datum}, to be overriden
    _height = 0     # height (C{meter})
    _lat    = 0     # latitude (C{degrees})
    _lon    = 0     # longitude (C{degrees})

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

    @deprecated_method
    def compassAngle(self, other, adjust=True, wrap=False):  # PYCHOK no cover
        '''DEPRECATED, use method C{compassAngleTo}.'''
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
           @kwarg wrap: Wrap and L{unroll180} longitudes and longitudinal
                       delta (C{bool}).

           @return: Compass angle from North (C{degrees360}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @note: Courtesy Martin Schultz.

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
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as the axes of
                    this point's datum ellipsoid).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{cosineAndoyerLambert} and methods
                 C{cosineForsytheAndoyerLambertTo}, C{cosineLawTo},
                 C{distanceTo*}, C{equirectangularTo}, C{euclideanTo},
                 C{flatLocalTo}/C{hubenyTo}, C{flatPolarTo},
                 C{haversineTo}, C{thomasTo} and C{vincentysTo}.
        '''
        return self._distanceTo_(cosineAndoyerLambert_, other, wrap=wrap)

    def cosineForsytheAndoyerLambertTo(self, other, wrap=False):
        '''Compute the distance between this and an other point using
           the U{Forsythe-Andoyer-Lambert correction
           <https://www2.UNB.CA/gge/Pubs/TR77.pdf>} of the U{Law of Cosines
           <https://www.Movable-Type.co.UK/scripts/latlong.html#cosine-law>}
           formula.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as the axes of
                    this point's datum ellipsoid).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{cosineForsytheAndoyerLambert} and methods
                 C{cosineAndoyerLambertTo}, C{cosineLawTo}, C{distanceTo*},
                 C{equirectangularTo}, C{euclideanTo},
                 C{flatLocalTo}/C{hubenyTo}, C{flatPolarTo},
                 C{haversineTo}, C{thomasTo} and C{vincentysTo}.
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
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{cosineLaw} and methods C{cosineAndoyerLambertTo},
                 C{cosineForsytheAndoyerLambertTo}, C{distanceTo*},
                 C{equirectangularTo}, C{euclideanTo}, C{flatLocalTo}/C{hubenyTo},
                 C{flatPolarTo}, C{haversineTo}, C{thomasTo} and C{vincentysTo}.
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
                               arguments, ignored is C{B{LatLon}=None}.

           @return: Destination as a C{B{LatLon}(lat, lon, **B{LatLon_kwds})}
                    instance or if C{B{LatLon}=None}, a L{LatLon3Tuple}C{(lat,
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
        from pygeodesy.ecef import EcefKarney
        return EcefKarney  # default

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
           <https://www.Movable-Type.co.UK/scripts/latlong.html>}.

           Suitable only for short, non-near-polar distances up to a
           few hundred Km or Miles.  Use method C{haversineTo} or
           C{distanceTo*} for more accurate and/or larger distances.

           See function L{equirectangular_} for more details, the
           available B{C{options}} and errors raised.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None} for
                          the mean radius of this point's datum ellipsoid.
           @kwarg options: Optional keyword arguments for function
                           L{equirectangular}.

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{equirectangular} and methods C{cosineAndoyerLambertTo},
                 C{cosineForsytheAndoyerLambertTo}, C{cosineLawTo},
                 C{distanceTo*}, C{euclideanTo}, C{flatLocalTo}/C{hubenyTo},
                 C{flatPolarTo}, C{haversineTo}, C{thomasTo} and C{vincentysTo}.
        '''
        return self._distanceTo(equirectangular, other, radius, **options)

    def euclideanTo(self, other, radius=None, **options):
        '''Approximate the C{Euclidian} distance between this and
           an other point.

           See function L{euclidean} for the available B{C{options}}.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None} for
                          the mean radius of this point's datum ellipsoid.
           @kwarg options: Optional keyword arguments for function
                             L{euclidean}.

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{euclidean} and methods C{cosineAndoyerLambertTo},
                 C{cosineForsytheAndoyerLambertTo}, C{cosineLawTo},
                 C{distanceTo*}, C{equirectangularTo},
                 C{flatLocalTo}/C{hubenyTo}, C{flatPolarTo},
                 C{haversineTo}, C{thomasTo} and C{vincentysTo}.
        '''
        return self._distanceTo(euclidean, other, radius, **options)

    def flatLocalTo(self, other, radius=None, wrap=False):
        '''Compute the distance between this and an other point using the
           U{ellipsoidal Earth to plane projection
           <https://WikiPedia.org/wiki/Geographical_distance#Ellipsoidal_Earth_projected_to_a_plane>}
           aka U{Hubeny<https://www.OVG.AT/de/vgi/files/pdf/3781/>} formula.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None} for the
                          major radius of this point's datum/ellipsoid.
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: Invalid B{C{radius}}.

           @see: Function L{flatLocal}/L{hubeny}, methods C{cosineAndoyerLambertTo},
                 C{cosineForsytheAndoyerLambertTo}, C{cosineLawTo},
                 C{distanceTo*}, C{equirectangularTo}, C{euclideanTo},
                 C{flatPolarTo}, C{haversineTo}, C{thomasTo} and
                 C{vincentysTo} and U{local, flat Earth approximation
                 <https://www.edwilliams.org/avform.htm#flat>}.
        '''
        E = self.datum.ellipsoid
        r = self._distanceTo_(flatLocal_, other, wrap=wrap) * E.a2_
        a = E.a if radius in (None, 1, _1_0) else Radius(radius)
        return r * a

    hubenyTo = flatLocalTo  # for Karl Hubeny

    def flatPolarTo(self, other, radius=None, wrap=False):
        '''Compute the distance between this and an other point using
           the U{polar coordinate flat-Earth
           <https://WikiPedia.org/wiki/Geographical_distance#Polar_coordinate_flat-Earth_formula>}
           formula.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None}
                          for the mean radius of this point's datum
                          ellipsoid.
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{flatPolar} and methods C{cosineAndoyerLambertTo},
                 C{cosineForsytheAndoyerLambertTo}, C{cosineLawTo},
                 C{distanceTo*}, C{equirectangularTo}, C{euclideanTo},
                 C{flatLocalTo}/C{hubenyTo}, C{haversineTo},
                 C{thomasTo} and C{vincentysTo}.
        '''
        return self._distanceTo(flatPolar, other, radius, wrap=wrap)

    def haversineTo(self, other, radius=None, wrap=False):
        '''Compute the distance between this and an other point using the
           U{Haversine<https://www.Movable-Type.co.UK/scripts/latlong.html>}
           formula.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None} for
                          the mean radius of this point's datum ellipsoid.
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{haversine} and methods C{cosineAndoyerLambertTo},
                 C{cosineForsytheAndoyerLambertTo}, C{cosineLawTo},
                 C{distanceTo*}, C{equirectangularTo}, C{euclideanTo},
                 C{flatLocalTo}/C{hubenyTo}, C{flatPolarTo},
                 C{thomasTo} and C{vincentysTo}.
        '''
        return self._distanceTo(haversine, other, radius, wrap=wrap)

    def _havg(self, other, f=_0_5):
        '''(INTERNAL) Weighted, average height.

           @arg other: An other point (C{LatLon}).
           @kwarg f: Optional fraction (C{float}).

           @return: Average, fractional height (C{float}).
        '''
        return favg(self.height, other.height, f=f)

    @property_doc_(''' the height (C{meter}).''')
    def height(self):
        '''Get the height (C{meter}).
        '''
        return self._height

    @height.setter  # PYCHOK setter!
    def height(self, height):
        '''Set the height.

           @arg height: New height (C{meter}).

           @raise TypeError: Invalid B{C{height}} C{type}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        h = Height(height)
        self._update(h != self.height)
        self._height = h

    def heightStr(self, prec=-2, m=_m_):
        '''Return a string for the height B{C{height}}.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg m: Optional unit of the height (C{str}).

           @see: Function L{hstr}.
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

           @example:

            >>> p = LatLon(52.205, 0.119)
            >>> q = LatLon(52.205, 0.119)
            >>> e = p.isequalTo(q)  # True
        '''
        self.others(other)

        if eps:
            return max(abs(self.lat - other.lat),
                       abs(self.lon - other.lon)) < Scalar_(eps=eps)
        else:
            return self.lat == other.lat and \
                   self.lon == other.lon

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

           @example:

            >>> p = LatLon(52.205, 0.119, 42)
            >>> q = LatLon(52.205, 0.119)
            >>> e = p.isequalTo3(q)  # False
        '''
        return self.height == other.height and self.isequalTo(other, eps=eps)

    @Property_RO
    def isSpherical(self):
        '''Check whether this point is spherical (C{bool} or C{None} if unknown).
        '''
        return self.datum.isSpherical if self._datum else None

    @Property_RO
    def lam(self):
        '''Get the longitude (B{C{radians}}).
        '''
        return self.philam.lam

    @property_doc_(''' the latitude (C{degrees90}).''')
    def lat(self):
        '''Get the latitude (C{degrees90}).
        '''
        return self._lat

    @lat.setter  # PYCHOK setter!
    def lat(self, lat):
        '''Set the latitude.

           @arg lat: New latitude (C{str[N|S]} or C{degrees}).

           @raise ValueError: Invalid B{C{lat}}.
        '''
        lat = Lat(lat)  # parseDMS(lat, suffix=_NS_, clip=90)
        self._update(lat != self._lat)
        self._lat = lat

    @Property
    def latlon(self):
        '''Get the lat- and longitude (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return LatLon2Tuple(self._lat, self._lon, name=self.name)

    @latlon.setter  # PYCHOK setter!
    def latlon(self, latlonh):
        '''Set the lat- and longitude and optionally the height.

           @arg latlonh: New lat-, longitude and height (2- or
                         3-tuple or comma- or space-separated C{str}
                         of C{degrees90}, C{degrees180} and C{meter}).

           @raise TypeError: Height of B{C{latlonh}} not C{scalar} or
                             B{C{latlonh}} not C{list} or C{tuple}.

           @raise ValueError: Invalid B{C{latlonh}} or M{len(latlonh)}.

           @see: Function L{parse3llh} to parse a B{C{latlonh}} string
                 into a 3-tuple (lat, lon, h).
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

        lat = Lat(latlonh[0])  # parseDMS2(latlonh[0], latlonh[1])
        lon = Lon(latlonh[1])
        self._update(lat != self._lat or
                     lon != self._lon or h != self.height)
        self._lat, self._lon, self._height = lat, lon, h

    def latlon2(self, ndigits=0):
        '''Return this point's lat- and longitude in C{degrees}, rounded.

           @kwarg ndigits: Number of decimal digits (C{int}).

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

    @property_doc_(''' the longitude (C{degrees180}).''')
    def lon(self):
        '''Get the longitude (C{degrees180}).
        '''
        return self._lon

    @lon.setter  # PYCHOK setter!
    def lon(self, lon):
        '''Set the longitude.

           @arg lon: New longitude (C{str[E|W]} or C{degrees}).

           @raise ValueError: Invalid B{C{lon}}.
        '''
        lon = Lon(lon)  # parseDMS(lon, suffix=_EW_, clip=180)
        self._update(lon != self._lon)
        self._lon = lon

    @Property_RO
    def _ltp(self):
        '''(INTERNAL) Cache for L{toLtp}.
        '''
        from pygeodesy.ltp import Ltp
        return Ltp(self, ecef=self.Ecef(self.datum), name=self.name)

    @Property_RO
    def _N_vector(self):
        '''(INTERNAL) Get the (C{nvectorBase._N_vector_})
        '''
        from pygeodesy.nvectorBase import _N_vector_
        return _N_vector_(*self.xyzh)

    @Property_RO
    def phi(self):
        '''Get the latitude (B{C{radians}}).
        '''
        return self.philam.phi

    @Property_RO
    def philam(self):
        '''Get the lat- and longitude (L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        return PhiLam2Tuple(radians(self.lat),
                            radians(self.lon), name=self.name)

    def philam2(self, ndigits=0):
        '''Return this point's lat- and longitude in C{radians}, rounded.

           @kwarg ndigits: Number of decimal digits (C{int}).

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

    def thomasTo(self, other, wrap=False):
        '''Compute the distance between this and an other point using
           U{Thomas'<https://apps.DTIC.mil/dtic/tr/fulltext/u2/703541.pdf>}
           formula.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: Wrap and L{unrollPI} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as the axes of
                    this point's datum ellipsoid).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{thomas} and methods C{cosineAndoyerLambertTo},
                 C{cosineForsytheAndoyerLambertTo}, C{cosineLawTo}, C{distanceTo*},
                 C{equirectangularTo}, C{euclideanTo}, C{flatLocalTo}/C{hubenyTo},
                 C{flatPolarTo}, C{haversineTo} and C{vincentysTo}.
        '''
        return self._distanceTo_(thomas_, other, wrap=wrap)

    @deprecated_method
    def to2ab(self):  # PYCHOK no cover
        '''DEPRECATED, use property L{philam}.'''
        return self.philam

    def toCartesian(self, Cartesian=None, **Cartesian_kwds):
        '''Convert this point to cartesian, I{geocentric} coordinates,
           also known as I{Earth-Centered, Earth-Fixed} (ECEF).

           @kwarg Cartesian: Optional class to return the geocentric
                             coordinates (C{Cartesian}) or C{None}.
           @kwarg Cartesian_kwds: Optional, additional B{C{Cartesian}}
                                  keyword arguments, ignored if
                                  C{B{Cartesian}=None}.

           @return: A B{C{Cartesian}} or if B{C{Cartesian}} is C{None},
                    an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M,
                    datum)} with C{C=0} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian}} or B{C{Cartesian_kwds}}.
        '''
        r = self._ecef9
        if Cartesian is not None:  # class or .classof
            r = self._xnamed(Cartesian(r, **Cartesian_kwds))
        _datum_datum(r.datum, self.datum)
        return r

    def toEcef(self, height=None, M=False):
        '''Convert this point to I{geocentric} coordinates,
           also known as I{Earth-Centered, Earth-Fixed}
           (U{ECEF<https://WikiPedia.org/wiki/ECEF>}).

           @kwarg height: Optional height, overriding this point's
                          height (C{meter}).
           @kwarg M: Optionally, include the rotation L{EcefMatrix}
                     (C{bool}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} 0 and C{M} if available.

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
                            arguments, ignored if C{B{Xyz}=None}.

           @return: An B{C{Xyz}} instance or if C{B{Xyz}=None},
                    a L{Local9Tuple}C{(x, y, z, lat, lon, height,
                    ltp, ecef, M)} with C{M=None}, always.

           @raise TypeError: Invalid B{C{ltp}}.
        '''
        p = self._ltp if ltp is None else self._xLtp(ltp)
        return p._ecef2local(self._ecef9, Xyz, Xyz_kwds)

    def toLtp(self, Ecef=None):
        '''Return the I{local tangent plane} (LTP) for this point.

           @kwarg Ecef: Optional ECEF I{class} (L{EcefKarney}, ...
                        L{EcefYou}), overriding this point's C{Ecef}.
        '''
        if Ecef in (None, self.Ecef):
            r = self._ltp
        else:
            from pygeodesy.ltp import Ltp
            r = Ltp(self, ecef=Ecef(self.datum), name=self.name)
        return r

    def toNvector(self, h=None, Nvector=None, **Nvector_kwds):
        '''Convert this point to C{n-vector} (normal to the earth's
           surface) components, I{including height}.

           @kwarg h: Optional height, overriding this point's
                     height (C{meter}).
           @kwarg Nvector: Optional class to return the C{n-vector}
                           components (C{Nvector}) or C{None}.
           @kwarg Nvector_kwds: Optional, additional B{C{Nvector}}
                                keyword arguments, ignored if
                                C{B{Nvector}=None}.

           @return: A B{C{Nvector}} or an L{Vector4Tuple}C{(x, y, z, h)}
                    if B{C{Nvector}} is C{None}.

           @raise TypeError: Invalid B{C{Nvector}} or B{C{Nvector_kwds}}.
        '''
        return self.toVector(Vector=Nvector, h=self.height if h is None else h,
                                            ll=self, **Nvector_kwds)

    def toStr(self, form=F_DMS, prec=None, m=_m_, sep=_COMMASPACE_):  # PYCHOK expected
        '''Convert this point to a "lat, lon [+/-height]" string,
           formatted in the given form.

           @kwarg form: Optional format, F_D, F_DM, F_DMS for
                        deg°, deg°min′, deg°min′sec″ (C{str}).
           @kwarg prec: Optional number of decimal digits (0..8 or C{None}).
           @kwarg m: Optional unit of the height (C{str}), use C{None} to
                     exclude height from the returned string.
           @kwarg sep: Optional separator to join (C{str}).

           @return: Point in the specified form (C{str}).

           @example:

            >>> LatLon(51.4778, -0.0016).toStr()  # 51°28′40″N, 000°00′06″W
            >>> LatLon(51.4778, -0.0016).toStr(F_D)  # 51.4778°N, 000.0016°W
            >>> LatLon(51.4778, -0.0016, 42).toStr()  # 51°28′40″N, 000°00′06″W, +42.00m

        '''
        t = [latDMS(self.lat, form=form, prec=prec),
             lonDMS(self.lon, form=form, prec=prec)]
        if self.height and m is not None:
            t += [self.heightStr(m=m)]
        return sep.join(t)

    def toVector(self, Vector=None, **Vector_kwds):
        '''Convert this point to C{n-vector} (normal to the earth's
           surface) components, I{ignoring height}.

           @kwarg Vector: Optional class to return the C{n-vector}
                          components (L{Vector3d}) or C{None}.
           @kwarg Vector_kwds: Optional, additional B{C{Vector}}
                               keyword arguments, ignored if
                               C{B{Vector}=None}.

           @return: A B{C{Vector}} or a L{Vector3Tuple}C{(x, y, z)}
                    if B{C{Vector}} is C{None}.

           @raise TypeError: Invalid B{C{Vector}} or B{C{kwds}}.

           @note: These are C{n-vector} x, y and z components,
                  I{NOT} geocentric (ECEF) x, y and z coordinates!
        '''
        r = self._vector3tuple
        if Vector is not None:
            r = self._xnamed(Vector(r.x, r.y, r.z, **Vector_kwds))
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
        return latlon2n_xyz(self.lat, self.lon, name=self.name)

    def vincentysTo(self, other, radius=None, wrap=False):
        '''Compute the distance between this and an other point using
           U{Vincenty's<https://WikiPedia.org/wiki/Great-circle_distance>}
           spherical formula.

           @arg other: The other point (C{LatLon}).
           @kwarg radius: Mean earth radius (C{meter}) or C{None}
                          for the mean radius of this point's datum
                          ellipsoid.
           @kwarg wrap: Wrap and L{unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{vincentys} and methods C{cosineAndoyerLambertTo},
                 C{cosineForsytheAndoyerLambertTo}, C{cosineLawTo},
                 C{distanceTo*}, C{equirectangularTo}, C{euclideanTo},
                 C{flatLocalTo}/C{hubenyTo}, C{flatPolarTo},
                 C{haversineTo} and C{thomasTo}.
        '''
        return self._distanceTo(vincentys, other, radius, wrap=wrap)

    @Property_RO
    def _xLtp(self):
        '''(INTERNAL) Import and cache function C{ltp._xLtp}.
        '''
        from pygeodesy.ltp import _xLtp
        return _xLtp

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


def _trilaterate5(p1, d1, p2, d2, p3, d3, area=True, eps=EPS1,
                                          radius=R_M, wrap=False):
    '''(INTERNAL) Trilaterate three points by area overlap or by
       perimeter intersection of three circles.

       @note: The B{C{radius}} is only needed for both n-vectorial and
              the sphericalTrigonometry C{LatLon.distanceTo} methods and
              silently ignored by the C{ellipsoidalExact} C{-GeodSolve},
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
            if _near_concentric_ in str(x):  # XXX ConcentricError?
                pc += 1

        p1, r1, p2, r2, p3, r3 = p2, r2, p3, r3, p1, r1  # rotate

    if t:  # get min, max, points and count ...
        t = tuple(sorted(t))
        n = len(t),  # as tuple
        # ... or for a single trilaterated result,
        # min *is* max, min- *is* maxPoint and n=1
        return Trilaterate5Tuple(*(t[0] + t[-1] + n))

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
