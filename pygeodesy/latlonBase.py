
# -*- coding: utf-8 -*-

u'''(INTERNAL) base classes for elliposiodal, spherical and vectorial C{LatLon}s.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**,
see U{https://www.Movable-Type.co.UK/scripts/latlong.html},
U{<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}
and U{https://www.Movable-Type.co.UK/scripts/latlong-vectors.html}.

@newfield example: Example, Examples
'''

from pygeodesy.basics import EPS, R_M, map1, property_doc_, property_RO, \
                            _xinstanceof
from pygeodesy.dms import F_D, F_DMS, latDMS, lonDMS  # parseDMS, parseDMS2
from pygeodesy.ecef import EcefKarney
from pygeodesy.errors import _datum_datum, _ValueError
from pygeodesy.fmath import favg
from pygeodesy.formy import antipode, compassAngle, cosineAndoyerLambert_, \
                            cosineForsytheAndoyerLambert_, cosineLaw, \
                            equirectangular, euclidean, flatLocal_, \
                            flatPolar, haversine, isantipode, \
                            latlon2n_xyz,points2, thomas_, vincentys
from pygeodesy.interns import _COMMA_SPACE_, NN, _item_sq, _m_
from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.named import Bounds2Tuple, LatLon2Tuple, _NamedBase, \
                            notOverloaded, PhiLam2Tuple, Vector3Tuple
from pygeodesy.streprs import hstr
from pygeodesy.units import Lat, Lon, Height, Radius, Radius_, Scalar_
from pygeodesy.utily import unrollPI
from pygeodesy.vector3d import Vector3d

from math import asin, cos, degrees, radians

__all__ = ()
__version__ = '20.07.12'


class LatLonBase(_NamedBase):
    '''(INTERNAL) Base class for C{LatLon} points on spherical or
       ellipsoidal earth models.
    '''
    _datum  = None        #: (INTERNAL) L{Datum}, to be overriden.
    _Ecef   = EcefKarney  #: (INTERNAL) Preferred C{Ecef...} class.
    _e9t    = None        #: (INTERNAL) Cached toEcef (L{Ecef9Tuple}).
    _height = 0           #: (INTERNAL) Height (C{meter}).
    _lat    = 0           #: (INTERNAL) Latitude (C{degrees}).
    _latlon = None        #: (INTERNAL) Cached (L{LatLon2Tuple}).
    _lon    = 0           #: (INTERNAL) Longitude (C{degrees}).
    _name   = NN          #: (INTERNAL) name (C{str}).
    _philam = None        #: (INTERNAL) Cached (L{PhiLam2Tuple}).
    _v3d    = None        #: (INTERNAL) Cached toVector3d (L{Vector3d}).
    _xyz    = None        #: (INTERNAL) Cached xyz (L{Vector3Tuple}).
    _xyzh   = None        #: (INTERNAL) Cached xyzh (L{Vector4Tuple}).

    def __init__(self, lat, lon, height=0, name=NN):
        '''New C{LatLon}.

           @arg lat: Latitude (C{degrees} or DMS C{str} with N or S suffix).
           @arg lon: Longitude (C{degrees} or DMS C{str} with E or W suffix).
           @kwarg height: Optional height (C{meter} above or below the earth surface).
           @kwarg name: Optional name (C{str}).

           @return: New instance (C{LatLon}).

           @raise RangeError: Value of B{C{lat}} or B{C{lon}} outside the valid
                              range and C{rangerrors} set to C{True}.

           @raise ValueError: Invalid B{C{lat}} or B{C{lon}}.

           @example:

           >>> p = LatLon(50.06632, -5.71475)
           >>> q = LatLon('50°03′59″N', """005°42'53"W""")
        '''
        self._lat = Lat(lat)  # parseDMS2(lat, lon)
        self._lon = Lon(lon)  # PYCHOK LatLon2Tuple
        if height:  # elevation
            self._height = Height(height)
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

           @arg other: An other point (C{LatLon}).
           @kwarg f: Optional fraction (C{float}).

           @return: Average, fractional height (C{float}).
        '''
        return favg(self.height, other.height, f=f)

    def _update(self, updated, *attrs):
        '''(INTERNAL) Zap cached attributes if updated.
        '''
        if updated:
            _NamedBase._update(self, updated, '_e9t', '_latlon', '_philam',
                                      '_v3d', '_xyz', '_xyzh', *attrs)

    def antipode(self, height=None):
        '''Return the antipode, the point diametrically opposite
           to this point.

           @kwarg height: Optional height of the antipode, height
                            of this point otherwise (C{meter}).

           @return: The antipodal point (C{LatLon}).
        '''
        a, b = antipode(self.lat, self.lon)  # PYCHOK LatLon2Tuple
        h = self.height if height is None else height
        return self.classof(a, b, height=h)

    def bounds(self, wide, high, radius=R_M):  # PYCHOK no cover
        '''DEPRECATED, use method C{boundsOf}.
        '''
        return self.boundsOf(wide, high, radius=radius)

    def boundsOf(self, wide, high, radius=R_M):
        '''Return the SE and NW lat-/longitude of a great circle
           bounding box centered at this location.

           @arg wide: Longitudinal box width (C{meter}, same units as
                      B{C{radius}} or C{degrees} if B{C{radius}} is C{None}).
           @arg high: Latitudinal box height (C{meter}, same units as
                      B{C{radius}} or C{degrees} if B{C{radius}} is C{None}).
           @kwarg radius: Mean earth radius (C{meter}).

           @return: A L{Bounds2Tuple}C{(latlonSW, latlonNE)}, the
                    lower-left and upper-right corner (C{LatLon}).

           @see: U{https://www.Movable-Type.co.UK/scripts/latlong-db.html}
        '''
        w = Scalar_(wide, name='wide') * 0.5
        h = Scalar_(high, name='high') * 0.5
        if radius is not None:
            r = Radius_(radius)
            c = cos(self.phi)
            if c > EPS:
                w = degrees(asin(w / r) / c)
            else:
                w = 0  # XXX
            h = degrees(h / r)
        w, h = abs(w), abs(h)

        r = Bounds2Tuple(self.classof(self.lat - h, self.lon - w, height=self.height),
                         self.classof(self.lat + h, self.lon + w, height=self.height))
        return self._xnamed(r)

    def compassAngle(self, other):  # PYCHOK no cover
        '''DEPRECATED, use method C{compassAngleTo}.
        '''
        return self.compassAngleTo(other)

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
        '''(INTERNAL) I{Must be overloaded}.
        '''
        notOverloaded(self, self.datum)

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

    @property_RO
    def Ecef(self):
        '''Get the ECEF C{class} (L{EcefKarney} or L{EcefVeness}).
        '''
        return self._Ecef

    def equals(self, other, eps=None):  # PYCHOK no cover
        '''DEPRECATED, use method C{isequalTo}.
        '''
        return self.isequalTo(other, eps=eps)

    def equals3(self, other, eps=None):  # PYCHOK no cover
        '''DEPRECATED, use method C{isequalTo3}.
        '''
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
        a = E.a if radius in (None, 1.0, 1) else Radius(radius)
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
        self._update(h != self._height)
        self._height = h

    def heightStr(self, prec=-2, m=_m_):
        '''Return a string for the height B{C{height}}.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg m: Optional unit of the height (C{str}).

           @see: Function L{hstr}.
        '''
        return hstr(self.height, prec=prec, m=m)

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

    def isantipode(self, other, eps=EPS):  # PYCHOK no cover
        '''DEPRECATED, use method C{isantipodeTo}.
        '''
        return self.isantipodeTo(other, eps=eps)

    @property_RO
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

           @raise ValueError: Invalid B{C{eps}}.

           @see: Method L{isequalTo3}.

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(52.205, 0.119)
           >>> e = p.isequalTo(q)  # True
        '''
        self.others(other)

        e = 0 if eps in (None, 0, 0.0) else Scalar_(eps, name='eps')
        if e > 0:
            return max(map1(abs, self.lat - other.lat,
                                 self.lon - other.lon)) < e
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
        return self.isequalTo(other, eps=eps) and self.height == other.height

    @property_RO
    def isSpherical(self):
        '''Check whether this point is spherical (C{bool} or C{None} if unknown).
        '''
        return self.datum.isSpherical if self._datum else None

    @property_RO
    def lam(self):
        '''Get the longitude (B{C{radians}}).
        '''
        return self.philam.lam if self._philam is None else self._philam.lam

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

    @property_doc_(''' the lat- and longitude, optionally height.''')
    def latlon(self):
        '''Get the lat- and longitude (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        if self._latlon is None:
            self._latlon = LatLon2Tuple(self._lat, self._lon)
        return self._xrenamed(self._latlon)

    @latlon.setter  # PYCHOK setter!
    def latlon(self, latlonh):
        '''Set the lat- and longitude and optionally the height.

           @arg latlonh: New lat-, longitude and height (2- or
                        3-tuple of C{degrees} and C{meter}).

           @raise TypeError: Height of B{C{latlonh}} not C{scalar} or
                             B{C{latlonh}} not C{list} or C{tuple}.

           @raise ValueError: Invalid B{C{latlonh}} or M{len(latlonh)}.

           @see: Function L{parse3llh} to parse a B{C{latlonh}} string
                 into a 3-tuple (lat, lon, h).
        '''
        _xinstanceof(list, tuple, latlonh=latlonh)

        if len(latlonh) == 3:
            h = Height(latlonh[2], name=_item_sq(latlonh=2))
        elif len(latlonh) != 2:
            raise _ValueError(latlonh=latlonh)
        else:
            h = self._height

        lat = Lat(latlonh[0])  # parseDMS2(latlonh[0], latlonh[1])
        lon = Lon(latlonh[1])
        self._update(lat != self._lat or
                     lon != self._lon or h != self._height)
        self._lat, self._lon, self._height = lat, lon, h

    def latlon_(self, ndigits=0):  # PYCHOK no cover
        '''DEPRECATED, use method C{latlon2}.
        '''
        return self.latlon2(ndigits)

    def latlon2(self, ndigits=0):
        '''Return this point's lat- and longitude in C{degrees}, rounded.

           @kwarg ndigits: Number of decimal digits (C{int}).

           @return: A L{LatLon2Tuple}C{(lat, lon)}, both C{float}
                    and rounded away from zero.

           @note: The C{round}ed values are always C{float}, also
                  if B{C{ndigits}} is omitted.
        '''
        r = LatLon2Tuple(round(self.lat, ndigits),
                         round(self.lon, ndigits))
        return self._xnamed(r)

    def latlon2round(self, ndigits=0):  # PYCHOK no cover
        '''DEPRECATED, use method C{latlon2}.
        '''
        return self.latlon2(ndigits)

    @property_RO
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

    @property_RO
    def _N_vector(self):
        '''(INTERNAL) Get the (C{nvectorBase._N_vector_})
        '''
        from pygeodesy.nvectorBase import _N_vector_
        r = self._xyz or self._v3d or self.toVector()
        return _N_vector_(r.x, r.y, r.z, h=self.height)

    @property_RO
    def phi(self):
        '''Get the latitude (B{C{radians}}).
        '''
        return self.philam.phi if self._philam is None else self._philam.phi

    @property_RO
    def philam(self):
        '''Get the lat- and longitude (L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        if self._philam is None:
            self._philam = PhiLam2Tuple(radians(self.lat),
                                        radians(self.lon))
        return self._xnamed(self._philam)

    def philam2(self, ndigits=0):
        '''Return this point's lat- and longitude in C{radians}, rounded.

           @kwarg ndigits: Number of decimal digits (C{int}).

           @return: A L{PhiLam2Tuple}C{(phi, lam)}, both C{float}
                    and rounded away from zero.

           @note: The C{round}ed values are always C{float}, also
                  if B{C{ndigits}} is omitted.
        '''
        r = PhiLam2Tuple(round(self.phi, ndigits),
                         round(self.lam, ndigits))
        return self._xnamed(r)

    @property_RO
    def philamheight(self):
        '''Get the lat-, longitude in C{radians} and height (L{PhiLam3Tuple}C{(phi, lam, height)}).
        '''
        return self.philam.to3Tuple(self.height)

    def points(self, points, closed=True):  # PYCHOK no cover
        '''DEPRECATED, use method C{points2}.
        '''
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

    def to2ab(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{philam}.

           @return: A L{PhiLam2Tuple}C{(phi, lam)}.
        '''
        return self.philam

    def to3llh(self, height=None):  # PYCHOK no cover
        '''DEPRECATED, use property C{latlonheight} or C{latlon.to3Tuple}C{(}B{C{height}}C{)}.

           @return: A L{LatLon3Tuple}C{(lat, lon, height)}.
        '''
        return self.latlonheight if height in (None, self.height) else \
               self.latlon.to3Tuple(height)

    def to3xyz(self):  # PYCHOK no cover
        '''DEPRECATED, use method C{toNvector}, C{toVector}, C{toVector3d}
           or property C{xyz} or perhaps (geocentric) C{toEcef}.

           @return: A L{Vector3Tuple}C{(x, y, z)}, see property C{xyz}.
        '''
        return self.xyz  # self.toVector()

    def toCartesian(self, Cartesian=None, **Cartesian_kwds):
        '''Convert this point to cartesian (ECEF) coordinates.

           @kwarg Cartesian: Optional class to return the geocentric
                             coordinates (C{Cartesian}) or C{None}.
           @kwarg Cartesian_kwds: Optional, additional B{C{Cartesian}}
                                  keyword arguments, ignored if
                                  B{C{Cartesian=None}}.

           @return: A B{C{Cartesian}} or if B{C{Cartesian}} is C{None},
                    an L{Ecef9Tuple}C{(x, y, z, lat, lon, height, C, M,
                    datum)} with C{C=0} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian}} or B{C{Cartesian_kwds}}.
        '''
        r = self.toEcef()
        if Cartesian is not None:  # class or .classof
            r = Cartesian(r, **Cartesian_kwds)
            r = self._xnamed(r)
        _datum_datum(r.datum, self.datum)
        return r

    def toEcef(self):
        '''Convert this point to geocentric coordinates, also
           known as I{Earth-Centered, Earth-Fixed} (U{ECEF
           <https://WikiPedia.org/wiki/ECEF>}).

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} 0 and C{M} if available.

           @raise EcefError: A C{.datum} or an ECEF issue.
        '''
        if self._e9t is None:
            e = self.Ecef(self.datum).forward(self, M=True)
            self._e9t = self._xnamed(e)
        return self._e9t

    def toNvector(self, h=None, Nvector=None, **Nvector_kwds):
        '''Convert this point to C{n-vector} (normal to the earth's
           surface) components, I{including height}.

           @kwarg h: Optional height, overriding this point's
                     height (C{meter}).
           @kwarg Nvector: Optional class to return the C{n-vector}
                           components (C{Nvector}) or C{None}.
           @kwarg Nvector_kwds: Optional, additional B{C{Nvector}}
                                keyword arguments, ignored if
                                B{C{Nvector=None}}.

           @return: A B{C{Nvector}} or an L{Vector4Tuple}C{(x, y, z, h)}
                    if B{C{Nvector}} is C{None}.

           @raise TypeError: Invalid B{C{Nvector}} or B{C{Nvector_kwds}}.
        '''
        return self.toVector(Vector=Nvector, h=self.height if h is None else h,
                                            ll=self, **Nvector_kwds)

    def toStr(self, form=F_DMS, prec=None, m=_m_, sep=_COMMA_SPACE_):  # PYCHOK expected
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
                               B{C{Vector=None}}.

           @return: A B{C{Vector}} or a L{Vector3Tuple}C{(x, y, z)}
                    if B{C{Vector}} is C{None}.

           @raise TypeError: Invalid B{C{Vector}} or B{C{kwds}}.

           @note: These are C{n-vector} x, y and z components,
                  I{NOT} geocentric (ECEF) x, y and z coordinates!
        '''
        r = latlon2n_xyz(self.lat, self.lon)
        if Vector is not None:
            r = Vector(r.x, r.y, r.z, **Vector_kwds)
        return self._xnamed(r)

    def toVector3d(self):
        '''Convert this point to C{n-vector} (normal to the earth's
           surface) components, I{ignoring height}.

           @return: Unit vector (L{Vector3d}).

           @note: These are C{n-vector} x, y and z components,
                  I{NOT} geocentric (ECEF) x, y and z coordinates!
        '''
        if self._v3d is None:
            self._v3d = self.toVector(Vector=Vector3d)  # XXX .unit()
        return self._xnamed(self._v3d)

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

    @property_RO
    def xyz(self):
        '''Get the C{n-vector} X, Y and Z components (L{Vector3Tuple}C{(x, y, z)})

           @note: These are C{n-vector} x, y and z components, I{NOT}
                  geocentric (ECEF) x, y and z coordinates!
        '''
        if self._xyz is None:
            self._xyz = self.toVector(Vector=Vector3Tuple)
        return self._xnamed(self._xyz)

    @property_RO
    def xyzh(self):
        '''Get the C{n-vector} X, Y, Z and H components (L{Vector4Tuple}C{(x, y, z, h)})

           @note: These are C{n-vector} x, y and z components, I{NOT}
                  geocentric (ECEF) x, y and z coordinates!
        '''
        if self._xyzh is None:
            self._xyzh = self.xyz.to4Tuple(self.height)
        return self._xnamed(self._xyzh)


__all__ += _ALL_DOCS(LatLonBase)

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
