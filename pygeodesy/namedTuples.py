
# -*- coding: utf-8 -*-

u'''(INTERNAL) Classes C{_Named}, C{_NamedDict} and C{_NamedTuple},
all with nameable instances and several subclasses thereof.

In addition, the items in a C{_NamedDict} are accessable as attributes
and the items in a C{_NamedTuple} can be named to be accessable as
attributes, similar to standard Python C{namedtuple}s.

Results previously returned as tuples by C{pygeodesy} functions and
class methods are now instances of some C{...Tuple} class, all
sub-classes of C{_NamedTuple} defined here.

@newfield example: Example, Examples
'''

# update imported names under if __name__ == '__main__':
from pygeodesy.basics import _xinstanceof
from pygeodesy.interns import _datum_, _distance_, _easting_, _h_, \
                              _height_, _lam_, _lat_, _lon_, _n_, \
                              _northing_, _number_, _phi_, _points_, \
                              _precision_, _radius_, _x_, _y_, _z_
from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.named import _NamedTuple, _Pass
from pygeodesy.units import Bearing, Degrees, Degrees2, Easting, Height, \
                            Lam, Lat, Lon, Meter, Northing, Number_, Phi, \
                            Precision_, Radius, Scalar

__all__ = ()
__version__ = '20.09.27'

# __DUNDER gets mangled in class
_final_   = 'final'
_initial_ = 'initial'


class Bearing2Tuple(_NamedTuple):
    '''2-Tuple C{(initial, final)} bearings, both in compass C{degrees360}.
    '''
    _Names_ = (_initial_, _final_)
    _Units_ = ( Bearing,   Bearing)


class Bounds2Tuple(_NamedTuple):  # .geohash.py, .latlonBase.py, .points.py
    '''2-Tuple C{(latlonSW, latlonNE)} with the bounds' lower-left and
       upper-right corner as C{LatLon} instance.
    '''
    _Names_ = ('latlonSW', 'latlonNE')
    _Units_ = (_Pass,      _Pass)


class Bounds4Tuple(_NamedTuple):  # .geohash.py, .points.py
    '''4-Tuple C{(latS, lonW, latN, lonE)} with the bounds' lower-left
       C{(LatS, LowW)} and upper-right C{(latN, lonE)} corner lat- and
       longitudes.
    '''
    _Names_ = ('latS', 'lonW', 'latN', 'lonE')
    _Units_ = ( Lat,    Lon,    Lat,    Lon)


class Destination2Tuple(_NamedTuple):  # .ellipsoidalKarney.py, -Vincenty.py
    '''2-Tuple C{(destination, final)}, C{destination} in C{LatLon}
       and C{final} bearing in compass C{degrees360}.
    '''
    _Names_ = ('destination', _final_)
    _Units_ = (_Pass,          Bearing)


class Destination3Tuple(_NamedTuple):  # .karney.py
    '''3-Tuple C{(lat, lon, final)}, destination C{lat}, C{lon} in
       and C{final} bearing in compass C{degrees360}.
    '''
    _Names_ = (_lat_, _lon_, _final_)
    _Units_ = ( Lat,   Lon,   Bearing)


class Distance2Tuple(_NamedTuple):  # .datum.py, .ellipsoidalBase.py
    '''2-Tuple C{(distance, initial)}, C{distance} in C{meter} and
       C{initial} bearing in compass C{degrees360}.
    '''
    _Names_ = (_distance_, _initial_)
    _Units_ = ( Meter,      Bearing)


class Distance3Tuple(_NamedTuple):  # .ellipsoidalKarney.py, -Vincenty.py
    '''3-Tuple C{(distance, initial, final)}, C{distance} in C{meter}
       and C{initial} and C{final} bearing, both in compass C{degrees360}.
    '''
    _Names_ = (_distance_, _initial_, _final_)
    _Units_ = ( Meter,      Bearing,   Bearing)


class Distance4Tuple(_NamedTuple):  # .formy.py, .points.py
    '''4-Tuple C{(distance2, delta_lat, delta_lon, unroll_lon2)} with
       the distance in C{degrees squared}, the latitudinal C{delta_lat}
       = B{C{lat2}}-B{C{lat1}}, the wrapped, unrolled and adjusted
       longitudinal C{delta_lon} = B{C{lon2}}-B{C{lon1}} and
       C{unroll_lon2}, the unrolled or original B{C{lon2}}.

       @note: Use Function L{degrees2m} to convert C{degrees squared}
              to C{meter} as M{degrees2m(sqrt(distance2), ...)} or
              M{degrees2m(hypot(delta_lat, delta_lon), ...)}.
    '''
    _Names_ = ('distance2', 'delta_lat', 'delta_lon', 'unroll_lon2')
    _Units_ = ( Degrees2,    Degrees,     Degrees,     Degrees)


class EasNor2Tuple(_NamedTuple):  # .css.py, .osgr.py, .ups.py, .utm.py, .utmupsBase.py
    '''2-Tuple C{(easting, northing)}, both in C{meter},
       conventionally.
    '''
    _Names_ = (_easting_, _northing_)
    _Units_ = ( Easting,   Northing)


class EasNor3Tuple(_NamedTuple):  # .css.py, .lcc.py
    '''3-Tuple C{(easting, northing, height)}, all in C{meter},
       conventionally.
    '''
    _Names_ = (_easting_, _northing_, _height_)
    _Units_ = ( Easting,   Northing,   Height)


class LatLon2Tuple(_NamedTuple):
    '''2-Tuple C{(lat, lon)} in C{degrees90} and C{degrees180}.
    '''
    _Names_ = (_lat_, _lon_)
    _Units_ = ( Lat,   Lon)

    def to3Tuple(self, height):
        '''Extend this L{LatLon2Tuple} to a L{LatLon3Tuple}.

           @arg height: The height to add (C{scalar}).

           @return: A L{LatLon3Tuple}C{(lat, lon, height)}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        return self._xtend(LatLon3Tuple, height)

    def to4Tuple(self, height, datum):
        '''Extend this L{LatLon2Tuple} to a L{LatLon4Tuple}.

           @arg height: The height to add (C{scalar}).
           @arg datum: The datum to add (C{Datum}).

           @return: A L{LatLon4Tuple}C{(lat, lon, height, datum)}.

           @raise TypeError: If B{C{datum}} not a C{Datum}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        return self.to3Tuple(height).to4Tuple(datum)


class LatLon3Tuple(_NamedTuple):
    '''3-Tuple C{(lat, lon, height)} in C{degrees90}, C{degrees180}
       and C{meter}, conventionally.
    '''
    _Names_ = (_lat_, _lon_, _height_)
    _Units_ = ( Lat,   Lon,   Height)

    def to4Tuple(self, datum):
        '''Extend this L{LatLon3Tuple} to a L{LatLon4Tuple}.

           @arg datum: The datum to add (C{Datum}).

           @return: A L{LatLon4Tuple}C{(lat, lon, height, datum)}.

           @raise TypeError: If B{C{datum}} not a C{Datum}.
        '''
        from pygeodesy.datums import Datum
        _xinstanceof(Datum, datum=datum)
        return self._xtend(LatLon4Tuple, datum)


class LatLon4Tuple(_NamedTuple):  # .cartesianBase.py, .css.py, .ecef.py, .lcc.py
    '''4-Tuple C{(lat, lon, height, datum)} in C{degrees90},
       C{degrees180}, C{meter} and L{Datum}.
    '''
    _Names_ = (_lat_, _lon_, _height_, _datum_)
    _Units_ = ( Lat,   Lon,   Height,  _Pass)


class LatLonDatum3Tuple(_NamedTuple):  # .lcc.py, .osgr.py
    '''3-Tuple C{(lat, lon, datum)} in C{degrees90}, C{degrees180}
       and L{Datum}.
    '''
    _Names_ = (_lat_, _lon_, _datum_)
    _Units_ = ( Lat,   Lon,  _Pass)


class LatLonPrec3Tuple(_NamedTuple):  # .gars.py, .wgrs.py
    '''3-Tuple C{(lat, lon, precision)} in C{degrees}, C{degrees}
       and C{int}.
    '''
    _Names_ = (_lat_, _lon_, _precision_)
    _Units_ = ( Lat,   Lon,   Precision_)

    def to5Tuple(self, height, radius):
        '''Extend this L{LatLonPrec3Tuple} to a L{LatLonPrec5Tuple}.

           @arg height: The height to add (C{float} or C{None}).
           @arg radius: The radius to add (C{float} or C{None}).

           @return: A L{LatLonPrec5Tuple}C{(lat, lon, precision,
                    height, radius)}.
        '''
        return self._xtend(LatLonPrec5Tuple, height, radius)


class LatLonPrec5Tuple(_NamedTuple):  # .wgrs.py
    '''5-Tuple C{(lat, lon, precision, height, radius)} in C{degrees},
       C{degrees}, C{int} and C{height} or C{radius} in C{meter} (or
       C{None} if missing).
    '''
    _Names_ = (_lat_, _lon_, _precision_, _height_, _radius_)
    _Units_ = ( Lat,   Lon,   Precision_,  Height,   Radius)


class NearestOn3Tuple(_NamedTuple):  # .points.py, .sphericalTrigonometry.py
    '''3-Tuple C{(closest, distance, angle)} of the C{closest}
       point on the polygon, either a C{LatLon} instance or a
       L{LatLon3Tuple}C{(lat, lon, height)} and the C{distance}
       and C{angle} to the C{closest} point are in C{meter}
       respectively compass C{degrees360}.
    '''
    _Names_ = ('closest', _distance_, 'angle')
    _Units_ = (_Pass,      Meter,      Degrees)


class PhiLam2Tuple(_NamedTuple):  # .frechet.py, .hausdorff.py, .latlonBase.py, .points.py, .vector3d.py
    '''2-Tuple C{(phi, lam)} with latitude C{phi} in C{radians[PI_2]}
       and longitude C{lam} in C{radians[PI]}.

       @note: Using C{phi/lambda} for lat-/longitude in C{radians}
              follows Chris Veness' U{convention
              <https://www.Movable-Type.co.UK/scripts/latlong.html>}.
    '''
    _Names_ = (_phi_, _lam_)
    _Units_ = ( Phi,   Lam)

    def to3Tuple(self, height):
        '''Extend this L{PhiLam2Tuple} to a L{PhiLam3Tuple}.

           @arg height: The height to add (C{scalar}).

           @return: A L{PhiLam3Tuple}C{(phi, lam, height)}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        return self._xtend(PhiLam3Tuple, height)

    def to4Tuple(self, height, datum):
        '''Extend this L{PhiLam2Tuple} to a L{PhiLam4Tuple}.

           @arg height: The height to add (C{scalar}).
           @arg datum: The datum to add (C{Datum}).

           @return: A L{PhiLam4Tuple}C{(phi, lam, height, datum)}.

           @raise TypeError: If B{C{datum}} not a C{Datum}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        return self.to3Tuple(height).to4Tuple(datum)


class PhiLam3Tuple(_NamedTuple):  # .nvector.py, extends -2Tuple
    '''3-Tuple C{(phi, lam, height)} with latitude C{phi} in
       C{radians[PI_2]}, longitude C{lam} in C{radians[PI]} and
       C{height} in C{meter}.

       @note: Using C{phi/lambda} for lat-/longitude in C{radians}
              follows Chris Veness' U{convention
              <https://www.Movable-Type.co.UK/scripts/latlong.html>}.
    '''
    _Names_ = (_phi_, _lam_, _height_)
    _Units_ = ( Phi,   Lam,   Height)

    def to4Tuple(self, datum):
        '''Extend this L{PhiLam3Tuple} to a L{PhiLam4Tuple}.

           @arg datum: The datum to add (C{Datum}).

           @return: A L{PhiLam4Tuple}C{(phi, lam, height, datum)}.

           @raise TypeError: If B{C{datum}} not a C{Datum}.
        '''
        from pygeodesy.datums import Datum
        _xinstanceof(Datum, datum=datum)
        return self._xtend(PhiLam4Tuple, datum)


class PhiLam4Tuple(_NamedTuple):  # extends -3Tuple
    '''4-Tuple C{(phi, lam, height, datum)} with latitude C{phi} in
       C{radians[PI_2]}, longitude C{lam} in C{radians[PI]}, C{height}
       in C{meter} and L{Datum}.

       @note: Using C{phi/lambda} for lat-/longitude in C{radians}
              follows Chris Veness' U{convention
              <https://www.Movable-Type.co.UK/scripts/latlong.html>}.
    '''
    _Names_ = (_phi_, _lam_, _height_, _datum_)
    _Units_ = ( Phi,   Lam,   Height,  _Pass)


class Points2Tuple(_NamedTuple):  # .formy.py, .latlonBase.py
    '''2-Tuple C{(number, points)} with the C{number} of points
       and -possible reduced- C{list} or C{tuple} of C{points}.
    '''
    _Names_ = (_number_, _points_)
    _Units_ = ( Number_, _Pass)


class Trilaterate5Tuple(_NamedTuple):  # .latlonBase.py, .nvector.py
    '''5-Tuple C{(min, minPoint, max, maxPoint, n)} with C{min} and C{max}
       in C{meter}, the corresponding trilaterated C{minPoint} and C{maxPoint}
       as C{LatLon} and the number C{n}.  For area overlap, C{min} and C{max}
       are the smallest respectively largest overlap found.  For perimeter
       intersection, C{min} and C{max} represent the closest respectively
       farthest intersection margin.  Count C{n} is the total number of
       trilaterated overlaps or intersections found, C{0, 1, 2...6} with
       C{0} meaning concentric.

       @see: The C{ellipsoidalKarney-}, C{ellipsoidalVincenty-} and
             C{sphericalTrigonometry.LatLon.trilaterate5} method for further
             details on corner cases, like concentric or single trilaterated
             results.
   '''
    _Names_ = (min.__name__, 'minPoint', max.__name__, 'maxPoint', _n_)
    _Units_ = (Meter,        _Pass,      Meter,        _Pass,       Number_)


class Vector3Tuple(_NamedTuple):
    '''3-Tuple C{(x, y, z)} of (geocentric) components, all in
       C{meter} or C{units}.
    '''
    _Names_ = (_x_,    _y_,    _z_)
    _Units_ = ( Scalar, Scalar, Scalar)

    def to4Tuple(self, h):
        '''Extend this L{Vector3Tuple} to a L{Vector4Tuple}.

           @arg h: The height to add (C{scalar}).

           @return: A L{Vector4Tuple}C{(x, y, z, h)}.

           @raise ValueError: Invalid B{C{h}}.
        '''
        return self._xtend(Vector4Tuple, h)


class Vector4Tuple(_NamedTuple):  # .nvector.py
    '''4-Tuple C{(x, y, z, h)} of (geocentric) components, all
       in C{meter} or C{units}.
    '''
    _Names_ = (_x_,    _y_,    _z_,    _h_)
    _Units_ = ( Scalar, Scalar, Scalar, Height)


__all__ += _ALL_DOCS(Bearing2Tuple, Bounds2Tuple, Bounds4Tuple,
                     Destination2Tuple, Destination3Tuple,
                     Distance2Tuple, Distance3Tuple, Distance4Tuple,
                     EasNor2Tuple, EasNor3Tuple,
                     LatLon2Tuple, LatLon3Tuple, LatLon4Tuple,
                     LatLonDatum3Tuple, LatLonPrec3Tuple, LatLonPrec5Tuple,
                     NearestOn3Tuple,
                     PhiLam2Tuple, PhiLam3Tuple, PhiLam4Tuple, Points2Tuple,
                     Trilaterate5Tuple,
                     Vector3Tuple, Vector4Tuple)

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

# % env PYGEODESY_FOR_DOCS=1 python -m pygeodesy.named
# all 71 locals OK
