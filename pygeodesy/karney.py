
# -*- coding: utf-8 -*-

u'''Wrapper around several C{geomath.Math} functions from I{Karney}'s Python package
U{geographiclib<https://PyPI.org/project/geographiclib>}, provided that package is installed.

The I{wrapped} class methods return a L{GDict} instance offering access to the C{dict} items
either by C{key} or by C{attribute} name.

With env variable C{PYGEODESY_GEOGRAPHICLIB} left undefined or set to C{"2"}, this module,
L{pygeodesy.geodesicw} and L{pygeodesy.geodesicx} will use U{GeographicLib 2.0
<https://GeographicLib.SourceForge.io/C++/doc/>} and newer transcoding, otherwise C{1.52}
or older.

Karney-based functionality
==========================

1. The following classes and functions in C{pygeodesy}

  - L{AlbersEqualArea}, L{AlbersEqualArea2}, L{AlbersEqualArea4},
    L{AlbersEqualAreaCylindrical}, L{AlbersEqualAreaNorth}, L{AlbersEqualAreaSouth} --
    U{AlbersEqualArea<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AlbersEqualArea.html>}

  - L{AuxAngle}, L{AuxDST}, L{AuxDLat}, L{AuxLat} -- U{AuxAngle
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AuxAngle.html>},
    U{DST<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1DST.html>},
    U{DAuxLatitude<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1DAuxLatitude.html>},
    U{AuxLatitude<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AuxLatitude.html>} in I{GeographicLib 2.2+}

  - L{CassiniSoldner} -- U{CassiniSoldner<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1CassiniSoldner.html>}

  - L{EcefKarney} -- U{Geocentric<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Geocentric.html>}

  - L{Elliptic} -- U{EllipticFunction<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1EllipticFunction.html>}

  - L{EquidistantExact}, L{EquidistantGeodSolve}, L{EquidistantKarney} -- U{AzimuthalEquidistant
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AzimuthalEquidistant.html>}

  - L{Etm}, L{ExactTransverseMercator} -- U{TransverseMercatorExact
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1TransverseMercatorExact.html>}

  - L{Geodesic}, L{GeodesicLine} -- I{wrapped} U{geodesic.Geodesic<https://PyPI.org/project/geographiclib>},
    I{wrapped} U{geodesicline.GeodesicLine<https://PyPI.org/project/geographiclib>}

  - L{GeodesicAreaExact}, L{PolygonArea} -- U{PolygonArea<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1PolygonAreaT.html>}

  - L{GeodesicExact}, L{GeodesicLineExact} -- U{GeodesicExact<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1GeodesicExact.html>}, U{GeodesicLineExact<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1GeodesicLineExact.html>}

  - L{GeoidKarney} -- U{Geoid<https://GeographicLib.SourceForge.io/C++/doc/geoid.html>}

  - L{Georef} -- U{Georef<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Georef.html>}

  - L{GnomonicExact}, L{GnomonicGeodSolve}, L{GnomonicKarney} -- U{Gnomonic
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Gnomonic.html>}

  - L{JacobiConformal} -- U{JacobiConformal
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1JacobiConformal.html#details>}

  - L{KTransverseMercator} - U{TransverseMercator
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1TransverseMercator.html>}

  - L{LocalCartesian}, L{Ltp} -- U{LocalCartesian<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1LocalCartesian.html>}

  - L{Osgr} -- U{OSGB<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1OSGB.html>}

  - L{rhumbaux}, L{RhumbAux}, L{RhumbLineAux} -- U{Rhumb
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Rhumb.html>} and U{RhumbLine
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1RhumbLine.html>} from I{GeographicLib 2.2+}

  - L{rhumbx}, L{Rhumb}, L{RhumbLine} -- U{Rhumb
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Rhumb.html>},
    U{RhumbLine<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1RhumbLine.html>},
    U{TransverseMercator<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1TransverseMercator.html>}
    from I{GeographicLib 2.0}

  - L{Ups} -- U{PolarStereographic<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1PolarStereographic.html>}

  - L{Utm} -- U{TransverseMercator<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1TransverseMercator.html>}

  - L{UtmUps}, L{Epsg} -- U{UTMUPS<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1UTMUPS.html>}

  - L{pygeodesy.atand}, L{pygeodesy.atan2d}, L{pygeodesy.sincos2}, L{pygeodesy.sincos2d}, L{pygeodesy.tand} -- U{geomath.Math
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Math.html>}

are I{transcoded} from C++ classes in I{Karney}'s U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/annotated.html>}.

2. These C{pygeodesy} modules and classes

  - L{ellipsoidalGeodSolve}, L{ellipsoidalKarney}, L{geodsolve}, L{karney}, L{rhumbsolve}
  - L{EquidistantKarney}, L{FrechetKarney}, L{GeodesicSolve}, L{GeodesicLineSolve}, L{GnomonicGeodSolve},
    L{GnomonicKarney}, L{HeightIDWkarney}

are or use I{wrappers} around I{Karney}'s Python U{geographiclib<https://PyPI.org/project/geographiclib>}
C{geodesic}, C++ utility U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>} or
C++ utility U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>}.

3. All C{pygeodesy} functions and methods to compute I{ellipsoidal} intersections, nearest points and trilaterations

  - L{ellipsoidalExact.intersection3}, L{ellipsoidalExact.intersections2}, L{ellipsoidalExact.nearestOn},
    L{ellipsoidalExact.LatLon.intersection3}, L{ellipsoidalExact.LatLon.intersections2},
    L{ellipsoidalExact.LatLon.nearestOn}, L{ellipsoidalExact.LatLon.trilaterate5}

  - L{ellipsoidalKarney.intersection3}, L{ellipsoidalKarney.intersections2}, L{ellipsoidalKarney.nearestOn},
    L{ellipsoidalKarney.LatLon.intersection3}, L{ellipsoidalKarney.LatLon.intersections2},
    L{ellipsoidalKarney.LatLon.nearestOn}, L{ellipsoidalKarney.LatLon.trilaterate5}

  - L{ellipsoidalVincenty.intersection3}, L{ellipsoidalVincenty.intersections2}, L{ellipsoidalVincenty.nearestOn},
    L{ellipsoidalVincenty.LatLon.intersection3}, L{ellipsoidalVincenty.LatLon.intersections2},
    L{ellipsoidalVincenty.LatLon.nearestOn}, L{ellipsoidalVincenty.LatLon.trilaterate5}

  - L{RhumbLineAux.intersection2} and L{RhumbLineAux.nearestOn4} C{(exact=None)} in L{rhumbaux} and
    L{RhumbLine.intersection2} and L{RhumbLine.nearestOn4} C{(exact=None)} in L{rhumbx}

are implementations of I{Karney}'s iterative solution posted under U{The B{ellipsoidal} case
<https://GIS.StackExchange.com/questions/48937/calculating-intersection-of-two-circles>} and in paper U{Geodesics
on an ellipsoid of revolution<https://ArXiv.org/pdf/1102.1215.pdf>} (pp 20-21, section B{14. MARITIME BOUNDARIES}).

4. Spherical functions

  - L{pygeodesy.excessKarney_}, L{sphericalTrigonometry.areaOf}

in C{pygeodesy} are based on I{Karney}'s post U{Area of a spherical polygon
<https://MathOverflow.net/questions/97711/the-area-of-spherical-polygons>}, 3rd Answer.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import _copysign, neg, unsigned0, _xgeographiclib, \
                             _xImportError, _xversion_info, _xinstanceof, \
                             _zip,  isodd  # PYCHOK shared
from pygeodesy.constants import NAN, _isfinite as _math_isfinite, _0_0, _1_16th, \
                               _1_0, _2_0, _180_0, _N_180_0, _360_0
from pygeodesy.datums import _a_ellipsoid, _EWGS84, _WGS84  # PYCHOK shared
# from pygeodesy.ellipsoids import _EWGS84  # from .datums
from pygeodesy.errors import GeodesicError, _ValueError, _xkwds, _xkwds_get
from pygeodesy.fmath import cbrt, fremainder, hypot as _hypot, norm2, \
                            Fsum, unstr  # PYCHOK shared
# from pygeodesy.fsums import Fsum  # from .fmath
from pygeodesy.interns import _2_, _a12_, _area_, _azi2_, _azi12_, _composite_, \
                              _lat1_, _lat2_, _lon1_, _lon2_, _m12_, _M12_, _M21_, \
                              _number_, _s12_, _S12_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS, _getenv
from pygeodesy.named import _Dict, _NamedBase, _NamedTuple, notImplemented, _Pass
from pygeodesy.props import deprecated_method, Property_RO, property_RO
# from pygeodesy.streps import unstr  # from .fmath
from pygeodesy.units import Bearing as _Azi, Degrees as _Deg, Lat, Lon, \
                            Meter as _M, Meter2 as _M2, Number_, \
                            Precision_, _1mm as _TOL_M  # PYCHOK shared
from pygeodesy.utily import atan2d, sincos2d, tand, _unrollon,  fabs

# from math import fabs  # from .utily

__all__ = _ALL_LAZY.karney
__version__ = '23.09.15'

_K_2_0      = _getenv('PYGEODESY_GEOGRAPHICLIB', _2_) == _2_
_perimeter_ = 'perimeter'


class _GTuple(_NamedTuple):  # in .testNamedTuples
    '''(INTERNAL) Helper.
    '''
    def toGDict(self, **updates):
        '''Convert this C{*Tuple} to a L{GDict}.

           @kwarg updates: Optional items to apply (C{nam=value} pairs)
        '''
        r = GDict(_zip(self._Names_, self))  # strict=True
        if updates:
            r.update(updates)
        if self._iteration is not None:
            r._iteration = self._iteration
        return r


class _Lat(Lat):
    '''(INTERNAL) Latitude B{C{lat}}.
    '''
    def __init__(self, *lat, **Error_name):
        kwds = _xkwds(Error_name, clip=0, Error=GeodesicError)
        Lat.__new__(_Lat, *lat, **kwds)


class _Lon(Lon):
    '''(INTERNAL) Longitude B{C{lon}}.
    '''
    def __init__(self, *lon, **Error_name):
        kwds = _xkwds(Error_name, clip=0, Error=GeodesicError)
        Lon.__new__(_Lon, *lon, **kwds)


class Area3Tuple(_NamedTuple):  # in .geodesicx.gxarea
    '''3-Tuple C{(number, perimeter, area)} with the C{number}
       of points on the polygon or polyline, the C{perimeter} in
       C{meter} and the C{area} in C{meter} I{squared}.
    '''
    _Names_ = (_number_, _perimeter_, _area_)
    _Units_ = ( Number_, _M,          _M2)


class Caps(object):  # PYCHOK
    '''(INTERNAL) Overriden by C{Caps} below.
    '''
    EMPTY          =  0        # formerly aka NONE
    _CAP_1         =  1 << 0   # for goedesicw only
    _CAP_1p        =  1 << 1   # for goedesicw only
#   _CAP_2         =  1 << 2
    _CAP_3         =  1 << 3   # for goedesicw only
#   _CAP_4         =  1 << 4
#   _CAP_ALL       =  0x1F
#   _CAP_MASK      = _CAP_ALL
    LATITUDE       =  1 << 7   # compute latitude C{lat2}
    LONGITUDE      =  1 << 8   # compute longitude C{lon2}  _CAP_3
    AZIMUTH        =  1 << 9   # azimuths C{azi1} and C{azi2}
    DISTANCE       =  1 << 10  # compute distance C{s12}  _CAP_1
    DISTANCE_IN    =  1 << 11  # allow distance C{s12} in Direct  _CAP_1 | _CAP_1p
    REDUCEDLENGTH  =  1 << 12  # compute reduced length C{m12}    _CAP_1 | _CAP_2
    GEODESICSCALE  =  1 << 13  # compute geodesic scales C{M12} and C{M21}  _CAP_1 | _CAP_2
    AREA           =  1 << 14  # compute area C{S12}  _CAP_4

    STANDARD       =  AZIMUTH | DISTANCE | DISTANCE_IN | LATITUDE | LONGITUDE
    ALL            =  0x7F80   # without LONG_UNROLL, LINE_OFF, REVERSE2 and _DEBUG_*

    _DIRECT3       =  AZIMUTH  | LATITUDE | LONGITUDE | _CAP_3   # for goedesicw only
    _INVERSE3      =  AZIMUTH  | DISTANCE | _CAP_1   # for goedesicw only
    _STD           =  STANDARD | _CAP_3   | _CAP_1   # for goedesicw only
    _STD_LINE      = _STD      | _CAP_1p   # for goedesicw only

    LINE_OFF       =  1 << 15  # Line without updates from parent geodesic or rhumb
    LONG_UNROLL    =  1 << 16  # unroll C{lon2} in .Direct and .Position
    REVERSE2       =  1 << 17  # reverse C{azi2}

    LATITUDE_LONGITUDE         = LATITUDE | LONGITUDE
    LATITUDE_LONGITUDE_AREA    = LATITUDE | LONGITUDE | AREA

    AZIMUTH_DISTANCE           = AZIMUTH | DISTANCE
    AZIMUTH_DISTANCE_AREA      = AZIMUTH | DISTANCE | AREA

    _ANGLE_ONLY    =  1 << 18  # angular distance C{a12} only
    _SALPs_CALPs   =  1 << 19  # (INTERNAL) GeodesicExact._GenInverse

    _DEBUG_AREA    =  1 << 20  # (INTERNAL) include Line details
    _DEBUG_DIRECT  =  1 << 21  # (INTERNAL) include Direct details
    _DEBUG_INVERSE =  1 << 22  # (INTERNAL) include Inverse details
    _DEBUG_LINE    =  1 << 23  # (INTERNAL) include Line details
    _DEBUG_ALL     = _DEBUG_AREA | _DEBUG_DIRECT | _DEBUG_INVERSE | \
                     _DEBUG_LINE | _ANGLE_ONLY | _SALPs_CALPs
    _OUT_ALL       =  ALL
    _OUT_MASK      =  ALL | LONG_UNROLL | REVERSE2 | _DEBUG_ALL

    _AZIMUTH_LATITUDE_LONGITUDE           =  AZIMUTH | LATITUDE | LONGITUDE
    _DEBUG_DIRECT_LINE                    = _DEBUG_DIRECT | _DEBUG_LINE
    _DISTANCE_IN_OUT                      =  DISTANCE_IN & _OUT_MASK
    _LINE                                 =  AZIMUTH | LATITUDE | LONG_UNROLL
    _REDUCEDLENGTH_GEODESICSCALE          =  REDUCEDLENGTH | GEODESICSCALE
#   _REDUCEDLENGTH_GEODESICSCALE_DISTANCE =  REDUCEDLENGTH | GEODESICSCALE | DISTANCE

Caps = Caps()  # PYCHOK singleton
'''I{Enum}-style masks to be bit-C{or}'ed to specify geodesic or
rhumb capabilities (C{caps}) and expected results (C{outmask}).

C{AREA} - compute area C{S12},

C{AZIMUTH} - include azimuths C{azi1} and C{azi2},

C{DISTANCE} - compute distance C{s12},

C{DISTANCE_IN} - allow distance C{s12} in C{.Direct},

C{EMPTY} - nothing, formerly aka C{NONE},

C{GEODESICSCALE} - compute geodesic scales C{M12} and C{M21},

C{LATITUDE} - compute latitude C{lat2},

C{LINE_OFF} - Line without updates from parent geodesic or rhumb.

C{LONGITUDE} - compute longitude C{lon2},

C{LONG_UNROLL} - unroll C{lon2} in C{.Direct},

C{REDUCEDLENGTH} - compute reduced length C{m12},

C{REVERSE2} - reverse C{azi2},

and C{ALL} - all of the above.

C{STANDARD} = C{AZIMUTH | DISTANCE | DISTANCE_IN | LATITUDE | LONGITUDE}'''


class _CapsBase(_NamedBase):  # in .auxilats, .geodesicx.gxbases
    '''(INTERNAL) Base class for C{[_]Geodesic*Exact}.
    '''
    ALL           = Caps.ALL
    AREA          = Caps.AREA
    AZIMUTH       = Caps.AZIMUTH
    DISTANCE      = Caps.DISTANCE
    DISTANCE_IN   = Caps.DISTANCE_IN
    EMPTY         = Caps.EMPTY  # aka NONE
    GEODESICSCALE = Caps.GEODESICSCALE
    LATITUDE      = Caps.LATITUDE
    LINE_OFF      = Caps.LINE_OFF
    LONGITUDE     = Caps.LONGITUDE
    LONG_UNROLL   = Caps.LONG_UNROLL
    REDUCEDLENGTH = Caps.REDUCEDLENGTH
    STANDARD      = Caps.STANDARD

    _caps         = 0  # None
    _debug        = 0  # or Caps._DEBUG_...

    @Property_RO
    def caps(self):
        '''Get the capabilities (bit-or'ed C{Caps}).
        '''
        return self._caps

    def caps_(self, caps):
        '''Check the available capabilities.

           @arg caps: Bit-or'ed combination of L{Caps} values
                      for all capabilities to be checked.

           @return: C{True} if I{all} B{C{caps}} are available,
                    C{False} otherwise (C{bool}).
        '''
        caps &= Caps._OUT_ALL
        return (self.caps & caps) == caps

    @property
    def debug(self):
        '''Get the C{debug} option (C{bool}).
        '''
        return bool(self._debug)

    @debug.setter  # PYCHOK setter!
    def debug(self, debug):
        '''Set the C{debug} option (C{bool}) to include
           more details in L{GDict} results.
        '''
        self._debug = Caps._DEBUG_ALL if debug else 0

    def _iter2tion(self, r, s):
        '''(INTERNAL) Copy C{C{s}.iter} into C{B{r}._iteration}.
        '''
        i = _xkwds_get(s, iter=None)
        if i is not None:
            self._iteration = r._iteration = i
        return r


class Direct9Tuple(_GTuple):
    '''9-Tuple C{(a12, lat2, lon2, azi2, s12, m12, M12, M21, S12)} with arc
       length C{a12}, angles C{lat2}, C{lon2} and azimuth C{azi2} in C{degrees},
       distance C{s12} and reduced length C{m12} in C{meter} and area C{S12} in
       C{meter} I{squared}.
    '''
    _Names_ = (_a12_, _lat2_, _lon2_, _azi2_, _s12_, _m12_, _M12_, _M21_, _S12_)
    _Units_ = (_Azi,  _Lat,   _Lon,   _Azi,   _M,    _Pass, _Pass, _Pass, _M2)


class GDict(_Dict):  # XXX _NamedDict
    '''Basic C{dict} with both key I{and} attribute access
       to the C{dict} items.

       Results of all C{geodesic} methods are returned as a
       L{GDict} instance.
    '''
    _iteration = None  # Iteration number (C{int}) or C{None}

    @property_RO
    def iteration(self):  # see .named._NamedBase
        '''Get the iteration number (C{int}) or
           C{None} if not available/applicable.
        '''
        return self._iteration

    def toDirect9Tuple(self, dflt=NAN):
        '''Convert this L{GDict} result to a 9-tuple, like I{Karney}'s
           method C{geographiclib.geodesic.Geodesic._GenDirect}.

           @kwarg dflt: Default value for missing items (C{any}).

           @return: L{Direct9Tuple}C{(a12, lat2, lon2, azi2,
                    s12, m12, M12, M21, S12)}
        '''
        return self._toTuple(Direct9Tuple, dflt)

    def toGeodSolve12Tuple(self, dflt=NAN):  # PYCHOK 12 args
        '''Convert this L{GDict} result to a 12-Tuple, compatible with I{Karney}'s
           U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>}
           result.

           @kwarg dflt: Default value for missing items (C{any}).

           @return: L{GeodSolve12Tuple}C{(lat1, lon1, azi1, lat2, lon2, azi2,
                    s12, a12, m12, M12, M21, S12)}.
        '''
        return self._toTuple(_MODS.geodsolve.GeodSolve12Tuple, dflt)

    def toInverse10Tuple(self, dflt=NAN):
        '''Convert this L{GDict} result to a 10-tuple, like I{Karney}'s
           method C{geographiclib.geodesic.Geodesic._GenInverse}.

           @kwarg dflt: Default value for missing items (C{any}).

           @return: L{Inverse10Tuple}C{(a12, s12, salp1, calp1,
                    salp2, calp2, m12, M12, M21, S12)}.
        '''
        return self._toTuple(Inverse10Tuple, dflt)

    @deprecated_method
    def toRhumb7Tuple(self, dflt=NAN):  # PYCHOK no cover
        '''DEPRECATED, used method C{toRhumb8Tuple}.

           @return: A I{DEPRECATED} L{Rhumb7Tuple}.
        '''
        return self._toTuple(_MODS.deprecated.Rhumb7Tuple, dflt)

    def toRhumb8Tuple(self, dflt=NAN):
        '''Convert this L{GDict} result to a 8-tuple.

           @kwarg dflt: Default value for missing items (C{any}).

           @return: L{Rhumb8Tuple}C{(lat1, lon1, lat2, lon2,
                    azi12, s12, S12, a12)}.
        '''
        return self._toTuple(Rhumb8Tuple, dflt)

    def toRhumbSolve7Tuple(self, dflt=NAN):
        '''Convert this L{GDict} result to a 8-tuple.

           @kwarg dflt: Default value for missing items (C{any}).

           @return: L{RhumbSolve7Tuple}C{(lat1, lon1, lat2, lon2,
                    azi12, s12, S12)}.
        '''
        return self._toTuple(_MODS.rhumbsolve.RhumbSolve7Tuple, dflt)

    def _toTuple(self, nTuple, dflt):
        '''(INTERNAL) Convert this C{GDict} to an B{C{nTuple}}.
        '''
        _g = getattr
        t  = tuple(_g(self, n, dflt) for n in nTuple._Names_)
        return nTuple(t, iteration=self._iteration)


class Inverse10Tuple(_GTuple):
    '''10-Tuple C{(a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12)} with arc length
       C{a12} in C{degrees}, distance C{s12} and reduced length C{m12} in C{meter}, area
       C{S12} in C{meter} I{squared} and the sines C{salp1}, C{salp2} and cosines C{calp1},
       C{calp2} of the initial C{1} and final C{2} foward azimuths.
    '''
    _Names_ = (_a12_, _s12_, 'salp1', 'calp1', 'salp2', 'calp2', _m12_, _M12_, _M21_, _S12_)
    _Units_ = (_Azi,  _M,    _Pass,   _Pass,   _Pass,   _Pass,   _Pass, _Pass, _Pass, _M2)

    def toGDict(self, **updates):
        '''Convert this C{Inverse10Tuple} to a L{GDict}.

           @kwarg updates: Optional items to apply (C{nam=value} pairs)
        '''
        return _GTuple.toGDict(self, azi1=atan2d(self.salp1, self.calp1),  # PYCHOK indent, namedTuple
                                     azi2=atan2d(self.salp2, self.calp2),  # PYCHOK namedTuple
                                   **updates)  # PYCHOK indent


class _kWrapped(object):  # in .geodesicw
    ''''(INTERNAL) Wrapper for some of I{Karney}'s U{geographiclib
        <https://PyPI.org/project/geographiclib>} classes.
    '''

    @Property_RO
    def geographiclib(self):
        '''Get the imported C{geographiclib}, provided the U{geographiclib
           <https://PyPI.org/project/geographiclib>} package is installed,
           otherwise an C{ImportError}.
        '''
        g = _xgeographiclib(self.__class__, 1, 49)
        from geographiclib.geodesic import Geodesic
        g.Geodesic = Geodesic
        from geographiclib.geodesicline import GeodesicLine
        g.GeodesicLine = GeodesicLine
        from geographiclib.geomath import Math
        g.Math = Math
        return g

    @Property_RO  # MCCABE 13
    def Math(self):
        '''Wrap the C{geomath.Math} class, provided the U{geographiclib
           <https://PyPI.org/project/geographiclib>} package is installed,
           otherwise C{None}.
        '''
        try:
            g = self.geographiclib
            M = g.Math
            if _xversion_info(g) < (2,):
                if _K_2_0:
                    M = None
#           elif not _K_2_0:  # XXX set 2.0?
#               _K_2_0 = True
        except (AttributeError, ImportError):
            M = None
        return M

_wrapped = _kWrapped()  # PYCHOK singleton, .datum, .test/base.py


class Rhumb8Tuple(_GTuple):
    '''8-Tuple C{(lat1, lon1, lat2, lon2, azi12, s12, S12, a12)} with lat- C{lat1},
       C{lat2} and longitudes C{lon1}, C{lon2} of both points, the azimuth of the
       rhumb line C{azi12}, the distance C{s12}, the area C{S12} under the rhumb
       line and the angular distance C{a12} between both points.
    '''
    _Names_ = (_lat1_, _lon1_, _lat2_, _lon2_, _azi12_, _s12_, _S12_,  _a12_)
    _Units_ = ( Lat,    Lon,    Lat,    Lon,   _Azi,    _M,    _M2,    _Deg)

    def toDirect9Tuple(self, dflt=NAN, **a12_azi1_azi2_m12_M12_M21):
        '''Convert this L{Rhumb8Tuple} result to a 9-tuple, like I{Karney}'s
           method C{geographiclib.geodesic.Geodesic._GenDirect}.

           @kwarg dflt: Default value for missing items (C{any}).
           @kwarg a12_azi1_azi2_m12_M12_M21: Optional keyword arguments
                     to specify or override L{Inverse10Tuple} items.

           @return: L{Direct9Tuple}C{(a12, lat2, lon2, azi2, s12,
                    m12, M12, M21, S12)}
        '''
        d = dict(azi1=self.azi12, M12=_1_0, m12=self.s12,  # PYCHOK attr
                 azi2=self.azi12, M21=_1_0)  # PYCHOK attr
        if a12_azi1_azi2_m12_M12_M21:
            d.update(a12_azi1_azi2_m12_M12_M21)
        return self._toTuple(Direct9Tuple, dflt, d)

    def toInverse10Tuple(self, dflt=NAN, **a12_m12_M12_M21_salp1_calp1_salp2_calp2):
        '''Convert this L{Rhumb8Tuple} to a 10-tuple, like I{Karney}'s
           method C{geographiclib.geodesic.Geodesic._GenInverse}.

           @kwarg dflt: Default value for missing items (C{any}).
           @kwarg a12_m12_M12_M21_salp1_calp1_salp2_calp2: Optional keyword
                     arguments to specify or override L{Inverse10Tuple} items.

           @return: L{Inverse10Tuple}C{(a12, s12, salp1, calp1, salp2, calp2,
                    m12, M12, M21, S12)}.
        '''
        s, c = sincos2d(self.azi12)  # PYCHOK attr
        d = dict(salp1=s, calp1=c, M12=_1_0, m12=self.s12,  # PYCHOK attr
                 salp2=s, calp2=c, M21=_1_0)
        if a12_m12_M12_M21_salp1_calp1_salp2_calp2:
            d.update(a12_m12_M12_M21_salp1_calp1_salp2_calp2)
        return self._toTuple(Inverse10Tuple, dflt, d)

    def _toTuple(self, nTuple, dflt, updates={}):
        '''(INTERNAL) Convert this C{Rhumb8Tuple} to an B{C{nTuple}}.
        '''
        _g = self.toGDict(**updates).get
        t  = tuple(_g(n, dflt) for n in nTuple._Names_)
        return nTuple(t, name=self.name)

    @deprecated_method
    def _to7Tuple(self):
        '''DEPRECATED, do not use!'''
        return _MODS.deprecated.Rhumb7Tuple(self[:-1])


def _around(x):  # in .utily.sincos2d
    '''I{Coarsen} a scalar by rounding small values to underflow to C{0.0}.

       @return: Coarsened value (C{float}).

       @see: I{Karney}'s U{geomath.Math.AngRound<https://SourceForge.net/p/
       geographiclib/code/ci/release/tree/python/geographiclib/geomath.py>}
    '''
    try:
        return _wrapped.Math.AngRound(x)
    except AttributeError:
        if x:
            y = _1_16th - fabs(x)
            if y > 0:  # fabs(x) < _1_16th
                x = _copysign(_1_16th - y, x)
        else:
            x = _0_0  # -0 to 0
        return x


def _atan2d(y, x):
    '''Return C{atan2(B{y}, B{x})} in C{degrees}.
    '''
    try:
        return _wrapped.Math.atan2d(y, x)
    except AttributeError:
        return atan2d(y, x)


def _cbrt(x):
    '''Return C{cubic root(B{x})}.
    '''
    try:
        return _wrapped.Math.cbrt(x)
    except AttributeError:
        return cbrt(x)


def _copyBit(x, y):
    '''Like C{copysign0(B{x}, B{y})}, with C{B{x} > 0}.
    '''
    return (-x) if _signBit(y) else x


def _2cos2x(cx, sx):  # in .auxDST, .auxLat, .gxbases
    '''Return M{2 * cos(2 * x)} from cos(x) and sin(x).
    '''
    r = cx - sx
    if r:
        r *= (cx + sx) * _2_0
    return r


def _diff182(deg0, deg, K_2_0=False):
    '''Compute C{deg - deg0}, reduced to C{[-180,180]} accurately.

       @return: 2-Tuple C{(delta_angle, residual)} in C{degrees}.
    '''
    try:
        return _wrapped.Math.AngDiff(deg0, deg)
    except AttributeError:
        if K_2_0 or _K_2_0:  # geographiclib 2.0
            _r, _360 = fremainder, _360_0
            d, t = _sum2(_r(-deg0, _360),
                         _r( deg,  _360))
            d, t = _sum2(_r( d,    _360), t)
            if d in (_0_0, _180_0, _N_180_0):
                d = _copysign(d, -t if t else (deg - deg0))
        else:
            _n = _norm180
            d, t = _sum2(_n(-deg0), _n(deg))
            d = _n(d)
            if t > 0 and d == _180_0:
                d = _N_180_0
            d, t = _sum2(d, t)
        return d, t


def _fix90(deg):  # mimick Math.LatFix
    '''Replace angle in C{degrees} outside [-90,90] by NAN.

       @return: Angle C{degrees} or NAN.
    '''
    try:
        return _wrapped.Math.LatFix(deg)
    except AttributeError:
        return NAN if fabs(deg) > 90 else deg


def _isfinite(x):  # mimick geomath.Math.isfinite
    '''Check finiteness of C{x}.

       @return: C{True} if finite.
    '''
    try:
        return _wrapped.Math.isfinite(x)
    except AttributeError:
        return _math_isfinite(x)  # and fabs(x) <= _MAX


def _norm2(x, y):  # mimick geomath.Math.norm
    '''Normalize C{B{x}} and C{B{y}}.

       @return: 2-Tuple of C{(B{x}, B{y})}, normalized.
    '''
    try:
        return _wrapped.Math.norm(x, y)
    except AttributeError:
        return norm2(x, y)


def _norm180(deg):  # mimick geomath.Math.AngNormalize
    '''Reduce angle in C{degrees} to (-180,180].

       @return: Reduced angle C{degrees}.
    '''
    try:
        return _wrapped.Math.AngNormalize(deg)
    except AttributeError:
        d = fremainder(deg, _360_0)
        if d in (_180_0, _N_180_0):
            d = _copysign(_180_0, deg) if _K_2_0 else _180_0
        return d


def _polygon(geodesic, points, closed, line, wrap):
    '''(INTERNAL) Compute the area or perimeter of a polygon,
        using a L{GeodesicExact}, L{GeodesicSolve} or (if the
        C{geographiclib} package is installed) a C{Geodesic}
        or C{geodesicw.Geodesic} instance.
    '''
    if not wrap:  # capability LONG_UNROLL can't be off
        notImplemented(None, wrap=wrap, up=3)

    if _MODS.booleans.isBoolean(points):
        # recursive call for each boolean clip

        def _a_p(clip, *args, **unused):
            return _polygon(geodesic, clip, *args)

        if not closed:  # closed only
            raise _ValueError(closed=closed, points=_composite_)

        return points._sum1(_a_p, closed, line, wrap)

    gP = geodesic.Polygon(line)
    _A = gP.AddPoint

    Ps = _MODS.iters.PointsIter(points, loop=1, wrap=wrap)  # base=LatLonEllipsoidalBase(0, 0)
    p1 =  p0 = Ps[0]

    # note, lon deltas are unrolled, by default
    _A(p1.lat, p1.lon)
    for p2 in Ps.iterate(closed=closed):
        if wrap and not Ps.looped:
            p2 = _unrollon(p1, p2)
        _A(p2.lat, p2.lon)
        p1 = p2
    if closed and line and p1 != p0:
        _A(p0.lat, p0.lon)

    # gP.Compute returns (number_of_points, perimeter, signed area)
    return gP.Compute(False, True)[1 if line else 2]


def _polynomial(x, cs, i, j):  # PYCHOK shared
    '''(INTERNAL) Like C++ C{GeographicLib.Math.hpp.polyval} but with a
       different signature and cascaded summation as C{karney._sum2_}.

       @return: M{sum(cs[k] * x**(j - k - 1) for k in range(i, j)}
    '''
    # assert 0 <= i <= j <= len(cs)
#   try:
#       return _wrapped.Math.polyval(j - i - 1, cs, i, x)
#   except AttributeError:
#       pass
    s, t = cs[i], _0_0
    for c in cs[i+1:j]:
        s, t = _sum2_(s * x, t * x, c)
    return s  # + t


def _remainder(x, y):
    '''Remainder of C{x / y}.

       @return: Remainder in the range M{[-y / 2, y / 2]}, preserving signed 0.0.
    '''
    try:
        return _wrapped.Math.remainder(x, y)
    except AttributeError:
        return fremainder(x, y)


if _K_2_0:
    from math import cos as _cos, sin as _sin

    def _sincos2(rad):
        return _sin(rad), _cos(rad)

    _signBit = _MODS.basics.signBit
else:
    _sincos2 = _MODS.utily.sincos2  # PYCHOK shared

    def _signBit(x):
        '''(INTERNAL) GeographicLib 1.52-.
        '''
        return x < 0


def _sincos2d(deg):
    '''Return sine and cosine of an angle in C{degrees}.

       @return: 2-Tuple C{(sin(B{deg}), cos(B{deg}))}.
    '''
    try:
        return _wrapped.Math.sincosd(deg)
    except AttributeError:
        return sincos2d(deg)


def _sincos2de(deg, t):
    '''Return sine and cosine of a corrected angle in C{degrees}.

       @return: 2-Tuple C{(sin(B{deg}), cos(B{deg}))}.
    '''
    try:
        return _wrapped.Math.sincosde(deg, t)
    except AttributeError:
        return sincos2d(deg, adeg=t)


def _sum2(u, v):  # mimick geomath.Math.sum, actually sum2
    '''Error-free summation like C{geomath.Math.sum}.

       @return: 2-Tuple C{(B{u} + B{v}, residual)}.

       @note: The C{residual} can be the same as B{C{u}} or B{C{v}}.

       @see: U{Algorithm 3.1<https://www.TUHH.De/ti3/paper/rump/OgRuOi05.pdf>}.
    '''
    try:
        return _wrapped.Math.sum(u, v)
    except AttributeError:
        s = u + v
        r = s - v
        t = s - r
        # if Algorithm_3_1:
        #   t = (u - t) + (v + r)
        # elif C_CPP:  # Math::sum C/C++
        #   r -= u
        #   t -= v
        #   t += r
        #   t = -t
        # else:
        t = (u - r) + (v - t)
        return s, t


def _sum2_(s, t, *vs):
    '''Accumulate any B{C{vs}} into a previous C{_sum2(s, t)}.

       @return: 2-Tuple C{(B{s} + B{t} + sum(B{vs}), residual)}.

       @see: I{Karney's} C++ U{Accumulator<https://GeographicLib.SourceForge.io/
             html/Accumulator_8hpp_source.html>} comments for more details and
             function C{_sum2} above.

       @note: NOT "error-free", see C{pygeodesy.test/testKarney.py}.
    '''
    _s2, _u0 = _sum2, unsigned0
    for v in vs:
        if v:
            t, u = _s2(t, v)  # start at the least-
            if s:
                s, t = _s2(s, t)  # significant end
                if s:
                    t += u  # accumulate u into t
#               elif t:  # s == 0 implies t == 0
#                   raise _AssertionError(t=t, txt=_not_(_0_))
                else:
                    s = _u0(u)  # result is u, t = 0
            else:
                s, t = _u0(t), u
    return s, t


def _tand(x):
    '''Return C{tan(B{x})} in C{degrees}.
    '''
    try:
        return _wrapped.Math.tand(x)
    except AttributeError:
        return tand(x)


def _unroll2(lon1, lon2, wrap=False):  # see .ellipsoidalBaseDI._intersects2
    '''Unroll B{C{lon2 - lon1}} like C{geodesic.Geodesic.Inverse}.

       @return: 2-Tuple C{(B{lon2} - B{lon1}, B{lon2})} with B{C{lon2}}
                unrolled if B{C{wrap}} is C{True}, normalized otherwise.
    '''
    if wrap:
        d, t = _diff182(lon1, lon2)
        lon2, _ = _sum2_(d, t, lon1)  # (lon1 + d) + t
    else:
        lon2 = _norm180(lon2)
    return (lon2 - lon1), lon2


def _unsigned2(x):
    '''(INTERNAL) Unsign B{C{x}}.
    '''
    return (neg(x), True) if _signBit(x) else (x, False)


__all__ += _ALL_DOCS(_CapsBase)

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
