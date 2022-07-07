
# -*- coding: utf-8 -*-

u'''I{Charles F.F. Karney}'s U{geographiclib<https://PyPI.org/project/geographiclib>} C{geodesic}, wrapped.

Wrapper around Python classes C{Geodesic} and C{GeodesicLine} and several C{Math} functions from
I{Karney}'s Python package U{geographiclib<https://PyPI.org/project/geographiclib>}, provided
that package is installed.

The I{wrapped} class methods return a L{GDict} instance offering access to the C{dict} items either
by C{key} or by C{attribute} name.

With env variable C{PYGEODESY_GEOGRAPHICLIB} left undefined or set to C{"2"}, this module and
L{pygeodesy.geodesicx} will use U{GeographicLib 2.0<https://GeographicLib.SourceForge.io/C++/doc/>}
transcoding, otherwise C{1.52} or older.

Karney-based functionality
==========================

1. The following classes and functions in C{pygeodesy}

  - L{AlbersEqualArea}, L{AlbersEqualArea2}, L{AlbersEqualArea4},
    L{AlbersEqualAreaCylindrical}, L{AlbersEqualAreaNorth}, L{AlbersEqualAreaSouth} --
    U{AlbersEqualArea<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1AlbersEqualArea.html>}

  - L{CassiniSoldner} -- U{CassiniSoldner<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1CassiniSoldner.html>}

  - L{EcefKarney} -- U{Geocentric<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1Geocentric.html>}

  - L{Elliptic} -- U{EllipticFunction<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1EllipticFunction.html>}

  - L{EquidistantExact}, L{EquidistantGeodSolve}, L{EquidistantKarney} -- U{AzimuthalEquidistant
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AzimuthalEquidistant.html>}

  - L{Etm}, L{ExactTransverseMercator} -- U{TransverseMercatorExact
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1TransverseMercatorExact.html>}

  - L{GeodesicAreaExact}, L{PolygonArea} -- U{PolygonArea<https://GeographicLib.SourceForge.io/
    html/classGeographicLib_1_1PolygonAreaT.html>}

  - L{GeodesicExact}, L{GeodesicLineExact} -- U{GeodesicExact<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1GeodesicExact.html>}, U{GeodesicLineExact<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1GeodesicLineExact.html>}

  - L{GeoidKarney} -- U{Geoid<https://GeographicLib.SourceForge.io/C++/doc/geoid.html>}

  - L{Georef} -- U{Georef<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1Georef.html>}

  - L{GnomonicExact}, L{GnomonicGeodSolve}, L{GnomonicKarney} -- U{Gnomonic
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Gnomonic.html>}

  - L{KTransverseMercator} - U{TransverseMercator
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1TransverseMercator.html>}

  - L{LocalCartesian}, L{Ltp} -- U{LocalCartesian<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1LocalCartesian.html>}

  - L{Rhumb}, L{RhumbLine} -- U{Rhumb<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Rhumb.html>},
    U{RhumbLine<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1RhumbLine.html>},
    U{TransverseMercator<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1TransverseMercator.html>}

  - L{Ups} -- U{PolarStereographic<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1PolarStereographic.html>}

  - L{Utm} -- U{TransverseMercator<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1TransverseMercator.html>}

  - L{UtmUps}, L{Epsg} -- U{UTMUPS<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1UTMUPS.html>}

  - L{pygeodesy.atand}, L{pygeodesy.atan2d}, L{pygeodesy.sincos2}, L{pygeodesy.sincos2d}, L{pygeodesy.tand} -- U{
    Math<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Math.html>}

are I{transcoded} from C++ classes in I{Karney}'s U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/annotated.html>}.

2. These C{pygeodesy} modules and classes

  - L{ellipsoidalGeodSolve}, L{ellipsoidalKarney}, L{geodsolve}, L{karney}, L{rhumbsolve}
  - L{EquidistantKarney}, L{FrechetKarney}, L{GeodesicSolve}, L{GeodesicLineSolve}, L{GnomonicGeodSolve},
    L{GnomonicKarney}, L{HeightIDWkarney}

are or use I{wrappers} around I{Karney}'s Python U{geographiclib<https://PyPI.org/project/geographiclib>}
C{geodesic}, C++ utility U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>} or
C++ utility U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>}.

3. All C{pygeodesy} functions and methods to compute I{ellipsoidal} intersections and trilaterations

  - L{ellipsoidalExact.intersection3}, L{ellipsoidalExact.intersections2}, L{ellipsoidalExact.nearestOn},
    L{ellipsoidalExact.LatLon.intersection3}, L{ellipsoidalExact.LatLon.intersections2},
    L{ellipsoidalExact.LatLon.nearestOn}, L{ellipsoidalExact.LatLon.trilaterate5}

  - L{ellipsoidalKarney.intersection3}, L{ellipsoidalKarney.intersections2}, L{ellipsoidalKarney.nearestOn},
    L{ellipsoidalKarney.LatLon.intersection3}, L{ellipsoidalKarney.LatLon.intersections2},
    L{ellipsoidalKarney.LatLon.nearestOn}, L{ellipsoidalKarney.LatLon.trilaterate5}

  - L{ellipsoidalVincenty.intersection3}, L{ellipsoidalVincenty.intersections2}, L{ellipsoidalVincenty.nearestOn},
    L{ellipsoidalVincenty.LatLon.intersection3}, L{ellipsoidalVincenty.LatLon.intersections2},
    L{ellipsoidalVincenty.LatLon.nearestOn}, L{ellipsoidalVincenty.LatLon.trilaterate5}

are implementations of I{Karney}'s solution posted under U{The B{ellipsoidal} case
<https://GIS.StackExchange.com/questions/48937/calculating-intersection-of-two-circles>} and in paper U{Geodesics
on an ellipsoid of revolution<https://ArXiv.org/pdf/1102.1215.pdf>} (pp 20-21, section B{14. MARITIME BOUNDARIES}).

4. Spherical functions

  - L{pygeodesy.excessKarney_}, L{sphericalTrigonometry.areaOf}

in C{pygeodesy} are based on I{Karney}'s post U{Area of a spherical polygon
<https://MathOverflow.net/questions/97711/the-area-of-spherical-polygons>}, 3rd Answer.
'''

from pygeodesy.basics import _copysign, _isfinite as _math_isfinite, neg, unsigned0, \
                             _xgeographiclib, _xImportError, _xversion_info, \
                             _xinstanceof, _zip, isodd  # PYCHOK shared
from pygeodesy.datums import Ellipsoid2, _ellipsoidal_datum, _WGS84
# from pygeodesy.ellipsoids import Ellipsoid2  # from .datums
from pygeodesy.errors import _AssertionError, _ValueError, _xkwds, _xkwds_get, \
                             _or  # PYCHOK shared
from pygeodesy.fmath import cbrt, fremainder, norm2, hypot as _hypot, unstr  # PYCHOK shared
from pygeodesy.interns import NAN, NN, _2_, _a12_, _area_, _azi2_, _DOT_, _lat2_, \
                             _lon2_, _m12_, _M12_, _M21_, _number_, _s12_, _S12_, \
                             _0_0, _180_0, _N_180_0, _360_0, _1_16th
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS, _getenv
from pygeodesy.named import callername, classname, _Dict, _NamedBase, \
                           _NamedTuple, _Pass
from pygeodesy.namedTuples import Destination3Tuple, Distance3Tuple
from pygeodesy.props import deprecated_method, Property, Property_RO, property_RO
# from pygeodesy.streps import unstr  # from .fmath
from pygeodesy.units import Bearing as _Azi, Degrees as _Deg, Lat, Lon, \
                            Meter as _M, Meter2 as _M2, Number_, \
                            Precision_, _1mm as _TOL_M  # PYCHOK shared
from pygeodesy.utily import atan2d, sincos2d, tand, unroll180, wrap360

__all__ = _ALL_LAZY.karney
__version__ = '22.07.03'

_EWGS84     = _WGS84.ellipsoid  # PYCHOK used!
_K_2_0      = _getenv('PYGEODESY_GEOGRAPHICLIB', _2_) == _2_
_perimeter_ = 'perimeter'


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


def _ellipsoid(a_ellipsoid, f=None, name=NN, raiser=True):  # in .geodesicx.gx, .geodsolve, .rhumbx
    '''(INTERNAL) Get an ellipsoid from C{(B{a_..}, B{f})} or C{B{.._ellipsoid}}.
    '''
    return Ellipsoid2(a_ellipsoid, f, name=name) if f is not None else \
          _ellipsoidal_datum(a_ellipsoid, name=name, raiser=raiser).ellipsoid


def _raiseX(inst, x, *args):  # PYCHOK no cover
    '''(INTERNAL) Throw a C{GeodesicError} for C{geographiclib} issue B{C{x}} .
    '''
    n = _DOT_(classname(inst), callername(up=2, underOK=True))
    raise GeodesicError(unstr(n, *args), txt=str(x))


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
    LATITUDE       =  1 << 7   # compute latitude C{lat2}
    LONGITUDE      =  1 << 8   # compute longitude C{lon2}
    AZIMUTH        =  1 << 9   # azimuths C{azi1} and C{azi2}
    DISTANCE       =  1 << 10  # compute distance C{s12}
    DISTANCE_IN    =  1 << 11  # allow distance C{s12} in Direct
    REDUCEDLENGTH  =  1 << 12  # compute reduced length C{m12}
    GEODESICSCALE  =  1 << 13  # compute geodesic scales C{M12} and C{M21}
    AREA           =  1 << 14  # compute area C{S12}

    STANDARD       =  AZIMUTH | DISTANCE | DISTANCE_IN | LATITUDE | LONGITUDE
    ALL            =  0x7F80   # without LONG_UNROLL, LINE_OFF, REVERSE2 and _DEBUG_*

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

C{LINE_OFF} - Line without updates from parent geodesic or rhumb.

C{LATITUDE} - compute latitude C{lat2},

C{LONGITUDE} - compute longitude C{lon2},

C{LONG_UNROLL} - unroll C{lon2} in C{.Direct},

C{REDUCEDLENGTH} - compute reduced length C{m12},

C{REVERSE2} - reverse C{azi2},

and C{ALL} - all of the above.

C{STANDARD} = C{AZIMUTH | DISTANCE | DISTANCE_IN | LATITUDE | LONGITUDE}'''


class _CapsBase(_NamedBase):  # in .geodesicx.gxbases, .rhumbx
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
    def toRhumb7Tuple(self, dflt=NAN):
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
        return self._toTuple(_MODS.rhumbx.Rhumb8Tuple, dflt)

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


class GeodesicError(_ValueError):
    '''Error raised for L{pygeodesy.geodesicx} lack of convergence
       or other L{pygeodesy.geodesicx} or L{pygeodesy.karney} issues.
    '''
    pass


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


class _Wrapped(object):
    ''''(INTERNAL) Wrapper for some of I{Karney}'s U{geographiclib
        <https://PyPI.org/project/geographiclib>} classes.
    '''

    @Property_RO  # MCCABE 24
    def Geodesic(self):
        '''Get the I{wrapped} C{Geodesic} class, provided the U{geographiclib
           <https://PyPI.org/project/geographiclib>} package is installed,
           otherwise an C{ImportError}.
        '''
        _Geodesic =  self.geographiclib.Geodesic
        _DIRECT3  = _Geodesic.AZIMUTH | _Geodesic.LATITUDE | _Geodesic.LONGITUDE
        _INVERSE3 = _Geodesic.AZIMUTH | _Geodesic.DISTANCE

        class Geodesic(_Geodesic):
            '''I{Karney}'s U{Geodesic<https://GeographicLib.SourceForge.io/C++/doc/
               python/code.html#geographiclib.geodesic.Geodesic>} wrapper.
            '''
            _debug   =  0  # like .geodesicx.bases._GeodesicBase
            _E       = _EWGS84
            LINE_OFF =  0  # in .azimuthal._GnomonicBase and .css.CassiniSoldner

            def __init__(self, a_ellipsoid=_EWGS84, f=None, name=NN):  # PYCHOK signature
                '''New C{Geodesic} instance.

                   @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}) or datum
                                     (L{Datum}) or the equatorial radius
                                     of the ellipsoid (C{meter}).
                   @arg f: The flattening of the ellipsoid (C{scalar}) if
                           B{C{a_ellipsoid}) is specified as C{meter}.
                   @kwarg name: Optional name (C{str}).
                '''
                if a_ellipsoid not in (Geodesic._E, None):  # spherical OK
                    self._E = _ellipsoid(a_ellipsoid, f, name=name, raiser=False)
                try:
                    _Geodesic.__init__(self, *self.ellipsoid.a_f)
                except (TypeError, ValueError) as x:
                    _raiseX(self, x, *self.ellipsoid.a_f)

            Area = _Geodesic.Polygon  # like GeodesicExact.Area

            @Property
            def debug(self):
                '''Get the C{debug} option (C{bool}).
                '''
                return bool(self._debug)

            @debug.setter  # PYCHOK setter!
            def debug(self, debug):
                '''Set the C{debug} option (C{bool}) to include
                   more details in L{GDict} results.
                '''
                self._debug = _MODS.geodesicx.Caps._DEBUG_ALL if debug else 0

            def Direct(self, lat1, lon1, azi1, s12, *outmask):
                '''Return the C{Direct} result.
                '''
                try:
                    d = _Geodesic.Direct(self, lat1, lon1, azi1, s12, *outmask)
                except (TypeError, ValueError) as x:
                    _raiseX(self, x, lat1, lon1, azi1, s12, *outmask)
                return GDict(d)

            def Direct3(self, lat1, lon1, azi1, s12):  # PYCHOK outmask
                '''Return the destination lat, lon and reverse azimuth
                   (final bearing) in C{degrees}.

                   @return: L{Destination3Tuple}C{(lat, lon, final)}.
                '''
                d = self.Direct(lat1, lon1, azi1, s12, _DIRECT3)
                return Destination3Tuple(d.lat2, d.lon2, d.azi2)

            @Property_RO
            def ellipsoid(self):
                '''Get this geodesic's ellipsoid (C{Ellipsoid[2]}).
                '''
                return self._E

            @Property_RO
            def f1(self):  # in .css.CassiniSoldner.reset
                '''Get the geodesic's ellipsoid I{1 - flattening} (C{float}).
                '''
                return getattr(self, '_f1', self.ellipsoid.f1)

            def _GDictDirect(self, lat, lon, azi, arcmode, s12_a12,
                                                  outmask=_Geodesic.STANDARD):
                '''(INTERNAL) Get C{._GenDirect} result as C{GDict}.
                '''
                try:
                    t = _Geodesic._GenDirect(self, lat, lon, azi, arcmode, s12_a12, outmask)
                except (TypeError, ValueError) as x:
                    _raiseX(self, x, lat, lon, azi, arcmode, s12_a12, outmask)
                return Direct9Tuple(t).toGDict()  # *t

            def _GDictInverse(self, lat1, lon1, lat2, lon2, outmask=_Geodesic.STANDARD):
                '''(INTERNAL) Get C{._GenInverse} result as C{GDict}.
                '''
                try:
                    t = _Geodesic._GenInverse(self, lat1, lon1, lat2, lon2, outmask)
                except (TypeError, ValueError) as x:
                    _raiseX(self, x, lat1, lon1, lat2, lon2, outmask)
                return Inverse10Tuple(t).toGDict(lon1=lon1, lon2=lon2)  # *t

            def Inverse(self, lat1, lon1, lat2, lon2, *outmask):
                '''Return the C{Inverse} result.
                '''
                try:
                    d = _Geodesic.Inverse(self, lat1, lon1, lat2, lon2, *outmask)
                except (TypeError, ValueError) as x:
                    _raiseX(self, x, lat1, lon1, lat2, lon2, *outmask)
                return GDict(d)

            def Inverse1(self, lat1, lon1, lat2, lon2, wrap=False):
                '''Return the non-negative, I{angular} distance in C{degrees}.
                '''
                # see .FrechetKarney.distance, .HausdorffKarney._distance
                # and .HeightIDWkarney._distances
                _, lon2 = unroll180(lon1, lon2, wrap=wrap)  # self.LONG_UNROLL
                r = self.Inverse(lat1, lon1, lat2, lon2)
                # XXX self.DISTANCE needed for 'a12'?
                return abs(r.a12)

            def Inverse3(self, lat1, lon1, lat2, lon2):  # PYCHOK outmask
                '''Return the distance in C{meter} and the forward and
                   reverse azimuths (initial and final bearing) in C{degrees}.

                   @return: L{Distance3Tuple}C{(distance, initial, final)}.
                '''
                r = self.Inverse(lat1, lon1, lat2, lon2, _INVERSE3)
                return Distance3Tuple(r.s12, wrap360(r.azi1), wrap360(r.azi2))

            def Line(self, lat1, lon1, azi1, *caps):
                '''Set up a L{GeodesicLine} to compute several points on a
                   single geodesic.
                '''
                return _wrapped.GeodesicLine(self, lat1, lon1, azi1, *caps)

        # Geodesic.Direct.__doc__  = _Geodesic.Direct.__doc__
        # Geodesic.Inverse.__doc__ = _Geodesic.Inverse.__doc__
        # Geodesic.Line.__doc__    = _Geodesic.Line.__doc__
        return Geodesic

    @Property_RO  # MCCABE 16
    def GeodesicLine(self):
        '''Get the I{wrapped} C{GeodesicLine} class, provided the U{geographiclib
           <https://PyPI.org/project/geographiclib>} package is installed,
           otherwise an C{ImportError}.
        '''
        _GeodesicLine = self.geographiclib.GeodesicLine

        class GeodesicLine(_GeodesicLine):
            '''I{Karney}'s U{GeodesicLine <https://GeographicLib.SourceForge.io/C++/doc/
               python/code.html#geographiclib.geodesicline.GeodesicLine>} wrapper.
            '''
            def __init__(self, lat1, lon1, azi1, *caps):
                try:
                    _GeodesicLine.__init__(self, lat1, lon1, azi1, *caps)
                except (TypeError, ValueError) as x:
                    _raiseX(self, x, lat1, lon1, azi1, *caps)

            @Property_RO
            def a1(self):
                '''Get the I{equatorial arc} (C{degrees}), the arc length between
                   the northward equatorial crossing and point C{(lat1, lon1)}.

                   @see: U{EquatorialArc<https://GeographicLib.SourceForge.io/
                         C++/doc/classGeographicLib_1_1GeodesicLine.html>}
                '''
                try:
                    return _atan2d(self._ssig1, self._csig1)
                except AttributeError:
                    return NAN  # see .geodesicx.gxline._GeodesicLineExact

            equatorarc = a1

            def ArcPosition(self, a12, *outmask):
                try:
                    d = _GeodesicLine.ArcPosition(self, a12, *outmask)
                except (TypeError, ValueError) as x:
                    _raiseX(self, x, a12, *outmask)
                return GDict(d)

            @Property_RO
            def azi0(self):  # see .css.CassiniSoldner.forward4
                '''Get the I{equatorial azimuth} (C{degrees}), the azimuth of the
                   geodesic line as it crosses the equator in a northward direction.

                   @see: U{EquatorialAzimuth<https://GeographicLib.SourceForge.io/
                         C++/doc/classGeographicLib_1_1GeodesicLine.html>}
                '''
                try:
                    return _atan2d(self._salp0, self._calp0)
                except AttributeError:
                    return NAN  # see .geodesicx.gxline._GeodesicLineExact

            equatorazimuth = azi0

            def Position(self, s12, *outmask):
                try:
                    d = _GeodesicLine.Position(self, s12, *outmask)
                except (TypeError, ValueError) as x:
                    _raiseX(self, x, s12, *outmask)
                return GDict(d)

        # GeodesicLine.ArcPosition.__doc__ = _GeodesicLine.ArcPosition.__doc__
        # GeodesicLine.Position.__doc__    = _GeodesicLine.Position.__doc__
        return GeodesicLine

    @Property_RO
    def Geodesic_WGS84(self):
        '''Get the I{wrapped} C{Geodesic.WGS84} I{instance} provided the
           U{geographiclib<https://PyPI.org/project/geographiclib>} package
           is installed, otherwise an C{ImportError}.
        '''
        return _EWGS84.geodesic

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
        '''Get the C{Math} class, provided the U{geographiclib
           <https://PyPI.org/project/geographiclib>} package is
           installed, otherwise C{None}.
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

_wrapped = _Wrapped()  # PYCHOK singleton, .datum, .test/base.py


def _around(x):  # in .utily.sincos2d
    '''I{Coarsen} a scalar by rounding small values to underflow to C{0.0}.

       @return: Coarsened value (C{float}).

       @see: I{Karney}'s U{Math.AngRound<https://SourceForge.net/p/
       geographiclib/code/ci/release/tree/python/geographiclib/geomath.py>}
    '''
    try:
        return _wrapped.Math.AngRound(x)
    except AttributeError:
        pass
    if x:
        y = _1_16th - abs(x)
        if y > 0:  # abs(x) < _1_16th
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


def _diff182(deg0, deg):
    '''Compute C{deg - deg0}, reduced to C{[-180,180]} accurately.

       @return: 2-Tuple C{(delta_angle, residual)} in C{degrees}.
    '''
    try:
        return _wrapped.Math.AngDiff(deg0, deg)
    except AttributeError:
        pass
    if _K_2_0:  # geographiclib 2.0
        d, t = _sum2(fremainder(-deg0, _360_0),
                     fremainder( deg,  _360_0))
        d, t = _sum2(fremainder( d,    _360_0), t)
        if d in (_0_0, _180_0, _N_180_0):
            d = _copysign(d, -t if t else (deg - deg0))
    else:
        d, t = _sum2(_norm180(-deg0), _norm180(deg))
        d = _norm180(d)
        if t > 0 and d == _180_0:
            d = _N_180_0
        d, t = _sum2(d, t)
    return d, t


# def _Equidistant(equidistant, exact=False, geodsolve=False):
#     # (INTERNAL) Get the C{EquidistantExact}, C{-GeodSolve} or
#     # C{-Karney} class if B{C{equidistant}} in not callable.
#     if equidistant is None or not callable(equidistant):
#         if exact:
#             equidistant = _MODS.azimuthal.EquidistantExact
#         elif geodsolve:
#             equidistant = _MODS.azimuthal.EquidistantGeodSolve
#         else:
#             equidistant = _MODS.azimuthal.EquidistantKarney
#     return equidistant


def _fix90(deg):  # mimick Math.LatFix
    '''Replace angle in C{degrees} outside [-90,90] by NAN.

       @return: Angle C{degrees} or NAN.
    '''
    try:
        return _wrapped.Math.LatFix(deg)
    except AttributeError:
        return NAN if abs(deg) > 90 else deg


def _isfinite(x):  # mimick Math.AngNormalize
    '''Check finiteness of C{x}.

       @return: C{True} if finite.
    '''
    try:
        return _wrapped.Math.isfinite(x)
    except AttributeError:
        return _math_isfinite(x)  # and abs(x) <= _MAX


def _norm2(x, y):  # mimick Math.norm
    '''Normalize C{B{x}} and C{B{y}}.

       @return: 2-Tuple of C{(B{x}, B{y})}, normalized.
    '''
    try:
        return _wrapped.Math.norm(x, y)
    except AttributeError:
        return norm2(x, y)


def _norm180(deg):  # mimick Math.AngNormalize
    '''Reduce angle in C{degrees} to (-180,180].

       @return: Reduced angle C{degrees}.
    '''
    try:
        return _wrapped.Math.AngNormalize(deg)
    except AttributeError:
        pass
    d = fremainder(deg, _360_0)
    if d in (_180_0, -_180_0):
        d = _copysign(_180_0, deg) if _K_2_0 else _180_0
    return d


def _polygon(geodesic, points, closed, line, wrap):
    '''(INTERNAL) Compute the area or perimeter of a polygon,
        using a L{GeodesicExact}, L{GeodesicSolve} or (if the
        C{geographiclib} package is installed) a C{Geodesic}
        or C{_wrapped.Geodesic} instance.
    '''
    if not wrap:  # capability LONG_UNROLL can't be off
        raise _ValueError(wrap=wrap)

    gP = geodesic.Polygon(line)
    pA = gP.AddPoint

    Ps = _MODS.iters.PointsIter(points, loop=1)  # base=LatLonEllipsoidalBase(0, 0)
    p0 =  Ps[0]

    # note, lon deltas are unrolled, by default
    pA(p0.lat, p0.lon)
    for p in Ps.iterate(closed=closed):
        pA(p.lat, p.lon)
    if closed and line and p != p0:
        pA(p0.lat, p0.lon)

    # gP.Compute returns (number_of_points, perimeter, signed area)
    return gP.Compute(False, True)[1 if line else 2]


def _polynomial(x, cs, i, j):  # PYCHOK shared
    '''(INTERNAL) Like C{GeographicLib.Math.hpp.polyval} but with a
       different signature and cascaded summation as C{karney._sum2_}.

       @return: M{sum(cs[k] * x**(j - k - 1) for k in range(i, j)}
    '''
    s, t = cs[i], _0_0
    _s2_ = _sum2_
    for c in cs[i+1:j]:
        s, t = _s2_(s * x, t * x, c)
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
    from pygeodesy.basics import signBit as _signBit
    from math import cos as _cos, sin as _sin

    def _sincos2(rad):
        return _sin(rad), _cos(rad)

else:
    from pygeodesy.utily import sincos2 as _sincos2  # PYCHOK shared

    def _signBit(x):
        '''(INTERNAL) GeographicLin 1.52-.
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


def _sum2(u, v):  # mimick Math::sum, actually sum2
    '''Error-free summation like C{Math::sum}.

       @return: 2-Tuple C{(B{u} + B{v}, residual)}.

       @note: The C{residual} can be the same as B{C{u}} or B{C{v}}.

       @see: U{Algorithm 3.1<https://www.TUHH.De/ti3/paper/rump/OgRuOi05.pdf>}.
    '''
    try:
        return _wrapped.Math.sum(u, v)
    except AttributeError:
        pass
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
