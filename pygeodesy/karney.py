
# -*- coding: utf-8 -*-

u'''Wrapper around several C{geomath.Math} functions from I{Karney}'s Python package U{geographiclib
<https://PyPI.org/project/geographiclib>}, provided that package is installed.

Methods of the I{wrapped} L{Geodesic<geodesicw.Geodesic>} and L{GeodesicLine<geodesicw.GeodesicLine>}
classes return a L{GDict} instance offering access to the C{dict} items either by C{key} or by
C{attribute} name.

With env variable C{PYGEODESY_GEOGRAPHICLIB} left undefined or set to C{"2"}, modules L{geodesicw},
L{geodesicx} and this module will use U{GeographicLib 2.0+<https://GeographicLib.SourceForge.io/C++/doc/>}
and newer transcoding, otherwise C{1.52} or older.  Set C{PYGEODESY_GEOGRAPHICLIB=2.4} to default to the
C{Jacobi amplitude} instead of C{Bulirsch}' function in methods L{ExactTransverseMercator.forward
<pygeodesy.ExactTransverseMercator.forward>} and L{reverse <pygeodesy.ExactTransverseMercator.reverse>}.

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

  - L{Intersector} -- U{Intersect
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Intersect.html>} from I{GeographicLib 2.3+}

  - L{JacobiConformal} -- U{JacobiConformal
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1experimental_1_1JacobiConformal.html>}

  - L{KTransverseMercator} - U{TransverseMercator
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1TransverseMercator.html>}

  - L{LocalCartesian}, L{Ltp} -- U{LocalCartesian<https://GeographicLib.SourceForge.io/C++/doc/
    classGeographicLib_1_1LocalCartesian.html>}

  - L{Osgr} -- U{OSGB<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1OSGB.html>}

  - L{rhumb.aux_}, L{RhumbAux}, L{RhumbLineAux} -- U{Rhumb
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Rhumb.html>} and U{RhumbLine
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1RhumbLine.html>} from I{GeographicLib 2.2+}

  - L{rhumb.ekx}, L{Rhumb}, L{RhumbLine} -- U{Rhumb
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Rhumb.html>},
    U{RhumbLine<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1RhumbLine.html>},
    U{TransverseMercator<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1TransverseMercator.html>}
    from I{GeographicLib 2.0}

  - L{Ups} -- U{PolarStereographic<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1PolarStereographic.html>}

  - L{Utm} -- U{TransverseMercator<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1TransverseMercator.html>}

  - L{UtmUps}, L{Epsg} -- U{UTMUPS<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1UTMUPS.html>}

  - L{atan1d}, L{atan2d}, L{sincos2}, L{sincos2d}, L{tand} -- U{geomath.Math
    <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Math.html>}

are I{transcoded} from C++ classes in I{Karney}'s U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/annotated.html>}.

2. These C{pygeodesy} modules and classes

  - L{ellipsoidalGeodSolve}, L{ellipsoidalKarney}, L{geodesici}, L{geodsolve}, L{karney}, L{rhumb.solve}
  - L{EquidistantKarney}, L{FrechetKarney}, L{GeodesicSolve}, L{GeodesicLineSolve}, L{GnomonicGeodSolve},
    L{GnomonicKarney}, L{HeightIDWkarney}, L{Intersectool}

are or use I{wrappers} around I{Karney}'s Python U{geographiclib<https://PyPI.org/project/geographiclib>} or
C++ utility U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>},
U{IntersectTool<https://GeographicLib.SourceForge.io/C++/doc/IntersectTool.1.html>} or
U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>}.

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

  - L{RhumbLineAux.Intersection} and L{RhumbLine.Intersection}

are implementations of I{Karney}'s iterative solution posted under U{The B{ellipsoidal} case
<https://GIS.StackExchange.com/questions/48937/calculating-intersection-of-two-circles>} and in paper U{Geodesics
on an ellipsoid of revolution<https://ArXiv.org/pdf/1102.1215.pdf>} (pp 20-21, section B{14. MARITIME BOUNDARIES}).

4. The C{pygeodesy} methods to compute I{ellipsoidal} intersections and nearest points

  - L{RhumbLineAux.Intersecant2}, L{RhumbLineAux.PlumbTo}, L{RhumbLine.Intersecant2} and L{RhumbLine.PlumbTo}

are I{transcoded} of I{Karney}'s iterative C++ function U{rhumb-intercept
<https://SourceForge.net/p/geographiclib/discussion/1026620/thread/2ddc295e/>}.

5. Spherical functions

  - L{pygeodesy.excessKarney_}, L{sphericalTrigonometry.areaOf}

in C{pygeodesy} are based on I{Karney}'s post U{Area of a spherical polygon
<https://MathOverflow.net/questions/97711/the-area-of-spherical-polygons>}, 3rd Answer.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import _copysign, isint, neg, unsigned0, \
                             _xgeographiclib, _zip
from pygeodesy.constants import NAN, _isfinite as _math_isfinite, _0_0, \
                               _1_16th, _1_0, _2_0, _180_0, _N_180_0, _360_0
from pygeodesy.errors import GeodesicError, _ValueError, _xkwds
from pygeodesy.fmath import cbrt, fremainder, norm2  # Fhorner, Fsum
from pygeodesy.internals import _getenv, _popen2, _PYGEODESY, _version_info
from pygeodesy.interns import NN, _a12_, _area_, _azi1_, _azi2_, _azi12_, \
                             _composite_, _lat1_, _lat2_, _lon1_, _lon2_, \
                             _m12_, _M12_, _M21_, _number_, _s12_, _S12_, \
                             _SPACE_, _UNDER_, _X_, _1_, _2_,  _BAR_  # PYCHOK used!
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import ADict, _NamedBase, _NamedTuple, notImplemented, _Pass
from pygeodesy.props import deprecated_method, Property_RO, property_RO, \
                                               property_ROnce
from pygeodesy.units import Azimuth as _Azi, Degrees as _Deg, Lat, Lon, \
                            Meter as _M, Meter2 as _M2, Number_
from pygeodesy.utily import atan2d, sincos2d, tand, _unrollon,  fabs

# from math import fabs  # from .utily

__all__ = _ALL_LAZY.karney
__version__ = '24.11.26'

_2_4_       = '2.4'
_K_2_0      = _getenv(_PYGEODESY(_xgeographiclib, 1), _2_)
_K_2_4      = _K_2_0 ==  _2_4_
_K_2_0      = _K_2_4 or (_K_2_0 == _2_)
_perimeter_ = 'perimeter'


class _GTuple(_NamedTuple):  # in .testNamedTuples
    '''(INTERNAL) Helper.
    '''
    def toGDict(self, **updates):  # NO name=NN
        '''Convert this C{*Tuple} to a L{GDict}.

           @kwarg updates: Optional items to apply (C{name=value} pairs)
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
       of points of the polygon or polyline, the C{perimeter} in
       C{meter} and the C{area} in C{meter} I{squared}.
    '''
    _Names_ = (_number_, _perimeter_, _area_)
    _Units_ = ( Number_, _M,          _M2)


class Caps(object):
    '''(INTERNAL) Overriden by C{Caps} below.
    '''
    EMPTY         =  0        # formerly aka NONE
    _CAP_1        =  1 << 0   # for goedesici/-w
    _CAP_1p       =  1 << 1   # for goedesici/-w
    _CAP_2        =  1 << 2
    _CAP_3        =  1 << 3   # for goedesici/-w
    _CAP_4        =  1 << 4   # for goedesicw
    _CAP_ALL      =  0x1F
#   _CAP_MASK     = _CAP_ALL
    LATITUDE      =  1 << 7   # compute latitude C{lat2}
    LONGITUDE     =  1 << 8  | _CAP_3  # compute longitude C{lon2}
    AZIMUTH       =  1 << 9   # azimuths C{azi1} and C{azi2}
    DISTANCE      =  1 << 10 | _CAP_1  # compute distance C{s12}
    DISTANCE_IN   =  1 << 11 | _CAP_1 | _CAP_1p  # allow distance C{s12} in Direct
    REDUCEDLENGTH =  1 << 12 | _CAP_1 | _CAP_2  # compute reduced length C{m12}
    GEODESICSCALE =  1 << 13 | _CAP_1 | _CAP_2  # compute geodesic scales C{M12} and C{M21}
    AREA          =  1 << 14 | _CAP_4  # compute area C{S12}

    STANDARD      =  AZIMUTH | DISTANCE | LATITUDE | LONGITUDE

    _DIRECT3      =  AZIMUTH  | LATITUDE | LONGITUDE  # for goedesicw only
    _INVERSE3     =  AZIMUTH  | DISTANCE  # for goedesicw only
    _STD          =  STANDARD  # for goedesicw only
    _STD_LINE     = _STD      | DISTANCE_IN  # for goedesici/-w

    LINE_CAPS     = _STD_LINE | REDUCEDLENGTH | GEODESICSCALE  # .geodesici only
    LONG_UNROLL   =  1 << 15  # unroll C{lon2} in .Direct and .Position
#                 =  1 << 16  # unused
    LINE_OFF      =  1 << 17  # Line without updates from parent geodesic or rhumb
    REVERSE2      =  1 << 18  # reverse C{azi2}
    ALL           =  0x7F80 | _CAP_ALL  # without LONG_UNROLL, LINE_OFF, REVERSE2 and _DEBUG_*

    LATITUDE_LONGITUDE         = LATITUDE | LONGITUDE
    LATITUDE_LONGITUDE_AREA    = LATITUDE | LONGITUDE | AREA

    AZIMUTH_DISTANCE           = AZIMUTH | DISTANCE
    AZIMUTH_DISTANCE_AREA      = AZIMUTH | DISTANCE | AREA

    _SALP_CALPs_   =  1 << 19  # (INTERNAL) GeodesicExact._GenInverse

    _DEBUG_AREA    =  1 << 20  # (INTERNAL) include Line details
    _DEBUG_DIRECT  =  1 << 21  # (INTERNAL) include Direct details
    _DEBUG_INVERSE =  1 << 22  # (INTERNAL) include Inverse details
    _DEBUG_LINE    =  1 << 23  # (INTERNAL) include Line details
    _DEBUG_ALL     = _DEBUG_AREA | _DEBUG_DIRECT | _DEBUG_INVERSE | \
                     _DEBUG_LINE | _SALP_CALPs_

    _OUT_ALL       =  ALL  # see geographiclib.geodesiccapabilities.py
    _OUT_MASK      =  ALL | LONG_UNROLL | REVERSE2 | _DEBUG_ALL

    _AZIMUTH_LATITUDE_LONGITUDE           =  AZIMUTH | LATITUDE | LONGITUDE
    _AZIMUTH_LATITUDE_LONG_UNROLL         =  AZIMUTH | LATITUDE | LONG_UNROLL
    _DEBUG_DIRECT_LINE                    = _DEBUG_DIRECT | _DEBUG_LINE
#   _DISTANCE_IN_OUT                      =  DISTANCE_IN & _OUT_MASK  # == DISTANCE_IN in .gx, .gxline
    _REDUCEDLENGTH_GEODESICSCALE          =  REDUCEDLENGTH | GEODESICSCALE
#   _REDUCEDLENGTH_GEODESICSCALE_DISTANCE =  REDUCEDLENGTH | GEODESICSCALE | DISTANCE

    def items(self):
        '''Yield all I{public} C{Caps} as 2-tuple C{(NAME, mask)}.
        '''
        for n, C in Caps.__class__.__dict__.items():
            if isint(C) and not n.startswith(_UNDER_) \
                        and n.replace(_UNDER_, NN).isupper():
                yield n, C

    def toStr(self, Cask, sep=_BAR_):
        '''Return C{Caps} or an C{outmask} as C{str} or tuple of C{str}s.
        '''
        s = (n for n, C in sorted(Caps.items())
                        if C and (Cask & C) == C)  # and int1s(C) <= 3
        return sep.join(s) if sep else tuple(s)

Caps = Caps()  # PYCHOK singleton
'''I{Enum}-style masks to be bit-C{or}'ed to specify geodesic or
rhumb capabilities (C{caps}) and results (C{outmask}).

C{AREA} - compute area C{S12},

C{AZIMUTH} - include azimuths C{azi1} and C{azi2},

C{DISTANCE} - compute distance C{s12} and angular distance C{a12},

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

_key2Caps = dict(a12 =Caps.DISTANCE,  # in GDict._unCaps
                 azi2=Caps.AZIMUTH,
                 lat2=Caps.LATITUDE,
                 lon2=Caps.LONGITUDE,
                 m12 =Caps.REDUCEDLENGTH,
                 M12 =Caps.GEODESICSCALE,
                 M21 =Caps.GEODESICSCALE,
                 s12 =Caps.DISTANCE,
                 S12 =Caps.AREA)


class _CapsBase(_NamedBase):  # in .auxilats, .geodesicx.gxbases
    '''(INTERNAL) Base class for C{Geodesic*}, C{Geodesic*Exact}, C{Intersectool} and C{Rhumb*Base}.
    '''
    ALL           = Caps.ALL
    AREA          = Caps.AREA
    AZIMUTH       = Caps.AZIMUTH
    DISTANCE      = Caps.DISTANCE
    DISTANCE_IN   = Caps.DISTANCE_IN
    EMPTY         = Caps.EMPTY  # aka NONE
    GEODESICSCALE = Caps.GEODESICSCALE
    LATITUDE      = Caps.LATITUDE
    LINE_CAPS     = Caps.LINE_CAPS
    LINE_OFF      = Caps.LINE_OFF
    LONGITUDE     = Caps.LONGITUDE
    LONG_UNROLL   = Caps.LONG_UNROLL
    REDUCEDLENGTH = Caps.REDUCEDLENGTH
    STANDARD      = Caps.STANDARD
    _STD_LINE     = Caps._STD_LINE  # for geodesici

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

    def _iter2tion(self, r, iter=None, **unused):
        '''(INTERNAL) Copy C{C{s}.iter} into C{B{r}._iteration}.
        '''
        if iter is not None:
            self._iteration = r._iteration = iter
        return r


class Direct9Tuple(_GTuple):
    '''9-Tuple C{(a12, lat2, lon2, azi2, s12, m12, M12, M21, S12)} with arc
       length C{a12}, angles C{lat2}, C{lon2} and azimuth C{azi2} in C{degrees},
       distance C{s12} and reduced length C{m12} in C{meter} and area C{S12} in
       C{meter} I{squared}.
    '''
    _Names_ = (_a12_, _lat2_, _lon2_, _azi2_, _s12_, _m12_, _M12_, _M21_, _S12_)
    _Units_ = (_Azi,  _Lat,   _Lon,   _Azi,   _M,    _Pass, _Pass, _Pass, _M2)


class GDict(ADict):  # XXX _NamedDict
    '''A C{dict} with both C{key} I{and} C{attribute} access to the C{dict} items.

       Results of all C{geodesic} and C{rhumb} methods (with capitalized named) are
       returned as L{GDict} instances, see for example L{GeodesicExact} and L{RhumbAux}.
    '''
    def toDirect9Tuple(self, dflt=NAN):
        '''Convert this L{GDict} result to a 9-tuple, like I{Karney}'s method
           C{geographiclib.geodesic.Geodesic._GenDirect}.

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

    def _toNAN(self, outmask):  # .GeodesicLineExact._GenPosition
        '''(INTERNAL) Convert this C{GDict} to all C{NAN}s.
        '''
        d = dict((n, NAN) for n in GeodSolve12Tuple._Names_)
        return self.set_(**d)._unCaps(outmask)

    @deprecated_method
    def toRhumb7Tuple(self, dflt=NAN):  # PYCHOK no cover
        '''DEPRECATED on 23.12.07, use method C{toRhumb8Tuple}.

           @return: A I{DEPRECATED} L{Rhumb7Tuple}.
        '''
        return self._toTuple(_MODS.deprecated.classes.Rhumb7Tuple, dflt)

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
        return self._toTuple(_MODS.rhumb.solve.RhumbSolve7Tuple, dflt)

    def _toTuple(self, nTuple, dflt):
        '''(INTERNAL) Convert this C{GDict} to an B{C{nTuple}}.
        '''
        _g = getattr
        t  = tuple(_g(self, n, dflt) for n in nTuple._Names_)
        return nTuple(t, iteration=self._iteration)

    def _2X(self, gl, _2X=_X_):  # .Intersectool, .Intersector
        '''(INTERNAL) Rename C{-2} attr to C{-X} or C{-M}.
        '''
        X = GDict(self)
        for n in (_lat2_, _lon2_, _azi2_, _s12_, _a12_):
            if n in X:  # X._X = X._2
                X[n[:-1] + _2X] = X.pop(n)
            v = getattr(gl, n, X)
            if v is not X:  # X._2 = gl._2
                X[n] = v
        return X

    def _unCaps(self, outmask):  # in .geodsolve
        '''(INTERNAL) Remove superfluous items.
        '''
        for k, C in _key2Caps.items():
            if k in self and (outmask & C) != C:
                self.pop(k)  # delattr(self, k)
        return self


class GeodSolve12Tuple(_GTuple):
    '''12-Tuple C{(lat1, lon1, azi1, lat2, lon2, azi2, s12, a12, m12, M12, M21, S12)} with
       angles C{lat1}, C{lon1}, C{azi1}, C{lat2}, C{lon2} and C{azi2} and arc C{a12} all in
       C{degrees}, initial C{azi1} and final C{azi2} forward azimuths, distance C{s12} and
       reduced length C{m12} in C{meter}, area C{S12} in C{meter} I{squared} and geodesic
       scale factors C{M12} and C{M21}, both C{scalar}, see U{GeodSolve
       <https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>}.
    '''
    # from GeodSolve --help option -f ... lat1 lon1 azi1 lat2 lon2 azi2 s12 a12 m12 M12 M21 S12
    _Names_ = (_lat1_, _lon1_, _azi1_, _lat2_, _lon2_, _azi2_, _s12_, _a12_, _m12_, _M12_, _M21_, _S12_)
    _Units_ = (_Lat,   _Lon,   _Azi,   _Lat,   _Lon,   _Azi,   _M,    _Deg,  _Pass, _Pass, _Pass, _M2)


class Inverse10Tuple(_GTuple):
    '''10-Tuple C{(a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12)} with arc length
       C{a12} in C{degrees}, distance C{s12} and reduced length C{m12} in C{meter}, area
       C{S12} in C{meter} I{squared} and the sines C{salp1}, C{salp2} and cosines C{calp1},
       C{calp2} of the initial C{1} and final C{2} (forward) azimuths.
    '''
    _Names_ = (_a12_, _s12_, 'salp1', 'calp1', 'salp2', 'calp2', _m12_, _M12_, _M21_, _S12_)
    _Units_ = (_Azi,  _M,    _Pass,   _Pass,   _Pass,   _Pass,   _Pass, _Pass, _Pass, _M2)

    def toGDict(self, **updates):
        '''Convert this C{Inverse10Tuple} to a L{GDict}.

           @kwarg updates: Optional items to apply (C{nam=value} pairs)
        '''
        return _GTuple.toGDict(self, azi1=atan2d(self.salp1, self.calp1),  # PYCHOK namedTuple
                                     azi2=atan2d(self.salp2, self.calp2),  # PYCHOK namedTuple
                                   **updates)  # PYCHOK indent


class _kWrapped(object):  # in .geodesicw
    ''''(INTERNAL) Wrapper for some of I{Karney}'s U{geographiclib
        <https://PyPI.org/project/geographiclib>} classes.
    '''

    @property_ROnce
    def geographiclib(self):
        '''Lazily import C{geographiclib}, provided the U{geographiclib
           <https://PyPI.org/project/geographiclib>} package is installed,
           otherwise raise a C{LazyImportError}.
        '''
        g = _xgeographiclib(self.__class__.__module__, 1, 49)
        from geographiclib.geodesic import Geodesic
        g.Geodesic = Geodesic
        from geographiclib.geodesicline import GeodesicLine
        g.GeodesicLine = GeodesicLine
        from geographiclib.geomath import Math
        g.Math = Math
        return g

    @property_ROnce
    def Math(self):
        '''Wrap the C{geomath.Math} class, provided the U{geographiclib
           <https://PyPI.org/project/geographiclib>} package is installed,
           otherwise C{None}.
        '''
        try:
            g = self.geographiclib
            M = g.Math
            if _version_info(g) < (2,):
                if _K_2_0:
                    M = None
#           elif not _K_2_0:  # XXX set 2.0?
#               _K_2_0 = True
        except (AttributeError, ImportError):
            M = None
        return M

    @property_RO
    def Math_K_2(self):
        return (_2_4_ if _K_2_4 else
               (_2_   if _K_2_0 else _1_)) if self.Math else NN

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
        return _MODS.deprecated.classes.Rhumb7Tuple(self[:-1])


class _Xables(object):
    '''(INTERNAL) Get I{Karney}'s executable paths from/and env vars.
    '''
    bin_ = '/opt/local/bin/'  # '/opt/local/Cellar/geographiclib/2.3/bin/'  # HomeBrew on macOS
    ENV  =  NN

    def GeoConvert(self, *dir_):
        return self._path(self.GeoConvert, *dir_)

    def GeodSolve(self, *dir_):
        return self._path(self.GeodSolve, *dir_)

    def IntersectTool(self, *dir_):
        return self._path(self.IntersectTool, *dir_)

    def RhumbSolve(self, *dir_):
        return self._path(self.RhumbSolve, *dir_)

    def name_version(self, path, base=True):
        # return C{(path + ' ' + version)} of an executable
        if path:
            try:
                r, s = _popen2((path, '--version'))
                if base:
                    path = _MODS.os.path.basename(path)
                    r    = _SPACE_(path, r.split()[-1])
                else:
                    r    = _MODS.streprs.Fmt.PARENSPACED(r, s)
                return r
            except (IndexError, IOError, OSError):
                pass
        return NN

    def _path(self, which, *dir_):
        self.ENV = E = _PYGEODESY(which)
        return _getenv(E, NN) or \
               (NN(dir_[0], which.__name__) if dir_ else E)

    def X_not(self, path):
        return 'env %s=%r not executable' % (self.ENV, path)

    def X_OK(self, path):  # is C{path} an executable?
        os = _MODS.os  # import os
        return os.access(path, os.X_OK) if path else False

_Xables = _Xables()  # PYCHOK singleton


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


def _llz2gl(gl, **llz2):  # see .geodesici._llz2G
    '''(INTERNAL) Set C{gl.lat2, .lon2, .azi2} from C{llz2}.
    '''
    if llz2:
        for n in (_lat2_, _lon2_, _azi2_):  # _lat1_, _lon1_, _azi1_
            v = llz2.get(n, None)
            if v is not None:
                setattr(gl, n, v)
    return gl


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


try:
    from math import fma as _fma  # since 3.13

    def _poly_fma(x, s, *cs):
        for c in cs:
            s = _fma(s, x, c)
        return s

except ImportError:  # Python 3.12-

    def _poly_fma(x, s, *cs):  # PYCHOK redef
        t = _0_0
        for c in cs:
            s, t, _ = _sum3(s * x, t * x, c)
        return s + t

#   def _poly_fma(x, *cs):
#       S = Fhorner(x, *cs, incx=False)
#       return float(S)

def _polynomial(x, cs, i, j):  # PYCHOK shared
    '''(INTERNAL) Like C++ C{GeographicLib.Math.hpp.polyval} but with a
       signature and cascaded summation different from C{karney._sum3}.

       @return: M{sum(x**(j - k - 1) * cs[k] for k in range(i, j)}
    '''
    if (i + 1) < j <= len(cs):  # load _Rtuple._tuple
        try:
            r = _wrapped.Math.polyval(j - i - 1, cs, i, x)
        except AttributeError:
            r = _poly_fma(x, *cs[i:j])
    else:
        r = cs[i]
    return float(r)


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


def _sum2(a, b):  # mimick geomath.Math.sum, actually sum2
    '''Error-free summation like C{geomath.Math.sum}.

       @return: 2-Tuple C{(B{a} + B{b}, residual)}.

       @note: The C{residual} can be the same as B{C{a}} or B{C{b}}.

       @see: U{TwoSum<https://accurate-algorithms.readthedocs.io/en/latest/ch04summation.html>}
             and I{Knuth}'s U{Algorithm 3.1<https://www.TUHH.De/ti3/paper/rump/OgRuOi05.pdf>}.
    '''
    try:
        return _wrapped.Math.sum(a, b)
    except AttributeError:
        # if Algorithm_3_1:
        s = a + b
        r = s - b
        t = s - r
        # elif C_CPP:  # Math::sum C/C++
        #   r -= a; t -= b; t += r; t = -t
        # else:
        t = (a - r) + (b - t)
        # assert fabs(s) >= fabs(t)
        return s, t


def _sum3(s, t, *xs):
    '''Accumulate any B{C{xs}} into a previous C{_sum2(s, t)}.

       @return: 3-Tuple C{(s, t, n)} where C{s} is the sum of B{s}, B{t} and all
                B{xs}, C{t} the residual and C{n} the number of zero C{xs}.

       @see: I{Karney's} C++ U{Accumulator<https://GeographicLib.SourceForge.io/
             C++/doc/Accumulator_8hpp_source.html>} comments for more details and
             function C{_sum2} above.

       @note: Not "error-free", see C{pygeodesy.test/testKarney.py}.
    '''
    z = 0
    for x in xs:
        if x:
            t, r = _sum2(t, x)  # start at the least-
            if s:
                s, t = _sum2(s, t)  # -significant end
                if s:
                    t += r  # accumulate r into t
                else:
                    # assert t == 0  # s == 0 implies t == 0
                    s = unsigned0(r)  # result is r, t = 0
            else:
                s, t = unsigned0(t), r
        else:
            z += 1
    # assert fabs(s) >= fabs(t)
    return s, t, z


def _tand(x):
    '''Return C{tan(B{x})} in C{degrees}.
    '''
    try:
        return _wrapped.Math.tand(x)
    except AttributeError:
        return tand(x)  # Error=None


def _unroll2(lon1, lon2, wrap=False):  # see .ellipsoidalBaseDI._intersects2
    '''Unroll B{C{lon2 - lon1}} like C{geodesic.Geodesic.Inverse}.

       @return: 2-Tuple C{(B{lon2} - B{lon1}, B{lon2})} with B{C{lon2}}
                unrolled if C{B{wrap} is True}, normalized otherwise.
    '''
    if wrap:
        d, t = _diff182(lon1, lon2)
        lon2, t, _ = _sum3(d, t, lon1)  # (lon1 + d) + t
        lon2 += t
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
# Copyright (C) 2016-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
