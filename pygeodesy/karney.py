
# -*- coding: utf-8 -*-

u'''I{Charles F.F. Karney}'s U{geographiclib<https://PyPI.org/project/geographiclib>} C{geodesic}, wrapped.

Wrapper around I{Karney}'s Python classes C{Geodesic} and C{GeodesicLine} and several C{Math}
functions from I{Karney}'s Python package U{geographiclib<https://PyPI.org/project/geographiclib>},
provided that package is installed.

The I{wrapped} class methods return a L{GDict} instance providing access to the C{dict} items
by C{key} or by C{attribute} name.

I{Karney}-based functionality
=============================

1. The following classes and functions in C{pygeodesy}

  - L{AlbersEqualArea}, L{AlbersEqualArea2}, L{AlbersEqualArea4},
    L{AlbersEqualAreaCylindrical}, L{AlbersEqualAreaNorth}, L{AlbersEqualAreaSouth} --
    U{AlbersEqualArea<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1AlbersEqualArea.html>}

  - L{CassiniSoldner} -- U{CassiniSoldner<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1CassiniSoldner.html>}

  - L{EcefKarney} -- U{Geocentric<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1Geocentric.html>}

  - L{Elliptic} -- U{EllipticFunction<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1EllipticFunction.html>}

  - L{EquidistantExact}, L{EquidistantGeodSolve}, L{EquidistantKarney} -- U{AzimuthalEquidistant
    <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1AzimuthalEquidistant.html>}

  - L{Etm}, L{ExactTransverseMercator} -- U{TransverseMercatorExact
    <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1TransverseMercatorExact.html>}

  - L{GeodesicAreaExact}, L{PolygonArea} -- U{PolygonArea<https://GeographicLib.SourceForge.io/
    html/classGeographicLib_1_1PolygonAreaT.html>}

  - L{GeodesicExact}, L{GeodesicLineExact} -- U{GeodesicExact<https://GeographicLib.SourceForge.io/
    html/classGeographicLib_1_1GeodesicExact.html>}, U{GeodesicLineExact<https://GeographicLib.SourceForge.io/
    html/classGeographicLib_1_1GeodesicLineExact.html>}

  - L{GeoidKarney} -- U{Geoid<https://GeographicLib.SourceForge.io/html/geoid.html>}

  - L{Georef} -- U{Georef<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1Georef.html>}

  - L{GnomonicExact}, L{GnomonicGeodSolve}, L{GnomonicKarney} -- U{Gnomonic
    <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Gnomonic.html>}

  - L{LocalCartesian}, L{Ltp} -- U{LocalCartesian<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1LocalCartesian.html>}

  - L{Ups} -- U{PolarStereographic<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1PolarStereographic.html>}

  - L{Utm} -- U{TransverseMercator<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1TransverseMercator.html>}

  - L{UtmUps}, L{Epsg} -- U{UTMUPS<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1UTMUPS.html>}

  - L{atand}, L{atan2d}, L{sincos2}, L{sincos2d}-- U{Math<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1Math.html>}

are I{transcoded} from C++ classes in I{Karney}'s U{GeographicLib<https://GeographicLib.SourceForge.io/html/annotated.html>}.

2. These C{pygeodesy} modules and classes

  - L{ellipsoidalGeodSolve}, L{ellipsoidalKarney}, L{geodsolve}, L{karney}
  - L{EquidistantKarney}, L{FrechetKarney}, L{GeodesicSolve}, L{GeodesicLineSolve}, L{GnomonicGeodSolve},
    L{GnomonicKarney}, L{HeightIDWkarney}

are or use I{wrappers} around I{Karney}'s Python U{geographiclib<https://PyPI.org/project/geographiclib>}
C{geodesic} or C++ utility U{GeodSolve<https://GeographicLib.SourceForge.io/html/utilities.html>}.

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

  - L{excessKarney_}, L{sphericalTrigonometry.areaOf}

in C{pygeodesy} are based on I{Karney}'s post U{Area of a spherical polygon
<http://OSGeo-org.1560.x6.Nabble.com/Area-of-a-spherical-polygon-td3841625.html>}.
'''

from pygeodesy.basics import copysign0, _xImportError, _xversion
from pygeodesy.datums import _ellipsoidal_datum, _WGS84  # PYCHOK used!
from pygeodesy.ellipsoids import Ellipsoid2
from pygeodesy.errors import _ValueError, _xkwds  # PYCHOK shared
from pygeodesy.interns import MAX as _MAX, NAN, NN, _lat1_, \
                             _lat2_, _lon2_, _0_0, _1_0, _2_0, \
                             _16_0, _180_0, _360_0
from pygeodesy.iters import PointsIter
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _Dict, _NamedTuple, _Pass
from pygeodesy.namedTuples import Destination3Tuple, Distance3Tuple
from pygeodesy.props import Property, Property_RO
from pygeodesy.units import Degrees as _Deg, Meter as _M, \
                            Meter2 as _M2, _1mm as _TOL_M  # PYCHOK shared
from pygeodesy.utily import atan2d, unroll180, wrap360

from math import fmod


__all__ = _ALL_LAZY.karney
__version__ = '21.08.21'

_16th = _1_0 / _16_0

_a12_  = 'a12'
_azi1_ = 'azi1'
_azi2_ = 'azi2'
_lon1_ = 'lon1'
_m12_  = 'm12'
_M12_  = 'M12'
_M21_  = 'M21'
_s12_  = 's12'
_S12_  = 'S12'


class _GTuple(_NamedTuple):  # in .testNamedTuples
    '''(INTERNAL) Helper.
    '''
    def toGDict(self, **updates):
        '''Convert this C{*Tuple} to a L{GDict}.

           @kwarg updates: Optional updates to apply
                           (C{nam=value} pairs)

        '''
        r = GDict(zip(self._Names_, self))
        if updates:
            r.update(updates)
        return r


class Direct9Tuple(_GTuple):
    '''9-Tuple C{(a12, lat2, lon2, azi2, s12, m12, M12, M21, S12)} with arc
       length C{a12}, angles C{lat2}, C{lon2} and azimuth C{azi2} in C{degrees},
       distance C{s12} and reduced length C{m12} in C{meter} and area C{S12} in
       C{meter} I{squared}.
    '''
    _Names_ = (_a12_, _lat2_, _lon2_, _azi2_, _s12_, _m12_, _M12_, _M21_, _S12_)
    _Units_ = (_Deg,  _Deg,   _Deg,   _Deg,   _M,    _Pass, _Pass, _Pass, _M2)


class GDict(_Dict):
    '''Basic C{dict} with both key I{and} attribute access
       to the C{dict} items.

       Results of all C{geodesic} methods are returned as a
       L{GDict} instance.
    '''
    def toDirect9Tuple(self, dflt=NAN):
        '''Convert this L{GDict} result to a 9-tuple, like I{Karney}'s
           method C{geographiclib.geodesic.Geodesic._GenDirect}.

           @return: L{Direct9Tuple}C{(a12, lat2, lon2, azi2,
                    s12, m12, M12, M21, S12)}
        '''

        def _9Tuple(a12=dflt, lat2=dflt, lon2=dflt, azi2=dflt,
                    s12=dflt, m12=dflt,  M12=dflt,  M21=dflt, S12=dflt,
                  **unused):
            return Direct9Tuple(a12, lat2, lon2, azi2,
                                s12, m12,  M12,  M21, S12)

        return _9Tuple(**self)

    def toGeodSolve12Tuple(self, dflt=NAN):  # PYCHOK 12 args
        '''Convert this L{GDict} result to a 12-Tuple, compatible with
           I{Karney}'s U{GeodSolve<https://GeographicLib.SourceForge.io/html/GeodSolve.1.html>}
           result.

           @return: L{GeodSolve12Tuple}C{(lat1, lon1, azi1, lat2, lon2, azi2,
                    s12, a12, m12, M12, M21, S12)}.
        '''

        def _12Tuple(lat1=dflt, lon1=dflt, azi1=dflt,
                     lat2=dflt, lon2=dflt, azi2=dflt, s12=dflt, a12=dflt,
                                                      m12=dflt, M12=dflt,
                                                      M21=dflt, S12=dflt,
                                                    **unused):
            return GeodSolve12Tuple(lat1, lon1, azi1, lat2, lon2, azi2,
                                    s12,  a12,  m12,  M12,  M21,  S12)

        return _12Tuple(**self)

    def toInverse10Tuple(self, dflt=NAN):
        '''Convert this L{GDict} result to a 10-tuple, like I{Karney}'s
           method C{geographiclib.geodesic.Geodesic._GenInverse}.

           @return: L{Inverse10Tuple}C{(a12, s12, salp1, calp1,
                    salp2, calp2, m12, M12, M21, S12)}.
        '''

        def _10Tuple(a12=dflt, s12=dflt, salp1=dflt, calp1=dflt,
                                         salp2=dflt, calp2=dflt,
                                         m12=dflt,   M12=dflt,
                                         M21=dflt,   S12=dflt,
                                       **unused):
            return Inverse10Tuple(a12, s12, salp1, calp1, salp2, calp2,
                                       m12,   M12,   M21,   S12)

        return _10Tuple(**self)


class GeodesicError(_ValueError):  # PYCHOK exported
    '''Error raised for I{geodesicx}'s methods for lack of convergence
       or other issues.
    '''
    pass


class GeodSolve12Tuple(_GTuple):
    '''12-Tuple C{(lat1, lon1, azi1, lat2, lon2, azi2, s12, a12, m12, M12, M21, S12)} with
       angles C{lat1}, C{lon1}, C{azi1}, C{lat2}, C{lon2} and C{azi2} and arc C{a12} all in
       C{degrees}, distance C{s12} and reduced length C{m12} in C{meter}, area C{S12} in
       C{meter} I{squared} and geodesic scales C{M12} and C{M21} C{scalar}, see U{GeodSolve
       <https://GeographicLib.SourceForge.io/html/GeodSolve.1.html>}.
    '''
    # from GeodSolve --help option -f ... lat1 lon1 azi1 lat2 lon2 azi2 s12 a12 m12 M12 M21 S12
    _Names_ = (_lat1_, _lon1_, _azi1_, _lat2_, _lon2_, _azi2_, _s12_, _a12_, _m12_, _M12_, _M21_, _S12_)
    _Units_ = (_Deg,   _Deg,   _Deg,   _Deg,   _Deg,   _Deg,   _M,    _Deg,  _Pass, _Pass, _Pass, _M2)


class Inverse10Tuple(_GTuple):
    '''10-Tuple C{(a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12)}
       with arc length C{a12} in C{degrees}, distance C{s12} and reduced
       length C{m12} in C{meter}, area C{S12} in C{meter} I{squared} and
       sines and cosines of initial and final (forward) azimuths.
    '''
    _Names_ = (_a12_, _s12_, 'salp1', 'calp1', 'salp2', 'calp2', _m12_, _M12_, _M21_, _S12_)
    _Units_ = (_Deg,  _M,    _Pass,   _Pass,   _Pass,   _Pass,   _Pass, _Pass, _Pass, _M2)


class _Wrapped(object):
    ''''(INTERNAL) Wrapper for some of I{Karney}'s U{geographiclib
        <https://PyPI.org/project/geographiclib>} classes.
    '''
    _geographiclib = None

    @Property_RO  # MCCABE 14
    def Geodesic(self):
        '''Get the I{wrapped} C{Geodesic} class, provided the U{geographiclib
           <https://PyPI.org/project/geographiclib>} package is installed,
           otherwise an C{ImportError}.
        '''
        self._xgeographiclib(_Wrapped.Geodesic)
        from geographiclib.geodesic import Geodesic as _Geodesic

        class Geodesic(_Geodesic):
            '''I{Karney}'s U{Geodesic<https://GeographicLib.SourceForge.io/html/
               python/code.html#geographiclib.geodesic.Geodesic>} wrapper.
            '''

            def __init__(self, a_ellipsoid, f=None, name=NN):  # PYCHOK signature
                '''New C{Geodesic} instance.

                   @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}) or datum
                                     (L{datum}) or the equatorial radius
                                     of the ellipsoid (C{meter}).
                   @arg f: The flattening of the ellipsoid (C{scalar}) if
                           B{C{a_ellipsoid}) is specified as C{meter}.
                   @kwarg name: Optional name (C{str}).
                '''
                if f is None:
                    self._E = _ellipsoidal_datum(a_ellipsoid, name=name).ellipsoid
                else:
                    self._E =  Ellipsoid2(a_ellipsoid, f, name=name)
                _Geodesic.__init__(self, self.ellipsoid.a, self.ellipsoid.f)

            Area   = _Geodesic.Polygon  # like GeodesicExact.Area

            _debug = 0  # like .geodesicx.bases._GeodesicBase

            @Property
            def debug(self):
                '''Get the C{debug} option (C{bool}).
                '''
                return bool(self._debug)

            @debug.setter  # PYCHOK setter!
            def debug(self, debug):
                '''Set the C{debug} option.

                   @arg debug: Include more details in results (C{bool}).
                '''
                from pygeodesy.geodesicx.gxbases import Caps
                self._debug = Caps._DEBUG_ALL if debug else 0

            def Direct(self, lat1, lon1, azi1, s12, *outmask):
                '''Return the C{Direct} result.
                '''
                d = _Geodesic.Direct(self, lat1, lon1, azi1, s12, *outmask)
                return GDict(d)

            def Direct3(self, lat1, lon1, azi1, s12):  # PYCHOK outmask
                '''Return the destination lat, lon and reverse azimuth
                   (final bearing) in C{degrees}.

                   @return: L{Destination3Tuple}C{(lat, lon, final)}.
                '''
                m = self.AZIMUTH | self.LATITUDE | self.LONGITUDE
                d = self.Direct(lat1, lon1, azi1, s12, m)
                return Destination3Tuple(d.lat2, d.lon2, d.azi2)

            @Property_RO
            def ellipsoid(self):
                '''Get this geodesic's ellipsoid (C{Ellipsoid[2]}).
                '''
                return self._E

            def _GDictDirect(self, lat, lon, azi, arcmode, s12_a12,
                                                  outmask=_Geodesic.STANDARD):
                '''(INTERNAL) Get C{._GenDirect} result as C{GDict}.
                '''
                def _toGdict(a12, lat2, lon2, azi2,
                             s12, m12,  M12,  M21, S12):
                    return GDict(a12=a12, lat2=lat2, lon2=lon2, azi2=azi2,
                                 s12=s12, m12=m12,   M12=M12,   M21=M21, S12=S12)

                t = _Geodesic._GenDirect(self, lat, lon, azi, arcmode, s12_a12, outmask)
                return _toGdict(*t)

            def _GDictInverse(self, lat1, lon1, lat2, lon2, outmask=_Geodesic.STANDARD):
                '''(INTERNAL) Get C{._GenInverse} result as C{GDict}.
                '''
                def _toGdict(a12, s12, salp1, calp1, salp2, calp2,
                                       m12,   M12,   M21,   S12):
                    return GDict(a12=a12, s12=s12,
                                 salp1=salp1, calp1=calp1, azi1=atan2d(salp1, calp1),
                                 salp2=salp2, calp2=calp2, azi2=atan2d(salp2, calp2),
                                 m12=m12, M12=M12, M21=M21, S12=S12, lon1=lon1)

                t =_Geodesic._GenInverse(self, lat1, lon1, lat2, lon2, outmask)
                return _toGdict(*t)

            def Inverse(self, lat1, lon1, lat2, lon2, *outmask):
                '''Return the C{Inverse} result.
                '''
                d = _Geodesic.Inverse(self, lat1, lon1, lat2, lon2, *outmask)
                return GDict(d)

            def Inverse1(self, lat1, lon1, lat2, lon2, wrap=False):
                '''Return the non-negative, I{angular} distance in C{degrees}.
                '''
                # see .FrechetKarney.distance, .HausdorffKarney._distance
                # and .HeightIDWkarney._distances
                _, lon2 = unroll180(lon1, lon2, wrap=wrap)  # self.LONG_UNROLL
                d = self.Inverse(lat1, lon1, lat2, lon2)
                # XXX self.DISTANCE needed for 'a12'?
                return abs(d.a12)

            def Inverse3(self, lat1, lon1, lat2, lon2):  # PYCHOK outmask
                '''Return the distance in C{meter} and the forward and
                   reverse azimuths (initial and final bearing) in C{degrees}.

                   @return: L{Distance3Tuple}C{(distance, initial, final)}.
                '''
                m = self.DISTANCE | self.AZIMUTH
                d = self.Inverse(lat1, lon1, lat2, lon2, m)
                return Distance3Tuple(d.s12, wrap360(d.azi1), wrap360(d.azi2))

            def Line(self, lat1, lon1, azi1, *caps):
                return _wrapped.GeodesicLine(self, lat1, lon1, azi1, *caps)

            _Line = Line  # for .azimuthal._GnomonicBase.reverse, matching
            # PYCHOK .geodesicx.GeodesictExact._Line and -._GDictDirect

        # Geodesic.Direct.__doc__  = _Geodesic.Direct.__doc__
        # Geodesic.Inverse.__doc__ = _Geodesic.Inverse.__doc__
        # Geodesic.Line.__doc__    = _Geodesic.Line.__doc__
        return Geodesic

    @Property_RO
    def GeodesicLine(self):
        '''Get the I{wrapped} C{GeodesicLine} class, provided the U{geographiclib
           <https://PyPI.org/project/geographiclib>} package is installed,
           otherwise an C{ImportError}.
        '''
        self._xgeographiclib(_Wrapped.GeodesicLine)
        from geographiclib.geodesicline import GeodesicLine as _GeodesicLine

        class GeodesicLine(_GeodesicLine):
            '''I{Karney}'s U{GeodesicLine <https://GeographicLib.SourceForge.io/html/
               python/code.html#geographiclib.geodesicline.GeodesicLine>} wrapper.
            '''
            def ArcPosition(self, a12, *outmask):
                d = _GeodesicLine.ArcPosition(self, a12, *outmask)
                return GDict(d)

            def Position(self, s12, *outmask):
                d = _GeodesicLine.Position(self, s12, *outmask)
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
        return _WGS84.ellipsoid.geodesic

    @Property_RO
    def geographiclib(self):
        '''Get the imported C{geographiclib}, provided the U{geographiclib
           <https://PyPI.org/project/geographiclib>} package is installed,
           otherwise an C{ImportError}.
        '''
        return self._xgeographiclib(_Wrapped.geographiclib)

    @Property_RO
    def Math(self):
        '''Get the C{Math} class, provided the U{geographiclib
           <https://PyPI.org/project/geographiclib>} package is
           installed, otherwise C{None}.
        '''
        try:
            self._xgeographiclib(_Wrapped.Math)
            from geographiclib.geomath import Math
            # replace karney. with Math. functions
            from pygeodesy import karney
            karney._around    = Math.AngRound
            karney._diff182   = Math.AngDiff
            karney._fix90     = Math.LatFix
            karney._isfinite  = Math.isfinite
            karney._norm180   = Math.AngNormalize
            karney._sum2      = Math.sum
            try:  # geographiclib 1.49
                karney._remainder = Math.remainder
            except AttributeError:
                pass
        except ImportError:
            Math = None
        return Math

    def _xgeographiclib(self, where):
        '''(INTERNAL) Import C{geographiclib}.
        '''
        if self._geographiclib is None:
            try:
                import geographiclib as g
            except ImportError as x:
                raise _xImportError(x, _Wrapped, name=where.name)
            self._geographiclib = _xversion(g, _Wrapped, 1, 49, name=where.name)
        return self._geographiclib

_wrapped = _Wrapped()  # PYCHOK singleton, .datum


def _around(x):
    '''Coarsen a scalar to zero.

       @return: Coarsened value (C{float}).
    '''
    M = _wrapped.Math
    if M:
        return M.AngRound(x)

    if x:
        y = abs(x)
        if y < _16th:
            y = _16th - y
            x = copysign0(_16th - y, x)
    else:
        x = _0_0  # -0 to 0
    return x


def _diff182(deg0, deg):  # mimick Math.AngDiff
    '''Compute C{deg - deg0}, reduced to C{[-180,180]} accurately.

       @return: 2-Tuple C{(delta_angle, residual)} in C{degrees}.
    '''
    M = _wrapped.Math
    if M:
        return M.AngDiff(deg0, deg)

    d, t = _sum2(_norm180(-deg0), _norm180(deg))
    d = _norm180(d)
    if t > 0 and d == _180_0:
        d = -_180_0
    return _sum2(d, t)


# def _Equidistant(equidistant, exact=False, geodsolve=False):
#     # (INTERNAL) Get the C{EquidistantExact}, C{-GeodSolve} or
#     # C{-Karney} class if B{C{equidistant}} in not callable.
#     if equidistant is None or not callable(equidistant):
#         if exact:
#             from pygeodesy.azimuthal import EquidistantExact as equidistant
#         elif geodsolve:
#             from pygeodesy.azimuthal import EquidistantGeodSolve as equidistant
#         else:
#             from pygeodesy.azimuthal import EquidistantKarney as equidistant
#     return equidistant


def _fix90(deg):  # mimick Math.LatFix
    '''Replace angle in C{degrees} outside [-90,90] by NAN.

       @return: Angle C{degrees} or NAN.
    '''
    M = _wrapped.Math
    if M:
        return M.LatFix(deg)

    return NAN if abs(deg) > 90 else deg


def _fsum2_(*vs):  # see .test/testKarney.py
    '''Cascaded summation, like C{.fmath.fsum_}.

       @arg vs: Values to be added (C{scalar}[]).

       @return: 2-Tuple C{(sum_of_vs, residual)}.

       @note: NOT "error-free", see .test/testKarney.py.

       @see: U{Algorithm 4.1<http://www.ti3.TUHH.De/paper/rump/OgRuOi05.pdf>}.
    '''
    s = t = _0_0
    for v in vs:
        s, t = _3sum2(s, t, v)
    return _sum2(s, t)


def _isfinite(x):  # mimick Math.AngNormalize
    '''Check finitenessof C{x}.

       @return: True if finite.
    '''
    M = _wrapped.Math
    if M:
        return M.isfinite(x)

    return abs(x) <= _MAX


def _norm180(deg):  # mimick Math.AngNormalize
    '''Reduce angle in C{degrees} to (-180,180].

       @return: Reduced angle C{degrees}.
    '''
    M = _wrapped.Math
    if M:
        return M.AngNormalize(deg)

    # with Python 2.7.16 and 3.7.3 on macOS 10.13.6
    #  fmod( 0,   360) ==  0.0
    #  fmod( 360, 360) ==  0.0
    #  fmod(-0,   360) ==  0.0
    #  fmod(-0.0, 360) == -0.0
    #  fmod(-360, 360) == -0.0
    # however, using the % operator ...
    #    0   % 360 == 0
    #  360   % 360 == 0
    #  360.0 % 360 == 0.0
    #   -0   % 360 == 0
    # -360   % 360 == 0
    #   -0.0 % 360 == 0.0
    # -360.0 % 360 == 0.0

    # On Windows 32-bit with Python 2.7, math.fmod(-0.0, 360)
    # == +0.0.  This fixes this bug.  See also Math::AngNormalize
    # in the C++ library, Math::sincosd has a similar fix.
    d = (fmod(deg, _360_0) if _isfinite(deg) else NAN) if deg else deg
    return (d + _360_0) if d < -_180_0 else (d  # XXX was <= twice
                        if d <  _180_0 else (d - _360_0))


def _polygon(geodesic, points, closed, line, wrap):
    '''(INTERNAL) Compute the area or perimeter of a polygon,
        using a L{GeodesicExact}, L{GeodesicSolve} or (if the
        C{geographiclib} package is installed) a C{Geodesic}
        or C{_wrapped.Geodesic} instance.
    '''
    if not wrap:  # capability LONG_UNROLL can't be off
        raise _ValueError(wrap=wrap)

    gP = geodesic.Polygon(line)
    p_ = gP.AddPoint

    Ps = PointsIter(points, loop=1)  # base=LatLonEllipsoidalBase(0, 0)
    p0 = Ps[0]

    # note, lon deltas are unrolled, by default
    p_(p0.lat, p0.lon)
    for p in Ps.iterate(closed=closed):
        p_(p.lat, p.lon)
    if closed and line and p != p0:
        p_(p0.lat, p0.lon)

    # gP.Compute returns (number_of_points, perimeter, signed area)
    return gP.Compute(False, True)[1 if line else 2]


def _remainder(x, y):
    '''Remainder of C{x / y}.

       @return: Remainder in the range M{[-y / 2, y / 2]}.
    '''
    M = _wrapped.Math
    if M:
        try:  # geographiclib 1.49
            return M.remainder(x, y)
        except AttributeError:
            pass

    z = (fmod(x, y) if _isfinite(x) else NAN) if x else x
    # On Windows 32-bit with python 2.7, math.fmod(-0.0, 360)
    # == +0.0.  This fixes this bug.  See also Math::AngNormalize
    # in the C++ library, Math.sincosd has a similar fix.
    h = y / _2_0
    return (z + y) if z < -h else (z
                   if z <  h else (z - y))


def _sum2(u, v):  # mimick Math::sum, actually sum2
    '''Error-free summation like C{Math::sum}.

       @return: 2-Tuple C{(sum_u_plus_v, residual)}.

       @note: The C{residual} can be the same as
              B{C{u}} or B{C{v}}.

       @see: U{Algorithm 3.1<http://www.ti3.TUHH.De/paper/rump/OgRuOi05.pdf>}.
    '''
    M = _wrapped.Math
    if M:
        return M.sum(u, v)

    s = u + v
    r = s - v
    t = s - r
    # if True:  # Algorithm 3.1
    t = (u - r) + (v - t)

    # else:  # in Math C/C++
    #   r -= u
    #   t -= v
    #   t = -(r + t)

    # u + v =       s      + t
    #       = round(u + v) + t
    return s, t


def _3sum2(s, t, x):
    '''Accumulate B{C{x}} into C{_sum2(s, t)}.

       @see: I{Karney's} C++ U{Accumulator<https://GeographicLib.SourceForge.io/
             html/Accumulator_8hpp_source.html>} comments for more details.
    '''
    if x:
        t, u = _sum2(t, x)  # start at least-
        s, t = _sum2(s, t)  # significant end
        if s:
            t += u  # accumulate u to t
        else:
            s  = u  # result is u
    return s, t


def _unroll2(lon1, lon2, wrap=False):  # see .ellipsoidalBaseDI._intersects2
    '''Unroll B{C{lon2 - lon1}} like C{geodesic.Geodesic.Inverse}.

       @return: 2-Tuple C{(lon2 - lon1, lon2)} with B{C{lon2}} unrolled
                if B{C{wrap}} is C{True}, normalized otherwise.
    '''
    if wrap:
        d, t = _diff182(lon1, lon2)
        lon2 = (lon1 + d) + t  # _fsum2_(lon1, d, t)
    else:
        lon2 = _norm180(lon2)
    return (lon2 - lon1), lon2

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
