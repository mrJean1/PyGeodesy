
# -*- coding: utf-8 -*-

u'''Wrapper around I{Charles Karney}'s Python classes C{Geodesic} and
C{GeodesicLine} and C{Math} functions C{AngDiff}, C{AngNormalize}, C{LatFix} and
C{sum} from I{Karney's} U{geographiclib<https://PyPI.org/project/geographiclib>},
Python package, provided that package is installed.

The I{wrapped} class methods return a C{named._NamedDict}-like instance providing
access to the C{dict} items by key or by attribute.

Following are U{PyGeodesy<https://PyPI.org/project/PyGeodesy>} classes
and functions I{transcribed} from I{Karney}'s original U{GeographicLib
<https://GeographicLib.SourceForge.io/html/annotated.html>} in C++:

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

  - L{EquidistantKarney} -- U{AzimuthalEquidistant<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1AzimuthalEquidistant.html>}

  - L{Etm}, L{ExactTransverseMercator} -- U{TransverseMercatorExact
    <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1TransverseMercatorExact.html>}

  - L{GeoidKarney} -- U{Geoid<https://GeographicLib.SourceForge.io/html/geoid.html>}

  - L{Georef} -- U{Georef<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1Georef.html>}

  - L{GnomonicKarney} -- U{Gnomonic<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1Gnomonic.html>}

  - L{LocalCartesian}, L{Ltp} -- U{LocalCartesian<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1LocalCartesian.html>}

  - L{Ups} -- U{PolarStereographic<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1PolarStereographic.html>}

  - L{Utm} -- U{TransverseMercator<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1TransverseMercator.html>}

  - L{UtmUps}, L{Epsg} -- U{UTMUPS<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1UTMUPS.html>}

  - L{atan2d}, L{sincos2}, L{sincos2d}-- U{Math<https://GeographicLib.sourceforge.io/html/
    classGeographicLib_1_1Math.html>}

The following U{PyGeodesy<https://PyPI.org/project/PyGeodesy>} module, classes and functions are I{wrappers}
around some of I{Karney}'s Python U{geographiclib<https://PyPI.org/project/geographiclib>}:

  - L{karney}, L{ellipsoidalKarney}, L{EquidistantKarney}, L{FrechetKarney}, L{GnomonicKarney}, L{HeightIDWkarney},
    L{ellipsoidalKarney.areaOf}, L{ellipsoidalKarney.isclockwise}, L{ellipsoidalKarney.perimeterOf}

Lastly, spherical functions:

  - L{excessKarney_}, L{sphericalTrigonometry.areaOf}

are based on I{Karney}'s post U{Area of a spherical polygon
<http://OSGeo-org.1560.x6.Nabble.com/Area-of-a-spherical-polygon-td3841625.html>} and ellipsoidal
functions and methods:

  - L{ellipsoidalKarney.intersections2}, L{ellipsoidalKarney.nearestOn}, L{ellipsoidalKarney.LatLon.intersections2}, L{ellipsoidalKarney.LatLon.nearestOn}

are implementations of I{Karney}'s post U{The B{ellipsoidal} case
<https://GIS.StackExchange.com/questions/48937/calculating-intersection-of-two-circles>} and paper
U{Geodesics on an ellipsoid of revolution<https://ArXiv.org/pdf/1102.1215.pdf>} (pp 20-21, section
B{14. MARITIME BOUNDARIES}).
'''

from pygeodesy.basics import _xversion
from pygeodesy.datums import _WGS84  # PYCHOK used!
from pygeodesy.interns import NAN, _0_0, _180_0, _360_0
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.namedTuples import Destination3Tuple, Distance3Tuple
from pygeodesy.props import Property_RO
from pygeodesy.utily import unroll180, wrap360

from math import fmod

__all__ = _ALL_LAZY.karney
__version__ = '21.04.15'


class _Adict(dict):
    '''(INTERNAL) Basic C{dict} with key I{and} attribute access
       to the C{dict} items, minimal version of C{named._NamedDict}.
    '''
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            return dict.__getattr__(self, name)


class _Wrapped(object):
    ''''(INTERNAL) Wrapper for some of I{Karney}'s U{geographiclib
        <https://PyPI.org/project/geographiclib>} classes.
    '''
    _geographiclib = None

    @Property_RO
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

            def Direct(self, lat1, lon1, azi1, s12, *outmask):
                '''Return the C{Direct} result.
                '''
                d = _Geodesic.Direct(self, lat1, lon1, azi1, s12, *outmask)
                return _Adict(d)

            def Direct3(self, lat1, lon1, azi1, s12):  # PYCHOK outmask
                '''Return the destination lat, lon and reverse azimuth
                   (final bearing) in C{degrees}.
                '''
                m = self.AZIMUTH | self.LATITUDE | self.LONGITUDE
                d = self.Direct(lat1, lon1, azi1, s12, m)
                return Destination3Tuple(d.lat2, d.lon2, d.azi2)

            def Inverse(self, lat1, lon1, lat2, lon2, *outmask):
                '''Return the C{Inverse} result.
                '''
                d = _Geodesic.Inverse(self, lat1, lon1, lat2, lon2, *outmask)
                return _Adict(d)

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
                '''
                m = self.DISTANCE | self.AZIMUTH
                d = self.Inverse(lat1, lon1, lat2, lon2, m)
                return Distance3Tuple(d.s12, wrap360(d.azi1), wrap360(d.azi2))

            def Line(self, lat1, lon1, azi1, *caps):
                return _wrapped.GeodesicLine(self, lat1, lon1, azi1, *caps)

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
                return _Adict(d)

            def Position(self, s12, *outmask):
                d = _GeodesicLine.Position(self, s12, *outmask)
                return _Adict(d)

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
            karney._diff182 = Math.AngDiff
            karney._fix90   = Math.LatFix
            karney._norm180 = Math.AngNormalize
            karney._sum2    = Math.sum
        except ImportError:
            Math = None
        return Math

    def _xgeographiclib(self, where):
        '''(INTERNAL) Import C{geographiclib}.
        '''
        if self._geographiclib is None:
            import geographiclib as g
            self._geographiclib = _xversion(g, _Wrapped, 1, 49, name=where.name)
        return self._geographiclib

_wrapped = _Wrapped()  # PYCHOK singleton, .datum


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
    s = r = t = _0_0
    for v in vs:
        s, t = _sum2(s, float(v))
        if t:
            r, t = _sum2(r, t)  # inlieu of r += t
    return (s + r), t


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
    # in the C++ library.  Math::sincosd has a similar fix.
    d = fmod(deg, _360_0) if deg else deg
    return (d + _360_0) if d <= -_180_0 else (d
                        if d <=  _180_0 else (d - _360_0))


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


def _unroll2(lon1, lon2, wrap=False):  # see .ellipsoidalBase._intersects2
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
