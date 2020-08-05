
# -*- coding: utf-8 -*-

u'''Wrapper around I{Charles Karney}'s Python implementation
of U{GeographicLib<https://PyPI.org/project/geographiclib>},
provided that package is installed.

Following are U{PyGeodesy<https://PyPI.org/project/PyGeodesy>} classes
and functions transcribed from I{Karney}'s original U{GeographicLib
<https://GeographicLib.SourceForge.io/html/annotated.html>} in C++:

  - L{CassiniSoldner} -- U{CassiniSoldner<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1CassiniSoldner.html>}

  - L{EcefKarney} -- U{Geocentric<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1Geocentric.html>}

  - L{EcefCartesian} -- U{LocalCartesian<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1LocalCartesian.html>}

  - L{EquidistantKarney} -- U{AzimuthalEquidistant<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1AzimuthalEquidistant.html>}

  - L{Elliptic} -- U{EllipticFunction<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1EllipticFunction.html>}

  - L{Etm}, L{ExactTransverseMercator} -- U{TransverseMercatorExact
    <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1TransverseMercatorExact.html>}

  - L{GeoidKarney} -- U{Geoid<https://GeographicLib.SourceForge.io/html/geoid.html>}

  - L{Georef} -- U{Georef<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1Georef.html>}

  - L{GnomonicKarney} -- U{Gnomonic<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1Gnomonic.html>}

  - L{Ups} -- U{PolarStereographic<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1PolarStereographic.html>}

  - L{Utm} -- U{TransverseMercator<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1TransverseMercator.html>}

  - L{UtmUps}, L{Epsg} -- U{UTMUPS<https://GeographicLib.SourceForge.io/html/
    classGeographicLib_1_1UTMUPS.html>}

  - L{atan2d}, L{sincos2}, L{sincos2d}-- U{Math<https://geographiclib.sourceforge.io/html/
    classGeographicLib_1_1Math.html>}

The following U{PyGeodesy<https://PyPI.org/project/PyGeodesy>} classes and
module are wrappers around some of I{Karney}'s Python U{GeographicLib
<https://PyPI.org/project/geographiclib>}:

  - L{ellipsoidalKarney}, L{FrechetKarney}, L{HeightIDWkarney}, L{karney}
'''

from pygeodesy.basics import NAN, property_RO
from pygeodesy.datum import Datums
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import Destination3Tuple, Distance3Tuple
from pygeodesy.utily import unroll180, wrap360

from math import fmod

__all__ = _ALL_LAZY.karney
__version__ = '20.08.04'


class _Adict(dict):
    '''(INTERNAL) Basic C{dict} with key I{and} attribute
       access to the items (minimal version of _NamedDict).
    '''
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            return dict.__getattr__(self, name)


class _Wrapped(object):
    ''''(INTERNAL) Wrapper for some of I{Karney}'s U{GeographicLib
        <https://PyPI.org/project/geographiclib>} classes.
    '''
    _Geodesic     = None
    _GeodesicLine = None
    _geoMath      = None
    _Math         = None

    @property_RO
    def Geodesic(self):
        '''(INTERNAL) Get the wrapped C{Geodesic} class, provided the
           U{GeographicLib<https://PyPI.org/project/geographiclib>}
           package is installed, otherwise throw an C{ImportError}.
        '''
        if _Wrapped._Geodesic is None:
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
                    '''Return the destination lat, lon and reverse azimth
                       in C{degrees}.
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
                    '''Return the non-negative I{angular} distance in C{degrees}.
                    '''
                    # see .FrechetKarney.distance, .HausdorffKarney._distance
                    # and .HeightIDWkarney._distances
                    _, lon2 = unroll180(lon1, lon2, wrap=wrap)  # self.LONG_UNROLL
                    d = self.Inverse(lat1, lon1, lat2, lon2)
                    # XXX self.DISTANCE needed for 'a12'?
                    return abs(d.a12)

                def Inverse3(self, lat1, lon1, lat2, lon2):  # PYCHOK outmask
                    '''Return the distance in C{meter} and the forward and
                       reverse azimuths in C{degrees}.
                    '''
                    m = self.DISTANCE | self.AZIMUTH
                    d = self.Inverse(lat1, lon1, lat2, lon2, m)
                    return Distance3Tuple(d.s12, wrap360(d.azi1), wrap360(d.azi2))

                def Line(self, lat1, lon1, azi1, *caps):
                    return _wrapped.GeodesicLine(self, lat1, lon1, azi1, *caps)

            # Geodesic.Direct.__doc__  = _Geodesic.Direct.__doc__
            # Geodesic.Inverse.__doc__ = _Geodesic.Inverse.__doc__
            # Geodesic.Line.__doc__    = _Geodesic.Line.__doc__

            _Wrapped._Geodesic = Geodesic
        return _Wrapped._Geodesic

    @property_RO
    def GeodesicLine(self):
        '''(INTERNAL) Get the wrapped C{GeodesicLine} class, provided
           the U{GeographicLib<https://PyPI.org/project/geographiclib>}
           package is installed, otherwise throw an C{ImportError}.
        '''
        if _Wrapped._GeodesicLine is None:
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

            _Wrapped._GeodesicLine = GeodesicLine
        return _Wrapped._GeodesicLine

    @property_RO
    def Geodesic_WGS84(self):
        '''(INTERNAL) Get the wrapped C{Geodesic.WGS84} I{instance} iff
           the U{GeographicLib<https://PyPI.org/project/geographiclib>}
           package is installed, otherwise throw an C{ImportError}.
        '''
        return Datums.WGS84.ellipsoid.geodesic

    @property_RO
    def geoMath(self):
        '''(INTERNAL) Get the C{Math} class if the U{GeographicLib
           <https://PyPI.org/project/geographiclib>} package is
           installed or C{False} otherwise.
        '''
        if self._geoMath is None:
            try:
                self._geoMath = self.Math
            except ImportError:
                self._geoMath = False
        return self._geoMath

    @property_RO
    def Math(self):
        '''(INTERNAL) Get the C{Math} class, provided the
           U{GeographicLib<https://PyPI.org/project/geographiclib>}
           package is installed, otherwise throw an C{ImportError}.
        '''
        if _Wrapped._Math is None:
            from geographiclib.geomath import Math
            _Wrapped._Math = Math
        return _Wrapped._Math


_wrapped = _Wrapped()  # imported by .datum


def _diff182(deg0, deg):  # mimick Math.AngDiff
    '''Compute C{deg - deg0}, reduced to C{[-180,180]} accurately.

       @return: 2-Tuple C{(delta_angle, residual)} in C{degrees}.
    '''
    M = _wrapped.geoMath
    if M:
        return M.AngDiff(deg0, deg)

    d, t = _sum2(_norm180(-deg0), _norm180(deg))
    d = _norm180(d)
    if d == 180 and t > 0:
        d = -180
    return _sum2(d, t)


def _fix90(deg):  # mimick Math.LatFix
    '''Replace angle in C{degrees} outside [-90,90] by NAN.

       @return: Angle C{degrees} or NAN.
    '''
    M = _wrapped.geoMath
    if M:
        return M.LatFix(deg)

    return NAN if abs(deg) > 90 else deg


def _fsum2_(*vs):  # see .test/testKarney.py
    '''Cascaded summation, like C{.fmath.fsum_}.

       @arg vs: Values to be added (C{scalar}[]).

       @return: 2-Tuple C{(sum_of_vs, residual)}.

       @note: NOT "error-free", see .test/testKarney.py.

       @see: U{Algorithm 4.1<https://www.ti3.TUHH.DE/paper/rump/OgRuOi05.pdf>}.
    '''
    s = r = 0.0
    for v in vs:
        s, t = _sum2(s, float(v))
        if t:
            r, t = _sum2(r, t)  # inlieu of r += t
    return (s + r), t


def _norm180(deg):  # mimick Math.AngNormalize
    '''Reduce angle in C{degrees} to (-180,180].

       @return: Reduced angle C{degrees}.
    '''
    M = _wrapped.geoMath
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
    d = fmod(deg, 360) if deg else deg
    return (d + 360) if d <= -180 else (d
                     if d <=  180 else (d - 360))


def _sum2(u, v):  # mimick Math::sum, actually sum2
    '''Error-free summation like C{Math::sum}.

       @return: 2-Tuple C{(sum_u_plus_v, residual)}.

       @note: The C{residual} can be the same as
              B{C{u}} or B{C{v}}.

       @see: U{Algorithm 3.1<https://www.ti3.TUHH.DE/paper/rump/OgRuOi05.pdf>}.
    '''
    M = _wrapped.geoMath
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
