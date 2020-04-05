
# -*- coding: utf-8 -*-

u'''Wrapper around I{Charles Karney's} Python implementation
of U{GeographicLib <https://PyPI.org/project/geographiclib>},
provided that package is installed.

'''
from pygeodesy.basics import NAN, _Adict, property_RO
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import Destination3Tuple, Distance3Tuple
from pygeodesy.utily import unroll180, wrap360

from math import fmod

__all__ = _ALL_LAZY.karney
__version__ = '20.04.05'


class _Wrapped(object):
    ''''(INTERNAL) Wrapper for some of Karney's U{GeographicLib
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
           package is installed, an C{ImportError} otherwise.
        '''
        if _Wrapped._Geodesic is None:

            from geographiclib.geodesic import Geodesic as _Geodesic

            class Geodesic(_Geodesic):
                '''Karney U{Geodesic <https://geographiclib.SourceForge.io/html/
                   python/code.html#geographiclib.geodesic.Geodesic>} wrapper.
                '''
                WGS84 = None

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
                    # see .FrechetKarney.distance, .HuasdorffKarney._distance
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

            Geodesic.WGS84 = Geodesic(_Geodesic.WGS84.a, _Geodesic.WGS84.f)

        return _Wrapped._Geodesic

    @property_RO
    def GeodesicLine(self):
        '''(INTERNAL) Get the wrapped C{GeodesicLine} class, provided
           the U{GeographicLib<https://PyPI.org/project/geographiclib>}
           package is installed, an C{ImportError} otherwise.
        '''
        if _Wrapped._GeodesicLine is None:

            from geographiclib.geodesicline import GeodesicLine as _GeodesicLine

            class GeodesicLine(_GeodesicLine):
                '''Karney U{GeodesicLine <https://geographiclib.SourceForge.io/html/
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
           package is installed, an C{ImportError} otherwise.
        '''
        if _Wrapped._Math is None:

            from geographiclib.geomath import Math

            _Wrapped._Math = Math

        return _Wrapped._Math


_wrapped = _Wrapped()


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
    '''Error-free transformation like C{Math::sum}.

       @return: 2-Tuple C{(sum_u_plus_v, residual)}.

       @note: The C{residual} can be the same as one
              of the first two arguments.
    '''
    M = _wrapped.geoMath
    if M:
        return M.sum(u, v)

    s = u + v
    r = s - v
    t = s - r
    r -= u
    t -= v
    t = -(r + t)
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
