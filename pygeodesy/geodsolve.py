
# -*- coding: utf-8 -*-

u'''Wrapper to invoke I{Karney}'s U{GeodSolve
<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>} utility
as an (exact) geodesic, but intended I{for testing purposes only}.

Set env variable C{PYGEODESY_GEODSOLVE} to the (fully qualified) path
of the C{GeodSolve} executable.
'''

# from pygeodesy.basics import _xinstanceof  # from .karney
from pygeodesy.interns import NN, _a12_, _azi1_, _azi2_, \
                             _lat1_, _lat2_, _lon1_, _lon2_, _m12_, \
                             _M12_, _M21_, _s12_, _S12_, _UNDER_
from pygeodesy.interns import _not_  # PYCHOK used!
from pygeodesy.karney import _Azi, Caps, _Deg, GeodesicError, _GTuple, \
                             _Pass, _Lat, _Lon, _M, _M2, _sincos2d, \
                             _xinstanceof, wrap360
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS, \
                             _getenv, printf
from pygeodesy.namedTuples import Destination3Tuple, Distance3Tuple
from pygeodesy.props import Property, Property_RO
from pygeodesy.solveBase import _LineSolveBase, _SolveBase
# from pygeodesy.utily import wrap360  # from .karney

__all__ = _ALL_LAZY.geodsolve
__version__ = '22.07.09'

_PYGEODESY_GEODSOLVE_ = 'PYGEODESY_GEODSOLVE'  # PYCHOK used!


class GeodSolve12Tuple(_GTuple):
    '''12-Tuple C{(lat1, lon1, azi1, lat2, lon2, azi2, s12, a12, m12, M12, M21, S12)} with
       angles C{lat1}, C{lon1}, C{azi1}, C{lat2}, C{lon2} and C{azi2} and arc C{a12} all in
       C{degrees}, initial C{azi1} and final C{azi2} forward azimuths, distance C{s12} and
       reduced length C{m12} in C{meter}, area C{S12} in C{meter} I{squared}  and geodesic
       scale factors C{M12} and C{M21}, both C{scalar}, see U{GeodSolve
       <https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>}.
    '''
    # from GeodSolve --help option -f ... lat1 lon1 azi1 lat2 lon2 azi2 s12 a12 m12 M12 M21 S12
    _Names_ = (_lat1_, _lon1_, _azi1_, _lat2_, _lon2_, _azi2_, _s12_, _a12_, _m12_, _M12_, _M21_, _S12_)
    _Units_ = (_Lat,   _Lon,   _Azi,   _Lat,   _Lon,   _Azi,   _M,    _Deg,  _Pass, _Pass, _Pass, _M2)


class _GeodesicSolveBase(_SolveBase):
    '''(INTERNAL) Base class for L{GeodesicSolve} and L{GeodesicLineSolve}.
    '''
    _Error         =  GeodesicError
    _Names_Direct  = \
    _Names_Inverse =  GeodSolve12Tuple._Names_
    _Solve_name    = 'GeodSolve'
    _Solve_path    = _getenv(_PYGEODESY_GEODSOLVE_, _PYGEODESY_GEODSOLVE_)

    @Property_RO
    def _b_option(self):
        return ('-b',) if self.reverse2 else ()

    @Property_RO
    def _cmdBasic(self):
        '''(INTERNAL) Get the basic C{GeodSolve} cmd (C{tuple}).
        '''
        return (self.GeodSolve,) + self._b_option \
                                 + self._e_option \
                                 + self._E_option \
                                 + self._p_option \
                                 + self._u_option + ('-f',)

    @Property_RO
    def _E_option(self):
        return ('-E',) if self.Exact else ()

    @Property
    def GeodSolve(self):
        '''Get the U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>}
           executable (C{filename}).
        '''
        return self._Solve_path

    @GeodSolve.setter  # PYCHOK setter!
    def GeodSolve(self, path):
        '''Set the U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>}
           executable (C{filename}), the (fully qualified) path to the C{GeodSolve} executable.

           @raise GeodesicError: Invalid B{C{path}}, B{C{path}} doesn't exist or
                                 isn't the C{GeodSolve} executable.
        '''
        self._setSolve(path)

    def toStr(self, **prec_sep):  # PYCHOK signature
        '''Return this C{GeodesicSolve} as string.

           @kwarg prec_sep: Keyword argumens C{B{prec}=6} and C{B{sep}=', '}
                      for the C{float} C{prec}ision, number of decimal digits
                      (0..9) and the C{sep}arator string to join.  Trailing
                      zero decimals are stripped for B{C{prec}} values of
                      1 and above, but kept for negative B{C{prec}} values.

           @return: GeodesicSolve items (C{str}).
        '''
        return _SolveBase._toStr(self, GeodSolve=self.GeodSolve, **prec_sep)

    @Property_RO
    def _u_option(self):
        return ('-u',) if self.unroll else ()


class GeodesicSolve(_GeodesicSolveBase):
    '''Wrapper to invoke I{Karney}'s U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>}
       as an C{Exact} version of I{Karney}'s Python class U{Geodesic<https://GeographicLib.SourceForge.io/C++/doc/
       python/code.html#geographiclib.geodesic.Geodesic>}.

       @note: Use property C{GeodSolve} or env variable C{PYGEODESY_GEODSOLVE} to specify the (fully
              qualified) path to the C{GeodSolve} executable.

       @note: This C{geodesic} is intended I{for testing purposes only}, it invokes the C{GeodSolve}
              executable for I{every} method call.
    '''

    def Area(self, polyline=False, name=NN):
        '''Set up an L{GeodesicAreaExact} to compute area and
           perimeter of a polygon.

           @kwarg polyline: If C{True} perimeter only, otherwise
                            area and perimeter (C{bool}).
           @kwarg name: Optional name (C{str}).

           @return: A L{GeodesicAreaExact} instance.

           @note: The B{C{debug}} setting is passed as C{verbose}
                  to the returned L{GeodesicAreaExact} instance.
        '''
        gaX = _MODS.geodesicx.GeodesicAreaExact(self, polyline=polyline,
                                                      name=name or self.name)
        if self.verbose or self.debug:  # PYCHOK no cover
            gaX.verbose = True
        return gaX

    Polygon = Area  # for C{geographiclib} compatibility

    def Direct3(self, lat1, lon1, azi1, s12):  # PYCHOK outmask
        '''Return the destination lat, lon and reverse azimuth
           (final bearing) in C{degrees}.

           @return: L{Destination3Tuple}C{(lat, lon, final)}.
        '''
        r = self._GDictDirect(lat1, lon1, azi1, False, s12, floats=False)
        return Destination3Tuple(float(r.lat2), float(r.lon2), wrap360(r.azi2),
                                 iteration=r._iteration)

    def Inverse3(self, lat1, lon1, lat2, lon2):  # PYCHOK outmask
        '''Return the distance in C{meter} and the forward and
           reverse azimuths (initial and final bearing) in C{degrees}.

           @return: L{Distance3Tuple}C{(distance, initial, final)}.
        '''
        r = self._GDictInverse(lat1, lon1, lat2, lon2, floats=False)
        return Distance3Tuple(float(r.s12), wrap360(r.azi1), wrap360(r.azi2),
                              iteration=r._iteration)

    def Line(self, lat1, lon1, azi1, caps=Caps.ALL):
        '''Set up an L{GeodesicLineSolve} to compute several points
           on a single geodesic.

           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi1: Azimuth at the first point (compass C{degrees}).
           @kwarg caps: Bit-or'ed combination of L{Caps} values specifying
                        the capabilities the L{GeodesicLineSolve} instance
                        should possess, always C{Caps.ALL}.

           @return: A L{GeodesicLineSolve} instance.

           @note: If the point is at a pole, the azimuth is defined by keeping
                  B{C{lon1}} fixed, writing C{B{lat1} = ±(90 − ε)}, and taking
                  the limit C{ε → 0+}.

           @see: C++ U{GeodesicExact.Line
                 <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1GeodesicExact.html>}
                 and Python U{Geodesic.Line<https://GeographicLib.SourceForge.io/C++/doc/python/code.html>}.
        '''
        return GeodesicLineSolve(self, lat1, lon1, azi1, caps=caps, name=self.name)


class GeodesicLineSolve(_GeodesicSolveBase, _LineSolveBase):
    '''Wrapper to invoke I{Karney}'s U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>}
       as an C{Exact} version of I{Karney}'s Python class U{GeodesicLine<https://GeographicLib.SourceForge.io/C++/doc/
       python/code.html#geographiclib.geodesicline.GeodesicLine>}.

       @note: Use property C{GeodSolve} or env variable C{PYGEODESY_GEODSOLVE} to specify the (fully
              qualified) path to the C{GeodSolve} executable.

       @note: This C{geodesic} is intended I{for testing purposes only}, it invokes the C{GeodSolve}
              executable for I{every} method call.
    '''

    def __init__(self, geodesic, lat1, lon1, azi1, caps=Caps.ALL, name=NN):
        '''New L{GeodesicLineSolve} instance, allowing points to be found along
           a geodesic starting at C{(B{lat1}, B{lon1})} with azimuth B{C{azi1}}.

           @arg geodesic: The geodesic to use (L{GeodesicSolve}).
           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi1: Azimuth at the first points (compass C{degrees}).
           @kwarg caps: Bit-or'ed combination of L{Caps} values specifying
                        the capabilities the L{GeodesicLineSolve} instance
                        should possess, always C{Caps.ALL}.  Use C{Caps.LINE_OFF}
                        if updates to the B{C{geodesic}} should I{not} be
                        reflected in this L{GeodesicLineSolve} instance.
           @kwarg name: Optional name (C{str}).

           @raise GeodesicError: Invalid path for the C{GeodSolve} executable or
                                 or isn't the C{GeodSolve} executable, see
                                 property C{geodesic.GeodSolve}.

           @raise TypeError: Invalid B{C{geodesic}}.
        '''
        _xinstanceof(GeodesicSolve, geodesic=geodesic)
        if (caps & Caps.LINE_OFF):  # copy to avoid updates
            geodesic = geodesic.copy(deep=False, name=NN(_UNDER_, geodesic.name))
        _LineSolveBase.__init__(self, geodesic, lat1, lon1, caps, name, azi1=azi1)
        try:
            self.GeodSolve = geodesic.GeodSolve  # geodesic or copy of geodesic
        except GeodesicError:
            pass

    def ArcPosition(self, a12, *unused):
        '''Find the position on the line given B{C{a12}}.

           @arg a12: Spherical arc length from the first point to the
                     second point (C{degrees}).

           @return: A C{GDict} with 12 items C{lat1, lon1, azi1, lat2, lon2,
                    azi2, m12, a12, s12, M12, M21, S12}.
        '''
        return self._GDictInvoke(self._cmdArc, True, self._Names_Direct, a12)

    @Property_RO
    def azi1(self):
        '''Get the azimuth at the first point (compass C{degrees}).
        '''
        return self._lla1.azi1

    azi12 = azi1  # like RhumbLineSolve

    @Property_RO
    def azi1_sincos2(self):
        '''Get the sine and cosine of the first point's azimuth (2-tuple C{(sin, cos)}).
        '''
        return _sincos2d(self.azi1)

    azi12_sincos2 = azi1_sincos2

    @Property_RO
    def _cmdArc(self):
        '''(INTERNAL) Get the C{GeodSolve} I{-a -L} cmd (C{tuple}).
        '''
        return self._cmdDistance + ('-a',)

    def Position(self, s12, *unused):
        '''Find the position on the line given B{C{s12}}.

           @arg s12: Distance from the first point to the second (C{meter}).

           @return: A C{GDict} with 12 items C{lat1, lon1, azi1, lat2, lon2,
                    azi2, m12, a12, s12, M12, M21, S12}, possibly C{a12=NAN}.
        '''
        return self._GDictInvoke(self._cmdDistance, True, self._Names_Direct, s12)

    def toStr(self, **prec_sep):  # PYCHOK signature
        '''Return this C{GeodesicLineSolve} as string.

           @kwarg prec_sep: Keyword argumens C{B{prec}=6} and C{B{sep}=', '}
                      for the C{float} C{prec}ision, number of decimal digits
                      (0..9) and the C{sep}arator string to join.  Trailing
                      zero decimals are stripped for B{C{prec}} values of
                      1 and above, but kept for negative B{C{prec}} values.

           @return: GeodesicLineSolve items (C{str}).
        '''
        return _LineSolveBase._toStr(self, azi1=self.azi1, geodesic=self._solve,
                                           GeodSolve=self.GeodSolve, **prec_sep)


__all__ += _ALL_DOCS(_GeodesicSolveBase)

if __name__ == '__main__':

    from sys import argv

    gS = GeodesicSolve(name='Test')
    gS.verbose = '--verbose' in argv  # or '-v' in argv

    if gS.GeodSolve in (_PYGEODESY_GEODSOLVE_, None):  # not set
        gS.GeodSolve = '/opt/local/bin/GeodSolve'  # '/opt/local/Cellar/geographiclib/1.51/bin/GeodSolve'  # HomeBrew
    printf('version: %s', gS.version)

    r = gS.Direct(40.6, -73.8, 51, 5.5e6)
    printf('Direct: %r', r, nl=1)
    printf('Direct3: %r', gS.Direct3(40.6, -73.8, 51, 5.5e6))

    printf('Inverse: %r',  gS.Inverse( 40.6, -73.8, 51.6, -0.5), nl=1)
    printf('Inverse1: %r', gS.Inverse1(40.6, -73.8, 51.6, -0.5))
    printf('Inverse3: %r', gS.Inverse3(40.6, -73.8, 51.6, -0.5))

    glS = GeodesicLineSolve(gS, 40.6, -73.8, 51, name='LineTest')
    p = glS.Position(5.5e6)
    printf('Position:    %s  %r', p == r, p, nl=1)
    p = glS.ArcPosition(49.475527)
    printf('ArcPosition: %s %r', p == r, p)

# % python3 -m pygeodesy.geodsolve

# version: /opt/local/bin/GeodSolve: GeographicLib version 1.51

# version: /opt/local/bin/GeodSolve: GeographicLib version 1.51

# Direct: GDict(M12=0.650911, M21=0.651229, S12=39735075134877.09375, a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, m12=4844148.703101, s12=5500000.0)
# Direct3: Destination3Tuple(lat=51.884565, lon=-1.141173, final=107.189397)

# Inverse: GDict(M12=0.64473, M21=0.645046, S12=40041368848742.53125, a12=49.94131, azi1=51.198883, azi2=107.821777, lat1=40.6, lat2=51.6, lon1=-73.8, lon2=-0.5, m12=4877684.602706, s12=5551759.400319)
# Inverse1: 49.94131021789904
# Inverse3: Distance3Tuple(distance=5551759.400319, initial=51.198883, final=107.821777)

# Position:    True  GDict(M12=0.650911, M21=0.651229, S12=39735075134877.09375, a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, m12=4844148.703101, s12=5500000.0)
# ArcPosition: False GDict(M12=0.650911, M21=0.651229, S12=39735074737272.734375, a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141174, m12=4844148.669561, s12=5499999.948497)


# % python3 -m pygeodesy.geodsolve --verbose

# GeodesicSolve 'Test' 1: /opt/local/bin/GeodSolve --version (invoke)
# GeodesicSolve 'Test' 1: /opt/local/bin/GeodSolve: GeographicLib version 1.51 (0)
# version: /opt/local/bin/GeodSolve: GeographicLib version 1.51
# GeodesicSolve 'Test' 2: /opt/local/bin/GeodSolve -E -p 10 -f \ 40.600000000000001 -73.799999999999997 51.0 5500000.0 (Direct)
# GeodesicSolve 'Test' 2: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.0, lat2=51.884564505606761, lon2=-1.141172861200829, azi2=107.189397162605886, s12=5500000.0, a12=49.475527463251467, m12=4844148.703101486, M12=0.65091056699808603, M21=0.65122865892196558, S12=39735075134877.094 (0)

# Direct: GDict(M12=0.650911, M21=0.651229, S12=39735075134877.09375, a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, m12=4844148.703101, s12=5500000.0)
# GeodesicSolve 'Test' 3: /opt/local/bin/GeodSolve -E -p 10 -f \ 40.600000000000001 -73.799999999999997 51.0 5500000.0 (Direct3)
# GeodesicSolve 'Test' 3: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.0, lat2=51.884564505606761, lon2=-1.141172861200829, azi2=107.189397162605886, s12=5500000.0, a12=49.475527463251467, m12=4844148.703101486, M12=0.65091056699808603, M21=0.65122865892196558, S12=39735075134877.094 (0)
# Direct3: Destination3Tuple(lat=51.884565, lon=-1.141173, final=107.189397)
# GeodesicSolve 'Test' 4: /opt/local/bin/GeodSolve -E -p 10 -f -i \ 40.600000000000001 -73.799999999999997 51.600000000000001 -0.5 (Inverse)
# GeodesicSolve 'Test' 4: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.198882845579824, lat2=51.600000000000001, lon2=-0.5, azi2=107.821776735514248, s12=5551759.4003186841, a12=49.941310217899037, m12=4877684.6027061976, M12=0.64472969205948238, M21=0.64504567852134398, S12=40041368848742.531 (0)

# Inverse: GDict(M12=0.64473, M21=0.645046, S12=40041368848742.53125, a12=49.94131, azi1=51.198883, azi2=107.821777, lat1=40.6, lat2=51.6, lon1=-73.8, lon2=-0.5, m12=4877684.602706, s12=5551759.400319)
# GeodesicSolve 'Test' 5: /opt/local/bin/GeodSolve -E -p 10 -f -i \ 40.600000000000001 -73.799999999999997 51.600000000000001 -0.5 (Inverse1)
# GeodesicSolve 'Test' 5: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.198882845579824, lat2=51.600000000000001, lon2=-0.5, azi2=107.821776735514248, s12=5551759.4003186841, a12=49.941310217899037, m12=4877684.6027061976, M12=0.64472969205948238, M21=0.64504567852134398, S12=40041368848742.531 (0)
# Inverse1: 49.94131021789904
# GeodesicSolve 'Test' 6: /opt/local/bin/GeodSolve -E -p 10 -f -i \ 40.600000000000001 -73.799999999999997 51.600000000000001 -0.5 (Inverse3)
# GeodesicSolve 'Test' 6: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.198882845579824, lat2=51.600000000000001, lon2=-0.5, azi2=107.821776735514248, s12=5551759.4003186841, a12=49.941310217899037, m12=4877684.6027061976, M12=0.64472969205948238, M21=0.64504567852134398, S12=40041368848742.531 (0)
# Inverse3: Distance3Tuple(distance=5551759.400319, initial=51.198883, final=107.821777)

# Position:    True  GDict(M12=0.650911, M21=0.651229, S12=39735075134877.09375, a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, m12=4844148.703101, s12=5500000.0)
# ArcPosition: False GDict(M12=0.650911, M21=0.651229, S12=39735074737272.734375, a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141174, m12=4844148.669561, s12=5499999.948497)


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
