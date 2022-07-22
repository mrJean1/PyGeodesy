
# -*- coding: utf-8 -*-

u'''Wrapper to invoke I{Karney}'s U{RhumbSolve
<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>} utility
as an (exact) rhumb or rhumb line, but intended I{for testing purposes only}.

Set env variable C{PYGEODESY_RHUMBSOLVE} to the (fully qualified) path
of the C{RhumbSolve} executable.
'''

# from pygeodesy.basics import _xinstanceof  # from .karney
# from pygeodesy.interns import NN  # from .karney
from pygeodesy.interns import NAN, NN, _azi12_, _lat2_, _lon2_, _s12_, _S12_, \
                             _UNDER_, _90_0, _180_0, _N_180_0, _not_  # PYCHOK used!
from pygeodesy.karney import GDict, _norm180, _sincos2d, _xinstanceof, wrap360
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS, _getenv, printf
from pygeodesy.namedTuples import Destination3Tuple, Distance3Tuple
from pygeodesy.props import deprecated_method, Property, Property_RO
from pygeodesy.rhumbx import Caps, RhumbError, Rhumb8Tuple  # PYCHOK used!
from pygeodesy.solveBase import _LineSolveBase, _SolveBase
# from pygeodesy.utily import wrap360  # from .karney

__all__ = _ALL_LAZY.rhumbsolve
__version__ = '22.07.09'

_PYGEODESY_RHUMBSOLVE_ = 'PYGEODESY_RHUMBSOLVE'  # PYCHOK used!


class _RhumbSolveBase(_SolveBase):
    '''(INTERNAL) Base class for L{RhumbSolve} and L{RhumbLineSolve}.
    '''
    _Error         =  RhumbError
    _Names_Direct  = _lat2_, _lon2_, _S12_
    _Names_Inverse = _azi12_, _s12_, _S12_
    _Solve_name    = 'RhumbSolve'
    _Solve_path    = _getenv(_PYGEODESY_RHUMBSOLVE_, _PYGEODESY_RHUMBSOLVE_)

    @Property_RO
    def _cmdBasic(self):
        '''(INTERNAL) Get the basic C{RhumbSolve} cmd (C{tuple}).
        '''
        return (self.RhumbSolve,) + self._e_option \
                                  + self._p_option \
                                  + self._s_option

    @Property
    def RhumbSolve(self):
        '''Get the U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>}
           executable (C{filename}).
        '''
        return self._Solve_path

    @RhumbSolve.setter  # PYCHOK setter!
    def RhumbSolve(self, path):
        '''Set the U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>}
           executable (C{filename}), the (fully qualified) path to the C{RhumbSolve} executable.

           @raise RhumbError: Invalid B{C{path}}, B{C{path}} doesn't exist or isn't
                              the C{RhumbSolve} executable.
        '''
        self._setSolve(path)

    @Property_RO
    def _s_option(self):  # == not -E for GeodSolve
        return () if self.Exact else ('-s', )

    def toStr(self, **prec_sep):  # PYCHOK signature
        '''Return this C{RhumbSolve} as string.

           @kwarg prec_sep: Keyword argumens C{B{prec}=6} and C{B{sep}=', '}
                      for the C{float} C{prec}ision, number of decimal digits
                      (0..9) and the C{sep}arator string to join.  Trailing
                      zero decimals are stripped for B{C{prec}} values of
                      1 and above, but kept for negative B{C{prec}} values.

           @return: RhumbSolve items (C{str}).
        '''
        return self._toStr(RhumbSolve=self.RhumbSolve, **prec_sep)

#   @Property_RO
#   def _u_option(self):
#       return '-u' if self.unroll else ()


class RhumbSolve(_RhumbSolveBase):
    '''Wrapper to invoke I{Karney}'s U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>},
       the C{Exact} version of class L{pygeodesy.Rhumb}.

       @note: Use property C{RhumbSolve} or env variable C{PYGEODESY_RHUMBSOLVE} to specify the (fully
              qualified) path to the C{RhumbSolve} executable.

       @note: This C{rhumb} is intended I{for testing purposes only}, it invokes the C{RhumbSolve}
              executable for I{every} method call.
    '''
#   def Area(self, polyline=False, name=NN):
#       '''Set up an L{RhumbArea} to compute area and
#          perimeter of a polygon.
#
#          @kwarg polyline: If C{True} perimeter only, otherwise
#                           area and perimeter (C{bool}).
#          @kwarg name: Optional name (C{str}).
#
#          @return: A L{RhumbArea} instance.
#
#          @note: The B{C{debug}} setting is passed as C{verbose}
#                 to the returned L{RhumbAreaExact} instance.
#       '''
#       rA = _MODS.rhumbx.RhumbArea(self, polyline=polyline,
#                                         name=name or self.name)
#       if self.verbose or self.debug:  # PYCHOK no cover
#           rA.verbose = True
#       return rA

#   Polygon = Area  # for C{geographiclib} compatibility

    def _azimuth_reverse(self, azimuth):
        '''(INTERNAL) Reverse final azimuth C{azimuth}.
        '''
        z = _norm180(float(azimuth))
        if self.reverse2:  # like .utils.atan2d
            z += _180_0 if z < 0 else _N_180_0
        return z

    def Direct3(self, lat1, lon1, azi1, s12):  # PYCHOK outmask
        '''Return the destination lat, lon and reverse azimuth
           (final bearing) in C{degrees}.

           @return: L{Destination3Tuple}C{(lat, lon, final)}.
        '''
        r = self._GDictDirect(lat1, lon1, azi1, False, s12, floats=False)
        z = self._azimuth_reverse(r.azi12)
        return Destination3Tuple(float(r.lat2), float(r.lon2), wrap360(z),
                                 iteration=r._iteration)

    def _GDictDirect(self, lat, lon, azi1, arcmode, s12_a12, *unused, **floats):
        '''(INTERNAL) Get C{_GenDirect}-like result as an 8-item C{GDict}.
        '''
        d = _RhumbSolveBase._GDictDirect(self, lat, lon, azi1, arcmode, s12_a12, **floats)
        r =  GDict(lat1=lat, lon1=lon, azi12=azi1, s12=s12_a12)  # a12=s12_a12 / self.ellipsoid._L_90
        r.update(d)
        return r

    def _GDictInverse(self, lat1, lon1, lat2, lon2, *unused, **floats):
        '''(INTERNAL) Get C{_GenInverse}-like result as an 8-item C{GDict}, but
           I{without} C{_SALPs_CALPs_}.
        '''
        i = _RhumbSolveBase._GDictInverse(self, lat1, lon1, lat2, lon2, **floats)
        a =  float(i.s12) / self.ellipsoid._L_90  # for .Inverse1
        r =  GDict(lat1=lat1, lon1=lon1, lat2=lat2, lon2=lon2, a12=a)
        r.update(i)
        return r

    def Inverse3(self, lat1, lon1, lat2, lon2):  # PYCHOK outmask
        '''Return the distance in C{meter} and the forward and
           reverse azimuths (initial and final bearing) in C{degrees}.

           @return: L{Distance3Tuple}C{(distance, initial, final)}.
        '''
        r = self._GDictInverse(lat1, lon1, lat2, lon2, floats=False)
        z = self._azimuth_reverse(r.azi12)
        return Distance3Tuple(float(r.s12), wrap360(r.azi12), wrap360(z),
                              iteration=r._iteration)

    def Line(self, lat1, lon1, azi1, caps=Caps.ALL):
        '''Set up an L{RhumbLineSolve} to compute several points
           on a single rhumb line.

           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi1: Azimuth at the first point (compass C{degrees}).
           @kwarg caps: Bit-or'ed combination of L{Caps} values specifying
                        the capabilities the L{RhumbLineSolve} instance
                        should possess, always C{Caps.ALL}.

           @return: A L{RhumbLineSolve} instance.

           @note: If the point is at a pole, the azimuth is defined by keeping
                  B{C{lon1}} fixed, writing C{B{lat1} = ±(90 − ε)}, and taking
                  the limit C{ε → 0+}.

           @see: C++ U{RhumbExact.Line
                 <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1RhumbExact.html>}
                 and Python U{Rhumb.Line<https://GeographicLib.SourceForge.io/C++/doc/python/code.html>}.
        '''
        return RhumbLineSolve(self, lat1, lon1, azi1, caps=caps, name=self.name)


class RhumbLineSolve(_RhumbSolveBase, _LineSolveBase):
    '''Wrapper to invoke I{Karney}'s U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>},
       the C{Exact} version of class L{pygeodesy.RhumbLine}.

       @note: Use property C{RhumbSolve} or env variable C{PYGEODESY_RHUMBSOLVE} to specify the (fully
              qualified) path to the C{RhumbSolve} executable.

       @note: This C{rhumb line} is intended I{for testing purposes only}, it invokes the C{RhumbSolve}
              executable for I{every} method call.
    '''
    def __init__(self, rhumb, lat1, lon1, azi12, caps=Caps.ALL, name=NN):
        '''New L{RhumbLineSolve} instance, allowing points to be found along
           a rhumb starting at C{(B{lat1}, B{lon1})} with azimuth B{C{azi12}}.

           @arg rhumb: The rhumb to use (L{RhumbSolve}).
           @arg lat1: Latitude of the first point (C{degrees90}).
           @arg lon1: Longitude of the first point (C{degrees180}).
           @arg azi12: Azimuth of the rhumb line (compass C{degrees180}).
           @kwarg caps: Bit-or'ed combination of L{Caps} values specifying
                        the capabilities the L{RhumbLineSolve} instance should
                        possess, always C{Caps.ALL}.  Use C{Caps.LINE_OFF}
                        if updates to the B{C{rhumb}} should I{not} be
                        reflected in this L{RhumbLineSolve} instance.

           @kwarg name: Optional name (C{str}).

           @raise RhumbError: Invalid path for C{RhumbSolve} executable or
                              isn't the C{RhumbSolve} executable, see
                              property C{B{rhumb}.RhumbSolve}.

           @raise TypeError: Invalid B{C{rhumb}}.
        '''
        _xinstanceof(RhumbSolve, rhumb=rhumb)
        if (caps & Caps.LINE_OFF):  # copy to avoid updates
            rhumb = rhumb.copy(deep=False, name=NN(_UNDER_, rhumb.name))
        _LineSolveBase.__init__(self, rhumb, lat1, lon1, caps, name, azi12=azi12)
        try:
            self.RhumbSolve = rhumb.RhumbSolve  # rhumb or copy of rhumb
        except RhumbError:
            pass

#   def ArcPosition(self, a12, *unused):
#       '''Find the position on the line given B{C{a12}}.
#
#          @arg a12: Spherical arc length from the first point to the
#                    second point (C{degrees}).
#
#          @return: A C{dict} with 8 items C{lat1, lon1, lat2, lon2,
#                   azi12, a12, s12, S12}.
#       '''
#       s = a12 * self.ellipsoid._L_90
#       a = self._GDictInvoke(self._cmdArc, True, self._Names_Direct, s)
#       r = GDict(a12=a12, s12=s, **self._lla1)
#       r.updated(a)
#       return r

    @Property_RO
    def azi12(self):
        '''Get this rhumb line's azimuth (compass C{degrees}).
        '''
        return self._lla1.azi12

    azi1 = azi12  # like GeodesicLineSolve

    @Property_RO
    def azi12_sincos2(self):  # PYCHOK no cover
        '''Get the sine and cosine of this rhumb line's azimuth (2-tuple C{(sin, cos)}).
        '''
        return _sincos2d(self.azi12)

    azi1_sincos2 = azi12_sincos2

#   @Property_RO
#   def _cmdArc(self):
#       '''(INTERNAL) Get the C{RhumbSolve} I{-a -L} cmd (C{tuple}).
#       '''
#       return self._cmdDistance + ('-a',)

    def Position(self, s12, *unused):
        '''Find the position on the line given B{C{s12}}.

           @arg s12: Distance from the first point to the second (C{meter}).

           @return: A L{GDict} with 7 items C{lat1, lon1, lat2, lon2,
                    azi12, s12, S12}.
        '''
        d = self._GDictInvoke(self._cmdDistance, True, self._Names_Direct, s12)
        r = GDict(s12=s12, **self._lla1)  # a12=s12 / self.ellipsoid._L_90
        r.update(d)
        return r

    def toStr(self, **prec_sep):  # PYCHOK signature
        '''Return this C{RhumbLineSolve} as string.

           @kwarg prec_sep: Keyword argumens C{B{prec}=6} and C{B{sep}=', '}
                      for the C{float} C{prec}ision, number of decimal digits
                      (0..9) and the C{sep}arator string to join.  Trailing
                      zero decimals are stripped for B{C{prec}} values of
                      1 and above, but kept for negative B{C{prec}} values.

           @return: RhumbLineSolve items (C{str}).
        '''
        return _LineSolveBase._toStr(self, azi12=self.azi12, rhumb=self._solve,
                                           RhumbSolve=self.RhumbSolve, **prec_sep)


class RhumbSolve7Tuple(Rhumb8Tuple):
    '''7-Tuple C{(lat1, lon1, lat2, lon2, azi12, s12, S12)} with lat- C{lat1},
       C{lat2} and longitudes C{lon1}, C{lon2} of both points, the azimuth of the
       rhumb line C{azi12}, the distance C{s12} and the area C{S12} under the
       rhumb line between both points.
    '''
    assert Rhumb8Tuple._Names_.index(_MODS.interns._a12_) == 7
    _Names_ = Rhumb8Tuple._Names_[:7]  # drop a12
    _Units_ = Rhumb8Tuple._Units_[:7]

    @deprecated_method
    def _to7Tuple(self):  # PYCHOK no cover
        '''DEPRECATED, I{don't use!}
        '''
        return _MODS.deprecated.Rhumb7Tuple(self[:7])


__all__ += _ALL_DOCS(_RhumbSolveBase)

if __name__ == '__main__':

    from sys import argv

    rS = RhumbSolve(name='Test')
    rS.verbose = '--verbose' in argv  # or '-v' in argv

    if rS.RhumbSolve in (_PYGEODESY_RHUMBSOLVE_, None):  # not set
        rS.RhumbSolve = '/opt/local/bin/RhumbSolve'  # '/opt/local/Cellar/geographiclib/1.51/bin/RhumbSolve'  # HomeBrew
    printf('version: %s', rS.version)

    r = rS.Direct(40.6, -73.8, 51, 5.5e6)
    printf('Direct: %r', r, nl=1)
    printf('Direct3: %r', rS.Direct3(40.6, -73.8, 51, 5.5e6))

    printf('Inverse: %r',  rS.Inverse( 40.6, -73.8, 51.6, -0.5), nl=1)
    printf('Inverse1: %r', rS.Inverse1(40.6, -73.8, 51.6, -0.5))
    printf('Inverse3: %r', rS.Inverse3(40.6, -73.8, 51.6, -0.5))

    printf('Inverse: %r',  rS.Inverse( 40.6, -73.8, 35.8, 140.3), nl=1)
    printf('Inverse1: %r', rS.Inverse1(40.6, -73.8, 35.8, 140.3))
    printf('Inverse3: %r', rS.Inverse3(40.6, -73.8, 35.8, 140.3))

    rlS = RhumbLineSolve(rS, 40.6, -73.8, 51, name='LineTest')
    p = rlS.Position(5.5e6)
    printf('Position:    %s  %r', p == r, p, nl=1)
#   p = rlS.ArcPosition(49.475527)
#   printf('ArcPosition: %s %r', p == r, p)

# % python3 -m pygeodesy.rhumbsolve

# version: /opt/local/bin/RhumbSolve: GeographicLib version 1.51
#
# Direct: GDict(S12=44095641862956.148438, azi12=51, lat1=40.6, lat2=71.6889, lon1=-73.8, lon2=0.25552, s12=5500000.0)
# Direct3: Destination3Tuple(lat=71.6889, lon=0.25552, final=51.0)
#
# Inverse: GDict(S12=37395209100030.367188, a12=51.929543, azi12=77.76839, lat1=40.6, lat2=51.6, lon1=-73.8, lon2=-0.5, s12=5771083.383328)
# Inverse1: 51.92954250756195
# Inverse3: Distance3Tuple(distance=5771083.383328, initial=77.76839, final=77.76839)
#
# Inverse: GDict(S12=-63760642939072.492188, a12=115.02062, azi12=-92.388888, lat1=40.6, lat2=35.8, lon1=-73.8, lon2=140.3, s12=12782581.067684)
# Inverse1: 115.02061966879258
# Inverse3: Distance3Tuple(distance=12782581.067684, initial=267.611112, final=267.611112)
#
# Position:    True  GDict(S12=44095641862956.148438, azi12=51, lat1=40.6, lat2=71.6889, lon1=-73.8, lon2=0.25552, s12=5500000.0)


# % python3 -m pygeodesy.rhumbsolve --verbose

# RhumbSolve 'Test' 1: /opt/local/bin/RhumbSolve --version (invoke)
# RhumbSolve 'Test' 1: /opt/local/bin/RhumbSolve: GeographicLib version 1.51 (0)
# version: /opt/local/bin/RhumbSolve: GeographicLib version 1.51
# RhumbSolve 'Test' 2: /opt/local/bin/RhumbSolve -p 10 \ 40.600000000000001 -73.799999999999997 51.0 5500000.0 (Direct)
# RhumbSolve 'Test' 2: lat2=71.688899882813047, lon2=0.255519824423445, S12=44095641862956.148 (0)

# Direct: GDict(S12=44095641862956.148438, azi12=51, lat1=40.6, lat2=71.6889, lon1=-73.8, lon2=0.25552, s12=5500000.0)
# RhumbSolve 'Test' 3: /opt/local/bin/RhumbSolve -p 10 \ 40.600000000000001 -73.799999999999997 51.0 5500000.0 (Direct3)
# RhumbSolve 'Test' 3: lat2=71.688899882813047, lon2=0.255519824423445, S12=44095641862956.148 (0)
# Direct3: Destination3Tuple(lat=71.6889, lon=0.25552, final=51.0)
# RhumbSolve 'Test' 4: /opt/local/bin/RhumbSolve -p 10 -i \ 40.600000000000001 -73.799999999999997 51.600000000000001 -0.5 (Inverse)
# RhumbSolve 'Test' 4: azi12=77.768389710255661, s12=5771083.3833280317, S12=37395209100030.367 (0)

# Inverse: GDict(S12=37395209100030.367188, a12=51.929543, azi12=77.76839, lat1=40.6, lat2=51.6, lon1=-73.8, lon2=-0.5, s12=5771083.383328)
# RhumbSolve 'Test' 5: /opt/local/bin/RhumbSolve -p 10 -i \ 40.600000000000001 -73.799999999999997 51.600000000000001 -0.5 (Inverse1)
# RhumbSolve 'Test' 5: azi12=77.768389710255661, s12=5771083.3833280317, S12=37395209100030.367 (0)
# Inverse1: 51.92954250756195
# RhumbSolve 'Test' 6: /opt/local/bin/RhumbSolve -p 10 -i \ 40.600000000000001 -73.799999999999997 51.600000000000001 -0.5 (Inverse3)
# RhumbSolve 'Test' 6: azi12=77.768389710255661, s12=5771083.3833280317, S12=37395209100030.367 (0)
# Inverse3: Distance3Tuple(distance=5771083.383328, initial=77.76839, final=77.76839)
# RhumbSolve 'Test' 7: /opt/local/bin/RhumbSolve -p 10 -i \ 40.600000000000001 -73.799999999999997 35.799999999999997 140.300000000000011 (Inverse)
# RhumbSolve 'Test' 7: azi12=-92.388887981699639, s12=12782581.0676841792, S12=-63760642939072.492 (0)

# Inverse: GDict(S12=-63760642939072.492188, a12=115.02062, azi12=-92.388888, lat1=40.6, lat2=35.8, lon1=-73.8, lon2=140.3, s12=12782581.067684)
# RhumbSolve 'Test' 8: /opt/local/bin/RhumbSolve -p 10 -i \ 40.600000000000001 -73.799999999999997 35.799999999999997 140.300000000000011 (Inverse1)
# RhumbSolve 'Test' 8: azi12=-92.388887981699639, s12=12782581.0676841792, S12=-63760642939072.492 (0)
# Inverse1: 115.02061966879258
# RhumbSolve 'Test' 9: /opt/local/bin/RhumbSolve -p 10 -i \ 40.600000000000001 -73.799999999999997 35.799999999999997 140.300000000000011 (Inverse3)
# RhumbSolve 'Test' 9: azi12=-92.388887981699639, s12=12782581.0676841792, S12=-63760642939072.492 (0)
# Inverse3: Distance3Tuple(distance=12782581.067684, initial=267.611112, final=267.611112)

# Position:    True  GDict(S12=44095641862956.148438, azi12=51, lat1=40.6, lat2=71.6889, lon1=-73.8, lon2=0.25552, s12=5500000.0)


# **) MIT License
#
# Copyright (C) 2022-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
