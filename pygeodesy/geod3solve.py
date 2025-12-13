
# -*- coding: utf-8 -*-

u'''Wrapper to invoke I{Karney}'s U{Geod3Solve
<https://GeographicLib.SourceForge.io/C++/doc/Geod3Solve.1.html>} utility
as a C{triaxial} geodesic, but intended I{mainly for testing purposes}.

Set env variable C{PYGEODESY_GEOD3SOLVE} to the (fully qualified) path
of the C{Geod3Solve} executable.
'''

from pygeodesy.angles import Ang, Deg, isAng,  hypot
from pygeodesy.basics import _xinstanceof  # typename
from pygeodesy.constants import _0_0, _0_5, _360_0
from pygeodesy.errors import GeodesicError, _xkwds_get
# from pygeodesy.fmath import hypot  # from .angles
# from pygeodesy.geodesicx import GeodesicAreaExact  # _MODS
from pygeodesy.interns import NN, _a12_, _DMAIN_, _s12_, _UNDER_, _UNUSED_
from pygeodesy.karney import Caps, _GTuple, _Xables,  sincos2d
# from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY  # from .solveBase
from pygeodesy.props import Property, Property_RO, property_RO
from pygeodesy.solveBase import _Solve3Base,  _ALL_DOCS, _ALL_LAZY
from pygeodesy.triaxials.triaxial3 import Triaxial3, Triaxial3s
from pygeodesy.units import Degrees, Meter
# from pygeodesy.utily import sincos2d  # from .karney

__all__ = _ALL_LAZY.geod3solve
__version__ = '25.12.12'

_Triaxial3_WGS84 = Triaxial3s.WGS84_3r  # a=6378172, b=6378102, c=6356752


class Geodesic3Error(GeodesicError):
    '''Error raised for issues in L{geod3solve<pygeodesy.geod3solve>}.
    '''
    pass


class Geod3Solve8Tuple(_GTuple):
    '''8-Tuple C{(bet1, omg1, alp1, bet2, omg2, alp2, s12, a12)} with C{ellipsoidal}
       latitudes C{bet1} and C{bet2}, C{ellipsoidal} longitudes C{omg1} and C{omg2},
       forward azimuths C{alp1} and C{alp2} in bearings from North, distanc C{s12} in
       C{meter}, conventionally and I{approximate} arc length {a12} in degrees, see
       U{Geod3Solve<https://GeographicLib.SourceForge.io/C++/doc/Geod3Solve.1.html>}.
    '''
    # from Geod3Solve --help option -f ... bet1 omg1 alp1 bet2 omg2 alp2 s12
    _Names_ = ('bet1', 'omg1', 'alp1', 'bet2', 'omg2', 'alp2', _s12_, _a12_)
    _Units_ = ( Deg,    Deg,    Deg,    Deg,    Deg,    Deg,    Meter, Deg)


class _Geodesic3SolveBase(_Solve3Base):
    '''(INTERNAL) Base class for L{Geodesic3Solve} and L{GeodesicLine3Solve}.
    '''
    _Error         =  Geodesic3Error
    _Names_Direct  = _Names_Distance = \
    _Names_Inverse =  Geod3Solve8Tuple._Names_[:7]  # 7 only, always
    _triaxial3     = _Triaxial3_WGS84
    _Xable_name    = _Xables.Geod3Solve.__name__  # typename
    _Xable_path    = _Xables.Geod3Solve()

    @Property_RO
    def a(self):
        '''Get the triaxial's I{major} radius, semi-axis (C{meter}).
        '''
        return self._triaxial3.a

    @Property_RO
    def b(self):
        '''Get the triaxial's I{middle} radius, semi-axis (C{meter}).
        '''
        return self._triaxial3.b

    @Property_RO
    def c(self):
        '''Get the triaxial's I{minor} radius, semi-axis (C{meter}).
        '''
        return self._triaxial3.c

    @Property_RO
    def _cmdBasic(self):
        '''(INTERNAL) Get the basic C{GeodSolve} cmd (C{tuple}).
        '''
        return (self.Geod3Solve, '-f') + (self._t_option +
                                          self._p_option +
                                          self._u_option)

    @property_RO
    def _e_option(self):
        return ()

    _mpd = _E_option = _e_option
    flattening = f = ellipsoid = datum = None

    @Property
    def Geod3Solve(self):
        '''Get the U{Geod3Solve<https://GeographicLib.SourceForge.io/C++/doc/Geod3Solve.1.html>}
           executable (C{filename}).
        '''
        return self._Xable_path

    @Geod3Solve.setter  # PYCHOK setter!
    def Geod3Solve(self, path):
        '''Set the U{Geod3Solve<https://GeographicLib.SourceForge.io/C++/doc/Geod3Solve.1.html>}
           executable (C{filename}), the (fully qualified) path to the C{Geod3Solve} executable.

           @raise GeodesicError: Invalid B{C{path}}, B{C{path}} doesn't exist or
                                 isn't the C{Geod3Solve} executable.
        '''
        self._setXable(path)

    @property_RO
    def _t_option(self):
        return ('-t',) + self._toStdin(self._triaxial3._abc3)

    def toStr(self, **prec_sep):  # PYCHOK signature
        '''Return this C{Geodesic3Solve} as string.

           @kwarg prec_sep: See L{toStr<pygeodesy.solveBase._SolveBase.toStr>}.

           @return: Geodesic3Solve items (C{str}).
        '''
        return _Solve3Base.toStr(self, Geod3Solve=self.Geod3Solve, **prec_sep)

    @Property_RO
    def _u_option(self):
        return ('-u',) if self.unroll else ()


class Geodesic3Solve(_Geodesic3SolveBase):
    '''Wrapper to invoke I{Karney}'s U{Geod3Solve<https://GeographicLib.SourceForge.io/C++/doc/Geod3Solve.1.html>}
       as a C{triaxial} version of I{Karney}'s Python class U{Geodesic<https://GeographicLib.SourceForge.io/C++/doc/
       python/code.html#geographiclib.geodesic.Geodesic>}.

       @note: Use property C{GeodSolve} or env variable C{PYGEODESY_GEOD3SOLVE} to specify the (fully
              qualified) path to the C{Geod3Solve} executable.

       @note: This C{geodesic} is intended I{for testing purposes only}, it invokes the C{GeodSolve}
              executable for I{every} method call.
    '''

    def __init__(self, triaxial3=_Triaxial3_WGS84, path=NN, **name):
        '''New C{Solve} instance.

           @arg triaxial3: A triaxial (L{Triaxial3}), default C{Triaxial3s.WGS84_3r}.
           @kwarg path: Optionally, the (fully qualified) path to the C{GeodSolve}
                        or C{RhumbSolve} executable (C{filename}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise TypeError: Invalid B{C{a_ellipsoid}} or B{C{f}}.
        '''
        _xinstanceof(Triaxial3, triaxial3=triaxial3)
        if self._triaxial3 != triaxial3:
            self._triaxial3 = triaxial3
        if name:
            self.name = name
        if path:
            self._setXable(path)

    def _a12d(self, r):
        '''(INTERNAL) Add arc C{a12} in degrees to C{GDict}.
        '''
        a = r.s12 or _0_0
        if a:
            z = _toAzi(r.alp1) + _toAzi(r.alp2)
            s, c = sincos2d(z * _0_5)
            t  = self.triaxial3
            a *= hypot(s / t._ab_elliperim,  # azimuth!
                       c / t._bc_elliperim) * _360_0
        r[_a12_] = a
        return r

    def Direct(self, bet1, omg1, alp1, s12, outmask=_UNUSED_, **unit):  # PYCHOK unused
        '''Return the C{Direct} result at distance C{s12}.
        '''
        bet1, omg1, alp1 = _toDegrees(bet1, omg1, alp1, **unit)
        r = self._GDictDirect(bet1, omg1, alp1, False, s12)
        return self._a12d(r)

    def DirectLine(self, bet1, omg1, alp1, caps=Caps.ALL, **unit_name):
        '''Set up a L{GeodesicLine3Solve} to compute several points
           on a single geodesic.

           @arg bet1: Ellipsoidal Latitude of the first point (C{Ang} or B{C{unit}}).
           @arg omg1: Ellipsoidal Longitude of the first point (C{Ang} or B{C{unit}}).
           @arg alp1: Azimuth at the first point (compass C{degrees}, C{Ang} or C{unit}).
           @kwarg caps: Desired capabilities for the L{GeodesicLine3Solve} instance.
           @kwarg unit_name: Optional C{B{name}=NN} (C{str}) and scalar C{B{unit}=}L{Degrees}
                             (or L{Degrees}).

           @return: A L{GeodesicLine3Solve} instance.

           @see: C++ U{GeodesicExact.Line
                 <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1GeodesicExact.html>}
                 and Python U{Geodesic.Line<https://GeographicLib.SourceForge.io/Python/doc/code.html>}.
        '''
        return GeodesicLine3Solve(self, bet1, omg1, alp1, caps=caps, **unit_name)

    Line = DirectLine  # ArcDirectLine

    def Inverse(self, bet1, omg1, bet2, omg2, outmask=_UNUSED_, **unit):  # PYCHOK unused
        '''Return the C{Inverse} result.
        '''
        r = self._GDictInverse(*_toDegrees(bet1, omg1, bet2, omg2, **unit))
        return self._a12d(r)

    def InverseLine(self, bet1, omg1, bet2, omg2, caps=Caps.ALL, **unit_name):
        '''Set up a L{GeodesicLine3Solve} to compute several points
           on a single geodesic.

           @arg bet1: Latitude of the first point (C{Ang} or B{C{unit}}).
           @arg omg1: Longitude of the first point (C{Ang} or B{C{unit}}).
           @arg bet2: Latitude of the second point  (C{Ang} or B{C{unit}}).
           @arg omg2: Longitude of the second point (C{Ang} or B{C{unit}}).
           @kwarg caps: Desired capabilities for the L{GeodesicLine3Solve} instance.
           @kwarg unit_name: Optional C{B{name}=NN} (C{str}) and scalar C{B{unit}=}L{Degrees}
                             (or L{Degrees}).

           @return: A L{GeodesicLine3Solve} instance.

           @note: Both B{C{bet1}} and B{C{bet2}} should in the range C{[-90, +90]}.

           @see: C++ U{GeodesicExact.InverseLine
                 <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1GeodesicExact.html>} and
                 Python U{Geodesic.InverseLine<https://GeographicLib.SourceForge.io/Python/doc/code.html>}.
        '''
        r  = self.Inverse(bet1, omg1, bet2, omg2, **unit_name)
        gl = GeodesicLine3Solve(self, r.bet1, r.omg1, r.alp1, caps=caps, **unit_name)
#       gl._a13 = r.a12  # gl.SetArc(r.a12)
#       gl._s13 = r.s12  # gl.SetDistance(r.s12)
        return gl

    def toStr(self, **prec_sep_other):  # PYCHOK signature
        '''Return this C{Geodesic3Solve} as string.

           @kwarg prec_sep: See L{toStr<pygeodesy.solveBase._Solve3Base.toStr>}.

           @return: Geodesic3Solve items (C{str}).
        '''
        return _Solve3Base.toStr(self, Geod3Solve=self.Geod3Solve, **prec_sep_other)

    @property_RO
    def triaxial3(self):
        '''Get the triaxial (C{Triaxial3}).
        '''
        return self._triaxial3


class GeodesicLine3Solve(_Geodesic3SolveBase):  # _SolveGDictLineBase):
    '''Wrapper to invoke I{Karney}'s U{Geod3Solve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.3.html>}
       as a C{triaxial} version of I{Karney}'s Python class U{GeodesicLine<https://GeographicLib.SourceForge.io/
       C++/doc/python/code.html#geographiclib.geodesicline.GeodesicLine>}.

       @note: Use property C{GeodSolve} or env variable C{PYGEODESY_GEODSOLVE} to specify the (fully
              qualified) path to the C{GeodSolve} executable.

       @note: This C{geodesic} is intended I{for testing purposes only}, it invokes the C{GeodSolve}
              executable for I{every} method call.
    '''
#   _a13   = \
#   _s13   = NAN  # see GeodesicSolve._InverseLine
    _solve = None  # L{Geodesic3Solve} instance

    def __init__(self, geodesic3, bet1, omg1, alp1, caps=Caps.ALL, **unit_name):
        '''New L{GeodesicLine3Solve} instance, allowing points to be found along
           a geodesic starting at C{(B{bet1}, B{omg1})} with azimuth B{C{alp1}}.

           @arg geodesic3: The geodesic to use (L{Geodesic3Solve}).
           @arg bet1: Ellipsoidal latitude of the first point (C{Ang} or B{C{unit}}).
           @arg omg1: Ellipsoidal longitude of the first point (C{Ang} or B{C{unit}}).
           @arg alp1: Azimuth at the first point (compass C{degrees}, C{Ang} or C{unit}).
           @kwarg caps: Bit-or'ed combination of L{Caps<pygeodesy.karney.Caps>} values
                        specifying the capabilities the L{GeodesicLine3Solve} instance
                        should possess, C{B{caps}=Caps.ALL} always.  Include C{Caps.LINE_OFF}
                        if updates to the B{C{geodesic3}} should I{not be reflected} in this
                        L{GeodesicLine3Solve} instance.
           @kwarg unit_name: Optional C{B{name}=NN} (C{str}) and scalar C{B{unit}=}L{Degrees}
                             (or L{Degrees}).

           @raise Geodesic3Error: Invalid path for the C{Geod3Solve} executable or isn't the
                                  C{Geod3Solve} executable, see property C{geodesic3.Geod3Solve}.

           @raise TypeError: Invalid B{C{geodesic3}}.
        '''
        _xinstanceof(Geodesic3Solve, geodesic3=geodesic3)
        if (caps & Caps.LINE_OFF):  # copy to avoid updates
            geodesic3 = geodesic3.copy(deep=False, name=_UNDER_(NN, geodesic3.name))  # NOT _under!
        self._bet1, \
        self._omg1, \
        self._alp1  = _toDegrees(bet1, omg1, alp1, **unit_name)
        self._caps  =  caps | Caps._AZIMUTH_LATITUDE_LONG_UNROLL
        self._debug =  geodesic3._debug & Caps._DEBUG_ALL
        self._solve =  geodesic3
        try:
            self.Geod3Solve = geodesic3.Geod3Solve  # geodesic or copy of geodesic
        except GeodesicError:
            pass
        if unit_name:
            name = _xkwds_get(unit_name, name=None)
            if name:
                self.name = name

    @Property_RO
    def alp1(self):
        '''Get the azimuth at the first point (compass C{degrees}).
        '''
        return self._alp1

#   azi12 = alp1  # like RhumbLineSolve

#   @Property_RO
#   def azi1_sincos2(self):
#       '''Get the sine and cosine of the first point's azimuth (2-tuple C{(sin, cos)}).
#       '''
#       return _sincos2d(self.alp1)
#
#   azi12_sincos2 = azi1_sincos2

    @Property_RO
    def bet1(self):
        '''Get the I{ellipsoidal} latitude of the first point (C{degrees}).
        '''
        return self._bet1

    @Property_RO
    def _cmdDistance(self):
        '''(INTERNAL) Get the C{Geod3Solve} I{-L} cmd (C{tuple}).
        '''
        t = self.bet1, self.omg1, self.alp1
        return self._cmdBasic + ('-L',) + self._toStdin(t)

    @property_RO
    def geodesic3(self):
        '''Get the C{triaxial} geodesic (L{Geodesic3Solve}).
        '''
        return self._solve  # see .solveBase._SolveLineBase

    def Intersecant2(self, bet0, omg0, radius, **kwds):  # PYCHOK no cover
        '''B{Not implemented}, throws a C{NotImplementedError} always.'''
        self._notImplemented(bet0, omg0, radius, **kwds)

    @Property_RO
    def omg1(self):
        '''Get the I{ellipsoidal} longitude of the first point (C{degrees}).
        '''
        return self._omg1

    def PlumbTo(self, bet0, omg0, **kwds):  # PYCHOK no cover
        '''B{Not implemented}, throws a C{NotImplementedError} always.'''
        self._notImplemented(bet0, omg0, **kwds)

    def Position(self, s12, outmask=Caps.ALL):  # PYCHOK usused
        '''Find the position on the line given B{C{s12}}.

           @arg s12: Distance from the first point to the second (C{meter}).

           @return: A C{GDict} with 7 items C{bet1, omg1, alp1, bet2, omg2, alp2, s12}.
        '''
        r = self._GDictInvoke(self._cmdDistance, self._Names_Distance, s12)  # ._unCaps(outmask)
        return self.geodesic3._a12d(r)

    def toStr(self, **prec_sep_other):  # PYCHOK signature
        '''Return this C{GeodesicLine3Solve} as string.

           @kwarg prec_sep: See L{toStr<pygeodesy.solveBase._Solve3Base.toStr>}.

           @return: GeodesicLine3Solve items (C{str}).
        '''
        return _Solve3Base.toStr(self, bet1=self.bet1, omg1=self.omg1, alp1=self.alp1,
                                       geodesic3=self.geodesic3, **prec_sep_other)

    @property_RO
    def triaxial3(self):
        '''Get the triaxial (C{Triaxial3}).
        '''
        return self.geodesic3.triaxial3


def _toAzi(alp):  # as degrees
    return alp.degrees0 if isAng(alp) else alp


def _toDegrees(*angs, **unit_name):
    unit = _xkwds_get(unit_name, unit=Degrees)
    for ang in angs:
        if not isAng(ang):
            ang = Ang.fromScalar(ang, unit=unit)
        yield ang.degrees


__all__ += _ALL_DOCS(_Geodesic3SolveBase)

if __name__ == _DMAIN_:

    def _main():
        from pygeodesy import printf
        from sys import argv

        gS = Geodesic3Solve(name='Test')
        gS.verbose = v = '--verbose' in argv  # or '-v' in argv

        if not _Xables.X_OK(gS.Geod3Solve):  # not set
            gS.Geod3Solve = _Xables.Geod3Solve(_Xables.bin_)
        printf('version: %s', gS.version, nt=v)

        r = gS.Direct(40.6, -73.8, 51, 5.5e6)
        printf('Direct: %r', r, nt=v)

        printf('Inverse: %r',  gS.Inverse( 40.6, -73.8, 51.6, -0.5), nt=v)

        glS = GeodesicLine3Solve(gS, 40.6, -73.8, 51, name='LineTest')
        printf('Line: %r', glS)
        p = glS.Position(5.5e6)
        printf('Position: %r %s', p, p == r)

    _main()

# % python3 -m pygeodesy.geod3solve

# version: /opt/local/bin/Geod3Solve: GeographicLib version 2.7
# Direct: GDict(a12=49.410276, alp1=51.0, alp2=107.340251, bet1=40.6, bet2=51.895816, omg1=-73.8, omg2=-1.038308, s12=5500000.0)
# Inverse: GDict(a12=49.816802, alp1=51.235527, alp2=107.918138, bet1=40.6, bet2=51.6, omg1=-73.8, omg2=-0.5, s12=5545275.651379)
# Line: GeodesicLine3Solve(alp1=51.0, bet1=40.6, geodesic3=Geodesic3Solve(Geod3Solve='/opt/local/bin/Geod3Solve', invokation=3, status=0), invokation=1, omg1=-73.8, status=0)
# Position: GDict(a12=49.410276, alp1=51.0, alp2=107.340251, bet1=40.6, bet2=51.895816, omg1=-73.8, omg2=-1.038308, s12=5500000.0) True


# % python3 -m pygeodesy.geod3solve --verbose

# Geodesic3Solve 'Test'@1: /opt/local/bin/Geod3Solve --version (invoke)
# Geodesic3Solve 'Test'@1: '/opt/local/bin/Geod3Solve: GeographicLib version 2.7' (0, stdout/-err)
# Geodesic3Solve 'Test'@1: /opt/local/bin/Geod3Solve: GeographicLib version 2.7 (0)
# version: /opt/local/bin/Geod3Solve: GeographicLib version 2.7
#
# Geodesic3Solve 'Test'@2: /opt/local/bin/Geod3Solve -f -t 6378172.0 6378102.0 6356752.0 -p 10 \ 40.6 -73.8 51.0 5500000.0 (Direct)
# Geodesic3Solve 'Test'@2: '40.600000000000001 -73.799999999999997 51.000000000000007 51.895816223972680 -1.038308043217667 107.340251322641734 5500000.0000000000' (0, stdout/-err)
# Geodesic3Solve 'Test'@2: bet1=40.600000000000001, omg1=-73.799999999999997, alp1=51.000000000000007, bet2=51.89581622397268, omg2=-1.038308043217667, alp2=107.340251322641734, s12=5500000.0 (0)
# Direct: GDict(a12=49.410276, alp1=51.0, alp2=107.340251, bet1=40.6, bet2=51.895816, omg1=-73.8, omg2=-1.038308, s12=5500000.0)
#
# Geodesic3Solve 'Test'@3: /opt/local/bin/Geod3Solve -f -t 6378172.0 6378102.0 6356752.0 -p 10 -i \ 40.6 -73.8 51.6 -0.5 (Inverse)
# Geodesic3Solve 'Test'@3: '40.600000000000001 -73.799999999999997 51.235527494379824 51.600000000000001 -0.500000000000000 107.918137616344865 5545275.6513788253' (0, stdout/-err)
# Geodesic3Solve 'Test'@3: bet1=40.600000000000001, omg1=-73.799999999999997, alp1=51.235527494379824, bet2=51.600000000000001, omg2=-0.5, alp2=107.918137616344865, s12=5545275.6513788253 (0)
# Inverse: GDict(a12=49.816802, alp1=51.235527, alp2=107.918138, bet1=40.6, bet2=51.6, omg1=-73.8, omg2=-0.5, s12=5545275.651379)
#
# Line: GeodesicLine3Solve(alp1=51.0, bet1=40.6, geodesic3=Geodesic3Solve(Geod3Solve='/opt/local/bin/Geod3Solve', invokation=3, status=0), invokation=1, omg1=-73.8, status=0)
# Position: GDict(a12=49.410276, alp1=51.0, alp2=107.340251, bet1=40.6, bet2=51.895816, omg1=-73.8, omg2=-1.038308, s12=5500000.0) True


# Examples <https://GeographicLib.SourceForge.io/C++/doc/Geod3Solve.1.html>

# % echo 40:38:23N 073:46:44W-19.43W X 01:21:33N 103:59:22E-19.43W | tr X '\n' | Cart3Convert -G | Cart3Convert -E -r | tr '\n' ' ' | Geod3Solve -i -p 0 -f
# 40.57193 -54.38111 3.20824 1.35529 123.41971 177.48319 15347602
# % echo 40:38:23N 073:46:44W-19.43W X 01:21:33N 103:59:22E-19.43W | tr X '\n' | Cart3Convert -G | Cart3Convert -E -r | tr '\n' ' ' | Geod3Solve -i -p 0
#                    3.20824                   177.48319 15347602
# % echo 40:38:23N 073:46:44W-19.43W X 01:21:33N 103:59:22E-19.43W | tr X '\n' | Cart3Convert -G | Cart3Convert -E -r | tr '\n' ' ' | Geod3Solve -i -p 9
#                    3.20824419242752          177.48319457984906 15347602.214888915

# % Geod3Solve -L 40.57193 -54.38111 3.20824
# 15347602
# 1.35529159 123.41971119 177.48319768
# 15347602
# 1.35529159 123.41971119 177.48319768
# % Geod3Solve -L 40.57193 -54.38111 3.20824 -f
# 15347602
# 40.57193000 -54.38111000 3.20824000 1.35529159 123.41971119 177.48319768 15347602.000
# % Geod3Solve -L 40.57193 -54.38111 3.20824 -f -p 9
# 15347602
# 40.57193000000000 -54.38111000000000 3.20824000000000 1.35529159010289 123.41971119123770 177.48319768195134 15347602.000000000

# **) MIT License
#
# Copyright (C) 2016-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
