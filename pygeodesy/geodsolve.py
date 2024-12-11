
# -*- coding: utf-8 -*-

u'''Wrapper to invoke I{Karney}'s U{GeodSolve
<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>} utility
as an (exact) geodesic, but intended I{for testing purposes only}.

Set env variable C{PYGEODESY_GEODSOLVE} to the (fully qualified) path
of the C{GeodSolve} executable.
'''

from pygeodesy.basics import _xinstanceof
# from pygeodesy.constants import NAN, _0_0  # from .karney
# from pygeodesy.geodesicx import GeodesicAreaExact  # _MODS
from pygeodesy.interns import NN, _UNDER_
from pygeodesy.karney import Caps, GeodesicError, GeodSolve12Tuple, \
                            _sincos2d, _Xables,  _0_0, NAN
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _name1__
from pygeodesy.namedTuples import Destination3Tuple, Distance3Tuple
from pygeodesy.props import Property, Property_RO, property_RO
from pygeodesy.solveBase import _SolveGDictBase, _SolveGDictLineBase
from pygeodesy.utily import _unrollon, _Wrap, wrap360

__all__ = _ALL_LAZY.geodsolve
__version__ = '24.11.02'


class _GeodesicSolveBase(_SolveGDictBase):
    '''(INTERNAL) Base class for L{GeodesicSolve} and L{GeodesicLineSolve}.
    '''
    _Error         =  GeodesicError
    _Names_Direct  = \
    _Names_Inverse =  GeodSolve12Tuple._Names_
    _Xable_name    = _Xables.GeodSolve.__name__
    _Xable_path    = _Xables.GeodSolve()

    @Property_RO
    def _b_option(self):
        return ('-b',) if self.reverse2 else ()

    @Property_RO
    def _cmdBasic(self):
        '''(INTERNAL) Get the basic C{GeodSolve} cmd (C{tuple}).
        '''
        return (self.GeodSolve, '-f') + (self._b_option +
                                         self._e_option +
                                         self._E_option +
                                         self._p_option +
                                         self._u_option)

    @Property
    def GeodSolve(self):
        '''Get the U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>}
           executable (C{filename}).
        '''
        return self._Xable_path

    @GeodSolve.setter  # PYCHOK setter!
    def GeodSolve(self, path):
        '''Set the U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>}
           executable (C{filename}), the (fully qualified) path to the C{GeodSolve} executable.

           @raise GeodesicError: Invalid B{C{path}}, B{C{path}} doesn't exist or
                                 isn't the C{GeodSolve} executable.
        '''
        self._setXable(path)

    def toStr(self, **prec_sep):  # PYCHOK signature
        '''Return this C{GeodesicSolve} as string.

           @kwarg prec_sep: Keyword argumens C{B{prec}=6} and C{B{sep}=", "}
                       for the C{float} C{prec}ision, number of decimal digits
                       (0..9) and the C{sep}arator string to join.  Trailing
                       zero decimals are stripped for B{C{prec}} values of 1
                       and above, but kept for negative B{C{prec}} values.

           @return: GeodesicSolve items (C{str}).
        '''
        return _SolveGDictBase._toStr(self, GeodSolve=self.GeodSolve, **prec_sep)

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

    def Area(self, polyline=False, **name):
        '''Set up a L{GeodesicAreaExact} to compute area and perimeter
           of a polygon.

           @kwarg polyline: If C{True}, compute the perimeter only, otherwise
                            perimeter and area (C{bool}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{GeodesicAreaExact} instance.

           @note: The B{C{debug}} setting is passed as C{verbose}
                  to the returned L{GeodesicAreaExact} instance.
        '''
        gaX = _MODS.geodesicx.GeodesicAreaExact(self, polyline=polyline, **name)
        if self.verbose or self.debug:  # PYCHOK no cover
            gaX.verbose = True
        return gaX

    Polygon = Area  # for C{geographiclib} compatibility

    def Direct3(self, lat1, lon1, azi1, s12):  # PYCHOK outmask
        '''Return the destination lat, lon and reverse azimuth (final bearing)
           in C{degrees}.

           @return: L{Destination3Tuple}C{(lat, lon, final)}.
        '''
        r = self._GDictDirect(lat1, lon1, azi1, False, s12, floats=False)
        return Destination3Tuple(float(r.lat2), float(r.lon2), wrap360(r.azi2),
                                 iteration=r._iteration)

    def _DirectLine(self, ll1, azi12, **caps_name):  # PYCHOK no cover
        '''(INTERNAL) Short-cut version.
        '''
        return self.DirectLine(ll1.lat, ll1.lon, azi12, **caps_name)

    def DirectLine(self, lat1, lon1, azi1, **caps_name):
        '''Set up a L{GeodesicLineSolve} to compute several points
           on a single geodesic.

           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi1: Azimuth at the first point (compass C{degrees}).
           @kwarg caps_name: Optional C{B{name}=NN} (C{str}) and keyword
                       argument C{B{caps}=Caps.ALL}, bit-or'ed combination
                       of L{Caps} values specifying the capabilities the
                       L{GeodesicLineSolve} instance should possess.

           @return: A L{GeodesicLineSolve} instance.

           @note: If the point is at a pole, the azimuth is defined by keeping
                  B{C{lon1}} fixed, writing C{B{lat1} = ±(90 − ε)}, and taking
                  the limit C{ε → 0+}.

           @see: C++ U{GeodesicExact.Line
                 <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1GeodesicExact.html>}
                 and Python U{Geodesic.Line<https://GeographicLib.SourceForge.io/Python/doc/code.html>}.
        '''
        return GeodesicLineSolve(self, lat1, lon1, azi1, **_name1__(caps_name, _or_nameof=self))

    Line = DirectLine

    def _Inverse(self, ll1, ll2, wrap, **outmask):  # PYCHOK no cover
        '''(INTERNAL) Short-cut version, see .ellipsoidalBaseDI.intersecant2.
        '''
        if wrap:
            ll2 = _unrollon(ll1, _Wrap.point(ll2))
        return self.Inverse(ll1.lat, ll1.lon, ll2.lat, ll2.lon, **outmask)

    def Inverse3(self, lat1, lon1, lat2, lon2):  # PYCHOK outmask
        '''Return the distance in C{meter} and the forward and
           reverse azimuths (initial and final bearing) in C{degrees}.

           @return: L{Distance3Tuple}C{(distance, initial, final)}.
        '''
        r = self._GDictInverse(lat1, lon1, lat2, lon2, floats=False)
        return Distance3Tuple(float(r.s12), wrap360(r.azi1), wrap360(r.azi2),
                              iteration=r._iteration)

    def _InverseLine(self, ll1, ll2, wrap, **caps_name):  # PYCHOK no cover
        '''(INTERNAL) Short-cut version.
        '''
        if wrap:
            ll2 = _unrollon(ll1, _Wrap.point(ll2))
        return self.InverseLine(ll1.lat, ll1.lon, ll2.lat, ll2.lon, **caps_name)

    def InverseLine(self, lat1, lon1, lat2, lon2, **caps_name):  # PYCHOK no cover
        '''Set up a L{GeodesicLineSolve} to compute several points
           on a single geodesic.

           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg lat2: Latitude of the second point (C{degrees}).
           @arg lon2: Longitude of the second point (C{degrees}).
           @kwarg caps_name: Optional C{B{name}=NN} (C{str}) and keyword
                       argument C{B{caps}=Caps.ALL}, bit-or'ed combination
                       of L{Caps} values specifying the capabilities the
                       L{GeodesicLineSolve} instance should possess.

           @return: A L{GeodesicLineSolve} instance.

           @note: Both B{C{lat1}} and B{C{lat2}} should in the range C{[-90, +90]}.

           @see: C++ U{GeodesicExact.InverseLine
                 <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1GeodesicExact.html>} and
                 Python U{Geodesic.InverseLine<https://GeographicLib.SourceForge.io/Python/doc/code.html>}.
        '''
        r  = self.Inverse(lat1, lon1, lat2, lon2)
        gl = GeodesicLineSolve(self, lat1, lon1, r.azi1, **_name1__(caps_name, _or_nameof=self))
        gl._a13 = r.a12  # gl.SetArc(r.a12)
        gl._s13 = r.s12  # gl.SetDistance(r.s12)
        return gl


class GeodesicLineSolve(_GeodesicSolveBase, _SolveGDictLineBase):
    '''Wrapper to invoke I{Karney}'s U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>}
       as an C{Exact} version of I{Karney}'s Python class U{GeodesicLine<https://GeographicLib.SourceForge.io/C++/doc/
       python/code.html#geographiclib.geodesicline.GeodesicLine>}.

       @note: Use property C{GeodSolve} or env variable C{PYGEODESY_GEODSOLVE} to specify the (fully
              qualified) path to the C{GeodSolve} executable.

       @note: This C{geodesic} is intended I{for testing purposes only}, it invokes the C{GeodSolve}
              executable for I{every} method call.
    '''
    _a13 = \
    _s13 = NAN  # see GeodesicSolve._InverseLine

    def __init__(self, geodesic, lat1, lon1, azi1, caps=Caps.ALL, **name):
        '''New L{GeodesicLineSolve} instance, allowing points to be found along
           a geodesic starting at C{(B{lat1}, B{lon1})} with azimuth B{C{azi1}}.

           @arg geodesic: The geodesic to use (L{GeodesicSolve}).
           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi1: Azimuth at the first points (compass C{degrees}).
           @kwarg caps: Bit-or'ed combination of L{Caps} values specifying the
                        capabilities the L{GeodesicLineSolve} instance should possess,
                        C{B{caps}=Caps.ALL} always.  Include C{Caps.LINE_OFF} if
                        updates to the B{C{geodesic}} should I{not} be reflected in
                        this L{GeodesicLineSolve} instance.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise GeodesicError: Invalid path for the C{GeodSolve} executable
                                 or isn't the C{GeodSolve} executable, see
                                 property C{geodesic.GeodSolve}.

           @raise TypeError: Invalid B{C{geodesic}}.
        '''
        _xinstanceof(GeodesicSolve, geodesic=geodesic)
        if (caps & Caps.LINE_OFF):  # copy to avoid updates
            geodesic = geodesic.copy(deep=False, name=_UNDER_(NN, geodesic.name))  # NOT _under!
        _SolveGDictLineBase.__init__(self, geodesic, lat1, lon1, caps, azi1=azi1, **name)
        try:
            self.GeodSolve = geodesic.GeodSolve  # geodesic or copy of geodesic
        except GeodesicError:
            pass

    @Property_RO
    def a13(self):
        '''Get the arc length to reference point 3 (C{degrees}).

           @see: Methods L{Arc} and L{SetArc}.
        '''
        return self._a13

    def Arc(self):  # PYCHOK no cover
        '''Return the arc length to reference point 3 (C{degrees} or C{NAN}).

           @see: Method L{SetArc} and property L{a13}.
        '''
        return self.a13

    def ArcPosition(self, a12, outmask=Caps.STANDARD):  # PYCHOK no cover
        '''Find the position on the line given B{C{a12}}.

           @arg a12: Spherical arc length from the first point to the
                     second point (C{degrees}).

           @return: A C{GDict} with 12 items C{lat1, lon1, azi1, lat2, lon2,
                    azi2, m12, a12, s12, M12, M21, S12}.
        '''
        return self._GDictInvoke(self._cmdArc, self._Names_Direct, a12)._unCaps(outmask)

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

    def Distance(self):
        '''Return the distance to reference point 3 (C{meter} or C{NAN}).
        '''
        return self.s13

    @property_RO
    def geodesic(self):
        '''Get the geodesic (L{GeodesicSolve}).
        '''
        return self._solve  # see .solveBase._SolveLineBase

    def Intersecant2(self, lat0, lon0, radius, **kwds):  # PYCHOK no cover
        '''B{Not implemented}, throws a C{NotImplementedError} always.'''
        self._notImplemented(lat0, lon0, radius, **kwds)

    def PlumbTo(self, lat0, lon0, **kwds):  # PYCHOK no cover
        '''B{Not implemented}, throws a C{NotImplementedError} always.'''
        self._notImplemented(lat0, lon0, **kwds)

    def Position(self, s12, outmask=Caps.STANDARD):
        '''Find the position on the line given B{C{s12}}.

           @arg s12: Distance from the first point to the second (C{meter}).

           @return: A C{GDict} with 12 items C{lat1, lon1, azi1, lat2, lon2,
                    azi2, m12, a12, s12, M12, M21, S12}, possibly C{a12=NAN}.
        '''
        return self._GDictInvoke(self._cmdDistance, self._Names_Direct, s12)._unCaps(outmask)

    @Property_RO
    def s13(self):
        '''Get the distance to reference point 3 (C{meter} or C{NAN}).

           @see: Methods L{Distance} and L{SetDistance}.
        '''
        return self._s13

    def SetArc(self, a13):  # PYCHOK no cover
        '''Set reference point 3 in terms relative to the first point.

           @arg a13: Spherical arc length from the first to the reference
                     point (C{degrees}).

           @return: The distance C{s13} (C{meter}) between the first and
                    the reference point or C{NAN}.
        '''
        if self._a13 != a13:
            self._a13 = a13
            self._s13 = self.ArcPosition(a13, outmask=Caps.DISTANCE).s12  # if a13 else _0_0
#           _update_all(self)
        return self._s13

    def SetDistance(self, s13):  # PYCHOK no cover
        '''Set reference point 3 in terms relative to the first point.

           @arg s13: Distance from the first to the reference point (C{meter}).

           @return: The arc length C{a13} (C{degrees}) between the first and
                    the reference point or C{NAN}.
        '''
        if self._s13 != s13:
            self._s13 = s13
            self._a13 = self.Position(s13, outmask=Caps.DISTANCE).a12 if s13 else _0_0
#           _update_all(self)
        return self._a13  # NAN for GeodesicLineExact without Cap.DISTANCE_IN

    def toStr(self, **prec_sep):  # PYCHOK signature
        '''Return this C{GeodesicLineSolve} as string.

           @kwarg prec_sep: Keyword argumens C{B{prec}=6} and C{B{sep}=", "}
                       for the C{float} C{prec}ision, number of decimal digits
                       (0..9) and the C{sep}arator string to join.  Trailing
                       zero decimals are stripped for B{C{prec}} values of 1
                       and above, but kept for negative B{C{prec}} values.

           @return: GeodesicLineSolve items (C{str}).
        '''
        return _SolveGDictLineBase._toStr(self, azi1=self.azi1, geodesic=self._solve,
                                                GeodSolve=self.GeodSolve, **prec_sep)


__all__ += _ALL_DOCS(_GeodesicSolveBase)

if __name__ == '__main__':

    def _main():
        from pygeodesy import printf
        from sys import argv

        gS = GeodesicSolve(name='Test')
        gS.verbose = '--verbose' in argv  # or '-v' in argv

        if not _Xables.X_OK(gS.GeodSolve):  # not set
            gS.GeodSolve = _Xables.GeodSolve(_Xables.bin_)
        printf('version: %s', gS.version)

        r = gS.Direct(40.6, -73.8, 51, 5.5e6)
        printf('Direct: %r', r, nl=1)
        printf('Direct3: %r', gS.Direct3(40.6, -73.8, 51, 5.5e6))

        printf('Inverse: %r',  gS.Inverse( 40.6, -73.8, 51.6, -0.5), nl=1)
        printf('Inverse1: %r', gS.Inverse1(40.6, -73.8, 51.6, -0.5))
        printf('Inverse3: %r', gS.Inverse3(40.6, -73.8, 51.6, -0.5))

        glS = GeodesicLineSolve(gS, 40.6, -73.8, 51, name='LineTest')
        p = glS.Position(5.5e6)
        printf('Position:    %5s %r', p == r, p, nl=1)
        p = glS.ArcPosition(49.475527)
        printf('ArcPosition: %5s %r', p == r, p)

    _main()

# % python3 -m pygeodesy.geodsolve

# version: /opt/local/bin/GeodSolve: GeographicLib version 2.2

# Direct: GDict(a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, m12=4844148.703101, M12=0.650911, M21=0.651229, s12=5500000.0, S12=39735075134877.09375)
# Direct3: Destination3Tuple(lat=51.884565, lon=-1.141173, final=107.189397)

# Inverse: GDict(a12=49.94131, azi1=51.198883, azi2=107.821777, lat1=40.6, lat2=51.6, lon1=-73.8, lon2=-0.5, m12=4877684.602706, M12=0.64473, M21=0.645046, s12=5551759.400319, S12=40041368848742.53125)
# Inverse1: 49.94131021789904
# Inverse3: Distance3Tuple(distance=5551759.400319, initial=51.198883, final=107.821777)

# Position:    True  GDict(a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, m12=4844148.703101, M12=0.650911, M21=0.651229, s12=5500000.0, S12=39735075134877.09375)
# ArcPosition: False GDict(a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141174, m12=4844148.669561, M12=0.650911, M21=0.651229, s12=5499999.948497, S12=39735074737272.734375)


# % python3 -m pygeodesy.geodsolve

# version: /opt/local/bin/GeodSolve: GeographicLib version 2.3

# Direct: GDict(a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, m12=4844148.703101, M12=0.650911, M21=0.651229, s12=5500000.0, S12=39735075134877.078125)
# Direct3: Destination3Tuple(lat=51.884565, lon=-1.141173, final=107.189397)

# Inverse: GDict(a12=49.94131, azi1=51.198883, azi2=107.821777, lat1=40.6, lat2=51.6, lon1=-73.8, lon2=-0.5, m12=4877684.602706, M12=0.64473, M21=0.645046, s12=5551759.400319, S12=40041368848742.53125)
# Inverse1: 49.94131021789904
# Inverse3: Distance3Tuple(distance=5551759.400319, initial=51.198883, final=107.821777)

# Position:    False GDict(a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, s12=5500000.0)
# ArcPosition: False GDict(a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141174, s12=5499999.948497)


# % python3 -m pygeodesy.geodsolve --verbose

# GeodesicSolve 'Test' 1: /opt/local/bin/GeodSolve --version (invoke)
# GeodesicSolve 'Test' 1: /opt/local/bin/GeodSolve: GeographicLib version 2.2 (0)
# version: /opt/local/bin/GeodSolve: GeographicLib version 2.2
# GeodesicSolve 'Test' 2: /opt/local/bin/GeodSolve -f -E -p 10 \ 40.600000000000001 -73.799999999999997 51.0 5500000.0 (Direct)
# GeodesicSolve 'Test' 2: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.0, lat2=51.884564505606761, lon2=-1.141172861200829, azi2=107.189397162605886, s12=5500000.0, a12=49.475527463251467, m12=4844148.703101486, M12=0.65091056699808603, M21=0.65122865892196558, S12=39735075134877.094 (0)

# Direct: GDict(a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, m12=4844148.703101, M12=0.650911, M21=0.651229, s12=5500000.0, S12=39735075134877.09375)
# GeodesicSolve 'Test' 3: /opt/local/bin/GeodSolve -f -E -p 10 \ 40.600000000000001 -73.799999999999997 51.0 5500000.0 (Direct3)
# GeodesicSolve 'Test' 3: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.0, lat2=51.884564505606761, lon2=-1.141172861200829, azi2=107.189397162605886, s12=5500000.0, a12=49.475527463251467, m12=4844148.703101486, M12=0.65091056699808603, M21=0.65122865892196558, S12=39735075134877.094 (0)
# Direct3: Destination3Tuple(lat=51.884565, lon=-1.141173, final=107.189397)
# GeodesicSolve 'Test' 4: /opt/local/bin/GeodSolve -f -E -p 10 -i \ 40.600000000000001 -73.799999999999997 51.600000000000001 -0.5 (Inverse)
# GeodesicSolve 'Test' 4: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.198882845579824, lat2=51.600000000000001, lon2=-0.5, azi2=107.821776735514248, s12=5551759.4003186841, a12=49.941310217899037, m12=4877684.6027061976, M12=0.64472969205948238, M21=0.64504567852134398, S12=40041368848742.531 (0)

# Inverse: GDict(a12=49.94131, azi1=51.198883, azi2=107.821777, lat1=40.6, lat2=51.6, lon1=-73.8, lon2=-0.5, m12=4877684.602706, M12=0.64473, M21=0.645046, s12=5551759.400319, S12=40041368848742.53125)
# GeodesicSolve 'Test' 5: /opt/local/bin/GeodSolve -f -E -p 10 -i \ 40.600000000000001 -73.799999999999997 51.600000000000001 -0.5 (Inverse1)
# GeodesicSolve 'Test' 5: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.198882845579824, lat2=51.600000000000001, lon2=-0.5, azi2=107.821776735514248, s12=5551759.4003186841, a12=49.941310217899037, m12=4877684.6027061976, M12=0.64472969205948238, M21=0.64504567852134398, S12=40041368848742.531 (0)
# Inverse1: 49.94131021789904
# GeodesicSolve 'Test' 6: /opt/local/bin/GeodSolve -f -E -p 10 -i \ 40.600000000000001 -73.799999999999997 51.600000000000001 -0.5 (Inverse3)
# GeodesicSolve 'Test' 6: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.198882845579824, lat2=51.600000000000001, lon2=-0.5, azi2=107.821776735514248, s12=5551759.4003186841, a12=49.941310217899037, m12=4877684.6027061976, M12=0.64472969205948238, M21=0.64504567852134398, S12=40041368848742.531 (0)
# Inverse3: Distance3Tuple(distance=5551759.400319, initial=51.198883, final=107.821777)

# Position:    True  GDict(a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, m12=4844148.703101, M12=0.650911, M21=0.651229, s12=5500000.0, S12=39735075134877.09375)
# ArcPosition: False GDict(a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141174, m12=4844148.669561, M12=0.650911, M21=0.651229, s12=5499999.948497, S12=39735074737272.734375)


# % python3 -m pygeodesy.geodsolve --verbose

# GeodesicSolve 'Test'@1: /opt/local/bin/GeodSolve --version (invoke)
# GeodesicSolve 'Test'@1: '/opt/local/bin/GeodSolve: GeographicLib version 2.3' (0, stdout/-err)
# GeodesicSolve 'Test'@1: /opt/local/bin/GeodSolve: GeographicLib version 2.3 (0)
# version: /opt/local/bin/GeodSolve: GeographicLib version 2.3
# GeodesicSolve 'Test'@2: /opt/local/bin/GeodSolve -f -E -p 10 \ 40.600000000000001 -73.799999999999997 51.0 5500000.0 (Direct)
# GeodesicSolve 'Test'@2: '40.600000000000001 -73.799999999999997 51.000000000000000 51.884564505606761 -1.141172861200843 107.189397162605871 5500000.0000000000 49.475527463251460 4844148.7031014860 0.65091056699808614 0.65122865892196569 39735075134877.078' (0, stdout/-err)
# GeodesicSolve 'Test'@2: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.0, lat2=51.884564505606761, lon2=-1.141172861200843, azi2=107.189397162605871, s12=5500000.0, a12=49.47552746325146, m12=4844148.703101486, M12=0.65091056699808614, M21=0.65122865892196569, S12=39735075134877.078 (0)

# Direct: GDict(a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, m12=4844148.703101, M12=0.650911, M21=0.651229, s12=5500000.0, S12=39735075134877.078125)
# GeodesicSolve 'Test'@3: /opt/local/bin/GeodSolve -f -E -p 10 \ 40.600000000000001 -73.799999999999997 51.0 5500000.0 (Direct3)
# GeodesicSolve 'Test'@3: '40.600000000000001 -73.799999999999997 51.000000000000000 51.884564505606761 -1.141172861200843 107.189397162605871 5500000.0000000000 49.475527463251460 4844148.7031014860 0.65091056699808614 0.65122865892196569 39735075134877.078' (0, stdout/-err)
# GeodesicSolve 'Test'@3: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.0, lat2=51.884564505606761, lon2=-1.141172861200843, azi2=107.189397162605871, s12=5500000.0, a12=49.47552746325146, m12=4844148.703101486, M12=0.65091056699808614, M21=0.65122865892196569, S12=39735075134877.078 (0)
# Direct3: Destination3Tuple(lat=51.884565, lon=-1.141173, final=107.189397)
# GeodesicSolve 'Test'@4: /opt/local/bin/GeodSolve -f -E -p 10 -i \ 40.600000000000001 -73.799999999999997 51.600000000000001 -0.5 (Inverse)
# GeodesicSolve 'Test'@4: '40.600000000000001 -73.799999999999997 51.198882845579824 51.600000000000001 -0.500000000000000 107.821776735514248 5551759.4003186813 49.941310217899037 4877684.6027061967 0.64472969205948238 0.64504567852134398 40041368848742.531' (0, stdout/-err)
# GeodesicSolve 'Test'@4: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.198882845579824, lat2=51.600000000000001, lon2=-0.5, azi2=107.821776735514248, s12=5551759.4003186813, a12=49.941310217899037, m12=4877684.6027061967, M12=0.64472969205948238, M21=0.64504567852134398, S12=40041368848742.531 (0)

# Inverse: GDict(a12=49.94131, azi1=51.198883, azi2=107.821777, lat1=40.6, lat2=51.6, lon1=-73.8, lon2=-0.5, m12=4877684.602706, M12=0.64473, M21=0.645046, s12=5551759.400319, S12=40041368848742.53125)
# GeodesicSolve 'Test'@5: /opt/local/bin/GeodSolve -f -E -p 10 -i \ 40.600000000000001 -73.799999999999997 51.600000000000001 -0.5 (Inverse1)
# GeodesicSolve 'Test'@5: '40.600000000000001 -73.799999999999997 51.198882845579824 51.600000000000001 -0.500000000000000 107.821776735514248 5551759.4003186813 49.941310217899037 4877684.6027061967 0.64472969205948238 0.64504567852134398 40041368848742.531' (0, stdout/-err)
# GeodesicSolve 'Test'@5: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.198882845579824, lat2=51.600000000000001, lon2=-0.5, azi2=107.821776735514248, s12=5551759.4003186813, a12=49.941310217899037, m12=4877684.6027061967, M12=0.64472969205948238, M21=0.64504567852134398, S12=40041368848742.531 (0)
# Inverse1: 49.94131021789904
# GeodesicSolve 'Test'@6: /opt/local/bin/GeodSolve -f -E -p 10 -i \ 40.600000000000001 -73.799999999999997 51.600000000000001 -0.5 (Inverse3)
# GeodesicSolve 'Test'@6: '40.600000000000001 -73.799999999999997 51.198882845579824 51.600000000000001 -0.500000000000000 107.821776735514248 5551759.4003186813 49.941310217899037 4877684.6027061967 0.64472969205948238 0.64504567852134398 40041368848742.531' (0, stdout/-err)
# GeodesicSolve 'Test'@6: lat1=40.600000000000001, lon1=-73.799999999999997, azi1=51.198882845579824, lat2=51.600000000000001, lon2=-0.5, azi2=107.821776735514248, s12=5551759.4003186813, a12=49.941310217899037, m12=4877684.6027061967, M12=0.64472969205948238, M21=0.64504567852134398, S12=40041368848742.531 (0)
# Inverse3: Distance3Tuple(distance=5551759.400319, initial=51.198883, final=107.821777)

# Position:    False GDict(a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, s12=5500000.0)
# ArcPosition: False GDict(a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141174, s12=5499999.948497)

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
