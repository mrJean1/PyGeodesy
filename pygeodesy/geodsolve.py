
# -*- coding: utf-8 -*-

u'''Wrapper to invoke I{Karney}'s U{GeodSolve
<https://GeographicLib.SourceForge.io/html/GeodSolve.1.htmlb>} utility
as an (exact) geodesic, but intended I{for testing purposes only}.

Set env variable C{PYGEODESY_GEODSOLVE} to the (fully qualified) path
of the C{GeodSolve} executable.
'''

from pygeodesy.basics import map2, _xinstanceof
from pygeodesy.datums import _ellipsoidal_datum
from pygeodesy.ellipsoids import Ellipsoids, Ellipsoid2
from pygeodesy.geodesicx.gxbases import _all_caps, Caps, _GeodesicBase
from pygeodesy.interns import DIG, NN, _0_, _COMMASPACE_, _SPACE_
from pygeodesy.interns import _not_  # PYCHOK used!
from pygeodesy.karney import GDict, GeodesicError, GeodSolve12Tuple
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, printf, _sys_version_info2
from pygeodesy.lazily import _env  # PYCHOK used!
from pygeodesy.named import callername
from pygeodesy.namedTuples import Destination3Tuple, Distance3Tuple
from pygeodesy.props import Property, Property_RO, property_RO
from pygeodesy.streprs import Fmt, fstr, fstrzs, pairs, strs
from pygeodesy.units import Precision_
from pygeodesy.utily import sincos2d, unroll180, wrap360

from subprocess import PIPE as _PIPE, Popen as _Popen, STDOUT as _STDOUT

__all__ = _ALL_LAZY.geodsolve
__version__ = '21.07.03'

_PYGEODESY_GEODSOLVE_ = 'PYGEODESY_GEODSOLVE'  # PYCHOK used!

_len_N  = len(GeodSolve12Tuple._Names_)
_SLASH_ = '/'
_stdin_ = 'stdin'


def _cmd_stdin_(cmd, stdin):  # PYCHOK no cover
    '''(INTERNAL) Cmd line, stdin and caller as sC{str}.
    '''
    c = Fmt.PAREN(callername(up=3))
    t = (c,) if stdin is None else (_SLASH_, str(stdin), c)
    return _SPACE_.join(cmd + t)


class _GeodesicSolveBase(_GeodesicBase):
    '''(NTERNAL) Base class for L{GeodesicSolve} and L{GeodesicLineSolve}.
    '''
    _E          =  Ellipsoids.WGS84
    _Exact      =  True
    _GeodSolve  = _env.get(_PYGEODESY_GEODSOLVE_, _PYGEODESY_GEODSOLVE_)
    _invokation =  0
    _prec       =  Precision_(prec=DIG)
    _reverse2   =  False
    _status     =  None
    _text_True  =  dict() if _sys_version_info2 < (3, 7) else dict(text=True)
    _unroll     =  False
    _verbose    =  False

    @Property_RO
    def a(self):
        '''Get the I{equatorial} radius, semi-axis (C{meter}).
        '''
        return self.ellipsoid.a

    @Property_RO
    def _b_option(self):
        return ('-b',) if self.reverse2 else ()

    @Property_RO
    def _cmdBasic(self):
        '''(INTERNAL) Get the C{GeodSolve} basic cmd (C{tuple}).
        '''
        return (self.GeodSolve,) + self._b_option \
                                 + self._e_option \
                                 + self._E_option \
                                 + self._p_option \
                                 + self._u_option + ('-f',)

    @Property_RO
    def _e_option(self):
        if self.ellipsoid is Ellipsoids.WGS84:
            return ()  # default
        a, f = strs((self.a, self.f), fmt=Fmt.F, prec=DIG + 3)  # not .G!
        return ('-e', a, f)

    @Property_RO
    def _E_option(self):
        return ('-E',) if self.Exact else ()

    @Property_RO
    def ellipsoid(self):
        '''Get the ellipsoid (C{Ellipsoid}).
        '''
        return self._E

    @Property
    def Exact(self):
        '''Get the C{GeodesicExact} usage (C{bool}).
        '''
        return self._Exact

    @Exact.setter  # PYCHOK setter!
    def Exact(self, Exact):
        '''Set the C{GeodesicExact} usage (C{bool}).

           @arg Exact: If C{True} use C{GeodesicExact},
                       otherwise use C{Geodesic} (C{bool}).
        '''
        Exact = bool(Exact)
        self._update(Exact != self.Exact)
        self._Exact = Exact

    @Property_RO
    def f(self):
        '''Get the ellipsoid's I{flattening} (C{float}), M{(a - b) / a}, C{0} for spherical, negative for prolate.
        '''
        return self.ellipsoid.f

    def _GDictInvoke(self, cmd, floats, *args):
        '''(INTERNAL) Invoke C{GeodSolve}, return C{GDict}.
        '''
        i = fstr(args, prec=DIG, fmt=Fmt.F, sep=_SPACE_) if args else None  # not Fmt.G!
        t = self._invoke(cmd, stdin=i).lstrip().split()  # 12-/+ tuple
        if len(t) > _len_N:  # instrumented?
            # unzip the name=value pairs to names and values
            n, v = zip(*(p.split('=') for p in t[:-_len_N]))
            v += tuple(t[-_len_N:])
            n += GeodSolve12Tuple._Names_
        else:
            n, v = GeodSolve12Tuple._Names_, t
        if self.verbose:  # PYCHOK no cover
            self._print(_COMMASPACE_.join(map(Fmt.EQUAL, n, map(fstrzs, v))))
        if floats:
            v = map(float, v)
        return GDict(zip(n, v))

    @Property
    def GeodSolve(self):
        '''Get the U{GeodSolve<https://GeographicLib.SourceForge.io/html/GeodSolve.1.html>}
           executable (C{filename}).
        '''
        return self._GeodSolve

    @GeodSolve.setter  # PYCHOK setter!
    def GeodSolve(self, path):
        '''Set the U{GeodSolve<https://GeographicLib.SourceForge.io/html/GeodSolve.1.html>}
           executable (C{filename}).

           @arg path: The (fully qualified) path to the C{GeodSolve} executable (C{str}).

           @raise GeodesicError: Invalid B{C{path}}, B{C{path}} doesn't exist or
                                 isn't the C{GeodSolve} executable.
        '''
        hold = self.GeodSolve
        self._update(path != hold)
        self._GeodSolve = path
        try:
            self.version  # test path and ...
            if self.status:  # ... return code
                raise GeodesicError(GeodSolve=path, status=self.status, txt=_not_(_0_))
            hold = path
        finally:  # restore in case of error
            self._update(hold != self.GeodSolve)
            self._GeodSolve = hold

    @property_RO
    def invokation(self):
        '''Get the most recent C{GeodSolve} invokation number (C{int}).
        '''
        return self._invokation

    def invoke(self, *options, **stdin):
        '''Invoke the C{GeodSolve} executable and return the result.

           @arg options: No, one or several C{GeodSolve} command line
                         options (C{str}s).
           @kwarg stdin: Optional input to pass to C{GeodSolve.stdin} (C{str}).

           @return: The C{GeodSolve.stdout} and C{.stderr} output (C{str}).

           @raise GeodesicError: On any error, including a non-zero return
                                 code from C{GeodSolve}.

           @note: The C{GeodSolve} return code is in L{status}.
        '''
        c = (self.GeodSolve,) + map2(str, options)
        i = stdin.get(_stdin_, None)
        r = self._invoke(c, stdin=i)
        s = self.status
        if s:
            raise GeodesicError(cmd=_cmd_stdin_(c, i), status=s,
                                txt=_not_(_0_))
        if self.verbose:  # PYCHOK no cover
            self._print(r)
        return r

    def _invoke(self, cmd, stdin=None):
        '''(INTERNAL) Invoke the C{GeodSolve} executable, with the
           given B{C{cmd}} line and optional input to B{C{stdin}}.
        '''
        self._invokation += 1
        self._status = t = None
        if self.verbose:  # PYCHOK no cover
            t = _cmd_stdin_(cmd, stdin)
            self._print(t)
        try:
            p = _Popen(cmd, creationflags=0,
                          # executable   =sys.executable,
                          # shell        =True,
                            stdin        =_PIPE,
                            stdout       =_PIPE,  # PYCHOK kwArgs
                            stderr       =_STDOUT,
                          **self._text_True)
            # invoke and write to stdin
            r = p.communicate(stdin)[0]
            if isinstance(r, bytes):  # Python 3+
                r = r.decode('utf-8')

            if len(r) < 6 or r[:5] in ('Error', 'ERROR'):
                raise ValueError(r)

            self._status = p.returncode
        except (IOError, OSError, TypeError, ValueError) as x:
            raise GeodesicError(cmd=t or _cmd_stdin_(cmd, stdin),
                                txt=str(x))
        return r.strip()

    @Property_RO
    def _p_option(self):
        return '-p', str(self.prec - 5)  # -p is distance prec

    @Property
    def prec(self):
        '''Get the precision, number of decimal digits (C{int}).
        '''
        return self._prec

    @prec.setter  # PYCHOK setter!
    def prec(self, prec):
        '''Set the precision for C{angles} in C{degrees}, like C{lat},
           C{lon}, C{azimuth} and C{arc}.

           @arg prec: Number of decimal digits (C{int}, C{0}..L{DIG}).

           @note: The precision for C{distance = B{prec} - 5} or up to
                  10 decimal digits for C{nanometer} and for C{area =
                  B{prec} - 12} or at most C{millimeter} I{squared}.
        '''
        prec = Precision_(prec=prec, high=DIG)
        self._update(prec != self.prec)
        self._prec = prec

    def _print(self, line):  # PYCHOK no cover
        '''(INTERNAL) Print a status line.
        '''
        if self.status is not None:
            line = _SPACE_(line, Fmt.PAREN(self.status))
        printf('%s %d: %s', self.named2, self.invokation, line)

    @Property
    def reverse2(self):
        '''Get the C{azi2} direction (C{bool}).
        '''
        return self._reverse2

    @reverse2.setter  # PYCHOK setter!
    def reverse2(self, reverse2):
        '''Set the direction for C{azi2} (C{bool}).

           @arg reverse2: If C{True} reverse C{azi2} (C{bool}).
        '''
        reverse2 = bool(reverse2)
        self._update(reverse2 != self.reverse2)
        self._reverse2 = reverse2

    @property_RO
    def status(self):
        '''Get the most recent C{GeodSolve} return code (C{int}, C{str})
           or C{None}.
        '''
        return self._status

    def toStr(self, prec=6, sep=_COMMASPACE_, **unused):  # PYCHOK signature
        '''Return this C{GeodesicSolve} as string.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Optional separator to join (C{str}).

           @return: GeodesicSolve items (C{str}).
        '''
        d = dict(ellipsoid=self.ellipsoid, GeodSolve=self.GeodSolve,
                 invokation=self.invokation, status=self.status)
        return sep.join(pairs(d, prec=prec))

    @Property
    def unroll(self):
        '''Get the C{lon2} unroll'ing (C{bool}).
        '''
        return self._unroll

    @unroll.setter  # PYCHOK setter!
    def unroll(self, unroll):
        '''Set unroll'ing for C{lon2} (C{bool}).

           @arg unroll: If C{True} unroll C{lon2},
                        otherwise don't (C{bool}).
        '''
        unroll = bool(unroll)
        self._update(unroll != self.unroll)
        self._unroll = unroll

    @Property_RO
    def _u_option(self):
        return '-u' if self.unroll else ()

    @Property
    def verbose(self):
        '''Get the C{verbose} option (C{bool}).
        '''
        return self._verbose

    @verbose.setter  # PYCHOK setter!
    def verbose(self, verbose):
        '''Set the C{verbose} option.

          @arg verbose: Print a message around each
                        C{GeodSolve} invokation (C{bool}).
        '''
        self._verbose = bool(verbose)

    @property_RO
    def version(self):
        '''Get the result of C{"GeodSolve --version"}.
        '''
        return self.invoke('--version')


class GeodesicSolve(_GeodesicSolveBase):
    '''Wrapper to invoke I{Karney}'s U{GeodSolve<https://GeographicLib.SourceForge.io/html/GeodSolve.1.html>}
       as an C{Exact} version of I{Karney}'s Python class U{Geodesic<https://GeographicLib.SourceForge.io/html/
       python/code.html#geographiclib.geodesic.Geodesic>}.

       @note: Use property C{GeodSolve} or env variable C{PYGEODESY_GEODSOLVE} to specify the (fully
              qualified) path to the C{GeodSolve} executable.

       @note: This C{geodesic} is intended I{for testing purposes only}, it invokes the C{GeodSolve}
              executable for I{every} method call.
    '''
    def __init__(self, a_ellipsoid=Ellipsoids.WGS84, f=None, name=NN):
        '''New L{GeodesicSolve} instance.

           @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}) or datum (L{datum}), or
                             the equatorial radius of the ellipsoid (C{meter}).
           @arg f: The flattening of the ellipsoid (C{scalar}) if B{C{a_ellipsoid}}
                   is specified as C{meter}.
           @kwarg name: Optional name (C{str}).
        '''
        if a_ellipsoid in (GeodesicSolve._E, None):
            pass  # ignore f, default WGS84
        elif f is None:
            self._E = _ellipsoidal_datum(a_ellipsoid, name=name, raiser=True).ellipsoid
        else:
            self._E =  Ellipsoid2(a_ellipsoid, f, name=name)

        if name:
            self.name = name

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
        from pygeodesy.geodesicx.gxarea import GeodesicAreaExact
        gaX = GeodesicAreaExact(self, polyline=polyline,
                                      name=name or self.name)
        if self.verbose or self.debug:  # PYCHOK no cover
            gaX.verbose = True
        return gaX

    Polygon = Area  # for C{geographiclib} compatibility

    @Property_RO
    def _cmdDirect(self):
        '''(INTERNAL) Get the C{GeodSolve} I{Direct} cmd (C{tuple}).
        '''
        return self._cmdBasic

    @Property_RO
    def _cmdInverse(self):
        '''(INTERNAL) Get the C{GeodSolve} I{Inverse} cmd (C{tuple}).
        '''
        return self._cmdBasic + ('-i',)

    def Direct(self, lat1, lon1, azi1, s12, *unused):
        '''Return the C{Direct} result.
        '''
        return self._GDictInvoke(self._cmdDirect, True, lat1, lon1, azi1, s12)

    def Direct3(self, lat1, lon1, azi1, s12):  # PYCHOK outmask
        '''Return the destination lat, lon and reverse azimuth
           (final bearing) in C{degrees}.

           @return: L{Destination3Tuple}C{(lat, lon, final)}.
        '''
        d = self._GDictInvoke(self._cmdDirect, False, lat1, lon1, azi1, s12)
        return Destination3Tuple(float(d.lat2), float(d.lon2), wrap360(d.azi2))

    def _GDictDirect(self, lat, lon, azi, arcmode, s12_a12, *unused):  # for .geodesicx.gxarea
        '''(INTERNAL) Get C{_GenDirect}-like result as C{GDict}.
        '''
        if arcmode:
            raise GeodesicError(arcmode=arcmode, txt=str(NotImplemented))
        return self._GDictInvoke(self._cmdDirect, True, lat, lon, azi, s12_a12)

    def _GDictInverse(self, lat1, lon1, lat2, lon2, *unused):  # for .geodesicx.gxarea
        '''(INTERNAL) Get C{_GenInverse}-like result as C{GDict}, but
           I{without} C{_SALPs_CALPs_}.
        '''
        return self._GDictInvoke(self._cmdInverse, True, lat1, lon1, lat2, lon2)

    def Inverse(self, lat1, lon1, lat2, lon2, *unused):
        '''Return the C{Inverse} result.
        '''
        return self._GDictInvoke(self._cmdInverse, True, lat1, lon1, lat2, lon2)

    def Inverse1(self, lat1, lon1, lat2, lon2, wrap=False):
        '''Return the non-negative, I{angular} distance in C{degrees}.
        '''
        # see .FrechetKarney.distance, .HausdorffKarney._distance
        # and .HeightIDWkarney._distances
        _, lon2 = unroll180(lon1, lon2, wrap=wrap)  # self.LONG_UNROLL
        d = self._GDictInvoke(self._cmdInverse, False, lat1, lon1, lat2, lon2)
        # XXX self.DISTANCE needed for 'a12'?
        return abs(float(d.a12))

    def Inverse3(self, lat1, lon1, lat2, lon2):  # PYCHOK outmask
        '''Return the distance in C{meter} and the forward and
           reverse azimuths (initial and final bearing) in C{degrees}.

           @return: L{Distance3Tuple}C{(distance, initial, final)}.
        '''
        d = self._GDictInvoke(self._cmdInverse, False, lat1, lon1, lat2, lon2)
        return Distance3Tuple(float(d.s12), wrap360(d.azi1), wrap360(d.azi2))

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
                 <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GeodesicExact.html>}
                 and Python U{Geodesic.Line<https://GeographicLib.SourceForge.io/html/python/code.html>}.
        '''
        return GeodesicLineSolve(self, lat1, lon1, azi1, caps=caps, name=self.name)

    _Line = Line


class GeodesicLineSolve(_GeodesicSolveBase):
    '''Wrapper to invoke I{Karney}'s U{GeodSolve<https://GeographicLib.SourceForge.io/html/GeodSolve.1.html>}
       as an C{Exact} version of I{Karney}'s Python class U{GeodesicLine<https://GeographicLib.SourceForge.io/html/
       python/code.html#geographiclib.geodesicline.GeodesicLine>}.

       @note: Use property C{GeodSolve} or env variable C{PYGEODESY_GEODSOLVE} to specify the (fully
              qualified) path to the C{GeodSolve} executable.

       @note: This C{geodesic} is intended I{for testing purposes only}, it invokes the C{GeodSolve}
              executable for I{every} method call.
    '''
    _caps =  0
    _gS   =  None  # Solve only
    _lla1 = ()

    def __init__(self, geodesic, lat1, lon1, azi1, caps=Caps.ALL, name=NN):
        '''New L{GeodesicLineSolve} instance, allowing points to be found along
           a geodesic starting at C{(B{lat1}, B{lon1})} with azimuth B{C{azi1}}.

           @arg geodesic: The geodesic to use (L{GeodesicSolve}).
           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi1: Azimuth at the first points (compass C{degrees}).
           @kwarg caps: Bit-or'ed combination of L{Caps} values specifying
                        the capabilities the L{GeodesicLineSolve} instance
                        should possess, always C{Caps.ALL}.
           @kwarg name: Optional name (C{str}).

           @raise TypeError: Invalid B{C{geodesic}}.
        '''
        _xinstanceof(GeodesicSolve, geodesic=geodesic)

        self._gS   = gS = geodesic
        self._lla1 = lat1, lon1, azi1
        self._caps = caps | Caps._LINE

        n = name or gS.name
        if n:
            self.name = n

        self._debug   = gS._debug

        self.Exact    = gS.Exact
        self.prec     = gS.prec
        self.reverse2 = gS.reverse2
        self.unroll   = gS.unroll
        self.verbose  = gS.verbose
        try:
            self.GeodSolve = gS.GeodSolve
        except GeodesicError:
            pass

    def ArcPosition(self, a12, *unused):
        '''Find the position on the line given B{C{a12}}.

           @arg a12: Spherical arc length from the first point to the
                     second point (C{degrees}).

           @return: A C{dict} with 12 items C{lat1, lon1, azi1, lat2, lon2,
                    azi2, m12, a12, s12, M12, M21, S12}.
        '''
        return self._GDictInvoke(self._cmdArc, True, a12)

    @Property_RO
    def azi1(self):
        '''Get the azimuth at the first point (compass C{degrees}).
        '''
        return self._lla1[2]

    @Property_RO
    def azi1_sincos2(self):
        '''Get the sine and cosine of the first point's azimuth (2-tuple C{(sin, cos)}).
        '''
        return sincos2d(self.azi1)

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
        return _all_caps(self.caps, caps)

    @Property_RO
    def _cmdArc(self):
        '''(INTERNAL) Get the C{GeodSolve} I{-a -L} cmd (C{tuple}).
        '''
        return self._cmdDistance + ('-a',)

    @Property_RO
    def _cmdDistance(self):
        '''(INTERNAL) Get the C{GeodSolve} I{-L} cmd (C{tuple}).
        '''
        return self._cmdBasic + ('-L',) + strs(self._lla1, prec=DIG, fmt=Fmt.F)

    @Property_RO
    def ellipsoid(self):
        '''Get the ellipsoid (C{Ellipsoid}).
        '''
        return self._gS.ellipsoid

    @Property_RO
    def lat1(self):
        '''Get the latitude of the first point (C{degrees}).
        '''
        return self._lla1[0]

    @Property_RO
    def lon1(self):
        '''Get the longitude of the first point (C{degrees}).
        '''
        return self._lla1[1]

    def Position(self, s12, *unused):
        '''Find the position on the line given B{C{s12}}.

           @arg s12: Distance from the first point to the second C({meter}).

           @return: A C{dict} with 12 items C{lat1, lon1, azi1, lat2, lon2,
                    azi2, m12, a12, s12, M12, M21, S12}, possibly C{a12=NAN}.
        '''
        return self._GDictInvoke(self._cmdDistance, True, s12)


__all__ += _ALL_DOCS(_GeodesicSolveBase)

if __name__ == '__main__':

    gS = GeodesicSolve(name='Test')
    if gS.GeodSolve in (_PYGEODESY_GEODSOLVE_, None):  # not set
        gS.GeodSolve = '/opt/local/Cellar/geographiclib/1.51/bin/GeodSolve'  # HomeBrew
    # gS.verbose = True
    printf('version: %s', gS.version)

    r = gS.Direct(40.6, -73.8, 51, 5.5e6)
    printf('Direct: %r', r, nl=1)  # GDict(M12=0.650911, M21=0.651229, S12=39735075134877.09375, a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, m12=4844148.703101, s12=5500000.0)
    printf('Direct3: %r', gS.Direct3(40.6, -73.8, 51, 5.5e6))  # Destination3Tuple(lat=51.884565, lon=-1.141173, final=107.189397)

    printf('Inverse: %r',  gS.Inverse( 40.6, -73.8, 51.6, -0.5), nl=1)  # GDict(M12=0.64473, M21=0.645046, S12=40041368848742.53125, a12=49.94131, azi1=51.198883, azi2=107.821777, lat1=40.6, lat2=51.6, lon1=-73.8, lon2=-0.5, m12=4877684.602706, s12=5551759.400319)
    printf('Inverse1: %r', gS.Inverse1(40.6, -73.8, 51.6, -0.5))  # 49.94131021789904
    printf('Inverse3: %r', gS.Inverse3(40.6, -73.8, 51.6, -0.5))  # Distance3Tuple(distance=5551759.400319, initial=51.198883, final=107.821777)

    glS = GeodesicLineSolve(gS, 40.6, -73.8, 51)
    p = glS.Position(5.5e6)
    printf('Position:    %s  %r', p == r, p, nl=1)
    p = glS.ArcPosition(49.475527)
    printf('ArcPosition: %s %r', p == r, p)

# % python3 -m pygeodesy.geodsolve
# version: /opt/local/Cellar/geographiclib/1.51/bin/GeodSolve: GeographicLib version 1.51
#
# Direct: GDict(M12=0.650911, M21=0.651229, S12=39735075134877.09375, a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, m12=4844148.703101, s12=5500000.0)
# Direct3: Destination3Tuple(lat=51.884565, lon=-1.141173, final=107.189397)
#
# Inverse: GDict(M12=0.64473, M21=0.645046, S12=40041368848742.53125, a12=49.94131, azi1=51.198883, azi2=107.821777, lat1=40.6, lat2=51.6, lon1=-73.8, lon2=-0.5, m12=4877684.602706, s12=5551759.400319)
# Inverse1: 49.94131021789904
# Inverse3: Distance3Tuple(distance=5551759.400319, initial=51.198883, final=107.821777)
#
# Position:    True  GDict(M12=0.650911, M21=0.651229, S12=39735075134877.09375,  a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141173, m12=4844148.703101, s12=5500000.0)
# ArcPosition: False GDict(M12=0.650911, M21=0.651229, S12=39735074737272.734375, a12=49.475527, azi1=51.0, azi2=107.189397, lat1=40.6, lat2=51.884565, lon1=-73.8, lon2=-1.141174, m12=4844148.669561, s12=5499999.948497)

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
