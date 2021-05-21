
# -*- coding: utf-8 -*-

u'''Wrapper to invoke I{Karney}'s U{GeodSolve
<https://GeographicLib.SourceForge.io/html/GeodSolve.1.htmlb>} utility,
intended I{for testing purposes only}.

Set env variable C{PYGEODESY_GEODSOLVE} to the I{fully qualified path}
of the C{GeodSolve} executable.
'''

from pygeodesy.basics import map1, map2
from pygeodesy.datums import _ellipsoidal_datum
from pygeodesy.ellipsoids import Ellipsoids, Ellipsoid2
from pygeodesy.interns import DIG, NN, _0_, _COMMASPACE_, _SPACE_
from pygeodesy.interns import _not_  # PYCHOK used!
from pygeodesy.karney import GDict, GeodesicError, GeodSolve12Tuple
from pygeodesy.lazily import _ALL_LAZY, printf, _sys
from pygeodesy.lazily import _env  # PYCHOK used!
from pygeodesy.named import callername, _NamedBase, notImplemented
from pygeodesy.namedTuples import Destination3Tuple, Distance3Tuple
from pygeodesy.props import Property, Property_RO, property_RO
from pygeodesy.streprs import Fmt, fstr, fstrzs, pairs, strs
from pygeodesy.units import Precision_
from pygeodesy.utily import unroll180

from subprocess import PIPE as _PIPE, Popen as _Popen, STDOUT as _STDOUT

__all__ = _ALL_LAZY.geodsolve
__version__ = '21.05.20'

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


class GeodesicSolve(_NamedBase):  # PYCHOK no cover
    '''Wrapper to invoke I{Karney}'s U{GeodSolve<https://geographiclib.sourceforge.io/html/GeodSolve.1.html>}
       as an C{Exact} version of I{Karney}'s Python class U{Geodesic<https://GeographicLib.SourceForge.io/html/
       python/code.html#geographiclib.geodesic.Geodesic>} I{wrapper}.

       @note: Use property C{GeodSolve} or env variable C{PYGEODESY_GEODSOLVE} to specify the path to the
              C{GeodSolve} executable.

       @note: This C{geodesic} is intended mainly for testing purposes, it invokes the C{GeodSolve} executable
              for I{each} method call.
    '''
    ALL         =  AREA = AZIMUTH = DISTANCE = DISTANCE_IN = GEODESICSCALE = \
    EMPTY       =  LATITUDE = LONGITUDE = LONG_UNROLL = REDUCEDLENGTH = STANDARD = 0
    _debug      =  0
    _E          =  Ellipsoids.WGS84
    _Exact      =  True
    _GeodSolve  = _env.get(_PYGEODESY_GEODSOLVE_, _PYGEODESY_GEODSOLVE_)
    _invokation =  0
    _map        =  True
    _prec       =  Precision_(prec=DIG)
    _reverse2   =  False
    _status     =  None
    _text_True  =  dict() if _sys.version_info[:2] < (3, 7) else dict(text=True)
    _unroll     =  False
    _verbose    =  False

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

    @Property_RO
    def a(self):
        '''Get the I{equatorial} radius, semi-axis (C{meter}).
        '''
        return self._E.a

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
        if self.debug:
            gaX.verbose = True
        return gaX

    Polygon = Area  # for C{geographiclib} compatibility

    @Property_RO
    def _b_option(self):
        return ('-b',) if self.reverse2 else ()

    @Property_RO
    def _cmdDirect(self):
        '''(INTERNAL) Get the C{GeodSolve} I{Direct} cmd (C{tuple}).
        '''
        return (self.GeodSolve,) + self._b_option \
                                 + self._e_option \
                                 + self._E_option \
                                 + self._p_option \
                                 + self._u_option + ('-f',)

    @Property_RO
    def _cmdInverse(self):
        '''(INTERNAL) Get the C{GeodSolve} I{Inverse} cmd (C{tuple}).
        '''
        return self._cmdDirect + ('-i',)

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
        from pygeodesy.geodesicx.bases import Caps
        self._debug = Caps._DEBUG_ALL if debug else 0

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
        return Destination3Tuple(*map1(float, d.lat2, d.lon2, d.azi2))

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
        return self._E.f

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

    def _GDictInvoke(self, cmd, floats, *args):
        '''(INTERNAL) Invoke C{GeodSolve}, get C{GDict}.
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
        if self.verbose:
            self._print(_COMMASPACE_.join(map(Fmt.EQUAL, n, map(fstrzs, v))))
        if floats:
            v = map(float, v)
        return GDict(zip(n, v))

    @Property
    def GeodSolve(self):
        '''Get the U{GeodSolve<https://geographiclib.sourceforge.io/html/GeodSolve.1.html>}
           executable (C{filename}).
        '''
        return self._GeodSolve

    @GeodSolve.setter  # PYCHOK setter!
    def GeodSolve(self, path):
        '''Set the U{GeodSolve<https://geographiclib.sourceforge.io/html/GeodSolve.1.html>}
           executable (C{filename}).

           @arg path: Fully qualified path to the C{GeodSolve} executable (C{str}).

           @raise GeodesicError: Invalid C{path} or C{geosolve} doesn't exist or isn't
                                 executable.
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
        return Distance3Tuple(*map1(float, d.s12, d.azi1, d.azi2))  # wrap360?

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
        if self.verbose:
            self._print(r)
        return r

    def _invoke(self, cmd, stdin=None):
        '''(INTERNAL) Invoke the C{GeodSolve} executable, with the
           given B{C{cmd}} line and optional input to B{C{stdin}}.
        '''
        self._invokation += 1
        self._status = t = None
        if self.verbose:
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

    def Line(self, *args, **kwds):
        '''Not implemented.

           @raise NotImplementedError: Always.
        '''
        notImplemented(self, self.Line.__name__, *args, **kwds)

    _Line = Line

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

           @note: The precision for C{distance=B{prec} - 5} or up to 10
                  decimal digits for C{nanometer} and for C{area=B{prec}
                  - 12} or at most C{millimeter} I{squared}.
        '''
        prec = Precision_(prec=prec, high=DIG)
        self._update(prec != self.prec)
        self._prec = prec

    def _print(self, line):
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

           @arg reverse2: Azimuth C{azi2} direction (C{bool}), C{True}
                          for I{reverse}, otherwise default I{forward}.

           @note: The precision for C{distance=B{prec} - 5}, up to 10
                  decimal digits for C{nanometer} and for C{area=B{prec}
                  - 12} or at most C{millimeter} I{squared}.
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


if __name__ == '__main__':

    gS = GeodesicSolve(name='Test')
    if gS.GeodSolve in (_PYGEODESY_GEODSOLVE_, None):  # not set
        gS.GeodSolve = '/opt/local/Cellar/geographiclib/1.51/bin/GeodSolve'  # HomeBrew
    # gS.verbose = True
    printf('version: %s', gS.version, nt=1)
    printf('Direct: %s', gS.Direct(40.6, -73.8, 51, 5.5e6))
    printf('Direct3: %s', gS.Direct3(40.6, -73.8, 51, 5.5e6))  # (51.884565, -1.141173, 107.189397)
    printf('Inverse: %s', gS.Inverse(40.6, -73.8, 51.6, -0.5))
    printf('Inverse1: %s', gS.Inverse1(40.6, -73.8, 51.6, -0.5))  # 49.94131021789904
    printf('Inverse3: %s', gS.Inverse3(40.6, -73.8, 51.6, -0.5), nt=1)  # (5551759.400319, 51.198883, 107.821777)

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
