
# -*- coding: utf-8 -*-

u'''(INTERNAL) Module with base classes for L{pygeodesy.geodsolve} and L{pygeodesy.rhumbsolve}.
'''

from pygeodesy.basics import map2, ub2str, _zip
from pygeodesy.errors import _AssertionError, _xkwds_get
from pygeodesy.interns import DIG, NN, _0_, _BACKSLASH_, _COMMASPACE_, \
                             _enquote, _EQUAL_, _not_, _SPACE_
from pygeodesy.karney import Caps, _CapsBase, _ellipsoid, _EWGS84, GDict, \
                             Precision_, unroll180
from pygeodesy.lazily import _ALL_DOCS, printf, _sys_version_info2
from pygeodesy.named import callername, notOverloaded
from pygeodesy.props import Property, Property_RO, property_RO, _update_all
from pygeodesy.streprs import Fmt, fstr, fstrzs, pairs, strs
# from pygeodesy.units import Precision_  # from .karney
# from pygeodesy.utily import unroll180  # from .karney

from subprocess import PIPE as _PIPE, Popen as _Popen, STDOUT as _STDOUT

__all__ = ()  # nothing public
__version__ = '22.08.02'

_Error_    = 'Error'
_ERROR_    = 'ERROR'
_text_True =  dict() if _sys_version_info2 < (3, 7) else dict(text=True)


def _cmd_stdin_(cmd, stdin):  # PYCHOK no cover
    '''(INTERNAL) Cmd line, stdin and caller as sC{str}.
    '''
    c = Fmt.PAREN(callername(up=3))
    t = (c,) if stdin is None else (_BACKSLASH_, str(stdin), c)
    return _SPACE_.join(cmd + t)


def _popen2(cmd, stdin=None):  # in .mgrs, .test.base, .test.testMgrs
    '''(INTERNAL) Invoke C{B{cmd} tuple} and return C{exitcode}
       and all output to C{stdout/-err}.
    '''
    p = _Popen(cmd, creationflags=0,
                  # executable=sys.executable, shell=True,
                    stdin=_PIPE, stdout=_PIPE, stderr=_STDOUT,
                 **_text_True)  # PYCHOK kwArgs
    r = p.communicate(stdin)[0]
    return p.returncode, ub2str(r).strip()


class _SolveLineSolveBase(_CapsBase):
    '''(NTERNAL) Base class for C{_Solve} and C{_LineSolve}.
    '''
    _E             = _EWGS84
    _Error         =  None
    _Exact         =  True
    _invokation    =  0
    _Names_Direct  = \
    _Names_Inverse = ()
    _prec          =  Precision_(prec=DIG)
    _reverse2      =  False
    _Solve_name    =  NN  # executable basename
    _Solve_path    =  NN  # executable path
    _status        =  None
    _unroll        =  False
    _verbose       =  False

    @Property_RO
    def a(self):
        '''Get the I{equatorial} radius, semi-axis (C{meter}).
        '''
        return self.ellipsoid.a

    @property_RO
    def _cmdBasic(self):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self)

    @Property_RO
    def ellipsoid(self):
        '''Get the ellipsoid (C{Ellipsoid}).
        '''
        return self._E

    @Property_RO
    def _e_option(self):
        E = self.ellipsoid
        if E is _EWGS84:
            return ()  # default
        a, f = strs(E.a_f, fmt=Fmt.F, prec=DIG + 3)  # not .G!
        return ('-e', a, f)

    @property
    def Exact(self):
        '''Get the Solve's C{exact} setting (C{bool}).
        '''
        return self._Exact

    @Exact.setter  # PYCHOK setter!
    def Exact(self, Exact):
        '''Set the Solve's C{exact} setting (C{bool}),
           if C{True} use I{exact} version.
        '''
        Exact = bool(Exact)
        if self._Exact != Exact:
            _update_all(self)
            self._Exact = Exact

    @Property_RO
    def f(self):
        '''Get the ellipsoid's I{flattening} (C{float}), M{(a - b) / a}, C{0} for spherical, negative for prolate.
        '''
        return self.ellipsoid.f

    def _GDictInvoke(self, cmd, floats, Names, *args):
        '''(INTERNAL) Invoke C{Solve}, return results as C{GDict}.
        '''
        N = len(Names)
        if N < 1:
            raise _AssertionError(cmd=cmd, Names=Names)
        i = fstr(args, prec=DIG, fmt=Fmt.F, sep=_SPACE_) if args else None  # not Fmt.G!
        t = self._invoke(cmd, stdin=i).lstrip().split()  # 12-/+ tuple
        if len(t) > N:  # PYCHOK no cover
            # unzip instrumented name=value pairs to names and values
            n, v = _zip(*(p.split(_EQUAL_) for p in t[:-N]))  # strict=True
            v += tuple(t[-N:])
            n += Names
        else:
            n, v = Names, t
        if self.verbose:  # PYCHOK no cover
            self._print(_COMMASPACE_.join(map(Fmt.EQUAL, n, map(fstrzs, v))))
        if floats:
            v = map(float, v)
        r = GDict(_zip(n, v))  # strict=True
        return self._iter2tion(r, r)

    @property_RO
    def invokation(self):
        '''Get the most recent C{Solve} invokation number (C{int}).
        '''
        return self._invokation

    def invoke(self, *options, **stdin):
        '''Invoke the C{Solve} executable and return the result.

           @arg options: No, one or several C{Solve} command line
                         options (C{str}s).
           @kwarg stdin: Optional input to pass to C{Solve.stdin} (C{str}).

           @return: The C{Solve.stdout} and C{.stderr} output (C{str}).

           @raise GeodesicError: On any error, including a non-zero return
                                 code from C{GeodSolve}.

           @raise RhumbError: On any error, including a non-zero return code
                              from C{RhumbSolve}.

           @note: The C{Solve} return code is in property L{status}.
        '''
        c = (self._Solve_path,) + map2(str, options)
        i = _xkwds_get(stdin, stdin=None)
        r =  self._invoke(c, stdin=i)
        s =  self.status
        if s:
            raise self._Error(cmd=_cmd_stdin_(c, i), status=s,
                              txt=_not_(_0_))
        if self.verbose:  # PYCHOK no cover
            self._print(r)
        return r

    def _invoke(self, cmd, stdin=None):
        '''(INTERNAL) Invoke the C{Solve} executable, with the
           given B{C{cmd}} line and optional input to B{C{stdin}}.
        '''
        self._invokation += 1
        self._status = t = None
        if self.verbose:  # PYCHOK no cover
            t = _cmd_stdin_(cmd, stdin)
            self._print(t)
        try:  # invoke and write to stdin
            s, r = _popen2(cmd, stdin)
            if len(r) < 6 or r[:5] in (_Error_, _ERROR_):
                raise ValueError(r)
        except (IOError, OSError, TypeError, ValueError) as x:
            raise self._Error(cmd=t or _cmd_stdin_(cmd, stdin),
                              txt=str(x))
        self._status = s
        return r

    @property_RO
    def _p_option(self):
        return '-p', str(self.prec - 5)  # -p is distance prec

    @Property
    def prec(self):
        '''Get the precision, number of (decimal) digits (C{int}).
        '''
        return self._prec

    @prec.setter  # PYCHOK setter!
    def prec(self, prec):
        '''Set the precision for C{angles} in C{degrees}, like C{lat}, C{lon},
           C{azimuth} and C{arc} in number of decimal digits (C{int}, C{0}..L{DIG}).

           @note: The precision for C{distance = B{prec} - 5} or up to
                  10 decimal digits for C{nanometer} and for C{area =
                  B{prec} - 12} or at most C{millimeter} I{squared}.
        '''
        prec = Precision_(prec=prec, high=DIG)
        if self._prec != prec:
            _update_all(self)
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
        '''Set the direction for C{azi2} (C{bool}), if C{True} reverse C{azi2}.
        '''
        reverse2 = bool(reverse2)
        if self._reverse2 != reverse2:
            _update_all(self)
            self._reverse2 = reverse2

    def _setSolve(self, path, **Solve_path):
        '''(INTERNAL) Set the executable C{path}.
        '''
        hold = self._Solve_path
        if hold != path:
            _update_all(self)
            self._Solve_path = path
        try:
            _ = self.version  # test path and ...
            if self.status:  # ... return code
                S_p = Solve_path or {self._Solve_name: _enquote(path)}
                raise self._Error(status=self.status, txt=_not_(_0_), **S_p)
            hold = path
        finally:  # restore in case of error
            if self._Solve_path != hold:
                _update_all(self)
                self._Solve_path = hold

    @property_RO
    def status(self):
        '''Get the most recent C{Solve} return code (C{int}, C{str})
           or C{None}.
        '''
        return self._status

    @Property
    def unroll(self):
        '''Get the C{lon2} unroll'ing (C{bool}).
        '''
        return self._unroll

    @unroll.setter  # PYCHOK setter!
    def unroll(self, unroll):
        '''Set unroll'ing for C{lon2} (C{bool}), if C{True} unroll C{lon2}, otherwise don't.
        '''
        unroll = bool(unroll)
        if self._unroll != unroll:
            _update_all(self)
            self._unroll = unroll

    @property
    def verbose(self):
        '''Get the C{verbose} option (C{bool}).
        '''
        return self._verbose

    @verbose.setter  # PYCHOK setter!
    def verbose(self, verbose):
        '''Set the C{verbose} option (C{bool}), C{True} prints
           a message around each C{RhumbSolve} invokation.
        '''
        self._verbose = bool(verbose)

    @Property_RO
    def version(self):
        '''Get the result of C{"GeodSolve --version"} or C{"RhumbSolve --version"}.
        '''
        return self.invoke('--version')


class _SolveBase(_SolveLineSolveBase):
    '''(NTERNAL) Base class for C{_GeodesicSolveBase} and C{_RhumbSolveBase}.
    '''
    def __init__(self, a_ellipsoid=_EWGS84, f=None, path=NN, name=NN):
        '''New C{Solve} instance.

           @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}) or datum (L{Datum}) or
                             the equatorial radius of the ellipsoid (C{scalar},
                             conventionally in C{meter}), see B{C{f}}.
           @arg f: The flattening of the ellipsoid (C{scalar}) if B{C{a_ellipsoid}}
                   is specified as C{scalar}.
           @kwarg path: Optionally, the (fully qualified) path to the C{GeodSolve}
                        or C{RhumbSolve} executable (C{filename}).
           @kwarg name: Optional name (C{str}).

           @raise TypeError: Invalid B{C{a_ellipsoid}} or B{C{f}}.
        '''
        if a_ellipsoid not in (self._E, None):  # NOT self.ellipsoid
            self._E = _ellipsoid(a_ellipsoid, f, name=name)
        if name:
            self.name = name
        if path:
            self._setSolve(path)

    @Property_RO
    def _cmdDirect(self):
        '''(INTERNAL) Get the C{Solve} I{Direct} cmd (C{tuple}).
        '''
        return self._cmdBasic

    @Property_RO
    def _cmdInverse(self):
        '''(INTERNAL) Get the C{Solve} I{Inverse} cmd (C{tuple}).
        '''
        return self._cmdBasic + ('-i',)

    def Direct(self, lat1, lon1, azi1, s12, *unused):
        '''Return the C{Direct} result.
        '''
        return self._GDictDirect(lat1, lon1, azi1, False, s12)

    def _GDictDirect(self, lat, lon, azi, arcmode, s12_a12, *unused, **floats):  # for .geodesicx.gxarea
        '''(INTERNAL) Get C{_GenDirect}-like result as C{GDict}.
        '''
        if arcmode:
            raise self._Error(arcmode=arcmode, txt=str(NotImplemented))
        floats = _xkwds_get(floats, floats=True)
        return self._GDictInvoke(self._cmdDirect, floats, self._Names_Direct,
                                                  lat, lon, azi, s12_a12)

    def _GDictInverse(self, lat1, lon1, lat2, lon2, *unused, **floats):  # for .geodesicx.gxarea
        '''(INTERNAL) Get C{_GenInverse}-like result as C{GDict}, but
           I{without} C{_SALPs_CALPs_}.
        '''
        floats = _xkwds_get(floats, floats=True)
        return self._GDictInvoke(self._cmdInverse, floats, self._Names_Inverse,
                                                   lat1, lon1, lat2, lon2)

    def Inverse(self, lat1, lon1, lat2, lon2, *unused):
        '''Return the C{Inverse} result.
        '''
        return self._GDictInverse(lat1, lon1, lat2, lon2)

    def Inverse1(self, lat1, lon1, lat2, lon2, wrap=False):
        '''Return the non-negative, I{angular} distance in C{degrees}.
        '''
        # see .FrechetKarney.distance, .HausdorffKarney._distance
        # and .HeightIDWkarney._distances
        _, lon2 = unroll180(lon1, lon2, wrap=wrap)  # self.LONG_UNROLL
        r = self._GDictInverse(lat1, lon1, lat2, lon2, floats=False)
        # XXX self.DISTANCE needed for 'a12'?
        return abs(float(r.a12))

    def _toStr(self, prec=6, sep=_COMMASPACE_, **Solve):  # PYCHOK signature
        '''(INTERNAL) Return this C{_Solve} as string..
        '''
        d = dict(ellipsoid=self.ellipsoid, invokation=self.invokation,
                                           status=self.status, **Solve)
        return sep.join(pairs(d, prec=prec))


class _LineSolveBase(_SolveLineSolveBase):
    '''(NTERNAL) Base class for C{GeodesicLineSolve} and C{RhumbLineSolve}.
    '''
#   _caps  =  0
#   _lla1  = {}
    _solve =  None  # L{GeodesicSolve} or L{RhumbSolve} instance

    def __init__(self, solve, lat1, lon1, caps, name, **azi):
        self._caps  = caps | Caps._LINE
        self._debug = solve._debug & Caps._DEBUG_ALL
        self._lla1  = GDict(lat1=lat1, lon1=lon1, **azi)
        self._solve = solve

        n = name or solve.name
        if n:
            self.name = n

    @Property_RO
    def _cmdDistance(self):
        '''(INTERNAL) Get the C{GeodSolve} I{-L} cmd (C{tuple}).
        '''
        def _lla3(lat1=0, lon1=0, **azi):
            _, azi = azi.popitem()
            return lat1, lon1, azi

        t = strs(_lla3(**self._lla1), prec=DIG, fmt=Fmt.F)  # self._solve.prec
        return self._cmdBasic + ('-L',) + t

    @property_RO
    def ellipsoid(self):
        '''Get the ellipsoid (C{Ellipsoid}).
        '''
        return self._solve.ellipsoid

    @Property_RO
    def lat1(self):
        '''Get the latitude of the first point (C{degrees}).
        '''
        return self._lla1.lat1

    @Property_RO
    def lon1(self):
        '''Get the longitude of the first point (C{degrees}).
        '''
        return self._lla1.lon1

    def _toStr(self, prec=6, sep=_COMMASPACE_, **solve):  # PYCHOK signature
        '''(INTERNAL) Return this C{_LineSolve} as string..
        '''
        d = dict(ellipsoid=self.ellipsoid, invokation=self._solve.invokation,
                                           lat1=self.lat1, lon1=self.lon1,
                                           status=self._solve.status, **solve)
        return sep.join(pairs(d, prec=prec))


__all__ += _ALL_DOCS(_SolveBase, _LineSolveBase, _SolveLineSolveBase)

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
