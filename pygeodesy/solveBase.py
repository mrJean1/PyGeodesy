
# -*- coding: utf-8 -*-

u'''(INTERNAL) Private base classes for L{pygeodesy.geodsolve} and L{pygeodesy.rhumb.solve}.
'''

from pygeodesy.basics import clips, map2, _zip
from pygeodesy.constants import DIG
from pygeodesy.datums import _earth_datum, _WGS84,  _EWGS84
# from pygeodesy.ellipsoids import _EWGS84  # from .datums
from pygeodesy.errors import _AssertionError, _xkwds_get, _xkwds_get1, \
                             _xkwds_item2
from pygeodesy.internals import _enquote, _popen2, printf
from pygeodesy.interns import NN, _0_, _AT_,_BACKSLASH_, _COLONSPACE_, \
                             _COMMASPACE_, _EQUAL_, _Error_, _SPACE_, \
                             _UNUSED_
from pygeodesy.karney import Caps, _CapsBase, GDict
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import callername, _name2__, notOverloaded
from pygeodesy.props import Property, Property_RO, property_RO, _update_all
from pygeodesy.streprs import Fmt, fstr, fstrzs, pairs, strs
from pygeodesy.units import Precision_
from pygeodesy.utily import unroll180,  wrap360  # PYCHOK shared

__all__ = _ALL_LAZY.solveBase
__version__ = '24.10.13'

_ERROR_ = 'ERROR'


def _cmd_stdin_(cmd, stdin):  # PYCHOK no cover
    '''(INTERNAL) Cmd line, stdin and caller as sC{str}.
    '''
    if stdin is not None:
        cmd += _BACKSLASH_, str(stdin)
    cmd += Fmt.PAREN(callername(up=3)),
    return _SPACE_.join(cmd)


# def _float_int(r):
#     '''(INTERNAL) Convert result into C{float} or C{int}.
#     '''
#     f = float(r)
#     i = int(f)
#     return i if float(i) == f else f  # PYCHOK inconsistent


class _SolveCapsBase(_CapsBase):
    '''(NTERNAL) Base class for C{_SolveBase} and C{_LineSolveBase}.
    '''
    _datum      = _WGS84
    _Error      =  None
    _Exact      =  True
    _invokat    = _AT_
    _invokation =  0
    _linelimit  =  0
    _prec       =  Precision_(prec=DIG)
    _prec2stdin =  DIG
    _Xable_name =  NN  # executable basename
    _Xable_path =  NN  # executable path
    _status     =  None
    _verbose    =  False

    @Property_RO
    def a(self):
        '''Get the I{equatorial} radius, semi-axis (C{meter}).
        '''
        return self.ellipsoid.a

    @property_RO
    def _cmdBasic(self):  # PYCHOK no covers        '''(INTERNAL) I{Must be overloaded}.'''
        notOverloaded(self, underOK=True)

    @property_RO
    def datum(self):
        '''Get the datum (C{Datum}).
        '''
        return self._datum

    def _Dict(self, Dict, n, v, floats=True, **unused):
        if self.verbose:  # PYCHOK no cover
            self._print(_COMMASPACE_.join(map(Fmt.EQUAL, n, map(fstrzs, v))))
        if floats:
            v = map(float, v)  # _float_int, see Intersectool._XDistInvoke
        return Dict(_zip(n, v))  # strict=True

    def _DictInvoke2(self, cmd, args, Names, Dict, **floats_R):
        '''(INTERNAL) Invoke C{Solve}, return results as C{Dict}.
        '''
        N = len(Names)
        if N < 1:
            raise _AssertionError(cmd=cmd, Names=Names)
        i = fstr(args, prec=self._prec2stdin, fmt=Fmt.F, sep=_SPACE_) if args else None  # NOT Fmt.G!
        t = self._invoke(cmd, stdin=i, **floats_R).lstrip().split()  # 12-/++ tuple
        if _xkwds_get(floats_R, _R=None):  # == '-R' in cmd
            return self._Dicts(Dict, Names, t, **floats_R), True
        elif len(t) > N:  # PYCHOK no cover
            # unzip instrumented name=value pairs to names and values
            n, v = _zip(*(p.split(_EQUAL_) for p in t[:-N]))  # strict=True
            v += tuple(t[-N:])
            n += Names
        else:
            n, v = Names, t
        r = self._Dict(Dict, n, t, **floats_R)
        return self._iter2tion(r, **r), None

    def _Dicts(self, Dict, Names, t, **floats_R):
        i, N = 0, len(Names)
        for x in range(0, len(t), N):
            if t[x] == 'nan':
                break
            X = self._Dict(Dict, Names, t[x:x + N], **floats_R)
            yield X.set_(iteration=i)
            i += 1

    @Property_RO
    def _E_option(self):
        return ('-E',) if self.Exact else ()

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
    def ellipsoid(self):
        '''Get the ellipsoid (C{Ellipsoid}).
        '''
        return self.datum.ellipsoid

    @Property_RO
    def _e_option(self):
        E = self.ellipsoid
        if E is _EWGS84:
            return ()  # default
        a, f = strs(E.a_f, fmt=Fmt.F, prec=DIG + 3)  # not .G!
        return ('-e', a, f)

    @Property_RO
    def flattening(self):
        '''Get the C{ellipsoid}'s I{flattening} (C{scalar}), M{(a - b) / a},
           C{0} for spherical, negative for prolate.
        '''
        return self.ellipsoid.f

    f = flattening

    def invokat(self, *prefix):
        '''Get and set the invokation number C{"@"} prefix (C{str}).

           @return: Previous prefix (C{str}).
        '''
        p = self._invokat
        if prefix:
            set._invokat = str(prefix[0])
        return p

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
        c = (self._Xable_path,) + map2(str, options)  # map2(_enquote, options)
        i = _xkwds_get1(stdin, stdin=None)
        r =  self._invoke(c, stdin=i)
        s =  self.status
        if s:
            raise self._Error(cmd=_cmd_stdin_(c, i), status=s,
                              txt_not_=_0_)
        if self.verbose:  # PYCHOK no cover
            self._print(r)
        return r

    def _invoke(self, cmd, stdin=None, **unused):  # _R=None
        '''(INTERNAL) Invoke the C{Solve} executable, with the
           given B{C{cmd}} line and optional input to B{C{stdin}}.
        '''
        self._invokation += 1
        self._status = t = None
        if self.verbose:  # PYCHOK no cover
            t = _cmd_stdin_(cmd, stdin)
            self._print(t)
        try:  # invoke and write to stdin
            r, s = _popen2(cmd, stdin)
            if len(r) < 6 or r[:5] in (_Error_, _ERROR_):
                raise ValueError(r)
        except (IOError, OSError, TypeError, ValueError) as x:
            raise self._Error(cmd=t or _cmd_stdin_(cmd, stdin), cause=x)
        self._status = s
        if self.verbose:  # and _R is None:  # PYCHOK no cover
            self._print(repr(r), 'stdout/-err')
        return r

    def linelimit(self, *limit):
        '''Set and get the print line length limit.

           @arg limit: New line limit (C{int}) or C{0}
                       or C{None} for unlimited.

           @return: Teh previous limit (C{int}).
        '''
        n = self._linelimit
        if limit:
            m = int(limit[0] or 0)
            self._linelimit = max(80, m) if m > 0 else (n if m < 0 else 0)
        return n

    @Property_RO
    def _mpd(self):  # meter per degree
        return self.ellipsoid._Lpd

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

    def _print(self, line, *suffix):  # PYCHOK no cover
        '''(INTERNAL) Print a status line.
        '''
        if self._linelimit:
            line =  clips(line, limit=self._linelimit, length=True)
        if self.status is not None:
            s    = _COMMASPACE_(self.status, *suffix)
            line = _SPACE_(line, Fmt.PAREN(s))
        p = NN(self.named2, self._invokat, self.invokation)
        printf(_COLONSPACE_(p, line))

    def _setXable(self, path, **Xable_path):
        '''(INTERNAL) Set the executable C{path}.
        '''
        hold = self._Xable_path
        if hold != path:
            _update_all(self)
            self._Xable_path = path
        try:
            _ = self.version  # test path and ...
            if self.status:  # ... return code
                S_p = Xable_path or {self._Xable_name: _enquote(path)}
                raise self._Error(status=self.status, txt_not_=_0_, **S_p)
            hold = path
        finally:  # restore in case of error
            if self._Xable_path != hold:
                _update_all(self)
                self._Xable_path = hold

    @property_RO
    def status(self):
        '''Get the most recent C{Solve} return code (C{int}, C{str})
           or C{None}.
        '''
        return self._status

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


class _SolveBase(_SolveCapsBase):
    '''(INTERNAL) Base class for C{_SolveBase} and C{_SolveLineBase}.
    '''
    _Names_Direct  = \
    _Names_Inverse = ()
    _reverse2      =  False
    _unroll        =  False

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


class _SolveGDictBase(_SolveBase):
    '''(NTERNAL) Base class for C{_GeodesicSolveBase} and C{_RhumbSolveBase}.
    '''

    def __init__(self, a_ellipsoid=_EWGS84, f=None, path=NN, **name):
        '''New C{Solve} instance.

           @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}) or datum (L{Datum}) or
                             the equatorial radius of the ellipsoid (C{scalar},
                             conventionally in C{meter}), see B{C{f}}.
           @arg f: The flattening of the ellipsoid (C{scalar}) if B{C{a_ellipsoid}}
                   is specified as C{scalar}.
           @kwarg path: Optionally, the (fully qualified) path to the C{GeodSolve}
                        or C{RhumbSolve} executable (C{filename}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise TypeError: Invalid B{C{a_ellipsoid}} or B{C{f}}.
        '''
        _earth_datum(self, a_ellipsoid, f=f, **name)
        if name:
            self.name = name
        if path:
            self._setXable(path)

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

    def Direct(self, lat1, lon1, azi1, s12, outmask=_UNUSED_):  # PYCHOK unused
        '''Return the C{Direct} result.
        '''
        return self._GDictDirect(lat1, lon1, azi1, False, s12)

    def _GDictDirect(self, lat, lon, azi, arcmode, s12_a12, outmask=_UNUSED_, **floats):  # PYCHOK for .geodesicx.gxarea
        '''(INTERNAL) Get C{_GenDirect}-like result as C{GDict}.
        '''
        if arcmode:
            raise self._Error(arcmode=arcmode, txt=str(NotImplemented))
        return self._GDictInvoke(self._cmdDirect, self._Names_Direct,
                                                  lat, lon, azi, s12_a12, **floats)

    def _GDictInverse(self, lat1, lon1, lat2, lon2, outmask=_UNUSED_, **floats):  # PYCHOK for .geodesicx.gxarea
        '''(INTERNAL) Get C{_GenInverse}-like result as C{GDict}, but I{without} C{_SALP_CALPs_}.
        '''
        return self._GDictInvoke(self._cmdInverse, self._Names_Inverse,
                                                   lat1, lon1, lat2, lon2, **floats)

    def _GDictInvoke(self, cmd,  Names, *args, **floats):
        '''(INTERNAL) Invoke C{Solve}, return results as C{Dict}.
        '''
        return self._DictInvoke2(cmd, args, Names, GDict, **floats)[0]  # _R

    def Inverse(self, lat1, lon1, lat2, lon2, outmask=_UNUSED_):  # PYCHOK unused
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


class _SolveGDictLineBase(_SolveGDictBase):
    '''(NTERNAL) Base class for C{GeodesicLineSolve} and C{RhumbLineSolve}.
    '''
#   _caps  =  0
#   _lla1  = {}
    _solve =  None  # L{GeodesicSolve} or L{RhumbSolve} instance

    def __init__(self, solve, lat1, lon1, caps, **azi_name):
        name, azi = _name2__(azi_name, _or_nameof=solve)
        if name:
            self.name = name

        self._caps  = caps | Caps._AZIMUTH_LATITUDE_LONG_UNROLL
        self._debug = solve._debug & Caps._DEBUG_ALL
        self._lla1  = GDict(lat1=lat1, lon1=lon1, **azi)
        self._solve = solve

    @Property_RO
    def _cmdDistance(self):
        '''(INTERNAL) Get the C{GeodSolve} I{-L} cmd (C{tuple}).
        '''
        def _lla3(lat1=0, lon1=0, **azi):
            _, azi = _xkwds_item2(azi)
            return lat1, lon1, azi

        t = strs(_lla3(**self._lla1), prec=DIG, fmt=Fmt.F)  # self._solve.prec
        return self._cmdBasic + ('-L',) + t

    @property_RO
    def datum(self):
        '''Get the datum (C{Datum}).
        '''
        return self._solve.datum

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


__all__ += _ALL_DOCS(_SolveBase, _SolveCapsBase, _SolveGDictBase, _SolveGDictLineBase)

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
