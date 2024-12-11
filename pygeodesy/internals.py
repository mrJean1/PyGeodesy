
# -*- coding: utf-8 -*-

u'''Mostly INTERNAL functions, except L{machine}, L{print_} and L{printf}.
'''
# from pygeodesy.basics import isiterablen, ubstr  # _MODS
# from pygeodesy.errors import _AttributeError, _error_init, _UnexpectedError, _xError2  # _MODS
from pygeodesy.interns import NN, _BAR_, _COLON_, _DASH_, _DOT_, _ELLIPSIS_, _EQUALSPACED_, \
                             _immutable_, _NL_, _pygeodesy_, _PyPy__, _python_, _QUOTE1_, \
                             _QUOTE2_, _s_, _SPACE_, _sys, _UNDER_
from pygeodesy.interns import _COMMA_, _Python_  # PYCHOK used!
# from pygeodesy.streprs import anstr, pairs, unstr  # _MODS

# import os  # _MODS
# import os.path  # _MODS
# import sys as _sys  # from .interns

_0_0      =  0.0  # PYCHOK in .basics, .constants
_100_0    =  100.0  # in .constants
_arm64_   = 'arm64'
_iOS_     = 'iOS'
_macOS_   = 'macOS'
_SIsecs   = 'fs', 'ps', 'ns', 'us', 'ms', 'sec'  # reversed
_Windows_ = 'Windows'


def _DUNDER_nameof(inst, *dflt):
    '''(INTERNAL) Get the DUNDER C{.__name__} attr.
    '''
    try:
        return inst.__name__
    except AttributeError:
        pass
    return dflt[0] if dflt else inst.__class__.__name__


def _DUNDER_nameof_(*names__):  # in .errors._IsnotError
    '''(INTERNAL) Yield the _DUNDER_nameof or name.
    '''
    return map(_DUNDER_nameof, names__, names__)


def _Property_RO(method):
    '''(INTERNAL) Can't import L{props.Property_RO}, I{recursively}.
    '''
    name = _DUNDER_nameof(method)

    def _del(inst, *unused):  # PYCHOK no cover
        inst.__dict__.pop(name, None)

    def _get(inst, *unused):  # PYCHOK 2 vs 3 args
        try:  # to get the cached value immediately
            v = inst.__dict__[name]
        except (AttributeError, KeyError):
            # cache the value in the instance' __dict__
            inst.__dict__[name] = v = method(inst)
        return v

    def _set(inst, val):  # PYCHOK no cover
        setattr(inst, name, val)  # force error

    return property(_get, _set, _del)


class _MODS_Base(object):
    '''(INTERNAL) Base-class for C{lazily._ALL_MODS}.
    '''
    def __delattr__(self, attr):  # PYCHOK no cover
        self.__dict__.pop(attr, None)

    def __setattr__(self, attr, value):  # PYCHOK no cover
        e = _MODS.errors
        n = _DOT_(self.name, attr)
        t = _EQUALSPACED_(n, repr(value))
        raise e._AttributeError(_immutable_, txt=t)

    @_Property_RO
    def basics(self):
        '''Get module C{pygeodesy.basics}, I{once}.
        '''
        from pygeodesy import basics as b  # DON'T _lazy_import2
        return b

    @_Property_RO
    def bits_machine2(self):
        '''Get platform 2-list C{[bits, machine]}, I{once}.
        '''
        import platform as p
        m = p.machine()  # ARM64, arm64, x86_64, iPhone13,2, etc.
        m = m.replace(_COMMA_, _UNDER_)
        if m.lower() == 'x86_64':  # PYCHOK on Intel or Rosetta2 ...
            v = p.mac_ver()[0]  # ... and only on macOS ...
            if v and _version2(v) > (10, 15):  # ... 11+ aka 10.16
                # <https://Developer.Apple.com/forums/thread/659846>
                #  _sysctl_uint('hw.optional.arm64') and \
                if _sysctl_uint('sysctl.proc_translated'):
                    m = _UNDER_(_arm64_, m)  # Apple Si emulating Intel x86-64
        return [p.architecture()[0],  # bits
                m]  # arm64, arm64_x86_64, x86_64, etc.

    @_Property_RO
    def ctypes3(self):
        '''Get C{ctypes.CDLL}, C{find_library} and C{dlopen}, I{once}.
        '''
        import ctypes as c
        from ctypes.util import find_library as f

        def dlopen(name):  # on macOS only
            return c._dlopen(name, c.DEFAULT_MODE)

        return c.CDLL, f, (dlopen if _ismacOS() else None)

    @_Property_RO
    def errors(self):
        '''Get module C{pygeodesy.errors}, I{once}.
        '''
        from pygeodesy import errors as e  # DON'T _lazy_import2
        return e

    @_Property_RO
    def inspect(self):  # in .basics
        '''Get module C{inspect}, I{once}.
        '''
        import inspect as i
        return i

    def ios_ver(self):
        '''Mimick C{platform.xxx_ver} for C{iOS}.
        '''
        try:  # Pythonista only
            from platform import iOS_ver
            t = iOS_ver()
        except (AttributeError, ImportError):
            t = NN, (NN, NN, NN), NN
        return t

    @_Property_RO
    def name(self):
        '''Get this name (C{str}).
        '''
        return _DUNDER_nameof(self.__class__)

    @_Property_RO
    def nix2(self):  # PYCHOK no cover
        '''Get Linux 2-tuple C{(distro, version)}, I{once}.
        '''
        from platform import uname
        v, n = NN, uname()[0]  # [0] == .system
        if n.lower() == 'linux':
            try:  # use distro only on Linux, not macOS, etc.
                import distro  # <https://PyPI.org/project/distro>
                _a = _MODS.streprs.anstr
                v  = _a(distro.version())  # first
                n  = _a(distro.id())  # .name()?
            except (AttributeError, ImportError):
                pass  # v = str(_0_0)
            n = n.capitalize()
        return n, v

    def nix_ver(self):  # PYCHOK no cover
        '''Mimick C{platform.xxx_ver} for C{*nix}.
        '''
        _, v = _MODS.nix2
        t = _version2(v, n=3) if v else (NN, NN, NN)
        return v, t, machine()

    @_Property_RO
    def os(self):
        '''Get module C{os}, I{once}.
        '''
        import os as o
        import os.path
        return o

    @_Property_RO
    def osversion2(self):
        '''Get 2-list C{[OS, release]}, I{once}.
        '''
        import platform as p
        _Nix, _ = _MODS.nix2
        # - mac_ver() returns ('10.12.5', ..., 'x86_64') on
        #   macOS and ('10.3.3', ..., 'iPad4,2') on iOS
        # - win32_ver is ('XP', ..., 'SP3', ...) on Windows XP SP3
        # - platform() returns 'Darwin-16.6.0-x86_64-i386-64bit'
        #   on macOS and 'Darwin-16.6.0-iPad4,2-64bit' on iOS
        # - sys.platform is 'darwin' on macOS, 'ios' on iOS,
        #   'win32' on Windows and 'cygwin' on Windows/Gygwin
        # - distro.id() and .name() return 'Darwin' on macOS
        for n, v in ((_iOS_, _MODS.ios_ver),
                     (_macOS_,   p.mac_ver),
                     (_Windows_, p.win32_ver),
                     (_Nix,  _MODS.nix_ver),
                     ('Java',    p.java_ver),
                     ('uname',   p.uname)):
            v = v()[0]
            if v and n:
                break
        else:
            n = v = NN  # XXX AssertionError?
        return [n, v]

    @_Property_RO
    def _Popen_kwds2(self):
        '''(INTERNAL) Get C{subprocess.Popen} and C{-kwds}.
        '''
        import subprocess as s
        kwds = dict(creationflags=0,  # executable=sys.executable, shell=True,
                    stdin=s.PIPE, stdout=s.PIPE, stderr=s.STDOUT)
        if _MODS.sys_version_info2 > (3, 6):
            kwds.update(text=True)
        return s.Popen, kwds

    @_Property_RO
    def Pythonarchine(self):
        '''Get 3- or 4-list C{[PyPy, Python, bits, machine]}, I{once}.
        '''
        v  =  _sys.version
        l3 = [_Python_(v)] + _MODS.bits_machine2
        pypy = _PyPy__(v)
        if pypy:  # PYCHOK no cover
            l3.insert(0, pypy)
        return l3

    @_Property_RO
    def streprs(self):
        '''Get module C{pygeodesy.streprs}, I{once}.
        '''
        from pygeodesy import streprs as s  # DON'T _lazy_import2
        return s

    @_Property_RO
    def sys_version_info2(self):
        '''Get C{sys.version_inf0[:2], I{once}.
        '''
        return _sys.version_info[:2]

    @_Property_RO
    def version(self):
        '''Get pygeodesy version, I{once}.
        '''
        from pygeodesy import version as v
        return v

_MODS = _MODS_Base()  # PYCHOK overwritten by .lazily


def _caller3(up, base=True):  # in .lazily, .named
    '''(INTERNAL) Get 3-tuple C{(caller name, file name, line number)}
       for the caller B{C{up}} frames back in the Python call stack.

       @kwarg base: Use C{B{base}=False} for the fully-qualified file
                    name, otherwise the base (module) name (C{bool}).
    '''
    f  =  None
    _b = _MODS.os.path.basename if base else _passarg
    try:
        f = _sys._getframe(up + 1)  # == inspect.stack()[up + 1][0]
        t = _MODS.inspect.getframeinfo(f)
        t =  t.function, _b(t.filename), t.lineno
# or ...
        # f = _sys._getframe(up + 1)
        # c =  f.f_code
        # t = (c.co_name,  # caller name
        #     _b(c.co_filename),  # file name .py
        #      f.f_lineno)  # line number
# or ...
        # t = _MODS.inspect.stack()[up + 1]  # (frame, filename, lineno, function, ...)
        # t =  t[3], _b(t[1]), t[2]
    except (AttributeError, IndexError, ValueError):
        # sys._getframe(1) ... 'importlib._bootstrap' line 1032,
        # may throw a ValueError('call stack not deep enough')
        t = NN, NN, 0
    finally:
        del f  # break ref cycle
    return t


def _enquote(strs, quote=_QUOTE2_, white=NN):  # in .basics, .solveBase
    '''(INTERNAL) Enquote a string containing whitespace or replace
       whitespace by C{white} if specified.
    '''
    if strs:
        t = strs.split()
        if len(t) > 1:
            strs = white.join(t if white else (quote, strs, quote))
    return strs


def _fper(p, q, per=_100_0, prec=1):
    '''Format a percentage C{B{p} * B{per} / B{q}} (C{str}).
    '''
    return '%.*f%%' % (prec, (float(p) * per / float(q)))


_getenv = _MODS.os.getenv  # PYCHOK in .lazily, ...


def _getPYGEODESY(which, dflt=NN):
    '''(INTERNAL) Return an C{PYGEODESY_...} ENV value or C{dflt}.
    '''
    return _getenv(_PYGEODESY(which), dflt)


def _headof(name):
    '''(INTERNAL) Get the head name of qualified C{name} or the C{name}.
    '''
    i = name.find(_DOT_)
    return name if i < 0 else name[:i]


# def _is(a, b):  # PYCHOK no cover
#     '''(INTERNAL) C{a is b}? in C{PyPy}
#     '''
#     return (a == b) if _isPyPy() else (a is b)


def _isAppleSi():  # PYCHOK no cover
    '''(INTERNAL) Is this C{macOS on Apple Silicon}? (C{bool})
    '''
    return _ismacOS() and machine().startswith(_arm64_)


def _is_DUNDER_main(name):
    '''(INTERNAL) Return C{bool(name == '__main__')}.
    '''
    return name == '__main__'


def _isiOS():  # in test/bases
    '''(INTERNAL) Is this C{iOS}? (C{bool})
    '''
    return _MODS.osversion2[0] is _iOS_


def _ismacOS():  # in test/bases
    '''(INTERNAL) Is this C{macOS}? (C{bool})
    '''
    return _sys.platform[:6]   == 'darwin' and \
           _MODS.osversion2[0] is _macOS_  # and _MODS.os.name == 'posix'


def _isNix():  # in test/bases
    '''(INTERNAL) Is this a C{Linux} distro? (C{str} or L{NN})
    '''
    return _MODS.nix2[0]


def _isPyChecker():  # PYCHOK no cover
    '''(INTERNAL) Is C{PyChecker} running? (C{bool}).
    '''
    # .../pychecker/checker.py --limit 0 --stdlib pygeodesy/<mod>/<name>.py
    return _sys.argv[0].endswith('/pychecker/checker.py')


def _isPyPy():  # in test/bases
    '''(INTERNAL) Is this C{PyPy}? (C{bool})
    '''
    # platform.python_implementation() == 'PyPy'
    return _MODS.Pythonarchine[0].startswith(_PyPy__)


def _isWindows():  # in test/bases
    '''(INTERNAL) Is this C{Windows}? (C{bool})
    '''
    return _sys.platform[:3]   == 'win' and \
           _MODS.osversion2[0] is _Windows_


def _load_lib(name):
    '''(INTERNAL) Load a C{dylib}, B{C{name}} must startwith('lib').
    '''
    CDLL, find_lib, dlopen = _MODS.ctypes3
    ns = find_lib(name), name
    if dlopen:
        # macOS 11+ (aka 10.16) no longer provides direct loading of
        # system libraries.  As a result, C{ctypes.util.find_library}
        # will not find any library, unless previously installed by a
        # low-level dlopen(name) call (with the library base C{name}).
        ns += (_DOT_(name, 'dylib'),
               _DOT_(name, 'framework'), _MODS.os.path.join(
               _DOT_(name, 'framework'),  name))
    else:  # not macOS
        dlopen = _passarg  # no-op

    for n in ns:
        try:
            if n and dlopen(n):  # pre-load handle
                lib = CDLL(n)  # == ctypes.cdll.LoadLibrary(n)
                if lib._name:  # has a qualified name
                    return lib
        except (AttributeError, OSError):
            pass

    return None  # raise OSError


def machine():
    '''Return standard C{platform.machine}, but distinguishing Intel I{native}
       from Intel I{emulation} on Apple Silicon (on macOS only).

       @return: Machine C{'arm64'} for Apple Silicon I{native}, C{'x86_64'}
                for Intel I{native}, C{"arm64_x86_64"} for Intel I{emulation},
                etc. (C{str} with C{comma}s replaced by C{underscore}s).
    '''
    return _MODS.bits_machine2[1]


def _name_version(pkg):
    '''(INTERNAL) Return C{pskg.__name__ + ' ' + .__version__}.
    '''
    return _SPACE_(pkg.__name__, pkg.__version__)


def _osversion2(sep=NN):  # in .lazily, test/bases.versions
    '''(INTERNAL) Get the O/S name and release as C{2-list} or C{str}.
    '''
    l2 = _MODS.osversion2
    return sep.join(l2) if sep else l2  # 2-list()


def _passarg(arg):
    '''(INTERNAL) Helper, no-op.
    '''
    return arg


def _passargs(*args):
    '''(INTERNAL) Helper, no-op.
    '''
    return args


def _plural(noun, n, nn=NN):
    '''(INTERNAL) Return C{noun}['s'] or C{NN}.
    '''
    return NN(noun, _s_) if n > 1 else (noun if n else nn)


def _popen2(cmd, stdin=None):  # in .mgrs, .solveBase, .testMgrs
    '''(INTERNAL) Invoke C{B{cmd} tuple} and return 2-tuple C{(std, status)}
       with all C{stdout/-err} output, I{stripped} and C{int} exit status.
    '''
    _Popen, kwds = _MODS._Popen_kwds2
    p = _Popen(cmd, **kwds)  # PYCHOK kwArgs
    r =  p.communicate(stdin)[0]  # stdout + NL + stderr
    return _MODS.basics.ub2str(r).strip(), p.returncode


def print_(*args, **nl_nt_prec_prefix__end_file_flush_sep__kwds):  # PYCHOK no cover
    '''Python 3+ C{print}-like formatting and printing.

       @arg args: Values to be converted to C{str} and joined by B{C{sep}},
                  all positional.

       @see: Function L{printf} for further details.
    '''
    return printf(NN, *args, **nl_nt_prec_prefix__end_file_flush_sep__kwds)


def printf(fmt, *args, **nl_nt_prec_prefix__end_file_flush_sep__kwds):
    '''C{Printf-style} and Python 3+ C{print}-like formatting and printing.

       @arg fmt: U{Printf-style<https://Docs.Python.org/3/library/stdtypes.html#
                 printf-style-string-formatting>} format specification (C{str}).
       @arg args: Arguments to be formatted (any C{type}, all positional).
       @kwarg nl_nt_prec_prefix__end_file_flush_sep__kwds: Optional keyword arguments
                 C{B{nl}=0} for the number of leading blank lines (C{int}), C{B{nt}=0}
                 the number of trailing blank lines (C{int}), C{B{prefix}=NN} to be
                 inserted before the formatted text (C{str}) and Python 3+ C{print}
                 keyword arguments C{B{end}}, C{B{sep}}, C{B{file}} and C{B{flush}}.
                 Any remaining C{B{kwds}} are C{printf-style} name-value pairs to be
                 formatted, I{iff no B{C{args}} are present} using C{B{prec}=6} for
                 the number of decimal digits (C{int}).

       @return: Number of bytes written.
    '''
    b, e, f, fl, p, s, kwds = _print7(**nl_nt_prec_prefix__end_file_flush_sep__kwds)
    try:
        if args:
            t = (fmt % args) if fmt else s.join(map(str, args))
        elif kwds:
            t = (fmt % kwds) if fmt else s.join(
                _MODS.streprs.pairs(kwds, prec=p))
        else:
            t =  fmt
    except Exception as x:
        _E, s = _MODS.errors._xError2(x)
        unstr = _MODS.streprs.unstr
        t = unstr(printf, fmt, *args, **nl_nt_prec_prefix__end_file_flush_sep__kwds)
        raise _E(s, txt=t, cause=x)
    try:
        n = f.write(NN(b, t, e))
    except UnicodeEncodeError:  # XXX only Windows
        t = t.replace('\u2032', _QUOTE1_).replace('\u2033', _QUOTE2_)
        n = f.write(NN(b, t, e))
    if fl:  # PYCHOK no cover
        f.flush()
    return n


def _print7(nl=0, nt=0, prec=6, prefix=NN, sep=_SPACE_, file=_sys.stdout,
                                           end=_NL_, flush=False, **kwds):
    '''(INTERNAL) Unravel the C{printf} and remaining keyword arguments.
    '''
    if nl > 0:
        prefix = NN(_NL_ * nl, prefix)
    if nt > 0:
        end = NN(end, _NL_ * nt)
    return prefix, end, file, flush, prec, sep, kwds


def _PYGEODESY(which, i=0):
    '''(INTERNAL) Return an ENV C{str} C{PYGEODESY_...}.
    '''
    try:
        w = which.__name__.lstrip(_UNDER_)[i:]
    except AttributeError:
        w = which
    return _UNDER_(_pygeodesy_, w).upper()


def _Pythonarchine(sep=NN):  # in .lazily, test/bases versions
    '''(INTERNAL) Get PyPy and Python versions, bits and machine as C{3- or 4-list} or C{str}.
    '''
    l3 = _MODS.Pythonarchine
    return sep.join(l3) if sep else l3  # 3- or 4-list


def _secs2str(secs):  # in .geoids, ../test/bases
    '''Convert a time in C{secs} to C{str}.
    '''
    if secs < _100_0:
        unit = len(_SIsecs) - 1
        while 0 < secs < 1 and unit > 0:
            secs *= 1e3  # _1000_0
            unit -= 1
        t = '%.3f %s' % (secs, _SIsecs[unit])
    else:
        m, s = divmod(secs, 60)
        if m < 60:
            t = '%d:%06.3f' % (int(m), s)
        else:
            h, m = divmod(int(m), 60)
            t = '%d:%02d:%06.3f' % (h, m, s)
    return t


def _sizeof(obj, deep=True):
    '''(INTERNAL) Recursively size an C{obj}ect.

       @kwarg deep: If C{True}, include the size of all
                    C{.__dict__.values()} (C{bool}).

       @return: The C{obj} size in bytes (C{int}), ignoring
                class attributes and counting instances only
                once or C{None}.

       @note: With C{PyPy}, the returned size is always C{None}.
    '''
    try:
        _zB = _sys.getsizeof
        _zD = _zB(None)  # some default
    except TypeError:  # PyPy3.10
        return None

    b = _MODS.basics
    _isiterablen = b.isiterablen
    _Str_Bytes   = b._Strs + b._Bytes  # + (range, map)

    def _zR(s, iterable):
        z, _s = 0, s.add
        for o in iterable:
            i = id(o)
            if i not in s:
                _s(i)
                z += _zB(o, _zD)
                if isinstance(o, dict):
                    z += _zR(s, o.keys())
                    z += _zR(s, o.values())
                elif _isiterablen(o) and not \
                      isinstance(o, _Str_Bytes):
                    z += _zR(s, o)
                elif deep:
                    try:  # size instance' attr values only
                        z += _zR(s, o.__dict__.values())
                    except AttributeError:  # None, int, etc.
                        pass
        return z

    return _zR(set(), (obj,))


def _sysctl_uint(name):
    '''(INTERNAL) Get an C{unsigned int sysctl} item by name, I{ONLY on macOS!}
    '''
    libc = _load_lib('libc') if _ismacOS() else None
    if libc:  # <https://StackOverflow.com/questions/759892/python-ctypes-and-sysctl>
        import ctypes as c
        n = c.c_char_p(_MODS.basics.str2ub(name))  # bytes(name, _utf_8_)
        u = c.c_uint(0)
        z = c.c_size_t(c.sizeof(u))
        r = libc.sysctlbyname(n, c.byref(u), c.byref(z), None, c.c_size_t(0))  # PYCHOK attr
    else:  # not macOS or couldn't find or load 'libc'=
        r = -2
    return int(r if r else u.value)  # -1 ENOENT error, -2 no libc or not macOS


def _tailof(name):
    '''(INTERNAL) Get the base name of qualified C{name} or the C{name}.
    '''
    i = name.rfind(_DOT_) + 1
    return name[i:] if i > 0 else name


def _under(name):  # PYCHOK in .datums, .auxilats, .ups, .utm, .utmupsBase, ...
    '''(INTERNAL) Prefix C{name} with an I{underscore}.
    '''
    return name if name.startswith(_UNDER_) else NN(_UNDER_, name)


def _usage(file_py, *args, **opts_help):  # in .etm, .geodesici  # PYCHOK no cover
    '''(INTERNAL) Build "usage: python -m ..." cmd line for module B{C{file_py}}.
    '''
    if opts_help:

        def _help(alts=(), help=NN, **unused):
            if alts and help:
                h = NN(help, _SPACE_).lstrip(_DASH_)
                for a in alts:
                    if a.startswith(h):
                        return NN(_DASH_, a),

        def _opts(opts=NN, alts=(), **unused):
            # opts='T--v-C-R meter-c|i|n|o'
            d, fmt = NN, _MODS.streprs.Fmt.SQUARE
            for o in (opts + _BAR_(*alts)).split(_DASH_):
                if o:
                    yield fmt(NN(d, _DASH_, o.replace(_BAR_, ' | -')))
                    d =  NN
                else:
                    d = _DASH_

        args = _help(**opts_help) or (tuple(_opts(**opts_help)) + args)

    u = _COLON_(_DUNDER_nameof(_usage)[1:], NN)
    return _SPACE_(u, *_usage_argv(file_py, *args))


def _usage_argv(argv0, *args):
    '''(INTERNAL) Return 3-tuple C{(python, '-m', module, *args)}.
    '''
    o = _MODS.os
    m =  o.path.dirname(argv0)
    m =  m.replace(o.getcwd(), _ELLIPSIS_) \
          .replace(o.sep,      _DOT_).strip()
    b =  o.path.basename(argv0)
    b, x = o.path.splitext(b)
    if x == '.py' and not _is_DUNDER_main(b):
        m = _DOT_(m or _pygeodesy_, b)
    p = NN(_python_, _MODS.sys_version_info2[0])
    return (p, '-m', _enquote(m)) + args


def _version2(version, n=2):
    '''(INTERNAL) Split C{B{version} str} into a C{1-, 2- or 3-tuple} of C{int}s.
    '''
    t = _version_ints(version.split(_DOT_, 2))
    if len(t) < n:
        t += (0,) * n
    return t[:n]


def _version_info(package):  # in .basics, .karney._kWrapped.Math
    '''(INTERNAL) Get the C{package.__version_info__} as a 2- or
       3-tuple C{(major, minor, revision)} if C{int}s.
    '''
    try:
        return _version_ints(package.__version_info__)
    except AttributeError:
        return _version2(package.__version__.strip(), n=3)


def _version_ints(vs):
    # helper for _version2 and _version_info above

    def _ints(vs):
        for v in vs:
            try:
                yield int(v.strip())
            except (TypeError, ValueError):
                pass

    return tuple(_ints(vs))


def _versions(sep=_SPACE_):
    '''(INTERNAL) Get pygeodesy, PyPy and Python versions, bits, machine and OS as C{8- or 9-list} or C{str}.
    '''
    l7 = [_pygeodesy_, _MODS.version] + _Pythonarchine() + _osversion2()
    return sep.join(l7) if sep else l7  # 5- or 6-list


__all__ = tuple(map(_DUNDER_nameof, (machine, print_, printf)))
__version__ = '24.11.06'

if _is_DUNDER_main(__name__):  # PYCHOK no cover

    def _main():
        from pygeodesy import _isfrozen, isLazy

        print_(*(_versions(sep=NN) + ['_isfrozen', _isfrozen,
                                      'isLazy',     isLazy]))

    _main()

# % python3 -m pygeodesy.internals
# pygeodesy 24.11.11 Python 3.13.0 64bit arm64 macOS 14.6.1 _isfrozen False isLazy 1

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
