# -*- coding: utf-8 -*-

u'''Mostly INTERNAL functions, except L{machine}, L{print_} and L{printf}.
'''
# from pygeodesy.basics import isiterablen  # _MODS
# from pygeodesy.errors import _AttributeError, _error_init, _UnexpectedError, _xError2  # _MODS
from pygeodesy.interns import NN, _COLON_, _DOT_, _ELLIPSIS_, _EQUALSPACED_, \
                             _immutable_, _NL_, _pygeodesy_, _PyPy__, _python_, \
                             _QUOTE1_, _QUOTE2_, _s_, _SPACE_, _sys, _UNDER_, _utf_8_
from pygeodesy.interns import _COMMA_, _Python_  # PYCHOK used!
# from pygeodesy.streprs import anstr, pairs, unstr  # _MODS

import os as _os  # in .lazily, ...
import os.path as _os_path
# import sys as _sys  # from .interns

_0_0      =  0.0  # PYCHOK in .basics, .constants
_arm64_   = 'arm64'
_iOS_     = 'iOS'
_macOS_   = 'macOS'
_Windows_ = 'Windows'


def _dunder_nameof(inst, *dflt):
    '''(INTERNAL) Get the double_underscore __name__ attr.
    '''
    try:
        return inst.__name__
    except AttributeError:
        pass
    return dflt[0] if dflt else inst.__class__.__name__


def _dunder_nameof_(*names__):  # in .errors._IsnotError
    '''(INTERNAL) Yield the _dunder_nameof or name.
    '''
    return map(_dunder_nameof, names__, names__)


def _Property_RO(method):
    '''(INTERNAL) Can't I{recursively} import L{props.property_RO}.
    '''
    name = _dunder_nameof(method)

    def _del(inst, attr):  # PYCHOK no cover
        delattr(inst, attr)  # force error

    def _get(inst, **unused):  # PYCHOK 2 vs 3 args
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
        m = _MODS.errors
        t = _EQUALSPACED_(self._DOT_(attr), repr(value))
        raise m._AttributeError(_immutable_, txt=t)

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
        '''Get 3-tuple C{(ctypes.CDLL, ._dlopen, .util.findlibrary)}, I{once}.
        '''
        if _ismacOS():
            from ctypes import CDLL, DEFAULT_MODE, _dlopen

            def dlopen(name):
                return _dlopen(name, DEFAULT_MODE)
        else:  # PYCHOK no cover
            from ctypes import CDLL
            dlopen = _passarg

        from ctypes.util import find_library
        return CDLL, dlopen, find_library

    @_Property_RO
    def ctypes5(self):
        '''Get 5-tuple C{(ctypes.byref, .c_char_p, .c_size_t, .c_uint, .sizeof)}, I{once}.
        '''
        from ctypes import byref, c_char_p, c_size_t, c_uint, sizeof  # get_errno
        return byref, c_char_p, c_size_t, c_uint, sizeof

    def _DOT_(self, name):  # PYCHOK no cover
        return _DOT_(self.name, name)

    @_Property_RO
    def errors(self):
        '''Get module C{pygeodesy.errors}, I{once}.
        '''
        from pygeodesy import errors  # DON'T _lazy_import2
        return errors

    def ios_ver(self):
        '''Mimick C{platform.xxx_ver} for C{iOS}.
        '''
        try:  # Pythonista only
            from platform import iOS_ver
            return iOS_ver()
        except (AttributeError, ImportError):
            return NN, (NN, NN, NN), NN

    @_Property_RO
    def libc(self):
        '''Load C{libc.dll|dylib}, I{once}.
        '''
        return _load_lib('libc')

    @_Property_RO
    def name(self):
        '''Get this name (C{str}).
        '''
        return _dunder_nameof(self.__class__)

    @_Property_RO
    def nix2(self):  # PYCHOK no cover
        '''Get Linux 2-list C{[distro, version]}, I{once}.
        '''
        import platform as p

        n, v = p.uname()[0], NN
        if n.lower() == 'linux':
            try:  # use distro only for Linux, not macOS, etc.
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
                n = v = NN  # XXX AssertioError?
        return [n, v]

    @_Property_RO
    def Pythonarchine(self):
        '''Get 3- or 4-list C{[PyPy, Python, bits, machine]}, I{once}.
        '''
        v  =  _sys.version
        l3 = [_Python_(v)] + self.bits_machine2
        pypy = _PyPy__(v)
        if pypy:  # PYCHOK no cover
            l3.insert(0, pypy)
        return l3

    @_Property_RO
    def streprs(self):
        '''Get module C{pygeodesy.streprs}, I{once}.
        '''
        from pygeodesy import streprs  # DON'T _lazy_import2
        return streprs

_MODS = _MODS_Base()  # PYCHOK overwritten by .lazily


def _caller3(up):  # in .lazily, .named
    '''(INTERNAL) Get 3-tuple C{(caller name, file name, line number)}
       for the caller B{C{up}} stack frames in the Python call stack.
    '''
    # sys._getframe(1) ... 'importlib._bootstrap' line 1032,
    # may throw a ValueError('call stack not deep enough')
    f = _sys._getframe(up + 1)
    c =  f.f_code
    return (c.co_name,  # caller name
           _os_path.basename(c.co_filename),  # file name .py
            f.f_lineno)  # line number


def _dunder_ismain(name):
    '''(INTERNAL) Return C{name == '__main__'}.
    '''
    return name == '__main__'


def _enquote(strs, quote=_QUOTE2_, white=NN):  # in .basics, .solveBase
    '''(INTERNAL) Enquote a string containing whitespace or replace
       whitespace by C{white} if specified.
    '''
    if strs:
        t = strs.split()
        if len(t) > 1:
            strs = white.join(t if white else (quote, strs, quote))
    return strs


def _headof(name):
    '''(INTERNAL) Get the head name of qualified C{name} or the C{name}.
    '''
    i = name.find(_DOT_)
    return name if i < 0 else name[:i]


# def _is(a, b):  # PYCHOK no cover
#     '''(INTERNAL) C{a is b}? in C{PyPy}
#     '''
#     return (a == b) if _isPyPy() else (a is b)


def _isAppleM():
    '''(INTERNAL) Is this C{Apple Silicon}? (C{bool})
    '''
    return _ismacOS() and machine().startswith(_arm64_)


def _isiOS():  # in test/bases.py
    '''(INTERNAL) Is this C{iOS}? (C{bool})
    '''
    return _MODS.osversion2[0] is _iOS_


def _ismacOS():  # in test/bases.py
    '''(INTERNAL) Is this C{macOS}? (C{bool})
    '''
    return _sys.platform[:6]   == 'darwin' and \
           _MODS.osversion2[0] is _macOS_  # and os.name == 'posix'


def _isNix():  # in test/bases.py
    '''(INTERNAL) Is this a C{Linux} distro? (C{str} or L{NN})
    '''
    return _MODS.nix2[0]


def _isPyPy():  # in test/bases.py
    '''(INTERNAL) Is this C{PyPy}? (C{bool})
    '''
    # platform.python_implementation() == 'PyPy'
    return _MODS.Pythonarchine[0].startswith(_PyPy__)


def _isWindows():  # in test/bases.py
    '''(INTERNAL) Is this C{Windows}? (C{bool})
    '''
    return _sys.platform[:3]   == 'win' and \
           _MODS.osversion2[0] is _Windows_


def _load_lib(name):
    '''(INTERNAL) Load a C{dylib}, B{C{name}} must startwith('lib').
    '''
    # macOS 11+ (aka 10.16) no longer provides direct loading of
    # system libraries.  As a result, C{ctypes.util.find_library}
    # will not find any library, unless previously installed by a
    # low-level dlopen(name) call (with the library base C{name}).
    CDLL, dlopen, find_lib = _MODS.ctypes3

    ns = find_lib(name), name
    if dlopen is not _passarg:  # _ismacOS()
        ns += (_DOT_(name, 'dylib'),
               _DOT_(name, 'framework'), _os_path.join(
               _DOT_(name, 'framework'),     name))
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


def _plural(noun, n):
    '''(INTERNAL) Return C{noun}['s'] or C{NN}.
    '''
    return NN(noun, _s_) if n > 1 else (noun if n else NN)


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


def _Pythonarchine(sep=NN):  # in .lazily, test/bases.py versions
    '''(INTERNAL) Get PyPy and Python versions, bit and machine as C{3- or 4-list} or C{str}.
    '''
    l3 = _MODS.Pythonarchine
    return sep.join(l3) if sep else l3  # 3- or 4-list


def _sizeof(obj):
    '''(INTERNAL) Recursively size an C{obj}ect.

       @return: The C{obj} size in bytes (C{int}),
                ignoring class attributes and
                counting duplicates only once or
                C{None}.

       @note: With C{PyPy}, the size is always C{None}.
    '''
    try:
        _zB = _sys.getsizeof
        _zD = _zB(None)  # some default
    except TypeError:  # PyPy3.10
        return None

    _isiterablen = _MODS.basics.isiterablen

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
                elif _isiterablen(o):  # not map, ...
                    z += _zR(s, o)
                else:
                    try:  # size instance' attr values only
                        z += _zR(s, o.__dict__.values())
                    except AttributeError:  # None, int, etc.
                        pass
        return z

    return _zR(set(), (obj,))


def _sysctl_uint(name):
    '''(INTERNAL) Get an unsigned int sysctl item by name, use on macOS ONLY!
    '''
    libc = _MODS.libc
    if libc:  # <https://StackOverflow.com/questions/759892/python-ctypes-and-sysctl>
        byref, char_p, size_t, uint, sizeof = _MODS.ctypes5
        n = name if str is bytes else bytes(name, _utf_8_)  # PYCHOK isPython2 = str is bytes
        u = uint(0)
        z = size_t(sizeof(u))
        r = libc.sysctlbyname(char_p(n), byref(u), byref(z), None, size_t(0))
    else:  # could find or load 'libc'
        r = -2
    return int(r if r else u.value)  # -1 ENOENT error, -2 no libc


def _tailof(name):
    '''(INTERNAL) Get the base name of qualified C{name} or the C{name}.
    '''
    i = name.rfind(_DOT_) + 1
    return name[i:] if i > 0 else name


def _under(name):  # PYCHOK in .datums, .auxilats, .ups, .utm, .utmupsBase, ...
    '''(INTERNAL) Prefix C{name} with an I{underscore}.
    '''
    return name if name.startswith(_UNDER_) else NN(_UNDER_, name)


def _usage(file_py, *args):  # in .etm
    '''(INTERNAL) Build "usage: python -m ..." cmd line for module B{C{file_py}}.
    '''
    m = _os_path.dirname(file_py).replace(_os.getcwd(), _ELLIPSIS_) \
                                 .replace(_os.sep, _DOT_).strip()
    b, x = _os_path.splitext(_os_path.basename(file_py))
    if x == '.py' and not _dunder_ismain(b):
        m = _DOT_(m or _pygeodesy_, b)
    p =  NN(_python_, _sys.version_info[0])
    u = _COLON_(_dunder_nameof(_usage)[1:], NN)
    return _SPACE_(u, p, '-m', _enquote(m), *args)


def _version2(version, n=2):
    '''(INTERNAL) Split C{B{version} str} into a C{1-, 2- or 3-tuple} of C{int}s.
    '''
    t = _version_ints(version.split(_DOT_, 2))
    if len(t) < n:
        t += (0,) * n
    return t[:n]


def _version_info(package):  # in .Base.karney, .basics
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


__all__ = tuple(map(_dunder_nameof, (machine, print_, printf)))
__version__ = '24.06.05'

if _dunder_ismain(__name__):  # PYCHOK no cover

    from pygeodesy import _isfrozen, isLazy, version as vs

    print_(_pygeodesy_, vs, *(_Pythonarchine() + _osversion2()
                                               + ['_isfrozen', _isfrozen,
                                                  'isLazy',     isLazy]))

# **) MIT License
#
# Copyright (C) 2016-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
