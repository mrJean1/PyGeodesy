
# -*- coding: utf-8 -*-

# Base classes and functions for PyGeodesy tests.

# After (C) Chris Veness 2011-2024 published under the same MIT Licence,
# see <https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>
# and <https://www.Movable-Type.co.UK/scripts/latlong.html>.

from glob import glob
from inspect import isclass, isfunction, ismethod, ismodule
import os
from os.path import abspath, basename, dirname, join as joined, splitext
from random import random, seed
import sys
from time import localtime, time

seed(localtime().tm_yday)
del localtime, seed

test_dir = dirname(abspath(__file__))
PyGeodesy_dir = dirname(test_dir)
# extend sys.path to include the ../.. directory,
# required for module .run.py to work
if PyGeodesy_dir not in sys.path:  # Python 3+ ModuleNotFoundError
    sys.path.insert(0, PyGeodesy_dir)

from pygeodesy import anstr, basics, clips, DeprecationWarnings, internals, interns, \
                      isint, isLazy, issubclassof, iterNumpy2over, LazyImportError, \
                      karney, map2, NN, normDMS, pairs, printf, property_RO, \
                      version as PyGeodesy_version  # PYCHOK expected
from pygeodesy.internals import _getenv, _name_version, _secs2str as secs2str

_DOT_     = interns._DOT_
_skipped_ = 'skipped'  # in .run
_SPACE_   = interns._SPACE_
_TILDE_   = interns._TILDE_
_HOME_dir = dirname(PyGeodesy_dir or _TILDE_) or _TILDE_

__all__ = ('bits_mach2', 'coverage', 'GeodSolve', 'geographiclib',
           'isiOS', 'ismacOS', 'isNix', 'isPyPy',  # 'isIntelPython'
           'isPython2', 'isPython3', 'isPython37', 'isPython39', 'isWindows',
           'numpy', 'PyGeodesy_dir', 'PythonX', 'scipy', 'test_dir',
           'RandomLatLon', 'TestsBase',  # classes
           'secs2str', 'tilde', 'type2str', 'versions')  # functions
__version__ = '24.12.12'

try:
    if float(_getenv('PYGEODESY_COVERAGE', '0')) > 0:
        import coverage
    else:
        coverage = None  # ignore coverage
except (ImportError, TypeError, ValueError):
    coverage = None
try:
    geographiclib = basics._xgeographiclib(basics, 1, 50)  # 1.50+
except ImportError:
    geographiclib = None
try:
    numpy = basics._xnumpy(basics, 1, 9)  # 1.9+
except ImportError:
    numpy = None
try:
    scipy = basics._xscipy(basics, 1, 0)  # 1.0+
except ImportError:
    scipy = None
_xcopy = basics._xcopy

try:
    _Ints = int, long
    _Strs = basestring
except NameError:  # Python 3+
    _Ints = int
    _Strs = str

_W_opts = sys.warnoptions or NN
if _W_opts:
    _W_opts = _SPACE_(*(_SPACE_('-W', _) for _ in _W_opts))

PythonX = sys.executable  # python or Pythonista path
# isIntelPython = 'intelpython' in PythonX

endswith   = str.endswith
startswith = str.startswith

v2 = sys.version_info[:2]
bits_mach2 = internals._MODS.bits_machine2
# isiOS is used by some tests known to fail on iOS only
isiOS      = internals._isiOS()  # public
ismacOS    = internals._ismacOS()  # public
isNix      = internals._isNix()
isPyPy     = internals._isPyPy()
isPython2  = v2[0] == 2
isPython3  = v2[0] == 3
isPython35 = v2 >= (3, 5)  # in .testCartesian
isPython37 = v2 >= (3, 7)  # in .run, .testLazy
isPython39 = v2 >= (3, 9)  # arm64 Apple Si
isWindows  = internals._isWindows()
del v2


class RandomLatLon(object):
    '''Random LatLon(lat, lon) generator.
    '''
    _ndigits = None

    def __init__(self, LatLon, lat_=170, lon_=350, ndigits=None):  # +/- ranges
        self._LatLon = LatLon
        self._lat = lat_
        self._lon = lon_
        if ndigits and isint(ndigits):
            self._ndigits = ndigits

    def __call__(self, **LatLon_kwds):
        lat = self._random_round(self._lat)
        lon = self._random_round(self._lon)
        return self._LatLon(lat, lon, **LatLon_kwds)

    def _random_round(self, scale):
        r = (random() - 0.5) * scale
        n =  self._ndigits
        if n is not None:
            r = round(r, n)  # ndigits=None
        return r


class TestsBase(object):
    '''Tests based on @examples from the original JavaScript code
       and examples in <https://www.EdWilliams.org/avform.htm> or
       elsewhere as indicated.
    '''
    _file     =  NN
    _name     =  NN
    _iterisk  =  NN
    _prefix   = _SPACE_ * 4
    _raiser   =  False
    _testX    =  False  # slow Exact test
    _time     =  None  # -prefix
    _time0    =  0
    _verbose  =  True  # print all tests, otherwise failures only
    _versions =  NN  # cached versions() string, set below

    failed  = 0
    known   = 0
    skipped = 0
    total   = 0

    def __init__(self, testfile, version, module=None, verbose=True,
                                          raiser=False, testX=False):
        self._file    = testfile
        self._name    = basename(testfile)
        self.title(self._name, version, module=module)
        self._verbose = verbose
        self._raiser  = raiser if raiser else '-raiser'.startswith(sys.argv[-1])
        self._testX   = testX  if testX else ('-testX' in sys.argv[1:] or
                                              '-Exact' in sys.argv[1:])
        self._time0   = time()
        try:
            sys.argv.remove('-prefix')
            self._time = self._time0
        except (IndexError, ValueError):
            pass

    def errors(self):
        '''Return the number of tests failures,
           excluding the KNOWN failures.
        '''
        return self.failed - self.known  # new failures

    def exit(self, errors=0):
        '''Exit with the number of test failures as exit status.
        '''
        sys.exit(min(errors + self.errors(), 99))

    def iterNumpy2over(self, n=None):
        '''Set the I{iterNumpy2} threshold.

           @keyword n: Optional, new threshold value (C{int}).

           @return: Previous threshold (C{int}).
        '''
        p = iterNumpy2over(n)
        self.test('iterNumpy2over', iterNumpy2over(), n)
        return p

    @property_RO
    def name(self):
        return self._name

    def printf(self, fmt, *args, **kwds):  # nl=0, nt=0
        '''Print a formatted line to sys.stdout.
        '''
        printf((fmt % args), prefix=self._prefix, **kwds)

    def pygeodesy_classes(self, Base=None, deprecated=False):
        '''Yield all PyGeodesy class definitions.
        '''
        B = Base if isclass(Base) else None
        for m in self.pygeodesy_modules(deprecated=deprecated):
            for n in dir(m):
                C = getattr(m, n)
                if isclass(C) and (issubclassof(C, B) if B else True) \
                              and C.__module__ == m.__name__:
                    yield C

    def pygeodesy_classes_(self, Base, *Classes):
        '''Yield all PyGeodesy class definitions, except any in C{Classes}.
        '''
        for C in self.pygeodesy_classes(Base):
            if C not in Classes:
                yield C

    def pygeodesy_modules(self, deprecated=False):
        '''Yield all PyGeodesy modules.
        '''
        import pygeodesy as _p  # PYCHOK expected
        for _, n in self.pygeodesy_names2(deprecated=deprecated):
            try:
                n = splitext(n)[0]
                if _DOT_ in n:
                    p, n = n.split(_DOT_)
                    m = getattr(_p, p)
                    if n != '__init__':
                        m = getattr(m, n)
                else:
                    m = getattr(_p, n)
            except (AttributeError, LazyImportError):
                continue
            if ismodule(m):
                yield m

    def pygeodesy_names2(self, deprecated=False):
        '''Yield all PyGeodesy module files and basenames.
        '''
        def _join(p, j, n):
            return (p + j + n) if p else n

        for p in (NN,) + interns._SUB_PACKAGES:
            m = _join(p, os.sep, '[_a-z]*.py')
            m =  joined(PyGeodesy_dir, interns._pygeodesy_, m)
            for n in sorted(glob(m)):
                m = _join(p, _DOT_, basename(n))
                if deprecated or not m.startswith('deprecated'):
                    yield n, m

    def results(self, passed='passed', nl=1):
        '''Summarize the test results.
        '''
        n = 'all'
        r = passed  # or _skipped_
        s = time() - self._time0
        t = self.total
        w = DeprecationWarnings()
        f = self.failed + (w or 0)
        if f:
            a = ', incl. '
            k = self.known or NN
            if k:
                if k == f:
                    k = ', ALL KNOWN'
                else:
                    k = '%s%s KNOWN' % (a, k)
                    a = ' plus '
            if w:
                w = '%s%s %s%s' % (a, w, DeprecationWarning.__name__,
                                         ('s' if w > 1 else ''))
            p = f * 100.0 / t
            r = '(%.1f%%) FAILED%s%s' % (p, k, w or '')
            n = '%s of %d' % (f, t)
        elif t:
            n = '%s %d' % (n, t)
        k = self.skipped or NN
        if k:
            k = ', %d %s' % (k, _skipped_)
        r = '%s%s (%s) %s' % (r, k, self.versions, secs2str(s))
        self.printf('%s %s tests %s', n, self._name, r, nl=nl)

    def skip(self, text=NN, n=1):
        '''Skip this test, leave a message.
        '''
        self.skipped += n
        if text and self._verbose:
            t = 'test' if n < 2 else '%s tests' % (n,)
            self.printf('%s %s (%d): %s', t, _skipped_, self.skipped, text, nl=1)

    def subtitle(self, module, testing='ing', **kwds):
        '''Print the subtitle of a test suite.
        '''
        t = (basename(module.__name__), module.__version__) + pairs(kwds.items())
        self.printf('test%s(%s)', testing, ', '.join(t), nl=1)

    def test(self, name, value, expect, error=0, **kwds):
        '''Compare a test value with the expected result.
        '''
        if self._time:
            self._prefix, self._time = prefix2(self._time)

        if self._iterisk:
            name += self._iterisk

        fmt, known, kwds = _get_kwds(**kwds)

        if not isinstance(expect, _Strs):
            expect = fmt % (expect,)  # expect as str

        f, r, v = NN, False, fmt % (value,)  # value as str
        if v != expect and v != normDMS(expect) and not (
           callable(known) and known(v, expect)):
            self.failed += 1  # failures
            if not known or callable(known):
                r = True
            else:  # failed before
                self.known += 1
                f = ', KNOWN'
            if error:
                f = '%s (%g)' % (f, error)
            f = '  FAILED%s, expected %s' % (f, expect)

        self.total += 1  # tests
        if f or self._verbose:
            self.printf('test %d %s: %s%s', self.total, name, v, f, **kwds)
        if r and self._raiser:
            raise TestError('test %d %s', self.total, name)

        if self._time:  # undo _prefix change
            self.__dict__.pop('_prefix')  # _xkwds_pop(self.__dict__, _prefix=None)

    def test_tol(self, name, value, expect, tol=1e-12, **kwds):
        e = 0 if value is None or expect is None else abs(value - expect)  # None iteration
        if e:
            m = max(abs(value), abs(expect))
            if m:
                e = min(e, m, e / m)
            return self.test(name, value, expect, error=e, known=e < tol, **kwds)
        else:
            return self.test(name, value, expect,          known=True, **kwds)

    def test__(self, fmt, *args, **kwds):
        '''Print subtotal test line.
        '''
        t = '-' * len(str(self.total))
        self.printf('test %s %s', t, (fmt % args), **kwds)

    def testCopy(self, inst, *attrs, **kwds):  # Clas=None
        C = kwds.get('Clas', inst.__class__)

        c = _xcopy(inst, **kwds)
        t = c.__class__, id(c) != id(inst)
        self.test('copy(%s)' % C.__name__, t, (C, True))
        for a in attrs:
            self.test('.' + a, getattr(c, a), getattr(inst, a))

        c = inst.copy(**kwds)
        t = c.__class__, id(c) != id(inst)
        self.test(C.__name__ + '.copy()', t, (C, True))
        for a in attrs:
            self.test('.' + a, getattr(c, a), getattr(inst, a))

    def testiter(self):
        '''Test with/-out I{iterNumpy2} threshold.
        '''
        yield iterNumpy2over(100000)
        self._iterisk = '*'
        yield iterNumpy2over(1)
        self._iterisk = NN

    def title(self, test, version, module=None):
        '''Print the title of the test suite.
        '''
        z = m = NN
        if module:
            m = ' (module %s %s)' % (basename(module.__name__),
                                     module.__version__)
        if isLazy:
            z = ' isLazy=%s' % (isLazy,)
        self.printf('testing %s %s%s%s', test, version, m, z, nl=1)

    @property
    def verbose(self):
        '''Get verbosity (C{bool}).
        '''
        return self._verbose

    @verbose.setter  # PYCHOK setter!
    def verbose(self, v):
        '''Set verbosity (C{bool}).
        '''
        self._verbose = bool(v)

    @property_RO
    def versions(self):
        '''Get versions (C{str}).
        '''
        return self._versions or versions()


class TestError(RuntimeError):  # ValueError's are often caught
    '''Error to show the line number of a test failure.
    '''
    def __init__(self, fmt, *args):
        RuntimeError.__init__(self, fmt % args)


def _env_c2(c, dev_null='>/dev/null', NUL_='>NUL:'):  # .testFrozen, .testLazily
    cmd = _SPACE_(PythonX, c)

    if ismacOS or isNix:
        env_cmd = _SPACE_('env %s', cmd, dev_null)
    elif isWindows:
        env_cmd = _SPACE_('set %s;', cmd, NUL_)
    else:
        env_cmd =  NN

    H = _getenv('HOME', NN)
    if H and cmd.startswith(H):
        cmd = NN(_TILDE_, cmd[len(H):])

    return env_cmd, cmd


def _fLate(f):
    t = 'proLate' if f < 0 else \
        ('obLate' if f > 0 else 'sphere')
    return 'f(%.1f)%s' % (f, t)


def _get_kwds(fmt='%s', prec=0, known=False, **kwds):
    '''(INTERNAL) Get C{fmt}, C{known} and other C{kwds}.
    '''
    if prec > 0:
        fmt = '%%.%df' % (prec,)
    elif prec < 0:
        fmt = '%%.%de' % (-prec,)
    return fmt, known, kwds


def prefix2(prev):  # in .run
    '''Get time prefix and time stamp.
    '''
    t = time()
    p = '%8.3f ' % ((t - prev),)
    p = p.replace('  0.', '   .')
    return p, t


def tilde(path):
    '''Return a shortened path, especially Pythonista.
    '''
    return path.replace(_HOME_dir, _TILDE_)


def type2str(obj, attr, **renamed):
    '''Return the type name of an object attribute.
    '''
    t = getattr(obj, attr, None)
    if isclass(t):
        t = '() class'
    elif isinstance(t, float):
        t = ' float'
    elif ismethod(t):
        t = '() method'
    elif isfunction(t):
        if isclass(obj):
            t = '() method'
        else:
            t = '() function'
    elif isinstance(t, _Ints):
        t = ' int'
    elif ismodule(t):
        t = ' module'
    elif isinstance(t, property):  # type(t) is property
        t = t.__class__.__name__
        if renamed:
            t = renamed.get(t, t)
        t = _SPACE_ + t
    elif isinstance(t, _Strs):
        t = ' str'
    else:
        t = str(type(t)).replace("'", NN)
        if t.startswith('<') and t.endswith('>'):
            t = _SPACE_ + t[1:-1]
        else:
            t = ' attribute'
    return t


def versions():
    '''Get pygeodesy, Python versions, bits, machine, OS,
       env variables and packages.
    '''
    vs = TestsBase._versions
    if not vs:

        vs =  internals._Pythonarchine(sep=_SPACE_)
        vs = 'PyGeodesy', PyGeodesy_version, vs
        for t in (coverage, geographiclib, numpy, scipy):
            if t:  # .internals._name_version
                vs += _name_version(t),

        if geographiclib:
            r = karney._wrapped
            t = r.Math_K_2
            if t:
                vs += r.Math.__name__, t

        for t in (GeoConvert, GeodSolve, IntersectTool, RhumbSolve):
            if t:
                vs += karney._Xables.name_version(t),

        t, r = internals._osversion2()
        if r:
            vs += t, r

        if isinstance(isLazy, int):
            vs += 'isLazy', str(isLazy)

        if _getenv('PYTHONDONTWRITEBYTECODE', None):
            vs += '-B',

        if _W_opts:
            vs += _W_opts,

        TestsBase._versions = vs = _SPACE_.join(vs)
    return vs

# versions()  # get versions once


def _X_OK(_Xables, which):
    '''(INTERNAL) Get a I{Karney}'s executable path or C{None}.
    '''
#   n = which.__name__
    p = which()  # sets _Xables.ENV
    if p is _Xables.ENV:  # or n.lower() in basics._XPACKAGES:
        p = None
    elif p and not _Xables.X_OK(p):
        # zap the envar to avoid double messages
        # when invoked as C{python -m test.bases}
        # environ[envar] = NN
        print(_Xables.X_not(p))
        p = None
    return p

_Xables       =  karney._Xables  # PYCHOK blank
GeoConvert    = _X_OK(_Xables, _Xables.GeoConvert)
GeodSolve     = _X_OK(_Xables, _Xables.GeodSolve)
IntersectTool = _X_OK(_Xables, _Xables.IntersectTool)
RhumbSolve    = _X_OK(_Xables, _Xables.RhumbSolve)
del basics, _Xables, _X_OK

if internals._is_DUNDER_main(__name__):
    try:
        import coverage  # PYCHOK re-imported
    except ImportError:
        pass
    print(versions())

# % python3.13 -m test.bases
# PyGeodesy 24.10.24 Python 3.13.0 64bit arm64 coverage 7.6.1 geographiclib 2.0 Math 2 macOS 14.6.1 isLazy 1

# % python3.12 -m test.bases
# PyGeodesy 24.10.24 Python 3.12.7 64bit arm64 coverage 7.6.1 geographiclib 2.0 numpy 2.1.0 scipy 1.14.1 Math 2 macOS 14.6.1 isLazy 1

# % python3.11 -m test.bases
# PyGeodesy 24.10.24 Python 3.11.5 64bit arm64 coverage 7.6.1 geographiclib 2.0 numpy 1.24.2 scipy 1.10.1 Math 2 macOS 14.4.1 isLazy 1 -W ignore

# % python3.10 -m test.bases
# PyGeodesy 24.10.24 Python 3.10.8 64bit arm64 coverage 7.6.1 geographiclib 2.0 numpy 1.23.3 scipy 1.9.1 Math 2 macOS 14.4.1 isLazy 1 -W ignore

# % python2 -m test.bases
# PyGeodesy 24.5.12 Python 2.7.18 64bit arm64_x86_64 coverage 5.5 geographiclib 1.50 numpy 1.16.6 scipy 1.2.2 macOS 10.16

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
