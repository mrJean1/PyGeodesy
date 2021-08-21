
# -*- coding: utf-8 -*-

# Base class and functions for PyGeodesy tests.

# After (C) Chris Veness 2011-2015 published under the same MIT Licence,
# see <https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>
# and <https://www.Movable-Type.co.UK/scripts/latlong.html>.

from copy import copy as _xcopy
from glob import glob
from inspect import isclass, isfunction, ismethod, ismodule
from os import getenv
from os.path import abspath, basename, dirname, join as joined, splitext
from platform import architecture, java_ver, mac_ver, win32_ver, uname
import sys
from time import time

try:
    if float(getenv('PYGEODESY_COVERAGE', '0')) > 0:
        import coverage
    else:
        coverage = None  # ignore coverage
except (ImportError, TypeError, ValueError):
    coverage = None

test_dir = dirname(abspath(__file__))
PyGeodesy_dir = dirname(test_dir)
# extend sys.path to include the ../.. directory,
# required for module .run.py to work
if PyGeodesy_dir not in sys.path:  # Python 3+ ModuleNotFoundError
    sys.path.insert(0, PyGeodesy_dir)

from pygeodesy import anstr, clips, DeprecationWarnings, isLazy, \
                      issubclassof, iterNumpy2over, LazyImportError, \
                      map2, NN, normDMS, pairs, printf, property_RO, \
                      version as PyGeodesy_version  # PYCHOK expected

__all__ = ('coverage', 'GeodSolve', 'geographiclib',  # constants
           'isIntelPython', 'isiOS', 'ismacOS', 'isNix', 'isPyPy',
           'isPython2', 'isPython3', 'isPython37', 'isWindows',
           'numpy', 'PyGeodesy_dir', 'PythonX', 'scipy', 'test_dir',
           'RandomLatLon', 'TestsBase',  # classes
           'ios_ver', 'nix_ver', 'secs2str',  # functions
           'tilde', 'type2str', 'versions')
__version__ = '21.08.16'

try:
    import geographiclib
except ImportError:
    geographiclib = None

# don't test with numpy and scypi older than 1.9 resp. 1.0
from pygeodesy import basics
try:
    numpy = basics._xnumpy(basics, 1, 9)
except ImportError:
    numpy = None
try:
    scipy = basics._xscipy(basics, 1, 0)
except ImportError:
    scipy = None
del basics

try:
    _Ints = int, long
    _Strs = basestring
except NameError:  # Python 3+
    _Ints = int
    _Strs = str

_os_bitstr = architecture()[0]  # XXX sys.maxsize
_pseudo_home_dir = dirname(PyGeodesy_dir or '~') or '~'
_SIsecs = 'fs', 'ps', 'ns', 'us', 'ms', 'sec'  # reversed
_SPACE_ = ' '

_W_opts = sys.warnoptions or ''
if _W_opts:
    _W_opts = '-W ' + ' -W '.join(_W_opts)

PythonX = sys.executable  # python or Pythonista path
isIntelPython = 'intelpython' in PythonX

# isiOS is used by some tests known to fail on iOS only
isiOS      = sys.platform == 'ios'  # public
ismacOS    = sys.platform == 'darwin'  # public
isNix      = uname()[0] in ('Linux', 'linux')
isPyPy     = '[PyPy ' in sys.version  # platform.python_implementation() == 'PyPy'
isPython2  = sys.version_info[0] == 2
isPython3  = sys.version_info[0] == 3
isPython37 = sys.version_info[:2] >= (3, 7)  # for testLazy
isWindows  = sys.platform[:3] == 'win'

GeodSolve  = getenv('PYGEODESY_GEODSOLVE', None)

try:
    # use distro only for Linux, not macOS, etc.
    if isNix:
        import distro  # <https://PyPI.org/project/distro>
    else:
        raise ImportError

    _Nix = anstr(distro.id()).capitalize()  # .name()?

    def nix_ver():  # *nix release
        try:  # no subprocess.check_output ...
            v = distro.version()
        except AttributeError:  # ... Python 2.6
            v = NN
        return anstr(v), _os_bitstr

except ImportError:
    _Nix = NN  # not linux?

    def nix_ver():  # PYCHOK expected
        return _Nix, _os_bitstr


class RandomLatLon(object):
    '''Random LatLon(lat, lon) generator.
    '''
    def __init__(self, LatLon, lat_=170, lon_=350):  # +/- ranges
        self._LatLon = LatLon
        self._lat_ = lat_
        self._lon_ = lon_

        from random import random
        self._random = random

    def __call__(self, **LatLon_kwds):
        return self._LatLon((self._random() - 0.5) * self._lat_,
                            (self._random() - 0.5) * self._lon_, **LatLon_kwds)


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
    _time     =  0
    _verbose  =  True  # print all tests, otherwise failures only
    _versions =  NN  # cached versions() string

    failed  = 0
    known   = 0
    skipped = 0
    total   = 0

    def __init__(self, testfile, version, module=None, verbose=True,
                                          raiser=False, testX=False):
        if not self._versions:  # get versions once
            TestsBase._versions = versions()
        self._file = testfile
        self._name = basename(testfile)
        self.title(self._name, version, module=module)
        self._time = time()
        self._verbose = verbose
        self._raiser = raiser if raiser else '-raiser'.startswith(sys.argv[-1])
        self._testX = testX if testX else ('-testX' in sys.argv[1:] or
                                           '-Exact' in sys.argv[1:])

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

    def pygeodesy_modules(self, deprecated=False):
        '''Yield all PyGeodesy modules.
        '''
        import pygeodesy as _p  # PYCHOK expected
        for _, n in self.pygeodesy_names2(deprecated=deprecated):
            try:
                m = getattr(_p, splitext(n)[0])
            except LazyImportError:
                continue
            if ismodule(m):
                yield m

    def pygeodesy_names2(self, deprecated=False):
        '''Yield all PyGeodesy module files and basenames.
        '''
        for n in sorted(glob(joined(PyGeodesy_dir, 'pygeodesy', '[a-z]*.py'))):
            m = basename(n)
            if deprecated or not m.startswith('deprecated'):
                yield n, m

    def printf(self, fmt, *args, **kwds):  # nl=0, nt=0
        '''Print a formatted line to sys.stdout.
        '''
        printf((fmt % args), prefix=self._prefix, **kwds)

    def results(self, passed='passed', nl=1):
        '''Summarize the test results.
        '''
        r = passed  # or 'skipped'
        s = time() - self._time
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
            r = '(%.1f%%) FAILED%s%s' % (100.0 * f / t, k, w or '')
            n = '%s of %d' % (f, t)
        elif t:
            n = 'all %d' % (t,)
        else:
            n = 'all'
        k = self.skipped or NN
        if k:
            k = ', %d skipped' % (k,)
        r = '%s%s (%s) %s' % (r, k, self._versions, secs2str(s))
        self.printf('%s %s tests %s', n, self._name, r, nl=nl)

    def skip(self, text=NN, n=1):
        '''Skip this test, leave a message.
        '''
        self.skipped += n
        if text and self._verbose:
            t = 'test' if n < 2 else '%s tests' % (n,)
            self.printf('%s skipped (%d): %s', t, self.skipped, text)

    def subtitle(self, module, testing='ing', **kwds):
        '''Print the subtitle of a test suite.
        '''
        t = (basename(module.__name__), module.__version__) + pairs(kwds.items())
        self.printf('test%s(%s)', testing, ', '.join(t), nl=1)

    def test(self, name, value, expect, **kwds):
        '''Compare a test value with the expected result.
        '''
        if self._iterisk:
            name += self._iterisk

        fmt, known, kwds = _fmt_known_kwds(**kwds)

        if not isinstance(expect, _Strs):
            expect = fmt % (expect,)  # expect as str

        f, r, v = NN, False, fmt % (value,)  # value as str
        if v != expect and v != normDMS(expect):
            self.failed += 1  # failures
            if known:  # failed before
                self.known += 1
                f = ', KNOWN'
            else:
                r = True
            f = '  FAILED%s, expected %s' % (f, expect)

        self.total += 1  # tests
        if f or self._verbose:
            self.printf('test %d %s: %s%s', self.total, name, v, f, **kwds)
        if r and self._raiser:
            raise TestError('test %d %s', self.total, name)

    def test_(self, name, value, *expects, **kwds):
        '''Compare a test value with one or several expected results.
        '''
        if self._iterisk:
            name += self._iterisk

        fmt, known, kwds = _fmt_known_kwds(**kwds)

        f, r, v = NN, False, fmt % (value,)  # value as str
        for x in expects:
            if v == x or v == normDMS(x):
                break
        else:
            self.failed += 1  # failures
            if known:  # failed before
                self.known += 1
                f = ', KNOWN'
            else:
                r = True
            f = '  FAILED%s, expected %s' % (f, ' or '.join(expects))

        self.total += 1  # tests
        if f or self._verbose:
            self.printf('test %d %s: %s%s', self.total, name, v, f, **kwds)
        if r and self._raiser:
            raise TestError('test %d %s', self.total, name)

    def test__(self, fmt, *args, **kwds):
        '''Print subtotal test line.
        '''
        t = '-' * len(str(self.total))
        self.printf('test %s %s', t, (fmt % args), **kwds)

    def testCopy(self, inst, *attrs, **kwds):  # Clas=None
        C = kwds.get('Clas', inst.__class__)

        c = _xcopy(inst)
        t = c.__class__, id(c) != id(inst)
        self.test('copy(%s)' % C.__name__, t, (C, True))
        for a in attrs:
            self.test('.' + a, getattr(c, a), getattr(inst, a))

        c = inst.copy()
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
        if module:
            m = ' (module %s %s)' % (basename(module.__name__),
                                     module.__version__)
        else:
            m = NN
        if isLazy:
            z = ' isLazy=%s' % (isLazy,)
        else:
            z = NN
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


class TestError(RuntimeError):  # ValueError's are often caught
    '''Error to show the line number of a test failure.
    '''
    def __init__(self, fmt, *args):
        RuntimeError.__init__(self, fmt % args)


def _fmt_known_kwds(fmt='%s', prec=0, known=False, **kwds):
    '''(INTERNAL) Get C{fmt}, C{known} and other C{kwds}.
    '''
    if prec > 0:
        fmt = '%%.%df' % (prec,)
    elif prec < 0:
        fmt = '%%.%de' % (-prec,)
    return fmt, known, kwds


if isiOS:

    def ios_ver(**unused):
        '''Get the iOS version information.
        '''
        # on iOS ('10.3.3', ..., 'iPad4,2')
        t = mac_ver()
        # append 'machine' iPad, iPhone to 'release'
        r = t[0] + _SPACE_ + t[2].split(',')[0]
        return (r,) + t[1:]

else:  # non-iOS

    def ios_ver(**unused):  # PYCHOK expected
        '''Get the iOS version information.
        '''
        return (NN, (NN, NN, NN), NN)


def secs2str(secs):
    '''Convert a seconds value to string.
    '''
    if secs < 100:
        unit = len(_SIsecs) - 1
        while 0 < secs < 1 and unit > 0:
            secs *= 1000.0
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


def tilde(path):
    '''Return a shortened path, especially Pythonista.
    '''
    return path.replace(_pseudo_home_dir, '~')


def type2str(obj, attr):
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
        t = _SPACE_ + t.__class__.__name__
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
    '''Get pygeodesy, Python versions, bits, machine, OS name and release.
    '''
    if not TestsBase._versions:

        from pygeodesy.interns import _pythonarchine

        vs = 'PyGeodesy', PyGeodesy_version, _pythonarchine(sep=_SPACE_)
        for t in (coverage, geographiclib, numpy, scipy):
            if t:
                vs += t.__name__, t.__version__

        # - mac_ver() returns ('10.12.5', ..., 'x86_64') on
        #   macOS and ('10.3.3', ..., 'iPad4,2') on iOS
        # - win32_ver is ('XP', ..., 'SP3', ...) on Windows XP SP3
        # - platform() returns 'Darwin-16.6.0-x86_64-i386-64bit'
        #   on macOS and 'Darwin-16.6.0-iPad4,2-64bit' on iOS
        # - sys.platform is 'darwin' on macOS, 'ios' on iOS,
        #   'win32' on Windows and 'cygwin' on Windows/Gygwin
        # - distro.id() and .name() return 'Darwin' on macOS
        for t, r in (('iOS',     ios_ver),
                     ('macOS',   mac_ver),
                     ('Windows', win32_ver),
                     (_Nix,      nix_ver),
                     ('Java',    java_ver),
                     ('uname',   uname)):
            r = r()[0]
            if r:
                vs += t, r
                break

        if isinstance(isLazy, int):
            vs += 'isLazy', str(isLazy)

        if getenv('PYTHONDONTWRITEBYTECODE', None):
            vs += '-B',

        if _W_opts:
            vs += _W_opts,

        TestsBase._versions = _SPACE_.join(vs)

    return TestsBase._versions


if __name__ == '__main__':

    print(versions())

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
