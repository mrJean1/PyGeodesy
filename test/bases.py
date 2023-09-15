
# -*- coding: utf-8 -*-

# Base classes and functions for PyGeodesy tests.

# After (C) Chris Veness 2011-2015 published under the same MIT Licence,
# see <https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>
# and <https://www.Movable-Type.co.UK/scripts/latlong.html>.

# from copy import copy as _xcopy
from glob import glob
from inspect import isclass, isfunction, ismethod, ismodule
from os import X_OK, access, getenv, sep as _SEP  # environ
from os.path import abspath, basename, dirname, join as joined, splitext
from platform import architecture, java_ver, mac_ver, win32_ver, uname
from random import gauss, random, seed, shuffle
import sys
from time import localtime, time

seed(localtime().tm_yday)
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

from pygeodesy import anstr, basics, clips, DeprecationWarnings, interns, isint, \
                      isLazy, issubclassof, iterNumpy2over, LazyImportError, \
                      lazily, map2, NN, normDMS, pairs, printf, property_RO, \
                      version as PyGeodesy_version  # PYCHOK expected

_DOT_     = interns._DOT_
_SIsecs   = 'fs', 'ps', 'ns', 'us', 'ms', 'sec'  # reversed
_skipped_ = 'skipped'  # in .run
_SPACE_   = interns._SPACE_
_TILDE_   = interns._TILDE_

__all__ = ('coverage', 'GeodSolve', 'geographiclib',  # constants
           'isIntelPython', 'isiOS', 'ismacOS', 'isNix', 'isPyPy',
           'isPython2', 'isPython3', 'isPython37', 'isWindows',
           'numpy', 'PyGeodesy_dir', 'PythonX', 'scipy', 'test_dir',
           'RandomLatLon', 'TestsBase',  # classes
           'ios_ver', 'nix_ver', 'secs2str',  # functions
           'tilde', 'type2str', 'versions')
__version__ = '23.09.14'

try:
    geographiclib = basics._xgeographiclib(basics, 1, 50)
except ImportError:
    geographiclib = None

# don't test with numpy and scypi older than 1.9 resp. 1.0
try:
    numpy = basics._xnumpy(basics, 1, 9)
except ImportError:
    numpy = None
try:
    scipy = basics._xscipy(basics, 1, 0)
except ImportError:
    scipy = None

_xcopy = basics._xcopy
del basics

try:
    _Ints = int, long
    _Strs = basestring
except NameError:  # Python 3+
    _Ints = int
    _Strs = str

_os_bitstr = architecture()[0]  # XXX sys.maxsize
_pseudo_home_dir = dirname(PyGeodesy_dir or _TILDE_) or _TILDE_

_W_opts = sys.warnoptions or NN
if _W_opts:
    _W_opts = _SPACE_(*(_SPACE_('-W', _) for _ in _W_opts))

PythonX = sys.executable  # python or Pythonista path
isIntelPython = 'intelpython' in PythonX

endswith   = str.endswith
startswith = str.startswith

# isiOS is used by some tests known to fail on iOS only
isiOS      = sys.platform[:3] == 'ios'  # public
ismacOS    = sys.platform[:6] == 'darwin'  # public
isNix      = uname()[0] in ('Linux', 'linux')
isPyPy     = interns._isPyPy()
isPython2  = sys.version_info[0] == 2
isPython3  = sys.version_info[0] == 3
isPython35 = sys.version_info[:2] >= (3, 5)  # in .testCartesian
isPython37 = sys.version_info[:2] >= (3, 7)  # in .run, .testLazy
isWindows  = sys.platform[:3] == 'win'

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

        for p in ('',) + interns._sub_packages:
            m = _join(p, _SEP, '[_a-z]*.py')
            for n in sorted(glob(joined(PyGeodesy_dir, 'pygeodesy', m))):
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

    def test(self, name, value, expect, **kwds):
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
            f = '  FAILED%s, expected %s' % (f, expect)

        self.total += 1  # tests
        if f or self._verbose:
            self.printf('test %d %s: %s%s', self.total, name, v, f, **kwds)
        if r and self._raiser:
            raise TestError('test %d %s', self.total, name)

        if self._time:  # undo _prefix change
            self.__dict__.pop('_prefix')

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


def _fLate(f):
    t = 'proLate' if f < 0 else \
        ('obLate' if f > 0 else 'sphere')
    return 'f(%.1f)%s' % (f, t)


def _getenv_path(envar):
    '''(INTERNAL) Get and validate the path of an executable.
    '''
    p = getenv(envar, None) or None
    if p and not access(p, X_OK):
        # zap the envar to avoid double messages
        # when invoked as C{python -m test.base}
        # environ[envar] = NN
        print('env %s=%r not executable' % (envar, p))
        p = None
    return p


def _get_kwds(fmt='%s', prec=0, known=False, **kwds):
    '''(INTERNAL) Get C{fmt}, C{known} and other C{kwds}.
    '''
    if prec > 0:
        fmt = '%%.%df' % (prec,)
    elif prec < 0:
        fmt = '%%.%de' % (-prec,)
    return fmt, known, kwds


try:  # Pythonista only
    from platform import iOS_ver as ios_ver
except (AttributeError, ImportError):

    def ios_ver(**unused):
        return (NN, (NN, NN, NN), NN)


def _name_version2(path):
    '''(INTERNAL) Get the C{(name, version)} of an executable.
    '''
    if path:
        from pygeodesy.solveBase import _popen2
        try:
            _, r = _popen2((path, '--version'))
            return basename(path), r.split()[-1]
        except (IndexError, IOError, OSError):
            pass
    return ()


def prefix2(prev):  # in .run
    '''Get time prefix and time stamp.
    '''
    t = time()
    p = '%7.3f ' % ((t - prev),)
    p = p.replace('  0.', '   .')
    return p, t


# <https://GitHub.com/ActiveState/code/blob/master/recipes/Python/
#        393090_Binary_floating_point_summatiaccurate_full/recipe-393090.py>
def randoms(n):
    '''Return a (lon) list of random floats.
    '''
    t  = [7, 1e100, -9e-20, -7, -1e100, 8e-20]
    t *= max(1, min(n, 10))
    s  = 0
    _a = t.append
    _g = gauss
    _r = random
    for _ in range(20):
        v  = _g(0, _r())**7 - s
        s += v
        _a(v)
    shuffle(t)
    return t


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
    return path.replace(_pseudo_home_dir, _TILDE_)


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
    '''Get pygeodesy, Python versions, bits, machine, OS name,
       env variables and packages.
    '''
    vs = TestsBase._versions
    if not vs:
        from pygeodesy.interns import _pythonarchine

        vs = 'PyGeodesy', PyGeodesy_version, _pythonarchine(sep=_SPACE_)
        for t in (coverage, numpy, scipy, geographiclib):
            if t:
                vs += t.__name__, t.__version__

        if geographiclib:
            from pygeodesy.karney import _K_2_0, _wrapped
            if _wrapped.Math:
                vs += 'Math', ('_K_2_0' if _K_2_0 else '_K_1_0')

        for t in (GeoConvert, GeodSolve, RhumbSolve):
            vs += _name_version2(t)

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

        TestsBase._versions = vs = _SPACE_.join(vs)
    return vs


GeoConvert = _getenv_path(lazily._PYGEODESY_GEOCONVERT_)
GeodSolve  = _getenv_path(lazily._PYGEODESY_GEODSOLVE_)
RhumbSolve = _getenv_path(lazily._PYGEODESY_RHUMBSOLVE_)
# versions()  # get versions once

if __name__ == '__main__':

    print(versions())

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
