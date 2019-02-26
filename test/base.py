
# -*- coding: utf-8 -*-

# Base class and functions for PyGeodesy tests.

# After (C) Chris Veness 2011-2015 published under the same MIT Licence,
# see <http://www.Movable-Type.co.UK/scripts/latlong-vectors.html>
# and <http://www.Movable-Type.co.UK/scripts/latlong.html>.

from inspect import isclass, isfunction, ismethod, ismodule
from os import getenv
from os.path import abspath, basename, dirname
from platform import architecture, java_ver, mac_ver, win32_ver, uname
import sys
from time import time

try:
    import geographiclib
except ImportError:
    geographiclib = None
try:
    import numpy
except ImportError:
    numpy = None
try:
    import scipy
except ImportError:
    scipy = None

test_dir = dirname(abspath(__file__))
PyGeodesy_dir = dirname(test_dir)
# extend sys.path to include the ../.. directory,
# required for module .run.py to work
if PyGeodesy_dir not in sys.path:  # Python 3+ ModuleNotFoundError
    sys.path.insert(0, PyGeodesy_dir)

from pygeodesy import anStr, isLazy, iterNumpy2over, normDMS, \
                      version as PyGeodesy_version  # PYCHOK expected

__all__ = ('geographiclib', 'numpy',  # constants
           'isIntelPython', 'isiOS', 'ismacOS', 'isNix', 'isPyPy',
           'isPython2', 'isPython3', 'isPython37', 'isWindows',
           'PyGeodesy_dir', 'PythonX', 'scipy',
           'TestsBase',  # classes
           'ios_ver', 'secs2str',  # functions
           'test_dir', 'tilde', 'type2str', 'versions')
__version__ = '19.02.23'

try:
    _Ints = int, long
    _Strs = basestring
except NameError:  # Python 3+
    _Ints = int
    _Strs = str

_os_bitstr = architecture()[0]  # XXX sys.maxsize
_pseudo_home_dir = dirname(PyGeodesy_dir or '~') or '~'
_SIsecs = 'fs', 'ps', 'ns', 'us', 'ms', 'sec'  # reversed

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
isWindows  = sys.platform.startswith('win')

try:
    # use distro only for Linux, not macOS, etc.
    if isNix:
        import distro  # <http://PyPI.org/project/distro>
    else:
        raise ImportError

    _Nix = anStr(distro.id()).capitalize()  # .name()?

    def nix_ver():  # *nix release
        try:  # no subprocess.check_output ...
            v = distro.version()
        except AttributeError:  # ... Python 2.6
            v = ''
        return anStr(v), _os_bitstr

except ImportError:
    _Nix = ''  # not linux?

    def nix_ver():  # PYCHOK expected
        return _Nix, _os_bitstr


class TestsBase(object):
    '''Tests based on @examples from the original JavaScript code
       and examples in <http://www.EdWilliams.org/avform.htm> or
       elsewhere as indicated.
    '''
    _file     = ''
    _name     = ''
    _iterisk  = ''
    _prefix   = '    '
    _time     = 0
    _verbose  = True  # print all tests, otherwise failures only
    _versions = ''  # cached versions() string

    failed  = 0
    known   = 0
    skipped = 0
    total   = 0

    def __init__(self, testfile, version, module=None, verbose=True):
        if not self._versions:  # get versions once
            TestsBase._versions = versions()
        self._file = testfile
        self._name = basename(testfile)
        self.title(self._name, version, module=module)
        self._time = time()
        self._verbose = verbose

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

    def printf(self, fmt, *args, **kwds):  # nl=0, nt=0
        '''Print a formatted line to sys.stdout.
        '''
        if kwds:
            nl = '\n' * kwds.get('nl', 0)
            nt = '\n' * kwds.get('nt', 0)
        else:
            nl = nt = ''
        print(''.join((nl, self._prefix, (fmt % args), nt)))

    def results(self, passed='passed', nl=1):
        '''Summarize the test results.
        '''
        s = time() - self._time
        r = passed  # or 'skipped'
        t = self.total
        f = self.failed
        if f:
            k = self.known or ''
            if k:
                if k == f:
                    k = ', ALL KNOWN'
                else:
                    k = ', incl. %s KNOWN' % (k,)
            r = '(%.1f%%) FAILED%s' % (100.0 * f / t, k)
            n = '%s of %d' % (f, t)
        elif t:
            n = 'all %d' % (t,)
        else:
            n = 'all'
        k = self.skipped or ''
        if k:
            k = ', %d skipped' % (k,)
        r = '%s%s (%s) %s' % (r, k, self._versions, secs2str(s))
        self.printf('%s %s tests %s', n, self._name, r, nl=nl)

    def skip(self, text='', n=1):
        '''Skip this test, leave a message.
        '''
        self.skipped += n
        if text and self._verbose:
            t = 'test' if n < 2 else '%s tests' % (n,)
            self.printf('%s skipped (%d): %s', t, self.skipped, text)

    def subtitle(self, module, testing='ing', **kwds):
        '''Print the subtitle of a test suite.
        '''
        t = (basename(module.__name__), module.__version__) + \
             tuple('%s=%s' % t for t in sorted(kwds.items()))
        self.printf('test%s(%s)', testing, ', '.join(t), nl=1)

    def test(self, name, value, expect, fmt='%s', known=False, **kwds):
        '''Compare a test value with the expected one.
        '''
        if self._iterisk:
            name += self._iterisk

        if not isinstance(expect, _Strs):
            expect = fmt % (expect,)  # expect as str

        f, v = '', fmt % (value,)  # value as str
        if v != expect and v != normDMS(expect):
            self.failed += 1  # failures
            if known:  # failed before
                self.known += 1
                f = ', KNOWN'
            f = '  FAILED%s, expected %s' % (f, expect)

        self.total += 1  # tests
        if f or self._verbose:
            self.printf('test %d %s: %s%s', self.total, name, v, f, **kwds)

    def test__(self, fmt, *args, **kwds):
        '''Print subtotal test line.
        '''
        t = '-' * len(str(self.total))
        self.printf('test %s %s', t, (fmt % args), **kwds)

    def testiter(self):
        '''Test with/-out I{iterNumpy2} threshold.
        '''
        yield iterNumpy2over(100000)
        self._iterisk = '*'
        yield iterNumpy2over(1)
        self._iterisk = ''

    def title(self, test, version, module=None):
        '''Print the title of the test suite.
        '''
        if module:
            m = ' (module %s %s)' % (basename(module.__name__),
                                     module.__version__)
        else:
            m = ''
        if isLazy:
            z = ' isLazy=%s' % (isLazy,)
        else:
            z = ''
        self.printf('testing %s %s%s%s', test, version, m, z, nl=1)


if isiOS:

    def ios_ver(**unused):
        '''Get the iOS version information.
        '''
        # on iOS ('10.3.3', ..., 'iPad4,2')
        t = mac_ver()
        # append 'machine' iPad, iPhone to 'release'
        r = t[0] + ' ' + t[2].split(',')[0]
        return (r,) + t[1:]

else:  # non-iOS

    def ios_ver(**unused):  # PYCHOK expected
        '''Get the iOS version information.
        '''
        return ('', ('', '', ''), '')


def secs2str(secs):
    '''Convert a seconds value to string.
    '''
    unit = len(_SIsecs) - 1
    while 0 < secs < 1 and unit > 0:
        secs *= 1000.0
        unit -= 1
    return '%.3f %s' % (secs, _SIsecs[unit])


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
    elif type(t) is property:
        t = ' property'
    elif isinstance(t, _Strs):
        t = ' str'
    else:
        t = str(type(t)).replace("'", '')
        if t.startswith('<') and t.endswith('>'):
            t = ' ' + t[1:-1]
        else:
            t = ' attribute'
    return t


def versions():
    '''Get pygeodesy, Python versions, size, OS name and release.
    '''
    vs = 'PyGeodesy', PyGeodesy_version
    if isPyPy:
        vs += 'PyPy', sys.version.split('[PyPy ')[1].split()[0]
    vs += 'Python', sys.version.split()[0], _os_bitstr
    if geographiclib:
        vs += 'geographiclib', geographiclib.__version__
    if numpy:
        vs += 'numpy', numpy.__version__
    if scipy:
        vs += 'scipy', scipy.__version__

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
        vs += 'B',

    return ' '.join(vs)
