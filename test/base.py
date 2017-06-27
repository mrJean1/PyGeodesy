
# -*- coding: utf-8 -*-

# Base class and functions for PyGeodesy tests.

# After (C) Chris Veness 2011-2015 published under the same MIT Licence,
# see <http://www.movable-type.co.uk/scripts/latlong-vectors.html>
# and <http://www.movable-type.co.uk/scripts/latlong.html>.

from inspect import isclass, isfunction, ismethod, ismodule
from os.path import abspath, basename, dirname
from platform import architecture, java_ver, mac_ver, win32_ver, uname
import sys
from time import time

PyGeodesy_dir = dirname(dirname(abspath(__file__)))
# extend sys.path to include the ../.. directory
if PyGeodesy_dir not in sys.path:  # Python 3+ ModuleNotFoundError
    sys.path.insert(0, PyGeodesy_dir)

from pygeodesy import version as PyGeodesy_version, \
                      normDMS  # PYCHOK expected

__all__ = ('isiOS', 'PyGeodesy_dir', 'Python_O',  # constants
           'TestsBase',
           'runner', 'secs2str', 'tilde', 'type2str', 'versions')
__version__ = '17.06.25'

try:
    _int = int, long
    _str = basestring
except NameError:  # Python 3+
    _int = int
    _str = str

_pseudo_home_dir = dirname(PyGeodesy_dir or '~') or '~'

# isiOS is used by some tests known to fail on iOS only
isiOS = True if sys.platform == 'ios' else False  # public

Python_O = sys.executable  # python or Pythonsta path
if not __debug__:
    Python_O += ' -OO'  # optimized


class TestsBase(object):
    '''Tests based on @examples from the original JavaScript code
       and examples in <http://www.edwilliams.org/avform.htm> or
       elsewhere as indicated.
    '''
    _file     = ''
    _name     = ''
    _prefix   = '    '
    _time     = 0
    _versions = ''  # cached versions() string

    failed = 0
    known  = 0
    total  = 0

    def __init__(self, testfile, version, module=None):
        if not self._versions:  # get versions once
            TestsBase._versions = versions()
        self._file = testfile
        self._name = basename(testfile)
        self.title(self._name, version, module=module)
        self._time = time()

    def errors(self):
        '''Return the number of tests failures,
           excluding the KNOWN failures.
        '''
        return self.failed - self.known  # new failures

    def exit(self, errors=0):
        '''Exit with the number of test failures as exist status.
        '''
        sys.exit(min(errors + self.errors(), 99))

    def printf(self, fmt, *args, **kwds):  # nl=0
        '''Print a formatted line to sys.stdout.
        '''
        nl = '\n' * kwds.get('nl', 0)
        print(nl + self._prefix + (fmt % args))

    def results(self, nl=1):
        '''Summarize the test results.
        '''
        s = time() - self._time
        n = self.failed
        if n:
            p = '' if n == 1 else 's'
            k = self.known or ''
            if k:
                k = ', incl. %s KNOWN' % (k,)
            r = '(%.1f%%) FAILED%s' % (100.0 * n / self.total, k)
        else:
            n, p, r = 'all','s', 'passed'
        r = '%s (%s) %s' % (r, self._versions, secs2str(s))
        self.printf('%s %s test%s %s', n, self._name, p, r, nl=nl)

    def subtitle(self, module, testing='ing', **kwds):
        '''Print the subtitle of a test suite.
        '''
        t = (basename(module.__name__), module.__version__) + \
             tuple('%s=%s' % t for t in sorted(kwds.items()))
        self.printf('test%s(%s)', testing, ', '.join(t), nl=1)

    def test(self, name, value, expect, fmt='%s', known=False):
        '''Compare a test value with the expected one.
        '''
        self.total += 1  # tests
        if not isinstance(expect, _str):
            expect = fmt % (expect,)  # expect as str
        v = fmt % (value,)  # value as str
        if v != expect and v != normDMS(expect):
            self.failed += 1  # failures
            f = '  FAILED'
            if known:  # failed before
                self.known += 1
                f += ', KNOWN'
            f = '%s, expected %s' % (f, expect)
        else:
            f = ''  # ' passed'
        self.printf('test %d %s: %s%s', self.total, name, v, f)

    def title(self, test, version, module=None):
        '''Print the title of the test suite.
        '''
        if module:
            m = ' (module %s %s)' % (basename(module.__name__),
                                     module.__version__)
        else:
            m = ''
        self.printf('testing %s %s%s', test, version, m, nl=1)


def secs2str(secs):
    '''Convert a time value to string.
    '''
    unit = ['sec', 'ms', 'us', 'ps']
    while secs < 1 and len(unit) > 1:
        secs *= 1000.0
        unit.pop(0)
    return '%.3f %s' % (secs, unit[0])


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
    elif isinstance(t, _int):
        t = ' int'
    elif ismodule(t):
        t = ' module'
    elif type(t) is property:
        t = ' property'
    elif isinstance(t, _str):
        t = ' str'
    else:
        t = ' attribute'
    return t


def versions():
    '''Get pygeodesy, Python versions, size, OS name and release.
    '''
    vs = 'PyGeodesy', PyGeodesy_version, \
         'Python', sys.version.split()[0], architecture()[0]

    xOS = 'iOS' if isiOS else 'macOS'
    # - mac_ver() returns ('10.12.5', ..., 'x86_64') on
    #   macOS and ('10.3.2', ..., 'iPad4,2') on iOS
    # - platform() returns 'Darwin-16.6.0-x86_64-i386-64bit'
    #   on macOS and 'Darwin-16.6.0-iPad4,2-64bit' on iOS
    # - sys.platform is 'darwin' on macOS and 'ios' on iOS
    for t, r in ((xOS,       mac_ver),
                 ('Windows', win32_ver),
                 ('Java',    java_ver),
                 ('uname',   uname)):
        r = r()[0]
        if r:
            vs += t, r
            break

    return ' '.join(vs)


if isiOS:

    try:  # prefer StringIO over io
        from StringIO import StringIO
    except ImportError:  # Python 3+
        from io import StringIO
    from runpy import run_path
    from traceback import format_exception

    def runner(test):
        '''Invoke one test module and return
           the exit status and console output.
        '''
        # mimick partial behavior of function runner
        # further below because subprocess.Popen is
        # not available on iOS/Pythonista/Python.

        x = None  # no exit, no exception

        sys3 = sys.argv, sys.stdout, sys.stderr
        sys.stdout = sys.stderr = std = StringIO()
        try:
            sys.argv = [test]
            run_path(test, run_name='__main__')
        except:  # PYCHOK have to on Pythonista
            x = sys.exc_info()
            if x[0] is SystemExit:
                x = x[1].code  # exit status
            else:  # append traceback
                x = [t for t in format_exception(*x)
                             if 'runpy.py", line' not in t]
                print(''.join(map(tilde, x)).rstrip())
                x = 1  # count as a failure
        sys.argv, sys.stdout, sys.stderr = sys3

        r = std.getvalue()
        if isinstance(r, bytes):  # Python 3+
            r = r.decode('utf-8')

        std.close()
        std = None  # del std

        if x is None:  # no exit or exception, count
            # failed tests excluding KNOWN ones
            x = r.count('FAILED, expected')
        return x, r

    Python_O = basename(Python_O)

else:  # non-iOS

    from subprocess import PIPE, STDOUT, Popen

    def runner(test):  # PYCHOK expected
        '''Invoke one test module and return
           the exit status and console output.
        '''
        c = [Python_O, test]
        p = Popen(c, creationflags=0,
                     executable   =sys.executable,
                   # shell        =True,
                     stdin        =None,
                     stdout       =PIPE,  # XXX
                     stderr       =STDOUT)  # XXX

        r = p.communicate()[0]
        if isinstance(r, bytes):  # Python 3+
            r = r.decode('utf-8')

        # the exit status reflects the number of
        # test failures in the tested module
        return p.returncode, r
