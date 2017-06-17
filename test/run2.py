
# -*- coding: utf-8 -*-

# Script to run some or all PyGeodesy tests with Python 2 or 3.

# Tested with 64-bit Python 3.6.1 on macOS 10.12.5 Sierra and
# with Pythonista 3.1 on iOS 10.3.2.

from glob import glob
from os import linesep as NL
from os.path import basename, dirname, join
from subprocess import PIPE, STDOUT, Popen
from time import time
import sys

_dir = dirname(__file__)
_dir_dir = dirname(_dir)

try:
    import base as _   # PYCHOK expected
except ImportError:  # Pythonista
    if _dir not in sys.path:
        sys.path.insert(0, _dir)
from base import isiOS, versions  # PYCHOK expected

__all__ = ('run',)
__version__ = '17.06.16'

_python_O = _python = sys.executable  # PYCHOK false

if isiOS:  # MCCABE 14

    from runpy import run_path
    try:
        from StringIO import StringIO
    except ImportError:  # Python 3+
        from io import StringIO
    from traceback import format_exception

    _python_O = basename(_python_O)

    def _tilda(path):  # shorten Pythonista paths
        return path.replace(_dir_dir, '~')

    def run(test):
        # mimick partial behavior of function run
        # further below because subprocess.Popen
        # in Pythonista on iOS is not available.

        x = None  # no exception
        t3 = sys.argv, sys.stdout, sys.stderr

        sys.stdout = sys.stderr = std = StringIO()
        try:
            sys.argv = [test]
            run_path(test, run_name='__main__')
        except:  # PYCHOK must on Pythonista
            x = sys.exc_info()
            if x[0] is SystemExit:
                x = x[1].code  # exit status
            else:
                x = [_ for _ in format_exception(*x)  # PYCHOK expected
                             if 'runpy.py", line' not in _]
                print(''.join(map(_tilda, x)).strip())
                x = 99
        sys.argv, sys.stdout, sys.stderr = t3

        r = std.getvalue()
        if isinstance(r, bytes):  # Python 3+
            r = r.decode('utf-8')

        std.close()
        std = None  # del std

        if x is None:  # no exception or exit,
            # get the number of failed tests
            # and exclude any KNOWN failures
            x = r.count('FAILED, expected')
        return x, r

else:  # non iOS

    def _tilda(path):  # PYCHOK expected
        return path

    def run(test):  # PYCHOK expected
        '''Invoke one test and return the
           exit status and console output.
        '''
        c = [_python_O, test]
        p = Popen(c, creationflags=0,
                     executable   =_python,
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


if not __debug__:
    _python_O += ' -OO'

# command line options
_failedonly = False
_raiser     = False
_results    = False
_verbose    = False

if __name__ == '__main__':  # MCCABE 26

    def _write(text):
        _results.write(text.encode('utf-8'))

    argv0, args = _tilda(sys.argv[0]), sys.argv[1:]

    if isiOS and not args:
        # allow this script to be used
        # with options inside Pythonista
        try:
            _input = raw_input  # old name
        except NameError:  # Python 3+
            _input = input
        args = _input('enter %s args: ' % (argv0,)).split()

    while args and args[0].startswith('-'):
        arg = args.pop(0)
        if '-help'.startswith(arg):
            print('usage: %s [-failedonly] [-raiser] [-results] [-verbose] [test/test...py ...]' % (argv0,))
            sys.exit(0)
        elif '-failedonly'.startswith(arg):
            _failedonly = True
        elif '-raiser'.startswith(arg):
            _raiser = True  # break on error
        elif '-results'.startswith(arg):
            _results = True
        elif '-verbose'.startswith(arg):
            _verbose = True
        else:
            print('%s invalid option: %s' % (argv0, arg))
            sys.exit(1)

    # shorten prompt
    p = _python_O
    if len(p) > 32:
        p = p[:16] + '...' + p[-16:]

    # PyGeodesy and Python versions, size, OS name and release
    v = versions()

#   import pygeodesy
#   v = ' '.join((v, _tilda(pygeodesy.__file__)))

    if not args:  # no tests specified, get all test*.py
        # scripts in the same directory as this one
        args = sorted(glob(join(_dir, 'test*.py')))

        if _results:  # save all test results
            t = '-'.join(['testresults'] + v.split()) + '.txt'
            t = join(_dir_dir, 'testresults', t)
            _results = open(t, 'wb')  # note, 'b' not 't'!
            _write('%s typical test results (%s)%s' % (argv0, v, NL))

    X, s = 0, time()
    for arg in args:

        t = _tilda('running %s %s' % (p, arg))
        print(t)
        x, r = run(arg)
        X += x  # failures, excl KNOWN ones

        if _results:
            _write(NL + t + NL)
            _write(r)

        if 'Traceback' in r:
            print('%s\n' % (r,))
            if _raiser:
                break
        elif _failedonly:
            for t in r.split('\n'):
                if 'FAILED' in t or 'passed' in t:
                    print(t.rstrip())
            print('')
        elif _verbose:
            print('%s\n' % (r,))

    if X:
        x = '%d FAILED' % (X,)
    else:
        x = 'all OK'

    s = time() - s
    t = '%s %s %s (%s) %.3f sec' % (argv0, p, x, v, s)
    print(t)
    if _results:
        _write(NL + t + NL)
        _results.close()
