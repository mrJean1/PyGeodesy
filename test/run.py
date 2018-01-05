
# -*- coding: utf-8 -*-

# Script to run some or all PyGeodesy tests with Python 2 or 3.

# Tested with 64-bit Python 2.6.9, 2.7,13, 3.5.3 and 3.6.4 on
# macOS 10.12.5 and 10.12.6 Sierra and with Pythonista 3.1 and
# 3.2 on iOS 10.3.2, 10.3.3, 11.0.3, 11.1.2 and 11.2.1.

from glob import glob
from os import environ, linesep as NL
from os.path import abspath, dirname, join
from time import time
import sys

_test_dir = dirname(abspath(__file__))
# extend sys.path to include the .. directory
if _test_dir not in sys.path:
    sys.path.insert(0, _test_dir)

from base import isiOS, PyGeodesy_dir, Python_O, \
          runner, secs2str, tilde, versions  # PYCHOK expected

__all__ = ()
__version__ = '18.01.04'

# command line options
_failedonly = False
_raiser     = False
_results    = False
_verbose    = False

if __name__ == '__main__':  # MCCABE 28

    def _write(text):
        _results.write(text.encode('utf-8'))

    argv0, args = tilde(sys.argv[0]), sys.argv[1:]

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

    # replace home dir with ~
    Python_O = Python_O.replace(environ.get('HOME', '~'), '~')

    # shorten Python path [-OO]
    if len(Python_O) > 32:
        Python_O = Python_O[:16] + '...' + Python_O[-16:]

    # PyGeodesy and Python versions, size, OS name and release
    v = versions()

#   import pygeodesy
#   v = ' '.join((v, tilde(pygeodesy.__file__)))

    if _results:  # save all test results
        t = '-'.join(['testresults'] + v.split()) + '.txt'
        t = join(PyGeodesy_dir, 'testresults', t)
        _results = open(t, 'wb')  # note, 'b' not 't'!
        _write('%s typical test results (%s)%s' % (argv0, v, NL))

    if not args:  # no tests specified, get all test*.py
        # scripts in the same directory as this one
        args = sorted(glob(join(_test_dir, 'test*.py')))

    T, X, s = 0, 0, time()
    for arg in args:

        t = 'running %s %s' % (Python_O, tilde(arg))
        print(t)

        x, r = runner(arg)
        X += x  # failures, excluding KNOWN ones

        if _results:
            _write(NL + t + NL)
            _write(r)

        if not X:  # count the tests
            T += r.count(NL + '    test ')

        if 'Traceback' in r:
            print(r + NL)
            if not x:  # count as failure
                X += 1
            if _raiser:
                break

        elif _failedonly:
            for t in r.split(NL):
                # print failures, KNOWN ones and totals
                if 'FAILED,' in t or 'passed' in t:
                    print(t.rstrip())
            print('')

        elif _verbose:
            print(r + NL)

        elif x:
            for t in r.split(NL):
                # print failures, KNOWN ones and totals
                if 'FAILED,' in t and 'KNOWN' not in t:
                    print(t.rstrip())

    if X:
        x = '%d FAILED' % (X,)
    elif T > 0:
        x = 'all %s tests OK' % (T,)
    else:
        x = 'all OK'

    s = secs2str(time() - s)
    t = '%s %s %s (%s) %s' % (argv0, Python_O, x, v, s)
    print(t)
    if _results:
        _write(NL + t + NL)
        _results.close()
