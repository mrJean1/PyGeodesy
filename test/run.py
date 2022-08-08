
# -*- coding: utf-8 -*-

# Script to run some or all PyGeodesy tests with Python 2 or 3.

from base import clips, coverage, isiOS, NN, PyGeodesy_dir, \
                 PythonX, secs2str, test_dir, _skipped_, _TILDE_, \
                 tilde, versions, _W_opts  # PYCHOK expected
from pygeodesy.basics import str2ub, ub2str

from os import access, environ, F_OK, linesep as NL
import sys

__all__ = ('run2',)
__version__ = '22.08.05'

if isiOS:  # MCCABE 14

    try:  # prefer StringIO over io
        from StringIO import StringIO
    except ImportError:  # Python 3+
        from io import StringIO
    from os.path import basename
    from runpy import run_path
    from traceback import format_exception

    def run2(test, *unused):
        '''Invoke one test module and return
           the exit status and console output.
        '''
        # Mimick partial behavior of function run2
        # further below because subprocess.Popen is
        # not available on iOS/Pythonista/Python.
        # One issue however, test scripts are all
        # imported and run in this same process.

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
                print(NN.join(map(tilde, x)).rstrip())
                x = 1  # count as a failure
        finally:
            sys.argv, sys.stdout, sys.stderr = sys3

        r = ub2str(std.getvalue())

        std.close()
        std = None  # del std

        if x is None:  # no exit status or exception:
            # count failed tests excluding KNOWN ones
            x = r.count('FAILED, expected')
        return x, r

    PythonX_ = basename(PythonX)

else:  # non-iOS
    from base import isPython37
    from subprocess import PIPE, Popen, STDOUT

    Popen_kwds = dict(creationflags=0, executable=sys.executable,
                    # shell=True,
                      stdin=None, stdout=PIPE, stderr=STDOUT)
    if isPython37:  # Python 3.7+
        Popen_kwds.update(text=True)

    pythonC_ = (PythonX,)  # python cmd tuple
    if not __debug__:
        pythonC_ += ('-O',)
    if _W_opts:  # include -W options
        pythonC_ += (_W_opts,)

    # replace home dir with _TILDE_
    PythonX_ = ' '.join(pythonC_).replace(environ.get('HOME', _TILDE_), _TILDE_)

    if coverage:
        pythonC_ += tuple('-m coverage run -a'.split())

    def run2(test, *opts):  # PYCHOK expected
        '''Invoke one test module and return the
           exit status and stdout/-err output.
        '''
        c = pythonC_ + (test,) + opts
        p = Popen(c, **Popen_kwds)
        r = ub2str(p.communicate()[0])  # no .strip()
        # the exit status reflects the number of
        # test failures in the tested module
        return p.returncode, r

# shorten Python path [-O]
PythonX_ = clips(PythonX_, 32)

# command line options
_failedonly = False
_prefix     = False
_raiser     = False
_results    = False  # or file
_verbose    = False

_Total = 0  # total tests
_FailX = 0  # failed tests


def _exit(last, text, exit):
    '''(INTERNAL) Close and exit.
    '''
    print(last)
    if _results:
        _write(NL + text + NL)
        _results.close()
    sys.exit(exit)


def _prefix2(prev):
    '''(INTERNAL) Get prefix and time stamp.
    '''
    t = time()
    if _prefix:
        p = '%7.3f ' % ((t - prev),)
        if p.startswith('  0.'):
            p = '   .' + p[4:]
    else:
        p = NN
    return p, t


def _run(prefix, test, *opts):  # MCCABE 13
    '''(INTERNAL) Run a test script and parse the result.
    '''
    global _Total, _FailX

    t = '%srunning %s %s' % (prefix, PythonX_, tilde(test))
    if access(test, F_OK):

        print(t)
        x, r = run2(test, *opts)
        r = r.replace(PyGeodesy_dir, '.')
        if _results:
            _write(NL + t + NL)
            _write(r)

        if 'Traceback' in r:
            print(r + NL)
            if not x:  # count as failure
                _FailX += 1
            if _raiser:
                raise SystemExit

        elif _failedonly:
            for t in _testlines(r, False):
                if ', KNOWN' not in t:
                    print(t)

        elif _verbose:
            print(r + NL)

        elif x:
            for t in _testlines(r, True):
                print(t)

    else:
        r = t + ' FAILED:  no such file' + NL
        x = 1
        if _results:
            _write(NL + r)
        print(r)

    _Total += r.count(NL + '    test ')  # number of tests
    _FailX += x  # failures, excluding KNOWN ones


def _testlines(r, skipped):
    '''(INTERNAL) Yield test lines.
    '''
    for t in r.split(NL):
        if 'FAILED,' in t or 'passed' in t or (skipped and _skipped_ in t):
            yield t.rstrip()
    yield NN


def _write(text):
    '''(INTERNAL) Write text to results.
    '''
    _results.write(str2ub(text))


if __name__ == '__main__':  # MCCABE 19

    from glob import glob
    from os.path import join
    from time import time

    argv0, args = tilde(sys.argv[0]), sys.argv[1:]

    while args and args[0].startswith('-'):
        arg = args.pop(0)
        if '-help'.startswith(arg):
            print('usage: %s [-B] [-failedonly] [-raiser] [-results] [-verbose] [-Z[0-9]] [test/test...py ...]' % (argv0,))
            sys.exit(0)
        elif arg.startswith('-B'):
            environ['PYTHONDONTWRITEBYTECODE'] = arg[2:]
        elif '-failedonly'.startswith(arg):
            _failedonly = True
        elif '-prefix'.startswith(arg):
            _prefix = True
        elif '-raiser'.startswith(arg):  # also .base.Test.__init__
            _raiser = True  # break on error
        elif '-results'.startswith(arg):
            _results = True
        elif '-verbose'.startswith(arg):
            _verbose = True
        elif arg.startswith('-Z'):
            environ['PYGEODESY_LAZY_IMPORT'] = arg[2:]
        else:
            print('%s invalid option: %s' % (argv0, arg))
            sys.exit(1)

    if not args:  # no tests specified, get all test*.py
        # scripts in the same directory as this one
        args = sorted(glob(join(test_dir, 'test[A-Z]*.py')))

    # PyGeodesy and Python versions, size, OS name and release
    v = versions()

    if _results:  # save all test results
        t = '-'.join(['testresults'] + v.split()) + '.txt'
        t = join(PyGeodesy_dir, 'testresults', t)
        _results = open(t, 'wb')  # note, 'b' not 't'!
        _write('%s typical test results (%s)%s' % (argv0, v, NL))

    s = t = time()
    try:
        for arg in args:
            p, t = _prefix2(t)
            _run(p, *arg.split())
    except KeyboardInterrupt:
        _exit(NN, '^C', 9)
    except SystemExit:
        pass
    p, t = _prefix2(t)
    s    =  t - s
    t    =  secs2str(s)
    if _Total > s > 1:
        t = '%s (%.3f tps)' % (t, _Total / s)

    if _FailX:
        s = NN if _FailX == 1 else 's'
        x = '%d (of %d) test%s FAILED' % (_FailX, _Total, s)
    elif _Total > 0:
        x = 'all %d tests OK' % (_Total,)
    else:
        x = 'all OK'

    t = '%s%s %s: %s (%s) %s' % (p, argv0, PythonX_, x, v, t)
    _exit(t, t, 2 if _FailX else 0)

# **) MIT License
#
# Copyright (C) 2016-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
