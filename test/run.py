
# -*- coding: utf-8 -*-
# ruff: noqa: E722

# Script to run some or all PyGeodesy tests with Python 2 or 3.

from bases import clips, coverage, _DOT_, isiOS, NN, prefix2, printf, \
                  PyGeodesy_dir, PythonX, secs2str, _skipped_, test_dir, \
                  tilde, _TILDE_, versions, _W_opts  # PYCHOK expected
from pygeodesy.basics import str2ub, ub2str

from os import access, environ, F_OK, linesep as LS, pathsep as PS
import sys
from traceback import format_exception

__all__ = ('run2',)
__version__ = '25.10.06'

NL = '\n'  # pygeodesy.interns._NL_
P  = None  # process

if isiOS:  # MCCABE 14

    try:  # prefer StringIO over io
        from StringIO import StringIO
    except ImportError:  # Python 3+
        from io import StringIO
    from os.path import basename
    from runpy import run_path

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
                printf(NN.join(map(tilde, x)).rstrip())
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
    _threads = False

else:  # non-iOS
    from bases import isAppleSi, isPython37
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
        global P
        c = pythonC_ + (test,) + opts
        P = Popen(c, **Popen_kwds)
        r = ub2str(P.communicate()[0])  # no .strip()
        # the exit status reflects the number of
        # test failures in the tested module
        return P.returncode, r

    _threads = isAppleSi and sys.stdout.isatty()

if _threads:
    from threading import Lock, Thread
else:  # non-threaded version
    class Lock(object):
        def __enter__(self):
            pass

        def __exit__(self, *unused):  # PYCHOK 4 args
            pass


class ThreadPool(object):
    '''Partially mimick standard U{ThreadPoolExecutor
       <https://docs.Python.org/3/library/concurrent.futures.html>}.
    '''
    def __init__(self, size=9):
        self.lock = Lock()
        self.size = min(size, 9) if _threads and size > 1 else 0
        self.work = 0

    def onedown(self):
        if self.work > 0:  # with self.lock
            self.work -= 1

    def shutdown(self, prefix=NN):
        w = self.work
        if w > 0:
            t = '%srunning %s' % (prefix, w)
            while w > 0:
                sleep(2)
                printf(t, end=NN, flush=True)
                v, w, t = w, self.work, _DOT_
                if w != v:
                    t += str(w)
            printf(NN)

    def submit(self, call, *args):
        z = self.size
        if z > 0:
            while self.work >= z:
                sleep(1)
            with self.lock:
                self.work += 1
            t = Thread(target=call, args=args)
            t.start()
        else:
            global P
            t = call(*args)
            P = None
        return t


p = environ.get('PYTHONPATH', '')
if p:  # prepend
    p = PS + p
environ['PYTHONPATH'] = PyGeodesy_dir + PS + test_dir + p
del p

# shorten Python path [-O]
PythonX_ = clips(PythonX_, 32)

# command line options
_failedonly = False
_prefix     = False
_raiser     = False
_results    = False  # or file
_verbose    = False

_Tpool = ThreadPool(0)
_Total = 0  # total tests
_FailX = 0  # failed tests


def _exit(last, text, exit):
    '''(INTERNAL) Close and exit.
    '''
    printf(last)
    if _results:
        _write(NL + text + NL)
        _results.close()
    if P:
        P.kill()
    sys.exit(exit)


def _prefix2(prev):
    '''(INTERNAL) get the prefix and current time.
    '''
    return prefix2(prev) if _prefix else (NN, time())


def _run(prefix, test, *opts):  # MCCABE 13
    '''(INTERNAL) Run a test script and parse the result.
    '''
    global _FailX, _Total  # _Tpool

    t = PythonX_
    if _Tpool.size > 0:
        t = '%s %s' % (_Tpool.work, t)
    t = '%srunning %s %s' % (prefix, t, tilde(test))
    if access(test, F_OK):
        printf(t)
        x,  r = run2(test, *opts)
        p = r = r.replace(PyGeodesy_dir, '.')
        if 'Traceback' in r:
            p += NL
            if not x:  # count as failure
                with _Tpool.lock:
                    _FailX += 1
            if _raiser:
                raise SystemExit
        elif _failedonly:
            p = NL.join(_ for _ in _testlines(r, False)
                                if ', KNOWN' not in _)
        elif _verbose:
            p += NL
        elif x:
            p = NL.join(_ for _ in _testlines(r, True))
        else:
            p = NN
    else:
        t += ' FAILED:  no such file'
        p = t + NL
        r = NN
        x = 1

    with _Tpool.lock:
        if p:
            printf(p)
        if _results:
            _write(t + NL + r)

        _Total += r.count(NL + '    test ')  # number of tests
        _FailX += x or 0  # failures, excluding KNOWN ones

        _Tpool.onedown()


def _testlines(r, skipped):
    '''(INTERNAL) Yield test lines.
    '''

    for t in r.split(LS if LS in r else NL):  # use NL on Windows, not LS
        if 'FAILED,' in t or 'passed' in t \
                          or (skipped and _skipped_ in t):
            yield t.rstrip()
    yield NN


def _write(text):
    '''(INTERNAL) Write text to results.
    '''
    _results.write(str2ub(text))


if __name__ == '__main__':  # MCCABE 19

    from glob import glob
    from os.path import join
    from time import sleep, time

    argv0, args = tilde(sys.argv[0]), sys.argv[1:]

    while args and args[0].startswith('-'):
        arg = args.pop(0)
        if '-help'.startswith(arg):
            printf('usage: %s [-B] [-failedonly] [-raiser] [-results]%s [-verbose] [-Z[0-9]] [test/test...py ...]',
                    argv0, (' [-threads <int>]' if _threads else NN))
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
        elif _threads and '-threads'.startswith(arg) \
                      and args and args[0].isdigit():
            _Tpool = ThreadPool(int(args.pop(0)))
        elif '-verbose'.startswith(arg):
            _verbose = True
        elif arg.startswith('-Z'):
            environ['PYGEODESY_LAZY_IMPORT'] = arg[2:]
        else:
            printf('%s invalid option: %s', argv0, arg)
            sys.exit(1)

    if not args:  # no tests specified, get all test*.py
        # scripts in the same directory as this one
        args = sorted(glob(join(test_dir, 'test[A-Z]*.py')))

    # PyGeodesy and Python versions, size, OS name and release
    v = versions()

    if _results:  # save all test results
        t = '-'.join(['testresults'] + v.split()) + '.txt'
        t = clips(t, 180)  # isWindows MAX?
        t = join(PyGeodesy_dir, 'testresults', t)
        _results = open(t, 'wb')  # note, 'b' not 't'!
        _write('%s typical test results (%s)%s' % (argv0, v, NL))

    s = t = time()
    try:
        for arg in args:
            p, t = _prefix2(t)
            _Tpool.submit(_run, p, *arg.split())
        p, t = _prefix2(t)
        _Tpool.shutdown(p)
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
# Copyright (C) 2016-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
