
# -*- coding: utf-8 -*-

# Script to run some or all PyGeodesy tests with Python 2 or 3.

# Tested with 64-bit Python 2.6.9, 2.7.13, 3.5.3 and 3.6.0 but only
# on MacOS 10.10 Mavericks, MacOS 10.11 El Capitan and MacOS 10.12.2
# and 10.12.3 Sierra.

from glob import glob
from os import linesep as NL
from os.path import dirname, join
from platform import architecture, mac_ver, win32_ver, uname
from subprocess import PIPE, STDOUT, Popen
import sys

__all__ = ('run',)
__version__ = '17.03.09'

_python_O = _python = sys.executable  # path
if not __debug__:
    _python_O += ' -OO'


def run(test):
    '''Invoke a test and return the
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


_failedonly = False
_raiser = False
_results = False
_verbose = False

if __name__ == '__main__':

    def _write(text):
        if _results:
            _results.write(text.encode('utf-8'))

    argv0, args = sys.argv[0], sys.argv[1:]
    while args and args[0].startswith('-'):
        arg = args.pop(0)
        if '-help'.startswith(arg):
            print('usage: %s [-failedonly] [-raiser] [-results] [-verbose] [tests/test...py ...]' % (argv0,))
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

    # get Python version and OS name and release
    v = 'Python', sys.version.split()[0], architecture()[0]
    for n, t in (('MacOS',   mac_ver),
                 ('Windows', win32_ver),
                 ('Uname',   uname)):
        t = t()[0]
        if t:
            v += n, t
            break
    v = ' '.join(v)

    if not args:  # no tests specified, get all test*.py
        # scripts in the same directory as this one
        args = sorted(glob(join(dirname(__file__), 'test*.py')))

        if _results:  # save all test results
            t = '-'.join(['testresults'] + v.split()) + '.txt'
            _results = open(join('testresults', t), 'wb')  # note, 'b' not 't'!
            _write('%s typical test results (%s)%s' % (argv0, v, NL))

    f = 0
    for arg in args:

        t = 'running %s %s' % (p, arg)
        print(t)
        x, r = run(arg)
        f += x  # failures, excl KNOWN ones
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

    if f:
        x = '%d FAILED' % (f,)
    else:
        x = 'all OK'
    t = '%s %s %s (%s)' % (argv0, p, x, v)
    print(t)
    if _results:
        _write(NL + t + NL)
        _results.close()
