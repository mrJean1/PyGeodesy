
# -*- coding: utf-8 -*-

# Script to run some or all PyGeodesy tests with Python 2 or 3.

# Tested only on MacOS 10.10 Mavericks and 10.11 El Capitan with
# Python 2.7.10 and Python 3.6.0rc1.

from glob import glob
from os.path import dirname, join
from subprocess import PIPE, STDOUT, Popen
import sys

__version__ = '16.12.09'

_python = sys.executable  # path

_failedonly = False
_verbose = False

argv0, args = sys.argv[0], sys.argv[1:]
while args and args[0].startswith('-'):
    arg = args.pop(0)
    if '-help'.startswith(arg):
        print('usage: %s [-failedonly] [-verbose] [tests/test...py ...]' % (argv0,))
        sys.exit(0)
    elif '-failedonly'.startswith(arg):
        _failedonly = True
    elif '-verbose'.startswith(arg):
        _verbose = True
    else:
        print('%s invalid option: %s' % (argv0, arg))
        sys.exit(1)

if not args:  # no tests specified, get all test*.py
    # scripts in the same directory as this one
    args = sorted(glob(join(dirname(__file__), 'test*.py')))

f = 0
for arg in args:
    print('%s %s %s' % (argv0, _python, arg))

    cmd = [_python, arg]
    p = Popen(cmd, creationflags=0,
                   executable   =_python,
                 # shell        =True,
                   stdin        =None,
                   stdout       =PIPE,  # XXX
                   stderr       =STDOUT)  # XXX

    r = p.communicate()[0]
    if isinstance(r, bytes):  # Python 3+
        r = r.decode('utf-8')

    f += p.returncode  # failures

    if _failedonly:
        for t in r.split('\n'):
            if 'FAILED' in t or 'passed' in t:
                print(t.rstrip())
    elif _verbose:
        print('%s' % (r,))
    else:
        continue

    print('')

v = '(Python %s)' % (sys.version.split()[0],)
if f:
    print('%s %s %s FAILED %s' % (argv0, _python, f, v))
else:
    print('%s %s all OK %s' % (argv0, _python, v))


# Typical test results (on MacOS X):

# % python tests/run.py
# tests/run.py /usr/bin/python tests/testBases.py
# tests/run.py /usr/bin/python tests/testDatum.py
# tests/run.py /usr/bin/python tests/testDms.py
# tests/run.py /usr/bin/python tests/testEllipsoidal.py
# tests/run.py /usr/bin/python tests/testGreatCircle.py
# tests/run.py /usr/bin/python tests/testLambert.py
# tests/run.py /usr/bin/python tests/testMgrs.py
# tests/run.py /usr/bin/python tests/testNavlabExamples.py
# tests/run.py /usr/bin/python tests/testOsgr.py
# tests/run.py /usr/bin/python tests/testSpherical.py
# tests/run.py /usr/bin/python tests/testUtm.py
# tests/run.py /usr/bin/python tests/tests.py
# tests/run.py /usr/bin/python all OK (Python 2.7.10)

# % python3 tests/run.py
# tests/run.py /usr/local/bin/python3 tests/testBases.py
# tests/run.py /usr/local/bin/python3 tests/testDatum.py
# tests/run.py /usr/local/bin/python3 tests/testDms.py
# tests/run.py /usr/local/bin/python3 tests/testEllipsoidal.py
# tests/run.py /usr/local/bin/python3 tests/testGreatCircle.py
# tests/run.py /usr/local/bin/python3 tests/testLambert.py
# tests/run.py /usr/local/bin/python3 tests/testMgrs.py
# tests/run.py /usr/local/bin/python3 tests/testNavlabExamples.py
# tests/run.py /usr/local/bin/python3 tests/testOsgr.py
# tests/run.py /usr/local/bin/python3 tests/testSpherical.py
# tests/run.py /usr/local/bin/python3 tests/testUtm.py
# tests/run.py /usr/bin/python tests/tests.py
# tests/run.py /usr/local/bin/python3 all OK (Python 3.6.0rc1)
