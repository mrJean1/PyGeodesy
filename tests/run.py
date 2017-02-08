
# -*- coding: utf-8 -*-

# Script to run some or all PyGeodesy tests with Python 2 or 3.

# Tested with 64-bit Python 2.6.9, 2.7.10, 2.7.13, 3.5.2, 3.5.2
# and 3.6.0 but only on MacOS 10.10 Mavericks, MacOS 10.11 El
# Capitan and MacOS 10.12.2 Sierra.

from glob import glob
from os.path import dirname, join
from platform import architecture
from subprocess import PIPE, STDOUT, Popen
import sys

__version__ = '17.02.05'

_python = sys.executable  # path

_failedonly = False
_raiser = False
_verbose = False

argv0, args = sys.argv[0], sys.argv[1:]
while args and args[0].startswith('-'):
    arg = args.pop(0)
    if '-help'.startswith(arg):
        print('usage: %s [-failedonly] [-raiser] [-verbose] [tests/test...py ...]' % (argv0,))
        sys.exit(0)
    elif '-failedonly'.startswith(arg):
        _failedonly = True
    elif '-raiser'.startswith(arg):
        _raiser = True  # break on error
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

    if 'Traceback' in r:
        print('%s' % (r,))
        if _raiser:
            break
    elif _failedonly:
        for t in r.split('\n'):
            if 'FAILED' in t or 'passed' in t:
                print(t.rstrip())
    elif _verbose:
        print('%s' % (r,))
    else:
        continue

    print('')

v = 'Python %s %s' % (sys.version.split()[0], architecture()[0])
if f:
    print('%s %s %s FAILED (%s)' % (argv0, _python, f, v))
else:
    print('%s %s all OK (%s)' % (argv0, _python, v))


# Typical test results (on MacOS X):

# % /usr/bin/python tests/run.py
# tests/run.py /usr/bin/python tests/testBases.py
# tests/run.py /usr/bin/python tests/testDatum.py
# tests/run.py /usr/bin/python tests/testDms.py
# tests/run.py /usr/bin/python tests/testEllipsoidal.py
# tests/run.py /usr/bin/python tests/testGreatCircle.py
# tests/run.py /usr/bin/python tests/testLcc.py
# tests/run.py /usr/bin/python tests/testMgrs.py
# tests/run.py /usr/bin/python tests/testNavlabExamples.py
# tests/run.py /usr/bin/python tests/testOsgr.py
# tests/run.py /usr/bin/python tests/testSpherical.py
# tests/run.py /usr/bin/python tests/testUtm.py
# tests/run.py /usr/bin/python tests/tests.py
# tests/run.py /usr/bin/python all OK (Python 2.7.10 64bit)


# % /usr/bin/python tests/run.py -h
# usage: tests/run.py [-failedonly] [-raiser] [-verbose] [tests/test...py ...]


# % /usr/bin/python tests/run.py -f
# tests/run.py /usr/bin/python tests/testBases.py
    # all geodesy.bases tests passed (Python 2.7.10 64bit)

# tests/run.py /usr/bin/python tests/testDatum.py
    # all geodesy.datum tests passed (Python 2.7.10 64bit)

# tests/run.py /usr/bin/python tests/testDms.py
    # all geodesy.dms tests passed (Python 2.7.10 64bit)

# tests/run.py /usr/bin/python tests/testEllipsoidal.py
    # all geodesy.ellipsoidalNvector tests passed (Python 2.7.10 64bit)
    # all geodesy.ellipsoidalVincenty tests passed (Python 2.7.10 64bit)

# tests/run.py /usr/bin/python tests/testGreatCircle.py
    # test 7 DistanceEiffelToVersailles: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 8 DistanceVersaillesToEiffel: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 21 MidpointEiffelToVersailles(m): 7042.15004788  FAILED, KNOWN, expected 7042.1597433
    # test 22 MidpointVersaillesToEiffel: 48.831495°N, 002.207536°E  FAILED, KNOWN, expected 48.831495°N, 002.207535°E
    # test 24 MidpointVersaillesToEiffel(m): 7042.15004788  FAILED, KNOWN, expected 7042.1597433
    # 5 geodesy.sphericalNvector tests (17.2%) FAILED, incl. 5 KNOWN (Python 2.7.10 64bit)
    # test 7 DistanceEiffelToVersailles: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 8 DistanceVersaillesToEiffel: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 21 MidpointEiffelToVersailles(m): 7042.15004788  FAILED, KNOWN, expected 7042.1597433
    # test 22 MidpointVersaillesToEiffel: 48.831495°N, 002.207536°E  FAILED, KNOWN, expected 48.831495°N, 002.207535°E
    # test 24 MidpointVersaillesToEiffel(m): 7042.15004788  FAILED, KNOWN, expected 7042.1597433
    # 5 geodesy.sphericalTrigonometry tests (17.2%) FAILED, incl. 5 KNOWN (Python 2.7.10 64bit)

# tests/run.py /usr/bin/python tests/testLcc.py
    # all testLcc.py tests passed (Python 2.7.10 64bit)

# tests/run.py /usr/bin/python tests/testMgrs.py
    # all geodesy.mgrs tests passed (Python 2.7.10 64bit)

# tests/run.py /usr/bin/python tests/testNavlabExamples.py
    # test 10 Example 2 destinationPoint: 53.327726°N, 063.464965°E, +301.02m  FAILED, KNOWN, expected 53.327726°N, 063.464965°E, +299.138m
    # 1 testNavlabExamples.py test (3.8%) FAILED, incl. 1 KNOWN (Python 2.7.10 64bit)

# tests/run.py /usr/bin/python tests/testOsgr.py
    # test 6 toLatLon1: 52°39′28.72″N, 001°43′00.63″E  FAILED, KNOWN, expected 52°39′28.72″N, 001°42′57.74″E
    # test 7 toLatLon1: 52.657979°N, 001.716843°E  FAILED, KNOWN, expected 52.657977°N, 001.716038°E
    # test 8 toOsgr1: 651463,313180  FAILED, KNOWN, expected 651409.903, 313177.270
    # test 9 toLatLon2: 52°39′27.25″N, 001°43′07.37″E  FAILED, KNOWN, expected 52°39′27.25″N, 001°43′04.47″E
    # test 10 toLatLon2: 52.65757°N, 001.718713°E  FAILED, KNOWN, expected 52.657568°N, 001.717908°E
    # test 11 toOsgr2: 651463,313180  FAILED, KNOWN, expected 651409,313177
    # 6 geodesy.osgr tests (25.0%) FAILED, incl. 6 KNOWN (Python 2.7.10 64bit)

# tests/run.py /usr/bin/python tests/testSpherical.py
    # all geodesy.sphericalNvector tests passed (Python 2.7.10 64bit)
    # all geodesy.sphericalTrigonometry tests passed (Python 2.7.10 64bit)

# tests/run.py /usr/bin/python tests/testUtm.py
    # all geodesy.utm tests passed (Python 2.7.10 64bit)

# tests/run.py /usr/bin/python tests/tests.py
    # all tests.py tests passed (Python 2.7.10 64bit)

# tests/run.py /usr/bin/python all OK (Python 2.7.10 64bit)


# % /usr/local/bin/python2.7 tests/run.py
# tests/run.py /usr/local/bin/python2.7 tests/testBases.py
# tests/run.py /usr/local/bin/python2.7 tests/testDatum.py
# tests/run.py /usr/local/bin/python2.7 tests/testDms.py
# tests/run.py /usr/local/bin/python2.7 tests/testEllipsoidal.py
# tests/run.py /usr/local/bin/python2.7 tests/testGreatCircle.py
# tests/run.py /usr/local/bin/python2.7 tests/testLcc.py
# tests/run.py /usr/local/bin/python2.7 tests/testMgrs.py
# tests/run.py /usr/local/bin/python2.7 tests/testNavlabExamples.py
# tests/run.py /usr/local/bin/python2.7 tests/testOsgr.py
# tests/run.py /usr/local/bin/python2.7 tests/testSpherical.py
# tests/run.py /usr/local/bin/python2.7 tests/testUtm.py
# tests/run.py /usr/local/bin/python2.7 tests/tests.py
# tests/run.py /usr/local/bin/python2.7 all OK (Python 2.7.13 64bit)


# % /usr/local/python3.5 tests/run.py
# tests/run.py /usr/local/bin/python3.5 tests/testBases.py
# tests/run.py /usr/local/bin/python3.5 tests/testDatum.py
# tests/run.py /usr/local/bin/python3.5 tests/testDms.py
# tests/run.py /usr/local/bin/python3.5 tests/testEllipsoidal.py
# tests/run.py /usr/local/bin/python3.5 tests/testGreatCircle.py
# tests/run.py /usr/local/bin/python3.5 tests/testLcc.py
# tests/run.py /usr/local/bin/python3.5 tests/testMgrs.py
# tests/run.py /usr/local/bin/python3.5 tests/testNavlabExamples.py
# tests/run.py /usr/local/bin/python3.5 tests/testOsgr.py
# tests/run.py /usr/local/bin/python3.5 tests/testSpherical.py
# tests/run.py /usr/local/bin/python3.5 tests/testUtm.py
# tests/run.py /usr/local/bin/python3.5 tests/tests.py
# tests/run.py /usr/local/bin/python3.5 all OK (Python 3.5.3 64bit)


# % python3.6 tests/run.py
# tests/run.py /usr/local/bin/python3.6 tests/testBases.py
# tests/run.py /usr/local/bin/python3.6 tests/testDatum.py
# tests/run.py /usr/local/bin/python3.6 tests/testDms.py
# tests/run.py /usr/local/bin/python3.6 tests/testEllipsoidal.py
# tests/run.py /usr/local/bin/python3.6 tests/testGreatCircle.py
# tests/run.py /usr/local/bin/python3.6 tests/testLcc.py
# tests/run.py /usr/local/bin/python3.6 tests/testMgrs.py
# tests/run.py /usr/local/bin/python3.6 tests/testNavlabExamples.py
# tests/run.py /usr/local/bin/python3.6 tests/testOsgr.py
# tests/run.py /usr/local/bin/python3.6 tests/testSpherical.py
# tests/run.py /usr/local/bin/python3.6 tests/testUtm.py
# tests/run.py /usr/local/bin/python3.6 tests/tests.py
# tests/run.py /usr/local/bin/python3.6 all OK (Python 3.6.0 64bit)


# % python3.6 tests/run.py -h
# usage: tests/run.py [-failedonly] [-raiser] [-verbose] [tests/test...py ...]


# % python3.6 tests/run.py -f
# tests/run.py /usr/local/bin/python3 tests/testBases.py
    # all geodesy.bases tests passed (Python 3.6.0 64bit)

# tests/run.py /usr/local/bin/python3 tests/testDatum.py
    # all geodesy.datum tests passed (Python 3.6.0 64bit)

# tests/run.py /usr/local/bin/python3 tests/testDms.py
    # all geodesy.dms tests passed (Python 3.6.0 64bit)

# tests/run.py /usr/local/bin/python3 tests/testEllipsoidal.py
    # all ellipsoidalNvector tests passed (Python 3.6.0 64bit)
    # all ellipsoidalVincenty tests passed (Python 3.6.0 64bit)

# tests/run.py /usr/local/bin/python3 tests/testGreatCircle.py
    # test 7 DistanceEiffelToVersailles: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 8 DistanceVersaillesToEiffel: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 21 MidpointEiffelToVersailles(m): 7042.150047881914  FAILED, KNOWN, expected 7042.159743304898
    # test 22 MidpointVersaillesToEiffel: 48.831495°N, 002.207536°E  FAILED, KNOWN, expected 48.831495°N, 002.207535°E
    # test 24 MidpointVersaillesToEiffel(m): 7042.150047882616  FAILED, KNOWN, expected 7042.159743304495
    # 5 sphericalNvector tests (17.2%) FAILED, incl. 5 KNOWN (Python 3.6.0 64bit)
    # test 7 DistanceEiffelToVersailles: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 8 DistanceVersaillesToEiffel: 14084.3001  FAILED, KNOWN, expected 14084.2807
    # test 21 MidpointEiffelToVersailles(m): 7042.150047881641  FAILED, KNOWN, expected 7042.159743304608
    # test 22 MidpointVersaillesToEiffel: 48.831495°N, 002.207536°E  FAILED, KNOWN, expected 48.831495°N, 002.207535°E
    # test 24 MidpointVersaillesToEiffel(m): 7042.150047882363  FAILED, KNOWN, expected 7042.159743304893
    # 5 sphericalTrigonometry tests (17.2%) FAILED, incl. 5 KNOWN (Python 3.6.0 64bit)

# tests/run.py /usr/local/bin/python3 tests/testLcc.py
    # all testLcc.py tests passed (Python 3.6.0 64bit)

# tests/run.py /usr/local/bin/python3 tests/testMgrs.py
    # all geodesy.mgrs tests passed (Python 3.6.0 64bit)

# tests/run.py /usr/local/bin/python3 tests/testNavlabExamples.py
    # test 10 Example 2 destinationPoint: 53.327726°N, 063.464965°E, +301.02m  FAILED, KNOWN, expected 53.327726°N, 063.464965°E, +299.138m
    # 1 testNavlabExamples.py test (3.8%) FAILED, incl. 1 KNOWN (Python 3.6.0 64bit)

# tests/run.py /usr/local/bin/python3 tests/testOsgr.py
    # test 6 toLatLon1: 52°39′28.72″N, 001°43′00.63″E  FAILED, KNOWN, expected 52°39′28.72″N, 001°42′57.74″E
    # test 7 toLatLon1: 52.657979°N, 001.716843°E  FAILED, KNOWN, expected 52.657977°N, 001.716038°E
    # test 8 toOsgr1: 651463,313180  FAILED, KNOWN, expected 651409.903, 313177.270
    # test 9 toLatLon2: 52°39′27.25″N, 001°43′07.37″E  FAILED, KNOWN, expected 52°39′27.25″N, 001°43′04.47″E
    # test 10 toLatLon2: 52.65757°N, 001.718713°E  FAILED, KNOWN, expected 52.657568°N, 001.717908°E
    # test 11 toOsgr2: 651463,313180  FAILED, KNOWN, expected 651409,313177
    # 6 geodesy.osgr tests (25.0%) FAILED, incl. 6 KNOWN (Python 3.6.0 64bit)

# tests/run.py /usr/local/bin/python3 tests/testSpherical.py
    # all sphericalNvector tests passed (Python 3.6.0 64bit)
    # all sphericalTrigonometry tests passed (Python 3.6.0 64bit)

# tests/run.py /usr/local/bin/python3 tests/testUtm.py
    # all geodesy.utm tests passed (Python 3.6.0 64bit)

# tests/run.py /usr/local/bin/python3 tests/tests.py
    # all tests.py tests passed (Python 3.6.0 64bit)

# tests/run.py /usr/local/bin/python3 all OK (Python 3.6.0 64bit)
