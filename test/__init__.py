
# -*- coding: utf-8 -*-

# Module to run all PyGeodesy tests as  python setup.py test

# Tested with 64-bit Python 2.7.13 and 3.6.1 on macOS 10.12.3,
# 10.12.4 and 10.12.5 Sierra.

from os.path import dirname, join
import unittest

_test_dir = dirname(__file__)

try:
    from run import run
except ImportError:  # Python 3+ ModuleNotFoundError
    import sys
    if _test_dir not in sys.path:
        sys.path.insert(0, _test_dir)
    from run import run  # PYCHOK expected

__all__ = ('TestSuite', 'run')
__version__ = '17.06.16'


class TestSuite(unittest.TestCase):
    '''Combine all test modules into a test suite/case
       and run each test module as a separate test.

       XXX rewrite the test modules to use the unittest.
    '''
    _runs = 1  # test_Ztotal is 1 run

    def _run(self, test):
        TestSuite._runs += 1  # pseudo global
        x, _ = run(join(_test_dir, test + '.py'))
        self.assertEqual(x, 0)

    def test_Bases(self):
        self._run('testBases')

    def test_Classes(self):
        self._run('testClasses')

    def test_Datum(self):
        self._run('testDatum')

    def test_Dms(self):
        self._run('testDms')

    def test_Ellipsoidal(self):
        self._run('testEllipsoidal')

    def test_Geohash(self):
        self._run('testGeohash')

    def test_GreatCircle(self):
        self._run('testGreatCircle')

    def test_LatLon(self):
        self._run('testLatLon')

    def test_Lcc(self):
        self._run('testLcc')

    def test_Mgrs(self):
        self._run('testMgrs')

    def test_Modules(self):
        self._run('testModules')

    def test_NavlabExamples(self):
        self._run('testNavlabExamples')

    def test_Osgr(self):
        self._run('testOsgr')

    def test_Simplify(self):
        self._run('testSimplify')

    def test_Spherical(self):
        self._run('testSpherical')

    def test_Utm(self):
        self._run('testUtm')

    def test_Vectorial(self):
        self._run('testVectorial')

    def test_Ztotal(self):
        # final test to make sure all tests were run
        t = sum(1 for t in dir(TestSuite) if t.startswith('test_'))
        self.assertEqual(TestSuite._runs, t)
