
# -*- coding: utf-8 -*-

# Module to run all PyGeodesy tests as  python setup.py test

# Tested with 64-bit Python 2.7.13 and 3.6.0 on macOS 10.12.3,
# 10.12.4 and 10.12.5 Sierra.

try:
    from run import run
except ImportError:  # Python 3+ ModuleNotFoundError
    from .run import run  # PYCHOK expected

from os.path import dirname, join
import unittest

__all__ = ('TestSuite', 'run')
__version__ = '17.05.25'

_tests_dir = dirname(__file__)


class TestSuite(unittest.TestCase):
    '''Combine all test modules into a test suite/case
       and run each test module as a separate test.

       XXX rewrite the test modules to use the unittest.
    '''

    def _run(self, test):
        x, _ = run(join(_tests_dir, test + '.py'))
        self.assertEqual(x, 0)

    def test_Bases(self):
        self._run('testBases')

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

    def test_Lcc(self):
        self._run('testLcc')

    def test_Mgrs(self):
        self._run('testMgrs')

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

    def test_Tests(self):
        self._run('tests')
