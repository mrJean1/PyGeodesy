
# -*- coding: utf-8 -*-

# Module to run all PyGeodesy tests as  python setup.py test

from glob import glob
from os.path import abspath, dirname, join
import sys
import unittest

_test_dir = dirname(abspath(__file__))
# extend sys.path to include the ../.. directory
if _test_dir not in sys.path:  # Python 3+ ModuleNotFoundError
    sys.path.insert(0, _test_dir)

from base import runner

__all__ = ('TestSuite',)
__version__ = '18.02.02'


class TestSuite(unittest.TestCase):
    '''Combine all test modules into a test suite/case
       and run each test module as a separate test.
    '''
    _runs = 0  # pseudo global

    def _run(self, test):
        TestSuite._runs += 1  # pseudo global
        x, _ = runner(join(_test_dir, test + '.py'))
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

    def test_Fmath(self):
        self._run('testFmath')

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

    def test_Points(self):
        self._run('testPoints')

    def test_Routes(self):
        self._run('testRoutes')

    def test_Simplify(self):
        self._run('testSimplify')

    def test_Spherical(self):
        self._run('testSpherical')

    def test_Utils(self):
        self._run('testUtils')

    def test_Utm(self):
        self._run('testUtm')

    def test_Vectorial(self):
        self._run('testVectorial')

    def test_WebMercator(self):
        self._run('testWebMercator')

    def test_Ztotal(self):
        # final test to make sure all tests were run
        t = len(glob(join(_test_dir, 'test[A-Z]*.py')))
        self.assertEqual(TestSuite._runs, t)
#       t = sum(1 for t in dir(TestSuite) if t.startswith('test_'))
#       self.assertEqual(TestSuite._runs, t)


if __name__ == '__main__':

    unittest.main(argv=sys.argv)  # catchbreak=None, failfast=None, verbosity=2
