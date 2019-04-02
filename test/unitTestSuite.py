
# -*- coding: utf-8 -*-

# Module to run all PyGeodesy tests as  python setup.py test

from base import test_dir
from run import run2

from glob import glob
from os.path import join
import unittest

__all__ = ('TestSuite',)
__version__ = '19.03.31'


class TestSuite(unittest.TestCase):
    '''Combine all test modules into a test suite/case
       and run each test module as a separate test.
    '''
    _runs = -1  # pseudo global, -1 for testGeoids

    def _run(self, test, *argv):
        TestSuite._runs += 1  # pseudo global
        x, _ = run2(join(test_dir, test + '.py'), *argv)
        self.assertEqual(x, 0)

    def test_Bases(self):
        self._run('testBases')

    def test_Classes(self):
        self._run('testClasses')

    def test_Clipy(self):
        self._run('testClipy')

    def test_Datum(self):
        self._run('testDatum')

    def test_Dms(self):
        self._run('testDms')

    def test_Elevations(self):
        self._run('testElevations')

    def test_Ellipsoidal(self):
        self._run('testEllipsoidal')

    def test_EllipsoidalGeodTest(self):
        self._run('testEllipsoidalGeodTest')

    def test_Fmath(self):
        self._run('testFmath')

    def test_Formy(self):
        self._run('testFormy')

    def test_Geohash(self):
        self._run('testGeohash')

    def test_Geoids(self):
        self._run('testGeoids', '-Karney')
        self._run('testGeoids', '-PGM', '-crop')  # -crop to cut time

    def test_GreatCircle(self):
        self._run('testGreatCircle')

    def test_Heights(self):
        self._run('testHeights')

    def test_LatLon(self):
        self._run('testLatLon')

    def test_Lazily(self):
        self._run('testLazily')

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

    def test_Utily(self):
        self._run('testUtily')

    def test_Utm(self):
        self._run('testUtm')

    def test_UtmTMcoords(self):
        self._run('testUtmTMcoords')

    def test_Vectorial(self):
        self._run('testVectorial')

    def test_WebMercator(self):
        self._run('testWebMercator')

    def test_Ztotal(self):
        # final test to make sure all tests were run
        t = len(glob(join(test_dir, 'test[A-Z]*.py')))
        self.assertEqual(TestSuite._runs, t)


if __name__ == '__main__':

    import sys

    unittest.main(argv=sys.argv)  # catchbreak=None, failfast=None, verbosity=2
