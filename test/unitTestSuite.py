
# -*- coding: utf-8 -*-

# Module to run all PyGeodesy tests as  python setup.py test

from base import test_dir
from run import run2

from glob import glob
from os.path import join
import unittest

__all__ = ('TestSuite',)
__version__ = '20.07.19'


class TestSuite(unittest.TestCase):
    '''Combine all test modules into a test suite/case
       and run each test module as a separate test.
    '''
    _runs = 0  # pseudo global, 0 for testGeoids

    def _run(self, test, *argv):
        TestSuite._runs += 1  # pseudo global
        x, _ = run2(join(test_dir, test + '.py'), *argv)
        self.assertEqual(x, 0)

    def test_Azimuthal(self):
        self._run('testAzimuthal')

    def test_Bases(self):
        self._run('testBases')

    def test_Basics(self):
        self._run('testBasics')

    def test_Cartesian(self):
        self._run('testCartesian')

    def test_Classes(self):
        self._run('testClasses')

    def test_Clipy(self):
        self._run('testClipy')

    def test_Css(self):
        self._run('testCss')

    def test_Datum(self):
        self._run('testDatum')

    def test_Deprecated(self):
        self._run('testDeprecated')

    def test_Dms(self):
        self._run('testDms')

    def test_Ecef(self):
        self._run('testEcef')

    def test_Elevations(self):
        self._run('testElevations')

    def test_Ellipsoidal(self):
        self._run('testEllipsoidal')

    def test_EllipsoidalGeodTest(self):
        self._run('testEllipsoidalGeodTest')

    def test_Elliptic(self):
        self._run('testElliptic')

    def test_Espg(self):
        self._run('testEpsg')

    def test_Errors(self):
        self._run('testErrors')

    def test_Etm(self):
        self._run('testEtm')

    def test_EtmTMcoords(self):
        self._run('testEtmTMcoords')

    def test_ExactTMcoords(self):
        self._run('testExactTMcoords')

    def test_Fmath(self):
        self._run('testFmath')

    def test_Formy(self):
        self._run('testFormy')

    def test_Frechet(self):
        self._run('testFrechet')

    def test_Gars(self):
        self._run('testGars')

    def test_Geohash(self):
        self._run('testGeohash')

    def test_Geoids(self):
        self._run('testGeoids')  # both Karney and PGM
        # self._run('testGeoids', '-Karney')
        # self._run('testGeoids', '-PGM', '-crop')  # -crop to cut time

    def test_GreatCircle(self):
        self._run('testGreatCircle')

    def test_Hausdorff(self):
        self._run('testHausdorff')

    def test_Heights(self):
        self._run('testHeights')

    def test_LatLon(self):
        self._run('testLatLon')

    def test_Karney(self):
        self._run('testKarney')

    def test_Lazily(self):
        self._run('testLazily')

    def test_Lcc(self):
        self._run('testLcc')

    def test_Mgrs(self):
        self._run('testMgrs')

    def test_Modules(self):
        self._run('testModules')

    def test_Named(self):
        self._run('testNamed')

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

    def test_Streprs(self):
        self._run('testStreprs')

    def test_TMcoords(self):
        self._run('testTMcoords')

    def test_Trf(self):
        self._run('testTrf')

    def test_Units(self):
        self._run('testUnits')

    def test_Ups(self):
        self._run('testUps')

    def test_Utily(self):
        self._run('testUtily')

    def test_Utm(self):
        self._run('testUtm')

    def test_UtmTMcoords(self):
        self._run('testUtmTMcoords')

    def test_UtmUps(self):
        self._run('testUtmUps')

    def test_UtmUpsTMcoords(self):
        self._run('testUtmUpsTMcoords')

    def test_Vectorial(self):
        self._run('testVectorial')

    def test_WebMercator(self):
        self._run('testWebMercator')

    def test_Wgrs(self):
        self._run('testWgrs')

    def test_Ztotal(self):
        # final test to make sure all tests were run
        t = len(glob(join(test_dir, 'test[A-Z]*.py')))
        self.assertEqual(TestSuite._runs, t)


if __name__ == '__main__':

    import sys

    unittest.main(argv=sys.argv)  # catchbreak=None, failfast=None, verbosity=2
