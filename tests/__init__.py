
# -*- coding: utf-8 -*-

# Script to run some or all PyGeodesy tests with Python 2 or 3.

# Tested with 64-bit Python 2.6.9, 2.7.10, 2.7.13, 3.5.2, 3.5.3
# and 3.6.0 but only on MacOS 10.10 Mavericks, MacOS 10.11 El
# Capitan and MacOS 10.12.2 and 10.12.3 Sierra.

from os.path import dirname, join
from subprocess import PIPE, STDOUT, Popen
import sys
import unittest

__version__ = '17.03.03'

_python = sys.executable  # path
_python_O = _python
if not __debug__:
    _python_O += ' -OO'

_tests_dir = dirname(__file__)


class TestSuite(unittest.TestCase):
    '''Combine all test modules into a single test suite.

       XXX The test modules should be rewriiten to use
       the standard unittest or soem other test module.
    '''

    def _run(self, name):
        arg = join(_tests_dir, name + '.py')
        cmd = [_python_O, arg]
        p = Popen(cmd, creationflags=0,
                       executable   =_python,
                     # shell        =True,
                       stdin        =None,
                       stdout       =PIPE,  # XXX
                       stderr       =STDOUT)  # XXX

        r = p.communicate()[0]
        if isinstance(r, bytes):  # Python 3+
            r = r.decode('utf-8')

        # the exit status reflects the number
        # of test failures in the test module
        self.assertEqual(p.returncode, 0)

    def test_Bases(self):
        self._run('testBases')

    def test_Datum(self):
        self._run('testDatum')

    def test_Dms(self):
        self._run('testDms')

    def test_Ellipsoidal(self):
        self._run('testEllipsoidal')

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

    def test_Spherical(self):
        self._run('testSpherical')

    def test_Utm(self):
        self._run('testUtm')

    def test_Tests(self):
        self._run('tests')
