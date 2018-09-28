
# -*- coding: utf-8 -*-

# Module to run all PyGeodesy tests as  python setup.py test

from os.path import abspath, dirname
import sys

_test_dir = dirname(abspath(__file__))
# setting __path__ is sufficient for importing
# modules internal to this test package
__path__ = [_test_dir, dirname(_test_dir)]
# extend sys.path to include the .. directory,
# required for module ..setup.py to work
if _test_dir not in sys.path:
    sys.path.insert(0, _test_dir)

from unitTestSuite import TestSuite  # PYCHOK for setup.py

__all__ = ('TestSuite',)
__version__ = '18.09.25'

del abspath, dirname, sys, _test_dir
