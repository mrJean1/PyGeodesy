
# -*- coding: utf-8 -*-

# Module to run all PyGeodesy tests as  python setup.py test

from os.path import abspath, dirname
import sys

_test_dir = dirname(abspath(__file__))
# extend sys.path to include the .. directory
if _test_dir not in sys.path:  # Python 3+ ModuleNotFoundError
    sys.path.insert(0, _test_dir)

from unitTestSuite import TestSuite  # for setup.py

__all__ = ('TestSuite',)
__version__ = '17.08.10'
