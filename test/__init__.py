
# -*- coding: utf-8 -*-

# Module to run all PyGeodesy tests as  python setup.py test

# Tested with 64-bit Python 2.7.13 and 3.6.1 on macOS 10.12.3,
# 10.12.4 and 10.12.5 Sierra and with Pythonista 3.1 using Python
# 2.7.12 and 3.5.1 on iOS 10.3.2.

from os.path import abspath, dirname
import sys

_test_dir = dirname(abspath(__file__))
# extend sys.path to include the .. directory
if _test_dir not in sys.path:  # Python 3+ ModuleNotFoundError
    sys.path.insert(0, _test_dir)

from unitTestSuite import TestSuite  # for setup.py

__all__ = ('TestSuite',)
__version__ = '17.06.19'
