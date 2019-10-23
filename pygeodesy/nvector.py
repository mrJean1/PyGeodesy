
# -*- coding: utf-8 -*-

u'''DEPRECATED, use module L{nvectorBase} instead.
'''

from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.nvectorBase import LatLonNvectorBase, \
                                  NorthPole, SouthPole, \
                                  NvectorBase as Nvector, sumOf  # PYCHOK exported
# all public constants, classes and functions
__all__ = _ALL_DOCS('LatLonNvectorBase') + (
          'NorthPole', 'SouthPole',  # constants
          'Nvector',  # classes
          'sumOf')  # functions
__version__ = '19.10.20'

# **) MIT License
#
# Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
