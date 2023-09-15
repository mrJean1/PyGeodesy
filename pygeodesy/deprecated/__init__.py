
# -*- coding: utf-8 -*-

u'''DEPRECATED classes, constants, functions, interns, methods, etc.

Kept and exported for backward compatibility, including deprecated modules
C{pygeodesy.bases}, C{pygeodesy.datum} and C{pygeodesy,nvector}, previously
inside the C{pygeodesy} package.

Use either C{from pygeodesy import bases} or C{from pygeodesy.deprecated import
bases}.  Likewise for C{datum} and C{nvector}.
'''

from pygeodesy.deprecated.classes import *  # PYCHOK expected
from pygeodesy.deprecated.consterns import *  # PYCHOK expected
from pygeodesy.deprecated.functions import *  # PYCHOK expected
from pygeodesy.lazily import _ALL_LAZY, isLazy as _isLazy

if _isLazy:  # XXX force import of all deprecated modules
    import pygeodesy.deprecated.bases as bases, \
           pygeodesy.deprecated.datum as datum, \
           pygeodesy.deprecated.nvector as nvector  # PYCHOK unused
    # XXX instead, use module_property or enhance .lazily

__all__ = _ALL_LAZY.deprecated
__version__ = '23.09.12'

# **) MIT License
#
# Copyright (C) 2018-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
