
# -*- coding: utf-8 -*-

u'''DEPRECATED classes, constants, functions, interns, methods, etc. and all
lazily imported.

Kept for backward compatibility, including DEPRECATED modules C{pygeodesy.bases},
C{pygeodesy.datum} and C{pygeodesy.nvector}.  Use either C{from pygeodesy import
bases} or C{from pygeodesy.deprecated import bases}.  Likewise for C{datum} and
C{nvector}.
'''

from pygeodesy.deprecated.bases import *  # PYCHOK not pygeodesy.__init__
from pygeodesy.deprecated.datum import *  # PYCHOK not pygeodesy.__init__
from pygeodesy.deprecated.nvector import *  # PYCHOK not pygeodesy.__init__

from pygeodesy.deprecated.classes import *  # PYCHOK expected
from pygeodesy.deprecated.consterns import *  # PYCHOK expected
from pygeodesy.deprecated.functions import *  # PYCHOK expected

from pygeodesy.lazily import _ALL_ATTRS, _ALL_DEPRECATED, _lazy_import_as, _unLazy0  # _lazy_import_star

__all__ = (_ALL_DEPRECATED.deprecated_bases +
           _ALL_DEPRECATED.deprecated_datum +
           _ALL_DEPRECATED.deprecated_nvector +

           _ALL_DEPRECATED.deprecated_classes +
           _ALL_DEPRECATED.deprecated_consterns +
           _ALL_DEPRECATED.deprecated_functions)
__version__ = '24.02.02'

if _unLazy0:
    from pygeodesy.deprecated import bases, datum, nvector, rhumbBase, \
                                     rhumbaux, rhumbsolve, rhumbx  # PYCHOK expected
    __all__ += _ALL_ATTRS(_ALL_DEPRECATED.deprecated)  # DEPRECATED modules

else:  # lazily import modules and exported attrs
    __getattr__ = _lazy_import_as(__name__)
#   __star__    = _lazy_import_star(__name__)

# **) MIT License
#
# Copyright (C) 2018-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
