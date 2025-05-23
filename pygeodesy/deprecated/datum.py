
# -*- coding: utf-8 -*-

u'''DEPRECATED on 2022.09.12, use module L{pygeodesy.datums} or L{pygeodesy.ellipsoids} instead.
'''

# XXX only the items previously public
from pygeodesy.constants import R_FM, R_KM, R_M, R_MA, R_MB, R_NM, R_SM, R_VM
from pygeodesy.datums import Datum, Datums, Transform, Transforms
from pygeodesy.ellipsoids import Ellipsoid, Ellipsoids, Curvature2Tuple
from pygeodesy.lazily import _ALL_DEPRECATED, _ALL_OTHER

__all__ = _ALL_DEPRECATED.deprecated_datum
__version__ = '24.12.31'

assert _ALL_OTHER(Curvature2Tuple, Datum,  Ellipsoid,  Transform) + tuple(_.name for _ in
                                  (Datums, Ellipsoids, Transforms,
                                   R_FM, R_KM, R_M, R_MA, R_MB, R_NM, R_SM, R_VM)) == __all__

# **) MIT License
#
# Copyright (C) 2016-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
