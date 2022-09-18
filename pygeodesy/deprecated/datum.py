
# -*- coding: utf-8 -*-

u'''DEPRECATED, use module L{datums} or L{ellipsoids} instead.
'''
# XXX only the items previously public
from pygeodesy.constants import R_M, R_MA, R_MB, R_KM, R_NM, R_SM, R_FM, R_VM
from pygeodesy.ellipsoids import Ellipsoid, Ellipsoids, Curvature2Tuple  # PYCHOK exported
from pygeodesy.datums import Datum, Datums, Transform, Transforms
# from pygeodesy.props import _deprecated_module

__all__ = (R_M.name,  R_MA.name, R_MB.name,
           R_KM.name, R_NM.name, R_SM.name, R_FM.name, R_VM.name,
           Datum.__name__, Ellipsoid.__name__, Transform.__name__, Curvature2Tuple.__name__,
           Datums.name,    Ellipsoids.name,    Transforms.name)
__version__ = '22.09.12'

# _deprecated_module(__name__)

# **) MIT License
#
# Copyright (C) 2016-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
