
# -*- coding: utf-8 -*-

u'''A pure Python version of I{Karney}'s C++ classes U{GeodesicExact
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GeodesicExact.html>}
and U{GeodesicLineExact
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GeodesicLineExact.html>}.

For more details, see the C++ U{GeographicLib<https://GeographicLib.SourceForge.io/html/index.html>}
documentation, especially the U{Class List<https://GeographicLib.SourceForge.io/html/annotated.html>}
and the background information on the page U{Geodesics on an ellipsoid of revolution
<https://GeographicLib.SourceForge.io/html/geodesic.html#geodseries>}.

Also, compare C{GeodesicExact} and C{GeodesicLineExact} to I{standard} classes C{Geodesic}
respectively C{GeodesicLine} from I{Karney}'s Python implementation U{geographiclib
<https://GeographicLib.SourceForge.io/html/other.html#python>}.
'''
from pygeodesy.geodesicx.gx import GeodesicExact, GeodesicLineExact  # PYCHOK exported
from pygeodesy.geodesicx.gxarea import GeodesicAreaExact, PolygonArea  # PYCHOK exported
from pygeodesy.geodesicx.gxbases import Caps, GeodesicError  # PYCHOK exported
# from pygeodesy.karney import GeodesicError  # from .geodesicx.bases
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY

__all__ = _ALL_LAZY.geodesicx + _ALL_DOCS(GeodesicError, GeodesicAreaExact,
                                          GeodesicExact, GeodesicLineExact,
                                          PolygonArea)
__version__ = '21.11.23'

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
