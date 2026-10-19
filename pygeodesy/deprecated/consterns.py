
# -*- coding: utf-8 -*-

u'''DEPRECATED constants, interns and singletons kept for backward compatibility.
'''

from pygeodesy.constants import EPS_2, MANT_DIG, _1_0
from pygeodesy.ellipses import Ellipse
from pygeodesy.lazily import _ALL_DEPRECATED, _FOR_DOCS
from pygeodesy.props import deprecated_method
from pygeodesy.units import Float, Int, Str

__all__ = _ALL_DEPRECATED.deprecated_consterns
__version__ = '26.02.12'


class _Deprecated_Float(Float):
    '''DEPRECATED on 2023.09.12, I{don't use}.'''
    pass  # PYLINT W0107


class _Deprecated_Int(Int):
    '''DEPRECATED on 2023.09.12, I{don't use}.'''
    pass  # PYLINT W0107


class _Deprecated_Str(Str):
    '''DEPRECATED on 2023.09.12, I{don't use}.'''
    pass  # PYLINT W0107


class Elliperim(object):
    '''DEPRECATED on 2026.02.06, use class L{Ellipse}.'''

    @deprecated_method
    def AGM(self, a, b, **unused):  # PYCHOK no cover
        '''DEPRECATED on 2026.02.12, use property L{Ellipse}{C(a, b).perimeterAGM}.'''
        return Ellipse(a, b).perimeterAGM

    @deprecated_method
    def Arc43(self, a, b):  # PYCHOK no cover
        '''DEPRECATED on 2026.02.12, use property L{Ellipse}{C(a, b).perimeter4Arc3}.'''
        return Ellipse(a, b).perimeter4Arc3

    @deprecated_method
    def E2k(self, a, b):  # PYCHOK no cover
        '''DEPRECATED on 2026.02.12, use property L{Ellipse}{C(a, b).perimeter2k}.'''
        return Ellipse(a, b).perimeter2k

    @deprecated_method
    def e2k(self, a, b, **unused):  # PYCHOK no cover
        '''DEPRECATED on 2026.02.12, use property L{Ellipse}{C(a, b).perimeter2k_}.'''
        return Ellipse(a, b).perimeter2k_

    @deprecated_method
    def GK(self, a, b):  # PYCHOK no cover
        '''DEPRECATED on 2026.02.12, use property L{Ellipse}{C(a, b).perimeterGK}.'''
        return Ellipse(a, b).perimeterGK

    @deprecated_method
    def HG(self, a, b, **unused):  # PYCHOK no cover
        '''DEPRECATED on 2026.02.12, use property L{Ellipse}{C(a, b).perimeterHGK}.'''
        return Ellipse(a, b).perimeterHGK

    @deprecated_method
    def R2(self, a, b):  # PYCHOK no cover
        '''DEPRECATED on 2026.02.12, use property L{Ellipse}{C(a, b).perimeter2R}.'''
        return Ellipse(a, b).perimeter2R

if not _FOR_DOCS:  # PYCHOK force epydoc
    Elliperim = Elliperim()  # singleton
del _FOR_DOCS

EPS1_2 = _Deprecated_Float(EPS1_2=_1_0 - EPS_2)
MANTIS = _Deprecated_Int(MANTIS=MANT_DIG)
OK     = _Deprecated_Str(OK='OK')

# **) MIT License
#
# Copyright (C) 2018-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
