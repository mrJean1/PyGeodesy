
# -*- coding: utf-8 -*-

u'''DEPRECATED constants and interns for export and backward compatibility.
'''

from pygeodesy.constants import EPS_2, MANT_DIG, _1_0
from pygeodesy.units import Float, Int, Str

__all__ = ()
__version__ = '23.09.12'


class _Deprecated_Float(Float):
    '''DEPRECATED, don't use.'''
    pass


class _Deprecated_Int(Int):
    '''DEPRECATED, don't use.'''
    pass


class _Deprecated_Str(Str):
    '''DEPRECATED, don't use.'''
    pass


EPS1_2 = _Deprecated_Float(EPS1_2=_1_0 - EPS_2)
MANTIS = _Deprecated_Int(MANTIS=MANT_DIG)
OK     = _Deprecated_Str(OK='OK')

__all__ += EPS1_2.name, MANTIS.name, OK.name

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
