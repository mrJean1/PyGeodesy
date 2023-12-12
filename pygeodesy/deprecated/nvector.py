
# -*- coding: utf-8 -*-

u'''DEPRECATED on 2021.05.20, use (INTERNAL) module L{pygeodesy.nvectorBase} instead.
'''

# from pygeodesy.interns import NN  # from .lazily
from pygeodesy.lazily import _ALL_DEPRECATED, _ALL_OTHER,  NN
from pygeodesy.nvectorBase import LatLonNvectorBase, NorthPole, NvectorBase, \
                                                     SouthPole, sumOf
from pygeodesy.props import deprecated_class

__all__ = _ALL_DEPRECATED.deprecated_nvector
__version__ = '23.11.26'


class Nvector(NvectorBase):  # PYCHOK no cover
    '''DEPRECATED on 2021.05.20, see (INTERNAL) class L{pygeodesy.nvectorBase.NvectorBase}.
    '''
    def __init__(self, x, y=None, z=None, h=0, ll=None, datum=None, name=NN):
        deprecated_class(self.__class__)
        NvectorBase.__init__(self, x, y=y, z=z, h=h, ll=ll, datum=datum, name=name)


assert (_ALL_OTHER(LatLonNvectorBase, Nvector, sumOf) +
                  (NorthPole.name, SouthPole.name)) == __all__

# **) MIT License
#
# Copyright (C) 2016-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
