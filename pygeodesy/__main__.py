
# -*- coding: utf-8 -*-

u'''Print L{pygeodesy} version, etc. using C{python -m pygeodesy}.
'''

import os
import os.path as os_path
import sys

__all__ = ()
__version__ = '20.11.03'

try:  # MCCABE 16
    from pygeodesy import interns, isLazy, pygeodesy_abspath, version, \
                                  _isfrozen
    from pygeodesy.interns import NN, _COMMASPACE_, _DOT_, \
                                 _pygeodesy_abspath_, _Python_, \
                                 _SPACE_, _version_
    from pygeodesy.streprs import Fmt

    def _dot_attr(name, value):
        return Fmt.DOT(Fmt.EQUAL(name, value))

    p = [_dot_attr(*t) for t in ((_version_,           version),
                                 (_pygeodesy_abspath_, pygeodesy_abspath),
                                 ('isLazy',            isLazy),
                                 ('_isfrozen',        _isfrozen),
                                 ('_floats',       len(interns._floats)))]

    def _name_version(pkg):
        return _SPACE_(pkg.__name__, pkg.__version__)

    v = []

    if '[PyPy ' in sys.version:  # see test/base.py
        v.append(_SPACE_('PyPy', sys.version.split('[PyPy ')[1].split()[0]))
    v.append(_SPACE_(_Python_, sys.version.split(None, 1)[0]))
    try:
        import platform
        v.append(platform.architecture()[0])  # bits
    except ImportError:
        pass
    try:
        import geographiclib
        v.append(_name_version(geographiclib))
    except ImportError:
        pass
    try:
        import numpy
        v.append(_name_version(numpy))
    except ImportError:
        pass
    try:
        import scipy
        v.append(_name_version(scipy))
    except ImportError:
        pass
    x = os_path.basename(pygeodesy_abspath)
    print('%s%s (%s)' % (x, _COMMASPACE_.join(p), _COMMASPACE_.join(v)))

except ImportError:
    m = os_path.dirname(__file__).replace(os.getcwd(), _DOT_).strip()
    if len(m.split()) > 1:
        m = Fmt.QUOTE2(m)
    v = sys.version_info[0]
    if v < 3:
        v = NN
    print('usage: python%s -m %s' % (v, m))

# **) MIT License
#
# Copyright (C) 2016-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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

# % python -m pygeodesy
# pygeodesy.version='19.10.9', .pygeodesy_abspath='.../PyGeodesy/pygeodesy', .isLazy=None, ._isfrozen=False (Python 2.7.16, geographiclib 1.50, numpy 1.16.4, scipy 1.2.2)

# % python3 -m pygeodesy
# pygeodesy.version='19.10.9', .pygeodesy_abspath='.../PyGeodesy/pygeodesy', .isLazy=1, ._isfrozen=False (Python 3.7.4, geographiclib 1.50, numpy 1.17.2, scipy 1.3.1)
