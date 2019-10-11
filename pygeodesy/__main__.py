
# -*- coding: utf-8 -*-

u'''Print L{pygeodesy} version, etc. using C{python -m pygeodesy}.
'''

import os
import os.path as os_path
import sys

__all__ = ()
__version__ = '19.10.09'

try:
    from pygeodesy import isLazy, pygeodesy_abspath, version, \
                         _isfrozen

    p = ['.%s=%r' % t for t in (('version',           version),
                                ('pygeodesy_abspath', pygeodesy_abspath),
                                ('isLazy',            isLazy),
                                ('_isfrozen',        _isfrozen))]
    v = []
    if '[PyPy ' in sys.version:
        v.append('PyPy ' + sys.version.split('[PyPy ')[1].split()[0])
    v.append('Python ' + sys.version.split(None, 1)[0])
    try:
        import geographiclib
        v.append('geographiclib ' + geographiclib.__version__)
    except ImportError:
        pass
    try:
        import numpy
        v.append('numpy ' + numpy.__version__)
    except ImportError:
        pass
    try:
        import scipy
        v.append('scipy ' + scipy.__version__)
    except ImportError:
        pass
    x = os_path.basename(pygeodesy_abspath)
    print('%s%s (%s)' % (x, ', '.join(p), ', '.join(v)))

except ImportError:
    m = os_path.dirname(__file__).replace(os.getcwd(), '.').strip()
    if len(m.split()) > 1:
        m = '"%s"' % (m,)
    v = sys.version_info[0]
    if v < 3:
        v = ''
    print('usage: python%s -m %s' % (v, m))

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

# % python -m pygeodesy
# pygeodesy.version='19.10.9', .pygeodesy_abspath='.../PyGeodesy/pygeodesy', .isLazy=None, ._isfrozen=False (Python 2.7.16, geographiclib 1.50, numpy 1.16.4, scipy 1.2.2)

# % python3 -m pygeodesy
# pygeodesy.version='19.10.9', .pygeodesy_abspath='.../PyGeodesy/pygeodesy', .isLazy=1, ._isfrozen=False (Python 3.7.4, geographiclib 1.50, numpy 1.17.2, scipy 1.3.1)
