
# -*- coding: utf-8 -*-

u'''Print L{pygeodesy} version, etc. using C{python -m pygeodesy}.
'''

__all__ = ()
__version__ = '21.11.28'


def _main():  # PYCHOK no cover

    import os.path as os_path

    try:
        from pygeodesy import _isfrozen, isLazy, pygeodesy_abspath, version
        from pygeodesy.interns import _COMMASPACE_, _floats, _pygeodesy_abspath_, \
                                      _pythonarchine, _SPACE_, _usage, _version_
        from pygeodesy.lazily import _a_l_l_, _all_imports, printf
        from pygeodesy.streprs import Fmt

        def _dot_attr(name, value):
            return Fmt.DOT(Fmt.EQUAL(name, value))

        p = [_dot_attr(*t) for t in ((_version_,            version),
                                     (_pygeodesy_abspath_,  pygeodesy_abspath),
                                     ('isLazy',             isLazy),
                                     ('_isfrozen',       _isfrozen),
                                     ('_floats',      len(_floats)),
                                     (_a_l_l_, len(_all_imports())))]

        def _name_version(pkg):
            return _SPACE_(pkg.__name__, pkg.__version__)

        v = _pythonarchine()
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
        printf('%s%s (%s)', x, _COMMASPACE_.join(p), _COMMASPACE_.join(v))

    except ImportError:
        printf(_usage(__file__))


_main()

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

# % python2 -m pygeodesy
# pygeodesy.version=21.6.30, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=None, ._isfrozen=False, ._floats=408 (Python 2.7.18, 64bit, geographiclib 1.50, numpy 1.16.6, scipy 1.2.2)

# % python3.8 -m pygeodesy
# pygeodesy.version=21.6.30, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=1, ._isfrozen=False, ._floats=53 (Python 3.8.6, 64bit, geographiclib 1.52, numpy 1.19.2, scipy 1.5.2)

# % python3.9 -m pygeodesy
# pygeodesy.version=21.6.30, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=1, ._isfrozen=False, ._floats=53 (Python 3.9.5, 64bit, geographiclib 1.52)
