
# -*- coding: utf-8 -*-

u'''Print L{pygeodesy} version, etc. using C{python -m pygeodesy}.
'''

__all__ = ()
__version__ = '23.03.29'

from os.path import basename, dirname


def _main():  # PYCHOK no cover

    try:
        from pygeodesy import _isfrozen, pygeodesy_abspath, version
        from pygeodesy.basics import _xgeographiclib, _xnumpy, _xscipy
        from pygeodesy.constants import _floats
        from pygeodesy.interns import _COMMASPACE_, _pygeodesy_abspath_, \
                                      _pythonarchine, _SPACE_, _usage, _version_
        from pygeodesy.lazily import _a_l_l_, _all_imports, isLazy, printf
        from pygeodesy.streprs import Fmt

        def _dot_attr(name, value):
            return Fmt.DOT(Fmt.EQUAL(name, value))

        p = [_dot_attr(*t) for t in ((_version_,            version),
                                     (_pygeodesy_abspath_,  pygeodesy_abspath),
                                     ('isLazy',             isLazy),
                                     ('_isfrozen',       _isfrozen),
                                     ('_floats',      len(_floats)),
                                     (_a_l_l_, len(_all_imports())))]

        def _nv(_xpkg, v):
            try:
                pkg = _xpkg(_main)
            except ImportError:
                pkg = None
            if pkg is not None:
                v.append(_SPACE_(pkg.__name__, pkg.__version__))

        v = _pythonarchine()
        _nv(_xgeographiclib, v)
        _nv(_xnumpy, v)
        _nv(_xscipy, v)

        x = basename(pygeodesy_abspath)
        printf('%s%s (%s)', x, _COMMASPACE_.join(p), _COMMASPACE_.join(v))

    except ImportError:
        from pygeodesy.interns import _usage
        from pygeodesy.lazily import printf
        printf(_usage(__file__))


try:
    _main()
except ImportError:
    from sys import executable as x
    print('%s: %s %s %s' % ('usage', basename(x), '-m', basename(dirname(__file__))))

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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

# python3 -m pygeodesy
# pygeodesy.version=23.1.6, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=1, ._isfrozen=False, ._floats=81, .__all__=908 (Python 3.11.0, 64bit, arm64, geographiclib 2.0)

# % python3.9 -m pygeodesy
# pygeodesy.version=23.1.6, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=1, ._isfrozen=False, ._floats=81, .__all__=908 (Python 3.9.6, 64bit, arm64)

# % python3.8 -m pygeodesy
# pygeodesy.version=23.1.6, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=1, ._isfrozen=False, ._floats=81, .__all__=908 (Python 3.8.10, 64bit, arm64_x86_64, geographiclib 1.52, numpy 1.19.2, scipy 1.5.2)

# % python2 -m pygeodesy
# pygeodesy.version=23.1.6, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=None, ._isfrozen=False, ._floats=560, .__all__=908 (Python 2.7.18, 64bit, arm64_x86_64, geographiclib 1.50, numpy 1.16.6, scipy 1.2.2)
