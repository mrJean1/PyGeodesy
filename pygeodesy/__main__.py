
# -*- coding: utf-8 -*-

u'''Print L{pygeodesy} version, etc. using C{python -m pygeodesy}.
'''

__all__ = ()
__version__ = '24.10.14'

from os.path import basename, dirname


def _main():  # PYCHOK no cover

    try:
        from pygeodesy import constants, _isfrozen, pygeodesy_abspath, version
        from pygeodesy.basics import _xcoverage,_xgeographiclib, _xnumpy, _xscipy
        from pygeodesy.internals import _name_version, printf, _usage, _versions
        from pygeodesy.interns import NN, _COMMASPACE_, _DEPRECATED_, _DOT_, _DUNDER_all_, \
                                     _EQUAL_, _pygeodesy_abspath_, _version_
        from pygeodesy.lazily import _all_deprecates, _all_imports, isLazy

        def _p(name_value):
            return _DOT_(NN, _EQUAL_(*name_value))

        p = list(_p(t) for t in ((_version_,                     version),
                                 (_pygeodesy_abspath_, pygeodesy_abspath),
                                 ('isLazy',                       isLazy),
                                 ('_isfrozen',                 _isfrozen),
                                 ('_floats',      len(constants._floats)),
                                 (_DUNDER_all_,      len(_all_imports())),
                                 (_DEPRECATED_,   len(_all_deprecates()))))

        def _nv(_xpkg, p):
            try:
                pkg = _xpkg(_main)
            except ImportError:
                pkg = None
            if pkg is not None:
                p.append(_name_version(pkg))

        _nv(_xcoverage, p)
        _nv(_xgeographiclib, p)
        _nv(_xnumpy, p)
        _nv(_xscipy, p)

        x = basename(pygeodesy_abspath)
        printf('%s%s (%s)', x, _COMMASPACE_.join(p), _versions())

    except ImportError:
        from pygeodesy.internals import printf, _usage
        printf(_usage(__file__))


try:
    _main()
except ImportError:
    from sys import executable as x
    t = 'usage:', basename(x), '-m', basename(dirname(__file__))
    print(' '.join(t))

# % python3.13 -m pygeodesy
# pygeodesy.version=24.9.6, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=1, ._isfrozen=False, ._floats=81, .__all__=936, .DEPRECATED=100, coverage 7.6.0 (pygeodesy 24.9.6 Python 3.13.0rc1 64bit arm64 macOS 14.6.1)

# % python3.12 -m pygeodesy
# pygeodesy.version=24.9.6, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=1, ._isfrozen=False, ._floats=81, .__all__=936, .DEPRECATED=100, coverage 7.6.0, geographiclib 2.0, numpy 2.1.0, scipy 1.14.1 (pygeodesy 24.9.6 Python 3.12.5 64bit arm64 macOS 14.6.1)

# % python3.11 -m pygeodesy
# pygeodesy.version=24.9.6, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=1, ._isfrozen=False, ._floats=81, .__all__=936, .DEPRECATED=100, coverage 7.6.0, geographiclib 2.0, numpy 1.24.2, scipy 1.10.1 (pygeodesy 24.9.6 Python 3.11.5 64bit arm64 macOS 14.6.1)

# % python3.10 -m pygeodesy
# pygeodesy.version=24.9.6, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=1, ._isfrozen=False, ._floats=81, .__all__=936, .DEPRECATED=100, coverage 7.6.0, geographiclib 2.0, numpy 1.23.3, scipy 1.9.1 (pygeodesy 24.9.6 Python 3.10.8 64bit arm64 macOS 14.6.1)

# % python3.9 -m pygeodesy
# pygeodesy.version=24.9.6, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=1, ._isfrozen=False, ._floats=81, .__all__=936, .DEPRECATED=100, coverage 7.2.2 (pygeodesy 24.9.6 Python 3.9.6 64bit arm64 macOS 14.6.1)

# % python3.8 -m pygeodesy
# pygeodesy.version=24.9.6, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=1, ._isfrozen=False, ._floats=81, .__all__=936, .DEPRECATED=100, coverage 7.2.2, geographiclib 1.52, numpy 1.19.2, scipy 1.5.2 (pygeodesy 24.9.6 Python 3.8.10 64bit arm64_x86_64 macOS 10.16)

# % python3.7 -m pygeodesy
# pygeodesy.version=24.9.6, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=1, ._isfrozen=False, ._floats=81, .__all__=936, .DEPRECATED=100, coverage 4.5.4, geographiclib 1.50, numpy 1.17.2, scipy 1.3.1 (pygeodesy 24.9.6 Python 3.7.6 64bit arm64_x86_64 macOS 10.16)

# % python2.7 -m pygeodesy
# pygeodesy.version=24.9.6, .pygeodesy_abspath=.../PyGeodesy/pygeodesy, .isLazy=None, ._isfrozen=False, ._floats=792, .__all__=936, .DEPRECATED=100, coverage 5.5, geographiclib 1.50, numpy 1.16.6, scipy 1.2.2 (pygeodesy 24.9.6 Python 2.7.18 64bit arm64_x86_64 macOS 10.16)

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
