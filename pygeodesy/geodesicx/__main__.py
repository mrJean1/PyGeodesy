
# -*- coding: utf-8 -*-

u'''Print L{geodesicx} version, etc. using C{python -m pygeodesy.geodesicx}.
'''

__all__ = ()
__version__ = '24.09.06'


def _main(**C4order):  # PYCHOK no cover

    try:
        from pygeodesy import GeodesicExact, geodesicx
        from pygeodesy.internals import _fper, _name_version, \
                                         printf, _sizeof, _versions
        from pygeodesy.interns import _COMMASPACE_, _EQUAL_
        try:
            import numpy
        except ImportError:
            numpy = None

        gX = GeodesicExact(**C4order)
        cs = geodesicx.gx._C4coeffs(gX.C4order)
        n  = len(cs)
        u  = n         if numpy else len(set(cs))
        z  = cs.nbytes if numpy else _sizeof(cs)
        p  = dict(C4order=gX.C4order, C4n=n, C4u=u,
                  C4u_n=_fper(u, n), C4x=len(gX._C4x),
                  C4t=type(cs).__name__, C4z=z)
        p  = list(_EQUAL_(*t) for t in p.items())
        if numpy:
            p.append(_name_version(numpy))
        try:
            import geographiclib
            p.append(_name_version(geographiclib))
        except ImportError:
            pass

        g = _name_version(geodesicx)
        printf('%s: %s (%s)', g, _COMMASPACE_.join(p), _versions())

    except ImportError:
        from pygeodesy.internals import _usage
        print(_usage(__file__))


from sys import argv  # .internals._isPyChecker
_main(C4order=int(argv[1])) if len(argv) == 2 and argv[1].isdigit() else _main()

# % python3.13 -m pygeodesy.geodesicx 30
# pygeodesy.geodesicx 24.09.06: C4order=30, C4n=5425, C4u=5107, C4u_n=94.1%, C4x=465, C4t=tuple, C4z=166008 (pygeodesy 24.9.6 Python 3.13.0rc1 64bit arm64 macOS 14.6.1)
# % python3.12 -m pygeodesy.geodesicx 30
# pygeodesy.geodesicx 24.09.06: C4order=30, C4n=5425, C4u=5425, C4u_n=100.0%, C4x=465, C4t=ndarray, C4z=43400, numpy 2.1.0, geographiclib 2.0 (pygeodesy 24.9.6 Python 3.12.5 64bit arm64 macOS 14.6.1)

# % python3.13 -m pygeodesy.geodesicx 27
# pygeodesy.geodesicx 24.09.06: C4order=27, C4n=4032, C4u=3764, C4u_n=93.4%, C4x=378, C4t=tuple, C4z=122632 (pygeodesy 24.9.6 Python 3.13.0rc1 64bit arm64 macOS 14.6.1)
# % python3.12 -m pygeodesy.geodesicx 27
# pygeodesy.geodesicx 24.09.06: C4order=27, C4n=4032, C4u=4032, C4u_n=100.0%, C4x=378, C4t=ndarray, C4z=32256, numpy 2.1.0, geographiclib 2.0 (pygeodesy 24.9.6 Python 3.12.5 64bit arm64 macOS 14.6.1)

# % python3.13 -m pygeodesy.geodesicx 24
# pygeodesy.geodesicx 24.09.06: C4order=24, C4n=2900, C4u=2708, C4u_n=93.4%, C4x=300, C4t=tuple, C4z=88232 (pygeodesy 24.9.6 Python 3.13.0rc1 64bit arm64 macOS 14.6.1)
# % python3.12 -m pygeodesy.geodesicx 24
# pygeodesy.geodesicx 24.09.06: C4order=24, C4n=2900, C4u=2900, C4u_n=100.0%, C4x=300, C4t=ndarray, C4z=23200, numpy 2.1.0, geographiclib 2.0 (pygeodesy 24.9.6 Python 3.12.5 64bit arm64 macOS 14.6.1)

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

# % python3 -m pygeodesy.geodesicx
# pygeodesy.geodesicx.version=24.05.31, .C4len=5425, .C4order=30, .C4set=5107, .C4set_len=94.1%, .C4x=465 (Python 3.12.3, 64bit, arm64, geographiclib 2.0)
