
# -*- coding: utf-8 -*-

u'''Print L{auxilats} version, etc. using C{python -m pygeodesy.auxilats}.
'''

__all__ = ()
__version__ = '24.05.31'


def _CXstats():  # PYCHOK no cover
    '''(INTERNAL) Get the C{CP} stats.
    '''
    from pygeodesy.auxilats import Aux, AuxLat, auxLat
    from pygeodesy.datums import _WGS84

    A  = AuxLat(_WGS84)
    ax = A._coeffs(Aux.XI, Aux.CHI)
    Cx = auxLat._CXcoeffs(A.ALorder)
    pc = '%.1f%%' % (Cx.u * 100.0 / Cx.n)
    return dict(ALorder=A.ALorder, CXlen=Cx.n, CXset=Cx.u, CXset_len=pc, CXx=len(ax))


def _main():  # PYCHOK no cover

    import os.path as _os_path

    try:
        from pygeodesy import auxilats, printf, pygeodesy_abspath
        from pygeodesy.internals import _name_version, _Pythonarchine, _usage
        from pygeodesy.interns import _COMMASPACE_, _DOT_, _version_
        from pygeodesy.streprs import Fmt

        def _dot_attr(name, value):
            return Fmt.DOT(Fmt.EQUAL(name, value))

        s = tuple(sorted(_CXstats().items()))
        p = [_dot_attr(*t) for t in (((_version_, auxilats.__version__),) + s)]

        v = _Pythonarchine()
        try:
            import geographiclib
            v.append(_name_version(geographiclib))
        except ImportError:
            pass

        a =  auxilats.__name__
        x = _os_path.basename(pygeodesy_abspath)
        if not a.startswith(x):
            a = _DOT_(x, a)
        printf('%s%s (%s)', a, _COMMASPACE_.join(p), _COMMASPACE_.join(v))

    except ImportError:
        raise
        printf(_usage(__file__))


_main()

# **) MIT License
#
# Copyright (C) 2023-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
# pygeodesy.auxilats.version=24.05.31, .ALorder=6, .CXlen=522, .CXset=418, .CXset_len=80.1%, .CXx=6 (Python 3.12.3, 64bit, arm64, geographiclib 2.0)
