
# -*- coding: utf-8 -*-

u'''Print L{auxilats} version, etc. using C{python -m pygeodesy.auxilats}.
'''

__all__ = ()
__version__ = '24.09.04'


def _main(**ALorder):  # PYCHOK no cover

    try:
        from pygeodesy import auxilats
        from pygeodesy.internals import _fper, _name_version, \
                                         printf, _versions
        from pygeodesy.interns import _COMMASPACE_, _EQUAL_

        A  = auxilats.AuxLat(**ALorder)
        Cx = A._CXcoeffs  # PropertyRO: Adict of _Rdicts
        b, n, u, z = Cx.bnuz4()
        p = dict(ALorder=A.ALorder,CXb=b, CXb_z=_fper(b, z),
                            CXn=n, CXu=u, CXu_n=_fper(u, n))
        p = list(_EQUAL_(*t) for t in p.items())
        try:
            import geographiclib
            p.append(_name_version(geographiclib))
        except ImportError:
            pass

        a = _name_version(auxilats)
        printf('%s: %s (%s)', a, _COMMASPACE_(*p), _versions())

    except ImportError:
        from pygeodesy.internals import _usage
        print(_usage(__file__))


from sys import argv  # .internals._isPyChecker
_main(ALorder=int(argv[1])) if len(argv) == 2 and argv[1].isdigit() else _main()

# % python3.12 -m pygeodesy.auxilats 8
# pygeodesy.auxilats 24.09.04: ALorder=8, CXb=20310, CXb_z=71.5%, CXn=888, CXu=780, CXu_n=87.8%, geographiclib 2.0 (pygeodesy 24.9.9 Python 3.12.5 64bit arm64 macOS 14.6.1)

# % python3.12 -m pygeodesy.auxilats 6
# pygeodesy.auxilats 24.09.04: ALorder=6, CXb=11099, CXb_z=64.1%, CXn=522, CXu=448, CXu_n=85.8%, geographiclib 2.0 (pygeodesy 24.9.9 Python 3.12.5 64bit arm64 macOS 14.6.1)

# % python3.12 -m pygeodesy.auxilats 4
# pygeodesy.auxilats 24.09.04: ALorder=4, CXb=5367, CXb_z=58.8%, CXn=252, CXu=203, CXu_n=80.6%, geographiclib 2.0 (pygeodesy 24.9.9 Python 3.12.5 64bit arm64 macOS 14.6.1)


# **) MIT License
#
# Copyright (C) 2023-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
