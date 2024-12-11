
# -*- coding: utf-8 -*-

u'''Coefficients for C{_AUXLATITUDE_ORDER} 4 from I{Karney}'s C++ class U{AuxLatitude
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AuxLatitude.html>}
transcoded to a Python C{_Rdict[auxout][auxin]} of C{_Rtuple}s.

Copyright (C) Charles Karney (2022-2024) Karney@Alum.MIT.edu> and licensed under the
MIT/X11 License.  For more information, see <https:#GeographicLib.SourceForge.io>.
'''
from pygeodesy.auxilats.auxily import Aux
from pygeodesy.auxilats._CX_Rs import _Rcoeffs, _Rdict, _Rtuple

__all__ = ()
__version__ = '24.09.03'

_coeffs_4 = _Rcoeffs(4, {  # GEOGRAPHICLIB_AUXLATITUDE_ORDER == 4
    Aux.PHI: _Rdict(38,
       # C[phi,phi] skipped
       _Rtuple(Aux.BETA, 6,  # C[phi,beta]; even coeffs only
                  '0, 1, 0, 1/2, 1/3, 1/4'),
       _Rtuple(Aux.THETA, 6,  # C[phi,theta]; even coeffs only
                  '-2, 2, -4, 2, 8/3, 4'),
       _Rtuple(Aux.MU, 6,  # C[phi,mu]; even coeffs only
                  '-27/32, 3/2, -55/32, 21/16, 151/96, 1097/512'),
       _Rtuple(Aux.CHI, 10,  # C[phi,chi]
                  '116/45, -2, -2/3, 2, -227/45, -8/5',
                  '7/3, -136/35, 56/15, 4279/630'),
       _Rtuple(Aux.XI, 10,  # C[phi,xi]
                  '-2582/14175, -16/35, 4/45, 4/3, -11966/14175',
                  '152/945, 46/45, 3802/14175, 3044/2835, 6059/4725')
    ),
    Aux.BETA: _Rdict(38,
       _Rtuple(Aux.PHI, 6,  # C[beta,phi]; even coeffs only
                  '0, -1, 0, 1/2, -1/3, 1/4'),
       # C[beta,beta] skipped
       _Rtuple(Aux.THETA, 6,  # C[beta,theta]; even coeffs only
                  '0, 1, 0, 1/2, 1/3, 1/4'),
       _Rtuple(Aux.MU, 6,  # C[beta,mu]; even coeffs only
                  '-9/32, 1/2, -37/96, 5/16, 29/96, 539/1536'),
       _Rtuple(Aux.CHI, 10,  # C[beta,chi]
                  '38/45, -1/3, -2/3, 1, -7/9, -14/15',
                  '5/6, -34/21, 16/15, 2069/1260'),
       _Rtuple(Aux.XI, 10,  # C[beta,xi]
                  '-1082/14175, -46/315, 4/45, 1/3, -338/2025',
                  '68/945, 17/90, 1102/14175, 461/2835, 3161/18900')
    ),
    Aux.THETA: _Rdict(38,
       _Rtuple(Aux.PHI, 6,  # C[theta,phi]; even coeffs only
                  '2, -2, -4, 2, -8/3, 4'),
       _Rtuple(Aux.BETA, 6,  # C[theta,beta]; even coeffs only
                  '0, -1, 0, 1/2, -1/3, 1/4'),
       # C[theta,theta] skipped
       _Rtuple(Aux.MU, 6,  # C[theta,mu]; even coeffs only
                  '-23/32, -1/2, -5/96, 5/16, 1/32, 283/1536'),
       _Rtuple(Aux.CHI, 10,  # C[theta,chi]
                  '4/9, -2/3, -2/3, 0, -23/45',
                  '-4/15, 1/3, -24/35, 2/5, 83/126'),
       _Rtuple(Aux.XI, 10,  # C[thet',a,xi]
                  '-2102/14175, -158/315, 4/45, -2/3, 934/14175',
                  '-16/945, 16/45, 922/14175, -232/2835, 719/4725')
    ),
    Aux.MU: _Rdict(38,
       _Rtuple(Aux.PHI, 6,  # C[mu,phi]; even coeffs only
                  '9/16, -3/2, -15/32, 15/16, -35/48, 315/512'),
       _Rtuple(Aux.BETA, 6,  # C[mu,beta]; even coeffs only
                  '3/16, -1/2, 1/32, -1/16, -1/48, -5/512'),
       _Rtuple(Aux.THETA, 6,  # C[mu,theta]; even coeffs only
                  '13/16, 1/2, 33/32, -1/16, -5/16, -261/512'),
       # C[mu,mu] skipped
       _Rtuple(Aux.CHI, 10,  # C[mu,chi]
                  '41/180, 5/16, -2/3, 1/2, 557/1440',
                  '-3/5, 13/48, -103/140, 61/240, 49561/161280'),
       _Rtuple(Aux.XI, 10,  # C[mu,xi]
                  '-1609/28350, 121/1680, 4/45, -1/6, 16463/453600',
                  '26/945, -29/720, 449/28350, -1003/45360, -40457/2419200')
    ),
    Aux.CHI: _Rdict(50,
       _Rtuple(Aux.PHI, 10,  # C[chi,phi]
                  '-82/45, 4/3, 2/3, -2, -13/9',
                  '-16/15, 5/3, 34/21, -26/15, 1237/630'),
       _Rtuple(Aux.BETA, 10,  # C[chi,beta]
                  '-16/45, 0, 2/3, -1, 19/45',
                  '-2/5, 1/6, 16/105, -1/15, 17/1260'),
       _Rtuple(Aux.THETA, 10,  # C[chi,theta]
                  '-2/9, 2/3, 2/3, 0, 43/45',
                  '4/15, -1/3, 2/105, -2/5, -55/126'),
       _Rtuple(Aux.MU, 10,  # C[chi,mu]
                  '1/360, -37/96, 2/3, -1/2, 437/1440',
                  '-1/15, -1/48, 37/840, -17/480, -4397/161280'),
       # C[chi,chi] skipped
       _Rtuple(Aux.XI, 10,  # C[chi,xi]
                  '-2312/14175, -88/315, 34/45, -2/3, 6079/14175',
                  '-184/945, 1/45, 772/14175, -106/2835, -167/9450')
    ),
    Aux.XI: _Rdict(50,
       _Rtuple(Aux.PHI, 10,  # C[xi,phi]
                  '538/4725, 88/315, -4/45, -4/3, -2482/14175',
                  '8/105, 34/45, -898/14175, -1532/2835, 6007/14175'),
       _Rtuple(Aux.BETA, 10,  # C[xi,beta]
                  '34/675, 32/315, -4/45, -1/3, 74/2025',
                  '-4/315, -7/90, 2/14175, -83/2835, -797/56700'),
       _Rtuple(Aux.THETA, 10,  # C[xi,theta]
                  '778/4725, 62/105, -4/45, 2/3, 12338/14175',
                  '-32/315, 4/45, -1618/14175, -524/2835, -5933/14175'),
       _Rtuple(Aux.MU, 10,  # C[xi,mu]
                  '1297/18900, -817/10080, -4/45, 1/6, -29609/453600',
                  '-2/35, 49/720, -2917/56700, 4463/90720, 331799/7257600'),
       _Rtuple(Aux.CHI, 10,  # C[xi,chi]
                  '2458/4725, 46/315, -34/45, 2/3, 3413/14175',
                  '-256/315, 19/45, -15958/14175, 248/567, 16049/28350')  # PYCHOK exported
       # C[xi,xi] skipped
    )
})
# _ptrs_4 = (0,   0,   6,  12,  18,  28,  38,  44,  44,  50,  56,  66,
#           76,  82,  88,  88,  94, 104, 114, 120, 126, 132, 132, 142,
#          152, 162, 172, 182, 192, 192, 202, 212, 222, 232, 242, 252,
#          252)  # PYCHOK exported
del _Rcoeffs, _Rdict, _Rtuple

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
