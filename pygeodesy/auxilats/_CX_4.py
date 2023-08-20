
# -*- coding: utf-8 -*-

u'''Coeficients for C{_AUXLATITUDE_ORDER} 4 from I{Karney}'s C++ class U{AuxLatitude
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AuxLatitude.html>}
trancoded to a double, uniquified Python C{dict[auxout][auxin]}.

Copyright (C) Charles Karney (2022-2023) Karney@Alum.MIT.edu> and licensed under the
MIT/X11 License.  For more information, see <https:#GeographicLib.SourceForge.io>.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.auxilats.auxily import Aux, _Ufloats
from pygeodesy.constants import _0_0, _0_25, _0_5, _1_0, _N_1_0, \
                                _1_5, _2_0, _N_2_0, _4_0

__all__ = ()
__version__ = '23.08.19'

_f, _u = float, _Ufloats()
_coeffs_4 = _u._Coeffs(4, {  # GEOGRAPHICLIB_AUXLATITUDE_ORDER == 4
    Aux.PHI: {
        # C[phi,phi] skipped
        Aux.BETA: _u(  # C[phi,beta]; even coeffs only
            _0_0, _1_0,
            _0_0, _0_5,
            1 / _f(3),
            _0_25,),
        Aux.THETA: _u(  # C[phi,theta]; even coeffs only
            _N_2_0, _2_0,
            -_4_0, _2_0,
            8 / _f(3),
            _4_0,),
        Aux.MU: _u(  # C[phi,mu]; even coeffs only
            -27 / _f(32), _1_5,
            -55 / _f(32), 21 / _f(16),
            151 / _f(96),
            1097 / _f(512),),
        Aux.CHI: _u(  # C[phi,chi]
            116 / _f(45), _N_2_0, -2 / _f(3), _2_0,
            -227 / _f(45), -8 / _f(5), 7 / _f(3),
            -136 / _f(35), 56 / _f(15),
            4279 / _f(630),),
        Aux.XI: _u(  # C[phi,xi]
            -2582 / _f(14175), -16 / _f(35), 4 / _f(45), 4 / _f(3),
            -11966 / _f(14175), 152 / _f(945), 46 / _f(45),
            3802 / _f(14175), 3044 / _f(2835),
            6059 / _f(4725),)
    },
    Aux.BETA: {
        Aux.PHI: _u(  # C[beta,phi]; even coeffs only
            _0_0, _N_1_0,
            _0_0, _0_5,
            -1 / _f(3),
            _0_25,),
        # C[beta,beta] skipped
        Aux.THETA: _u(  # C[beta,theta]; even coeffs only
            _0_0, _1_0,
            _0_0, _0_5,
            1 / _f(3),
            _0_25,),
        Aux.MU: _u(  # C[beta,mu]; even coeffs only
            -9 / _f(32), _0_5,
            -37 / _f(96), 5 / _f(16),
            29 / _f(96),
            539 / _f(1536),),
        Aux.CHI: _u(  # C[beta,chi]
            38 / _f(45), -1 / _f(3), -2 / _f(3), _1_0,
            -7 / _f(9), -14 / _f(15), 5 / _f(6),
            -34 / _f(21), 16 / _f(15),
            2069 / _f(1260),),
        Aux.XI: _u(  # C[beta,xi]
            -1082 / _f(14175), -46 / _f(315), 4 / _f(45), 1 / _f(3),
            -338 / _f(2025), 68 / _f(945), 17 / _f(90),
            1102 / _f(14175), 461 / _f(2835),
            3161 / _f(18900),)
    },
    Aux.THETA: {
        Aux.PHI: _u(  # C[theta,phi]; even coeffs only
            _2_0, _N_2_0,
            -_4_0, _2_0,
            -8 / _f(3),
            _4_0,),
        Aux.BETA: _u(  # C[theta,beta]; even coeffs only
            _0_0, _N_1_0,
            _0_0, _0_5,
            -1 / _f(3),
            _0_25,),
        # C[theta,theta] skipped
        Aux.MU: _u(  # C[theta,mu]; even coeffs only
            -23 / _f(32), -1 / _f(2),
            -5 / _f(96), 5 / _f(16),
            1 / _f(32),
            283 / _f(1536),),
        Aux.CHI: _u(  # C[theta,chi]
            4 / _f(9), -2 / _f(3), -2 / _f(3), _0_0,
            -23 / _f(45), -4 / _f(15), 1 / _f(3),
            -24 / _f(35), 2 / _f(5),
            83 / _f(126),),
        Aux.XI: _u(  # C[thet),a,xi]
            -2102 / _f(14175), -158 / _f(315), 4 / _f(45), -2 / _f(3),
            934 / _f(14175), -16 / _f(945), 16 / _f(45),
            922 / _f(14175), -232 / _f(2835),
            719 / _f(4725),)
    },
    Aux.MU: {
        Aux.PHI: _u(  # C[mu,phi]; even coeffs only
            9 / _f(16), -3 / _f(2),
            -15 / _f(32), 15 / _f(16),
            -35 / _f(48),
            315 / _f(512),),
        Aux.BETA: _u(  # C[mu,beta]; even coeffs only
            3 / _f(16), -1 / _f(2),
            1 / _f(32), -1 / _f(16),
            -1 / _f(48),
            -5 / _f(512),),
        Aux.THETA: _u(  # C[mu,theta]; even coeffs only
            13 / _f(16), _0_5,
            33 / _f(32), -1 / _f(16),
            -5 / _f(16),
            -261 / _f(512),),
        # C[mu,mu] skipped
        Aux.CHI: _u(  # C[mu,chi]
            41 / _f(180), 5 / _f(16), -2 / _f(3), _0_5,
            557 / _f(1440), -3 / _f(5), 13 / _f(48),
            -103 / _f(140), 61 / _f(240),
            49561 / _f(161280),),
        Aux.XI: _u(  # C[mu,xi]
            -1609 / _f(28350), 121 / _f(1680), 4 / _f(45), -1 / _f(6),
            16463 / _f(453600), 26 / _f(945), -29 / _f(720),
            449 / _f(28350), -1003 / _f(45360),
            -40457 / _f(2419200),)
    },
    Aux.CHI: {
        Aux.PHI: _u(  # C[chi,phi]
            -82 / _f(45), 4 / _f(3), 2 / _f(3), _N_2_0,
            -13 / _f(9), -16 / _f(15), 5 / _f(3),
            34 / _f(21), -26 / _f(15),
            1237 / _f(630),),
        Aux.BETA: _u(  # C[chi,beta]
            -16 / _f(45), _0_0, 2 / _f(3), _N_1_0,
            19 / _f(45), -2 / _f(5), 1 / _f(6),
            16 / _f(105), -1 / _f(15),
            17 / _f(1260),),
        Aux.THETA: _u(  # C[chi,theta]
            -2 / _f(9), 2 / _f(3), 2 / _f(3), _0_0,
            43 / _f(45), 4 / _f(15), -1 / _f(3),
            2 / _f(105), -2 / _f(5),
            -55 / _f(126),),
        Aux.MU: _u(  # C[chi,mu]
            1 / _f(360), -37 / _f(96), 2 / _f(3), -1 / _f(2),
            437 / _f(1440), -1 / _f(15), -1 / _f(48),
            37 / _f(840), -17 / _f(480),
            -4397 / _f(161280),),
        # C[chi,chi] skipped
        Aux.XI: _u(  # C[chi,xi]
            -2312 / _f(14175), -88 / _f(315), 34 / _f(45), -2 / _f(3),
            6079 / _f(14175), -184 / _f(945), 1 / _f(45),
            772 / _f(14175), -106 / _f(2835),
            -167 / _f(9450),)
    },
    Aux.XI: {
        Aux.PHI: _u(  # C[xi,phi]
            538 / _f(4725), 88 / _f(315), -4 / _f(45), -4 / _f(3),
            -2482 / _f(14175), 8 / _f(105), 34 / _f(45),
            -898 / _f(14175), -1532 / _f(2835),
            6007 / _f(14175),),
        Aux.BETA: _u(  # C[xi,beta]
            34 / _f(675), 32 / _f(315), -4 / _f(45), -1 / _f(3),
            74 / _f(2025), -4 / _f(315), -7 / _f(90),
            2 / _f(14175), -83 / _f(2835),
            -797 / _f(56700),),
        Aux.THETA: _u(  # C[xi,theta]
            778 / _f(4725), 62 / _f(105), -4 / _f(45), 2 / _f(3),
            12338 / _f(14175), -32 / _f(315), 4 / _f(45),
            -1618 / _f(14175), -524 / _f(2835),
            -5933 / _f(14175),),
        Aux.MU: _u(  # C[xi,mu]
            1297 / _f(18900), -817 / _f(10080), -4 / _f(45), 1 / _f(6),
            -29609 / _f(453600), -2 / _f(35), 49 / _f(720),
            -2917 / _f(56700), 4463 / _f(90720),
            331799 / _f(7257600),),
        Aux.CHI: _u(  # C[xi,chi]
            2458 / _f(4725), 46 / _f(315), -34 / _f(45), 2 / _f(3),
            3413 / _f(14175), -256 / _f(315), 19 / _f(45),
            -15958 / _f(14175), 248 / _f(567),
            16049 / _f(28350),)  # PYCHOK exported
        # C[xi,xi] skipped
    }
})
# _ptrs_4 = (0,   0,   6,  12,  18,  28,  38,  44,  44,  50,  56,  66,
#           76,  82,  88,  88,  94, 104, 114, 120, 126, 132, 132, 142,
#          152, 162, 172, 182, 192, 192, 202, 212, 222, 232, 242, 252,
#          252)  # PYCHOK exported
del _f, _u

# **) MIT License
#
# Copyright (C) 2023-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
