
# -*- coding: utf-8 -*-

u'''Coefficients for C{_AUXLATITUDE_ORDER} 6 from I{Karney}'s C++ class U{AuxLatitude
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AuxLatitude.html>}
transcoded to a Python C{_Rdict[auxout][auxin]} of C{_Rtuple}s.

Copyright (C) Charles Karney (2022-2024) Karney@Alum.MIT.edu> and licensed under the
MIT/X11 License.  For more information, see <https:#GeographicLib.SourceForge.io>.
'''
from pygeodesy.auxilats.auxily import Aux
from pygeodesy.auxilats._CX_Rs import _Rcoeffs, _Rdict, _Rtuple

__all__ = ()
__version__ = '24.09.02'

_coeffs_6 = _Rcoeffs(6, {  # GEOGRAPHICLIB_AUXLATITUDE_ORDER == 6
    Aux.PHI: _Rdict(78,
       # C[phi,phi] skipped
       _Rtuple(Aux.BETA, 12,  # C[phi,beta]; even coeffs only
                  '0, 0, 1, 0, 0, 1/2, 0, 1/3, 0, 1/4, 1/5, 1/6'),
       _Rtuple(Aux.THETA, 12,  # C[phi,theta]; even coeffs only
                  '2, -2, 2, 6, -4, 2, -8, 8/3, -16, 4, 32/5, 32/3'),
       _Rtuple(Aux.MU, 12,  # C[phi,mu]; even coeffs only
                  '269/512, -27/32, 3/2, 6759/4096, -55/32, 21/16, -417/128',
                  '151/96, -15543/2560, 1097/512, 8011/2560, 293393/61440'),
       _Rtuple(Aux.CHI, 21,  # C[phi,chi]
                  '-2854/675, 26/45, 116/45, -2, -2/3, 2, 2323/945, 2704/315, -227/45',
                  '-8/5, 7/3, 73814/2835, -1262/105, -136/35, 56/15, -399572/14175',
                  '-332/35, 4279/630, -144838/6237, 4174/315, 601676/22275'),
       _Rtuple(Aux.XI, 21,  # C[phi,xi]
                  '28112932/212837625, 60136/467775, -2582/14175, -16/35, 4/45',
                  '4/3, 251310128/638512875, -21016/51975, -11966/14175, 152/945',
                  '46/45, -8797648/10945935, -94388/66825, 3802/14175, 3044/2835',
                  '-1472637812/638512875, 41072/93555, 6059/4725, 455935736/638512875',
                  '768272/467775, 4210684958/1915538625')
    ),
    Aux.BETA: _Rdict(78,
       _Rtuple(Aux.PHI, 12,  # C[beta,phi]; even coeffs only
                  '0, 0, -1, 0, 0, 1/2, 0, -1/3, 0, 1/4, -1/5, 1/6'),
       # C[beta,beta] skipped
       _Rtuple(Aux.THETA, 12,  # C[beta,theta]; even coeffs only
                  '0, 0, 1, 0, 0, 1/2, 0, 1/3, 0, 1/4, 1/5, 1/6'),
       _Rtuple(Aux.MU, 12,  # C[beta,mu]; even coeffs only
                  '205/1536, -9/32, 1/2, 1335/4096, -37/96, 5/16, -75/128',
                  '29/96, -2391/2560, 539/1536, 3467/7680, 38081/61440'),
       _Rtuple(Aux.CHI, 21,  # C[beta,chi]
                  '-3118/4725, -1/3, 38/45, -1/3, -2/3, 1, -247/270, 50/21, -7/9',
                  '-14/15, 5/6, 17564/2835, -5/3, -34/21, 16/15, -49877/14175',
                  '-28/9, 2069/1260, -28244/4455, 883/315, 797222/155925'),
       _Rtuple(Aux.XI, 21,  # C[beta,xi]
                  '7947332/212837625, 11824/467775, -1082/14175, -46/315, 4/45',
                  '1/3, 39946703/638512875, -16672/155925, -338/2025, 68/945',
                  '17/90, -255454/1563705, -101069/467775, 1102/14175, 461/2835',
                  '-189032762/638512875, 1786/18711, 3161/18900',
                  '80274086/638512875, 88868/467775, 880980241/3831077250')
    ),
    Aux.THETA: _Rdict(78,
       _Rtuple(Aux.PHI, 12,  # C[theta,phi]; even coeffs only
                  '-2, 2, -2, 6, -4, 2, 8, -8/3, -16, 4, -32/5, 32/3'),
       _Rtuple(Aux.BETA, 12,  # C[theta,beta]; even coeffs only
                  '0, 0, -1, 0, 0, 1/2, 0, -1/3, 0, 1/4, -1/5, 1/6'),
       # C[theta,theta] skipped
       _Rtuple(Aux.MU, 12,  # C[theta,mu]; even coeffs only
                  '499/1536, -23/32, -1/2, 6565/12288, -5/96, 5/16, -77/128',
                  '1/32, -4037/7680, 283/1536, 1301/7680, 17089/61440'),
       _Rtuple(Aux.CHI, 21,  # C[theta,chi]
                  '-3658/4725, 2/9, 4/9, -2/3, -2/3, 0, 61/135, 68/45, -23/45',
                  '-4/15, 1/3, 9446/2835, -46/35, -24/35, 2/5, -34712/14175',
                  '-80/63, 83/126, -2362/891, 52/45, 335882/155925'),
       _Rtuple(Aux.XI, 21,  # C[theta,xi]
                  '216932/2627625, 109042/467775, -2102/14175, -158/315, 4/45',
                  '-2/3, 117952358/638512875, -7256/155925, 934/14175, -16/945',
                  '16/45, -7391576/54729675, -25286/66825, 922/14175, -232/2835',
                  '-67048172/638512875, 268/18711, 719/4725',
                  '46774256/638512875, 14354/467775, 253129538/1915538625')
    ),
    Aux.MU: _Rdict(78,
       _Rtuple(Aux.PHI, 12,  # C[mu,phi]; even coeffs only
                  '-3/32, 9/16, -3/2, 135/2048, -15/32, 15/16, 105/256, -35/48',
                  '-189/512, 315/512, -693/1280, 1001/2048'),
       _Rtuple(Aux.BETA, 12,  # C[mu,beta]; even coeffs only
                  '-1/32, 3/16, -1/2, -9/2048, 1/32, -1/16, 3/256, -1/48',
                  '3/512, -5/512, -7/1280, -7/2048'),
       _Rtuple(Aux.THETA, 12,  # C[mu,theta]; even coeffs only
                  '-15/32, 13/16, 1/2, -1673/2048, 33/32, -1/16, 349/256, -5/16',
                  '963/512, -261/512, -921/1280, -6037/6144'),
       # C[mu,mu] skipped
       _Rtuple(Aux.CHI, 21,  # C[mu,chi]
                  '7891/37800, -127/288, 41/180, 5/16, -2/3, 1/2',
                  '-1983433/1935360, 281/630, 557/1440, -3/5, 13/48',
                  '167603/181440, 15061/26880, -103/140, 61/240',
                  '6601661/7257600, -179/168, 49561/161280',
                  '-3418889/1995840, 34729/80640, 212378941/319334400'),
       _Rtuple(Aux.XI, 21,  # C[mu,xi]
                  '12674323/851350500, -384229/14968800, -1609/28350, 121/1680',
                  '4/45, -1/6, -31621753811/1307674368000, -431/17325, 16463/453600',
                  '26/945, -29/720, -32844781/1751349600, 3746047/119750400, 449/28350',
                  '-1003/45360, 10650637121/326918592000, 629/53460, -40457/2419200',
                  '205072597/20432412000, -1800439/119750400, -59109051671/3923023104000')
    ),
    Aux.CHI: _Rdict(105,
       _Rtuple(Aux.PHI, 21,  # C[chi,phi]
                  '4642/4725, 32/45, -82/45, 4/3, 2/3, -2, -1522/945, 904/315',
                  '-13/9, -16/15, 5/3, -12686/2835, 8/5, 34/21, -26/15, -24832/14175',
                  '-12/5, 1237/630, 109598/31185, -734/315, 444337/155925'),
       _Rtuple(Aux.BETA, 21,  # C[chi,beta]
                  '-998/4725, 2/5, -16/45, 0, 2/3, -1, -2/27, -22/105, 19/45',
                  '-2/5, 1/6, 116/567, -22/105, 16/105, -1/15, 2123/14175, -8/105',
                  '17/1260, 128/4455, -1/105, 149/311850'),
       _Rtuple(Aux.THETA, 21,  # C[chi,theta]
                  '1042/4725, -14/45, -2/9, 2/3, 2/3, 0, -712/945, -4/45, 43/45',
                  '4/15, -1/3, 274/2835, 124/105, 2/105, -2/5, 21068/14175, -16/105',
                  '-55/126, -9202/31185, -22/45, -90263/155925'),
       _Rtuple(Aux.MU, 21,  # C[chi,mu]
                  '-96199/604800, 81/512, 1/360, -37/96, 2/3, -1/2',
                  '1118711/3870720, -46/105, 437/1440, -1/15, -1/48',
                  '-5569/90720, 209/4480, 37/840, -17/480',
                  '830251/7257600, 11/504, -4397/161280',
                  '108847/3991680, -4583/161280, -20648693/638668800'),
       # C[chi,chi] skipped
       _Rtuple(Aux.XI, 21,  # C[chi,xi]
                  '-55271278/212837625, 27128/93555, -2312/14175, -88/315, 34/45',
                  '-2/3, 106691108/638512875, -65864/155925, 6079/14175, -184/945',
                  '1/45, 5921152/54729675, -14246/467775, 772/14175, -106/2835',
                  '75594328/638512875, -5312/467775, -167/9450',
                  '2837636/638512875, -248/13365, -34761247/1915538625')
    ),
    Aux.XI: _Rdict(105,
       _Rtuple(Aux.PHI, 21,  # C[xi,phi]
                  '-44732/2837835, 20824/467775, 538/4725, 88/315, -4/45, -4/3',
                  '-12467764/212837625, -37192/467775, -2482/14175, 8/105, 34/45',
                  '100320856/1915538625, 54968/467775, -898/14175, -1532/2835',
                  '-5884124/70945875, 24496/467775, 6007/14175',
                  '-839792/19348875, -23356/66825, 570284222/1915538625'),
       _Rtuple(Aux.BETA, 21,  # C[xi,beta]
                  '-70496/8513505, 2476/467775, 34/675, 32/315, -4/45, -1/3',
                  '53836/212837625, 3992/467775, 74/2025, -4/315, -7/90',
                  '-661844/1915538625, 7052/467775, 2/14175, -83/2835',
                  '1425778/212837625, 934/467775, -797/56700',
                  '390088/212837625, -3673/467775, -18623681/3831077250'),
       _Rtuple(Aux.THETA, 21,  # C[xi,theta]
                  '-4286228/42567525, -193082/467775, 778/4725, 62/105, -4/45, 2/3',
                  '-61623938/70945875, 92696/467775, 12338/14175, -32/315, 4/45',
                  '427003576/1915538625, 612536/467775, -1618/14175, -524/2835',
                  '427770788/212837625, -8324/66825, -5933/14175',
                  '-9153184/70945875, -320044/467775, -1978771378/1915538625'),
       _Rtuple(Aux.MU, 21,  # C[xi,mu]
                  '-9292991/302702400, 7764059/239500800, 1297/18900, -817/10080, -4/45, 1/6',
                  '36019108271/871782912000, 35474/467775, -29609/453600, -2/35, 49/720',
                  '3026004511/30648618000, -4306823/59875200, -2917/56700, 4463/90720',
                  '-368661577/4036032000, -102293/1871100, 331799/7257600',
                  '-875457073/13621608000, 11744233/239500800, 453002260127/7846046208000'),
       _Rtuple(Aux.CHI, 21,  # C[xi,chi]
                  '2706758/42567525, -55222/93555, 2458/4725, 46/315, -34/45, 2/3',
                  '-340492279/212837625, 516944/467775, 3413/14175, -256/315, 19/45',
                  '4430783356/1915538625, 206834/467775, -15958/14175, 248/567',
                  '62016436/70945875, -832976/467775, 16049/28350',
                  '-651151712/212837625, 15602/18711, 2561772812/1915538625')  # PYCHOK exported
       # C[xi,xi] skipped
    )
})
# _ptrs_6 = (0,   0,  12,  24,  36,  57,  78,  90,  90, 102, 114, 135,
#          156, 168, 180, 180, 192, 213, 234, 246, 258, 270, 270, 291,
#          312, 333, 354, 375, 396, 396, 417, 438, 459, 480, 501, 522,
#          522)  # PYCHOK exported
del _Rcoeffs, _Rdict, _Rtuple

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
