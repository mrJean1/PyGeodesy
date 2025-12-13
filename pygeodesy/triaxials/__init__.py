
# -*- coding: utf-8 -*-

u'''Package of lazily imported modules L{triaxials.conformal3}, L{triaxials.triaxial3}
and L{triaxials.triaxial5} for triaxial ellipsoids.

Mostly transcoded to pure Python from I{Karney}'s U{GeographicLib 2.7 Triaxial<https://
GeographicLib.SourceForge.io/C++/doc/namespaceGeographicLib_1_1Triaxial.html} classes
and the I{experimental} U{GeographicLib 2.52 Jacobi<https://GeographicLib.SourceForge.io/
C++/2.5.2/jacobi.html>} class.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2008-2024, 2024-2025) and
licensed under the MIT/X11 License.  For more information, see the U{GeographicLib 2.5.2
and 2.7<https://GeographicLib.SourceForge.io>} documentation.
'''
from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER, _lazy_import_as, _unLazy0
# from pygeodesy.triaxials.triaxial5 import *  # PYCHOK for backward compatibility
# from pygeodesy.triaxials.bases import TriaxialError  # likewise

__all__ = _ALL_LAZY.triaxials  # _triaxial5  # likewise
__version__ = '25.12.04'

if _unLazy0:  # or _isfrozen
    from pygeodesy.triaxials.bases import LLK, TriaxialError
    from pygeodesy.triaxials.conformal3 import BetOmgGam5Tuple, Conformal3, Conformal3B, \
                                               Conformal3Sphere, Conformal5Tuple
    from pygeodesy.triaxials.triaxial3 import BetOmgAlp5Tuple, Cartesian5Tuple, PhiLamZet5Tuple, \
                                              Triaxial3, Triaxial3B, Triaxial3s
    from pygeodesy.triaxials.triaxial5 import BetaOmega2Tuple, BetaOmega3Tuple, \
                                              Conformal, ConformalSphere, Conformal2Tuple, \
                                              Triaxial, Triaxial_, Triaxials, \
                                              hartzell4, height4

    __all__ += _ALL_OTHER(LLK, TriaxialError,
                          BetOmgGam5Tuple, Conformal3, Conformal3B, ConformalSphere,
                          Conformal5Tuple,
                          BetOmgAlp5Tuple, Cartesian5Tuple, PhiLamZet5Tuple,
                          Triaxial3, Triaxial3B, Triaxial3s,
                          BetaOmega2Tuple, BetaOmega3Tuple,
                          Conformal, Conformal3Sphere, Conformal2Tuple, Triaxial,
                          Triaxial_, Triaxials, hartzell4, height4)
#   assert set(_ALL_LAZY.triaxials + _ALL_LAZY.triaxials_bases
#                                  + _ALL_LAZY.triaxials_conformal3
#                                  + _ALL_LAZY.triaxials_triaxial3
#                                  + _ALL_LAZY.triaxials_triaxial5) == set(__all__)

else:  # lazily import modules and exported attrs
    __getattr__ = _lazy_import_as(__name__)

# **) MIT License
#
# Copyright (C) 2025-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
