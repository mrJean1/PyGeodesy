
# -*- coding: utf-8 -*-

from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER, _lazy_import_as, _unLazy0

__all__ = _ALL_LAZY.rhumb
__version__ = '23.11.30'

if _unLazy0:
    from pygeodesy.rhumb.aux import RhumbAux, RhumbLineAux
    from pygeodesy.rhumb.ekx import Rhumb, RhumbLine
    from pygeodesy.rhumb.solve import RhumbSolve, RhumbLineSolve, RhumbSolve7Tuple

    __all__ += _ALL_OTHER(RhumbAux, RhumbLineAux, Rhumb, RhumbLine,
                          RhumbSolve, RhumbLineSolve, RhumbSolve7Tuple)
    assert _ALL_LAZY.rhumb_aux + _ALL_LAZY.rhumb_ekx + _ALL_LAZY.rhumb_solve == __all__

else:  # lazily import modules only
    __getattr__ = _lazy_import_as(__name__)

# **) MIT License
#
# Copyright (C) 2018-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
