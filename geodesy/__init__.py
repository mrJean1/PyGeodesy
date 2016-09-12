
# -*- coding: utf-8 -*-

# Python implementation of geodesy tools for various ellipsoidal and
# spherical earth models using trigonometric and vectorial mathematical
# methods.  Transcribed from JavaScript originals by (C) Chris Veness
# 2005-2016 and published under the same MIT Licence.  More details at
# <https://github.com/chrisveness/geodesy/>,
# <http://www.movable-type.co.uk/scripts/latlong.html>,
# <http://www.movable-type.co.uk/scripts/latlong-vincenty.html>,
# <http://www.movable-type.co.uk/scripts/latlong-vectors.html> and
# <http://www.movable-type.co.uk/scripts/latlong-os-gridref.html>

# all public contants, classes and functions
from datum import *  # __all__
from dms   import *  # __all__
from utils import *  # __all__
import ellipsoidalNvector
import ellipsoidalVincenty
import sphericalNvector
import sphericalTrigonometry

VincentyError = ellipsoidalVincenty.VincentyError

import datum as _datum, dms as _dms, utils as _utils
# all public contants, classes and functions
__all__ = _datum.__all__ + _dms.__all__ + (
          'ellipsoidalNvector', 'ellipsoidalVincenty',
          'sphericalNvector', 'sphericalTrigonometry',
          'VincentyError') + _utils.__all__
__version__ = '16.09.05'
del _datum, _dms, _utils
