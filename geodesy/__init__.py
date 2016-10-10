
# -*- coding: utf-8 -*-

# A Python implementation of geodesy tools for various ellipsoidal and
# spherical earth models using trigonometric and vector-based methods.

# Transcribed from JavaScript originals by (C) Chris Veness 2005-2016
# and published under the same MIT Licence.

# For more information and details see:
#
# <https://github.com/chrisveness/geodesy/>
# <http://www.movable-type.co.uk/scripts/latlong.html>
# <http://www.movable-type.co.uk/scripts/latlong-vincenty.html>
# <http://www.movable-type.co.uk/scripts/latlong-vectors.html>
# <http://www.movable-type.co.uk/scripts/latlong-utm-mgrs.html>

try:
    import datum as _  # PYCHOK expected
except ImportError:
    # extend sys.path for Python 3+
    import os, sys  # PYCHOK expected
    sys.path.insert(0, os.path.dirname(__file__))
    del os, sys

from datum import *  # PYCHOK __all__
from dms   import *  # PYCHOK __all__
from mgrs  import *  # PYCHOK __all__
from utils import *  # PYCHOK __all__
from utm   import *  # PYCHOK __all__
import ellipsoidalNvector  # PYCHOK false
import ellipsoidalVincenty  # PYCHOK false
import sphericalNvector  # PYCHOK false
import sphericalTrigonometry  # PYCHOK false

VincentyError = ellipsoidalVincenty.VincentyError

import datum as _datum, dms as _dms, mgrs as _mgrs, \
       utils as _utils, utm as _utm  # PYCHOK expected

# all public contants, classes and functions
__all__ = _datum.__all__ + _dms.__all__ + _mgrs.__all__ + (
          'ellipsoidalNvector', 'ellipsoidalVincenty',
          'sphericalNvector', 'sphericalTrigonometry',
          'VincentyError') + _utils.__all__ + _utm.__all__
__version__ = '16.10.10'

del _datum, _dms, _mgrs, _utils, _utm
