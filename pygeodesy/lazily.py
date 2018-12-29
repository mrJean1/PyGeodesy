
# -*- coding: utf-8 -*-

u'''Lazily import C{pygeodesy} modules and attributes.

I{Lazy import is supported only for U{Python 3.7+
<http://Snarky.CA/lazy-importing-in-python-3-7>} and is enabled by default
in U{PyGeodesy 18.11.10+<http://PyPI.org/project/PyGeodesy>}}.

To disable lazy import set environment variable C{PYGEODESY_LAZY_IMPORT}
to C{0} or an empty string.  Use C{2} or higher to print a message for
each lazily imported module and attribute, similar to environment variable
C{PYTHONVERBOSE} showing all imports.

Lazy import only applies to top-level modules of C{pygeodesy}.  A lazy
import of a top-level module also loads all sub-modules imported by
that top-level module.

Based on U{modutil.lazy_import<http://modutil.ReadTheDocs.io/en/latest/>}
by Brett Cannon's U{modutil<http://PyPI.org/project/modutil>}.

@var isLazy: Lazy import setting (C{int}) from C{PYGEODESY_LAZY_IMPORT},
             or C{None} if lazy import is not supported or not enabled,
             or C{False} if initializing lazy import failed.
'''

__all__ = ('LazyImportError',
           'isLazy')
__version__ = '18.12.26'

# @module_property[_RO?] <http://GitHub.com/jtushman/proxy_tools/
isLazy = None  # or 0..1+


class LazyImportError(ImportError):
    '''Lazy import is not supported, disabled or failed some other way.
    '''
    def __init__(self, fmt, *args):
        ImportError.__init__(self, (fmt % args) if args else fmt)


def _all_imports(*more):
    '''(INTERNAL) Build C{dict} of all lazy imports.
    '''
    imports = {}
    # imports naming conventions below:
    #  <module> == import <module>
    #  <module.name> == from <module> import <name>
    #  <module.name> as <othername> == from <module> import <name> as <as_name>
    for name in ('ellipsoidalKarney', 'ellipsoidalNvector', 'nvector', 'sphericalNvector', 'sphericalTrigonometry',
                 'points.areaOf as areaof', 'points.perimeterOf as perimeterof',  # DEPRECATED
                 'geohash', 'geohash.Geohash',
                 'ellipsoidalVincenty', 'ellipsoidalVincenty.VincentyError',
                 'vector3d', 'vector3d.CrossError', 'vector3d.crosserrors',
                 'bases', 'bases.LatLonHeightBase', 'bases.classname', 'bases.classnaming', 'bases.inStr',
                 'clipy', 'clipy.clipCS3', 'clipy.clipSH', 'clipy.clipSH3',
                 'datum', 'datum.R_M', 'datum.R_MA', 'datum.R_MB', 'datum.R_KM', 'datum.R_NM', 'datum.R_SM', 'datum.R_FM', 'datum.R_VM',
                          'datum.Datum', 'datum.Ellipsoid', 'datum.Transform', 'datum.Datums', 'datum.Ellipsoids', 'datum.Transforms',
                 'dms', 'dms.F_D', 'dms.F_DM', 'dms.F_DMS', 'dms.F_DEG', 'dms.F_MIN', 'dms.F_SEC', 'dms.F_RAD', 'dms.S_DEG', 'dms.S_MIN', 'dms.S_SEC', 'dms.S_RAD',
                        'dms.S_SEP', 'dms.RangeError', 'dms.bearingDMS', 'dms.clipDMS', 'dms.compassDMS', 'dms.compassPoint', 'dms.latDMS', 'dms.lonDMS',
                        'dms.normDMS', 'dms.parseDMS', 'dms.parseDMS2', 'dms.parse3llh', 'dms.precision', 'dms.rangerrors', 'dms.toDMS',
                 'elevations', 'elevations.elevation2', 'elevations.geoidHeight2',
                 'fmath', 'fmath.EPS', 'fmath.EPS1', 'fmath.Fsum', 'fmath.acos1', 'fmath.cbrt', 'fmath.cbrt2', 'fmath.favg', 'fmath.fdot', 'fmath.fdot3',
                          'fmath.fmean', 'fmath.fhorner', 'fmath.fpolynomial', 'fmath.fpowers', 'fmath.fStr', 'fmath.fStrzs', 'fmath.fsum', 'fmath.fsum_',
                          'fmath.hypot', 'fmath.hypot1', 'fmath.hypot3', 'fmath.isfinite', 'fmath.isint', 'fmath.isscalar', 'fmath.len2',
                          'fmath.map1', 'fmath.map2', 'fmath.scalar', 'fmath.sqrt3',
                 'formy', 'formy.antipode', 'formy.bearing', 'formy.bearing_', 'formy.compassAngle', 'formy.equirectangular', 'formy.equirectangular3',  # DEPRECATED
                          'formy.equirectangular_', 'formy.haversine', 'formy.haversine_', 'formy.heightOf', 'formy.horizon', 'formy.isantipode',
                 'lazily', 'lazily.isLazy', 'lazily.LazyImportError',
                 'lcc', 'lcc.Conic', 'lcc.Conics', 'lcc.Lcc', 'lcc.toLcc',
                 'mgrs', 'mgrs.Mgrs', 'mgrs.parseMGRS', 'mgrs.toMgrs',
                 'osgr', 'osgr.Osgr', 'osgr.parseOSGR', 'osgr.toOsgr',
                 'points', 'points.LatLon_', 'points.LatLon2psxy', 'points.Numpy2LatLon', 'points.Tuple2LatLon', 'points.areaOf', 'points.bounds',
                           'points.isclockwise', 'points.isconvex', 'points.isconvex_', 'points.isenclosedBy', 'points.isenclosedby', 'points.ispolar',
                           'points.nearestOn3', 'points.nearestOn4', 'points.perimeterOf',
                 'simplify', 'simplify.simplify1', 'simplify.simplify2', 'simplify.simplifyRDP', 'simplify.simplifyRDPm',
                             'simplify.simplifyRW', 'simplify.simplifyVW', 'simplify.simplifyVWm',
                 'utily', 'utily.PI', 'utily.PI2', 'utily.PI_2', 'utily.PI_4', 'utily.R_M', 'utily.LimitError', 'utily.anStr',
                          'utily.degrees', 'utily.degrees90', 'utily.degrees180', 'utily.degrees360', 'utily.degrees2m',
                          'utily.enStr2', 'utily.false2f', 'utily.ft2m', 'utily.halfs2',
                          'utily.issequence', 'utily.isNumpy2', 'utily.isPoints2', 'utily.isTuple2', 'utily.iterNumpy2', 'utily.iterNumpy2over',
                          'utily.limiterrors', 'utily.m2degrees', 'utily.m2ft', 'utily.m2km', 'utily.m2NM', 'utily.m2SM',
                          'utily.points2', 'utily.polygon', 'utily.property_RO',
                          'utily.radians', 'utily.radiansPI_2', 'utily.radiansPI', 'utily.radiansPI2', 'utily.tan_2', 'utily.tanPI_2_2',
                          'utily.unroll180', 'utily.unrollPI', 'utily.unStr',
                          'utily.wrap90', 'utily.wrap180', 'utily.wrap360', 'utily.wrapPI_2','utily.wrapPI', 'utily.wrapPI2',
                 'utm', 'utm.Utm', 'utm.UTMError', 'utm.parseUTM', 'utm.toUtm', 'utm.utmZoneBand2',
                 'webmercator', 'webmercator.Wm', 'webmercator.parseWM', 'webmercator.toWm') + more:
        import_from, _, as_name = name.partition(' as ')
        if not as_name:
            _, _, as_name = import_from.rpartition('.')
        imports[as_name] = import_from
    return imports


def _all_missing2(_all_):
    '''(INTERNAL) Get deltas between pygeodesy.__all__ and lazily._all_imports.
    '''
    _alzy = _all_imports('pygeodesy_abspath', 'version')  # dict
    return (('lazily._all_imports', ', '.join(a for a in _all_ if a not in _alzy)),
            ('pygeodesy.__all__',   ', '.join(a for a in _alzy if a not in _all_)))


def _lazy_import2(package_name):  # MCCABE 20
    '''Check for and set up lazy importing.

       @param package_name: The name of the package (C{str}) performing
                            the imports, to help facilitate resolving
                            relative imports.

       @return: 2-Tuple (package, getattr) of the importing package for
                easy reference within itself the a callable to be set
                to `__getattr__`.

       @raise LazyImportError: Lazy import not supported, an import
                               failed or a module name or attribute
                               name is invalid or does not exist.

       @note: This is the original function U{modutil.lazy_import
              <http://GitHub.com/brettcannon/modutil/blob/master/modutil.py>}
              modified to handle the C{__all__} and C{__dir__} attributes
              and call C{importlib.import_module(<module>.<name>, ...)}
              without causing a C{ModuleNotFoundError}.

       @see: The original U{modutil<http://PyPi.org/project/modutil>},
             U{PEP 562<http://www.Python.org/dev/peps/pep-0562>} and
             U{Werkzeug<http://GitHub.com/pallets/werkzeug/blob/master/werkzeug/__init__.py>}.
    '''
    global isLazy

    import sys
    if sys.version_info[:2] < (3, 7):  # not supported
        raise LazyImportError('no %s.%s for Python %s', package_name,
                             _lazy_import2.__name__, sys.version.split()[0])

    import os
    z = os.environ.get('PYGEODESY_LAZY_IMPORT', None)
    if z is None:  # PYGEODESY_LAZY_IMPORT not set
        isLazy = 1  # on by default on 3.7
    else:
        z = z.strip()  # like PYTHONVERBOSE et.al.
        isLazy = int(z) if z.isdigit() else (1 if z else 0)
    if isLazy < 1:  # not enabled
        raise LazyImportError('env %s=%r', 'PYGEODESY_LAZY_IMPORT', z)
    if os.environ.get('PYTHONVERBOSE', None):
        isLazy += 1

    try:  # to initialize
        from importlib import import_module

        package = import_module(package_name)
        parent = package.__spec__.parent  # __spec__ only in Python 3.7+
        if parent != package_name:  # assertion
            raise ImportError('parent %r vs %r' % (parent, package_name))
    except (AttributeError, ImportError) as x:
        isLazy = False  # failed
        raise LazyImportError('init failed: %s', x)

    if isLazy > 2:  # trim import path names
        cwdir = os.getcwd()
        cwdir = cwdir[:-len(os.path.basename(cwdir))]
    else:  # no import path names
        cwdir = ''

    imports = _all_imports()

    def __getattr__(as_name):  # __getattr__ only for Python 3.7+
        # only called once for each undefined pygeodesy attribute
        if as_name in imports:
            # importlib.import_module() implicitly sets sub-modules
            # on this module as appropriate for direct imports (see
            # note in the _lazy_import.__doc__ above).
            module_name, _, attr_name = imports[as_name].partition('.')
            if module_name not in imports:
                raise LazyImportError('no %s %s.%s', 'module', package_name, module_name)
            imported = import_module(module_name, parent)
            if attr_name:
                if not hasattr(imported, attr_name):
                    raise LazyImportError('no %s %s.%s', 'attribute', module_name, attr_name)
                imported = getattr(imported, attr_name)
        elif as_name in ('__all__',):  # '__dir__', '__members__'
            imported = imports.keys()
            module_name = ''
        else:
            raise LazyImportError('no %s %s.%s', 'attribute', package_name, as_name)

        setattr(package, as_name, imported)
        if isLazy > 1:
            m = n = ''
            if module_name and module_name != as_name:
                m = ' from .%s' % (module_name,)
            if isLazy > 2:
                # sys._getframe(1) ... 'importlib._bootstrap' line 1032
                f = sys._getframe(2)  # import line ...
                n = f.f_code.co_filename
                if cwdir and n.startswith(cwdir):
                    n = n[len(cwdir):]
                n = ' by %r line %d' % (n, f.f_lineno)
            print('# lazily imported %s.%s%s%s' % (package_name, as_name, m, n))
        return imported

    return package, __getattr__

# **) MIT License
#
# Copyright (C) 2018-2019 -- mrJean1 at Gmail dot com
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
