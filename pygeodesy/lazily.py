
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

@var isLazy: Lazy import setting (C{int} 0, 1, 2 or 3+) from environment
             variable C{PYGEODESY_LAZY_IMPORT}, or C{None} if lazy import
             is not supported or not enabled, or C{False} if initializing
             lazy import failed.
'''

# @module_property[_RO?] <http://GitHub.com/jtushman/proxy_tools/>
isLazy = None  # see @var isLazy above


class LazyImportError(ImportError):
    '''Lazy import is not supported, disabled or failed some other way.
    '''
    def __init__(self, fmt, *args):
        ImportError.__init__(self, (fmt % args) if args else fmt)


class _Enum_RO(dict):
    '''(INTERNAL) C{Read_Only} enum-like C{dict} sub-class.
    '''
    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:
            raise AttributeError("%s.%s doesn't exist" % (self._name, attr))  # PYCHOK expected

    def __setattr__(self, attr, value):
        raise TypeError('Read_Only %s.%s = %r' % (self._name, attr, value))  # PYCHOK expected


# DEPRECATED __all__ values for backward compatibility
_ALL_DEPRECATED = _Enum_RO(_name='_ALL_DEPRECATED',
                           formy=('equirectangular3',),  # DEPRECATED
                          points=('areaOf as areaof', 'boundsOf as bounds',  # DEPRECATED
                                  'isenclosedBy as isenclosedby', 'perimeterOf as perimeterof'))

# __all__ value for most modules, accessible as _ALL_LAZY.<module>
_ALL_LAZY = _Enum_RO(_name='_ALL_LAZY',
                     bases=('LatLonHeightBase', 'classname', 'classnaming', 'inStr'),
                     clipy=('clipCS3', 'clipSH', 'clipSH3'),
                     datum=('R_M', 'R_MA', 'R_MB', 'R_KM', 'R_NM', 'R_SM', 'R_FM', 'R_VM',
                            'Datum', 'Ellipsoid', 'Transform', 'Datums', 'Ellipsoids', 'Transforms'),
                       dms=('F_D', 'F_DM', 'F_DMS', 'F_DEG', 'F_MIN', 'F_SEC', 'F_RAD', 'S_DEG', 'S_MIN', 'S_SEC', 'S_RAD', 'S_SEP',
                            'RangeError', 'bearingDMS', 'clipDMS', 'compassDMS', 'compassPoint', 'latDMS', 'lonDMS',
                            'normDMS', 'parseDMS', 'parseDMS2', 'parse3llh', 'precision', 'rangerrors', 'toDMS'),
                elevations=('elevation2', 'geoidHeight2'),
         ellipsoidalKarney=(),  # module only
        ellipsoidalNvector=(),  # module only
       ellipsoidalVincenty=('VincentyError',),  # nothing else
                     fmath=('EPS', 'EPS1', 'Fdot', 'Fhorner', 'Fpolynomial', 'Fsum',
                            'acos1', 'cbrt', 'cbrt2',
                            'favg', 'fdot', 'fdot3', 'fmean', 'fhorner', 'fpolynomial', 'fpowers', 'fprod', 'frange', 'freduce', 'fStr', 'fStrzs', 'fsum', 'fsum_',
                            'hypot', 'hypot1', 'hypot3', 'isfinite', 'isint', 'isscalar', 'len2', 'map1', 'map2', 'scalar', 'sqrt3'),
                     formy=('antipode', 'bearing', 'bearing_', 'compassAngle', 'euclidean', 'euclidean_', 'equirectangular', 'equirectangular_',
                            'haversine', 'haversine_', 'heightOf', 'horizon', 'isantipode', 'vincentys', 'vincentys_'),
                   geohash=('Geohash',),  # nothing else
                    geoids=('GeoidError', 'GeoidG2012B', 'GeoidKarney', 'GeoidPGM', 'egmGeoidHeights', 'PGMError'),
                   heights=('HeightError', 'SciPyError', 'SciPyWarning',
                            'HeightCubic', 'HeightIDW', 'HeightIDW2', 'HeightIDW3', 'HeightLinear', 'HeightLSQBiSpline', 'HeightSmoothBiSpline'),
                    lazily=('LazyImportError', 'isLazy'),
                       lcc=('Conic', 'Conics', 'Lcc', 'toLcc'),
                      mgrs=('Mgrs', 'parseMGRS', 'toMgrs'),
                   nvector=(),  # module only
                      osgr=('Osgr', 'parseOSGR', 'toOsgr'),
                    points=('LatLon_', 'LatLon2psxy', 'Numpy2LatLon', 'Tuple2LatLon',
                            'areaOf', 'boundsOf', 'centroidOf',
                            'isclockwise', 'isconvex', 'isconvex_', 'isenclosedBy', 'ispolar',
                            'nearestOn3', 'nearestOn4', 'nearestOn5', 'perimeterOf'),
          sphericalNvector=(),  # module only
     sphericalTrigonometry=(),  # module only
                  simplify=('simplify1', 'simplify2', 'simplifyRDP', 'simplifyRDPm', 'simplifyRW', 'simplifyVW', 'simplifyVWm'),
                     utily=('PI', 'PI2', 'PI_2', 'PI_4', 'R_M', 'LimitError',
                            'anStr', 'clipStr',
                            'degrees', 'degrees90', 'degrees180', 'degrees360', 'degrees2m',
                            'enStr2', 'false2f', 'ft2m', 'halfs2',
                            'issequence', 'isNumpy2', 'isPoints2', 'isTuple2', 'iterNumpy2', 'iterNumpy2over',
                            'limiterrors', 'm2degrees', 'm2ft', 'm2km', 'm2NM', 'm2SM',
                            'points2', 'polygon', 'property_RO',
                            'radians', 'radiansPI_2', 'radiansPI', 'radiansPI2',
                            'sincos2', 'sincos2d', 'splice', 'tan_2', 'tanPI_2_2',
                            'unroll180', 'unrollPI', 'unStr',
                            'wrap90', 'wrap180', 'wrap360', 'wrapPI_2','wrapPI', 'wrapPI2'),
                       utm=('Utm', 'UTMError', 'parseUTM', 'toUtm', 'utmZoneBand2'),
                  vector3d=('CrossError', 'crosserrors'),  # nothing else
               webmercator=('Wm', 'parseWM', 'toWm'))

__all__ = _ALL_LAZY.lazily
__version__ = '19.04.05'


def _all_imports(**more):
    '''(INTERNAL) Build C{dict} of all lazy imports.
    '''
    # imports naming conventions stored below - [<key>] = <from>:
    #  import <module>                        - [<module>] = <module>
    #  from <module> import <attr>            - [<attr>] = <module>
    #  from pygeodesy import <attr>           - [<attr>] = <attr>
    #  from <module> import <attr> as <name>  - [<name>] = <module>.<attr>
    imports = {}
    for _all_ in (_ALL_LAZY, _ALL_DEPRECATED, more):
        for mod_name, names in _all_.items():
            if isinstance(names, tuple):  # and not mod_name.startswith('_'):
                imports[mod_name] = mod_name
                for name in names:
                    name, _, _as_ = name.partition(' as ')
                    if _as_:
                        imports[_as_] = mod_name + '.' + name
                    else:
                        imports[name] = mod_name
    return imports


def _all_missing2(_all_):
    '''(INTERNAL) Get deltas between pygeodesy.__all__ and lazily._all_imports.
    '''
    _alzy = _all_imports(pygeodesy_abspath=(), version=())
    return (('lazily._all_imports', ', '.join(a for a in _all_ if a not in _alzy)),
            ('pygeodesy.__all__',   ', '.join(a for a in _alzy if a not in _all_)))


def _lazy_import2(pack_name):  # MCCABE 23
    '''Check for and set up lazy importing.

       @param pack_name: The name of the package (C{str}) performing
                         the imports, to help facilitate resolving
                         relative imports.

       @return: 2-Tuple (package, getattr) of the importing package for
                easy reference within itself and the callable to be set
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
        raise LazyImportError('no %s.%s for Python %s', pack_name,
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

        package = import_module(pack_name)
        parent = package.__spec__.parent  # __spec__ only in Python 3.7+
        if parent != pack_name:  # assertion
            raise ImportError('parent %r vs %r' % (parent, pack_name))
    except (AttributeError, ImportError) as x:
        isLazy = False  # failed
        raise LazyImportError('init failed: %s', x)

    if isLazy > 2:  # trim import path names
        cwdir = os.getcwd()
        cwdir = cwdir[:-len(os.path.basename(cwdir))]
    else:  # no import path names
        cwdir = ''

    del os, z

    imports = _all_imports()

    def __getattr__(name):  # __getattr__ only for Python 3.7+
        # only called once for each undefined pygeodesy attribute
        if name in imports:
            # importlib.import_module() implicitly sets sub-modules
            # on this module as appropriate for direct imports (see
            # note in the _lazy_import.__doc__ above).
            mod_name, _, attr = imports[name].partition('.')
            if mod_name not in imports:
                raise LazyImportError('no %s %s.%s', 'module', pack_name, mod_name)
            imported = import_module(mod_name, parent)
            try:  # import the module attribute
                if attr:
                    imported = getattr(imported, attr)
                elif name != mod_name:
                    imported = getattr(imported, name)
            except AttributeError:
                raise LazyImportError('no %s %s.%s', 'attribute', mod_name, attr or name)

        elif name in ('__all__',):  # XXX '__dir__', '__members__'?
            imported = imports.keys()
            mod_name = ''
        else:
            raise LazyImportError('no %s %s.%s', 'module or attribute', pack_name, name)

        setattr(package, name, imported)
        if isLazy > 1:
            z = ''
            if mod_name and mod_name != name:
                z = ' from .%s' % (mod_name,)
            if isLazy > 2:
                # sys._getframe(1) ... 'importlib._bootstrap' line 1032,
                # may throw a ValueError('call stack not deep enough')
                try:
                    f = sys._getframe(2)  # importing line ...
                    n = f.f_code.co_filename
                    if cwdir and n.startswith(cwdir):
                        n = n[len(cwdir):]
                    z = '%s by %s line %d' % (z, n, f.f_lineno)
                except ValueError:
                    pass
            print('# lazily imported %s.%s%s' % (pack_name, name, z))

        return imported  # __getattr__

    return package, __getattr__  # _lazy_import2

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
