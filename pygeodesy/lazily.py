
# -*- coding: utf-8 -*-

u'''Lazily import C{pygeodesy} modules and attributes, based on
U{lazy_import<https://modutil.ReadTheDocs.io/en/latest/#lazy_import>}
from Brett Cannon's U{modutil<https://PyPI.org/project/modutil>}.

C{Lazy import} is I{supported only for U{Python 3.7+
<https://Snarky.CA/lazy-importing-in-python-3-7>}} and is I{enabled by
default in U{PyGeodesy 18.11.10+<https://PyPI.org/project/PyGeodesy>}
and later}.

To disable C{lazy import}, set environment variable C{PYGEODESY_LAZY_IMPORT}
to C{0} or an empty string.  Use C{2} or higher to print a message for
each lazily imported module and attribute, similar to environment variable
C{PYTHONVERBOSE} showing imports.  Using C{3} or higher also shows the
importing file name and line number.

@note: C{Lazy import} applies only to top-level modules of C{pygeodesy}.
A C{lazy import} of a top-level module also loads all sub-modules
imported by that top-level module.

@var isLazy: Lazy import setting (C{int} 0, 1, 2 or 3+) from environment
             variable C{PYGEODESY_LAZY_IMPORT}, or C{None} if C{lazy import}
             is not supported or not enabled, or C{False} if initializing
             C{lazy import} failed.
'''
from pygeodesy.interns import _areaOf_, _COMMA_SPACE_, _doesn_t_exist_, \
                              _DOT_, _dot_, _dunder_name, _isclockwise_, \
                              _ispolar_, _item_sq, _Missing, NN, \
                              _perimeterOf_, _UNDERSCORE_

from os import environ as _environ

_FOR_DOCS = _environ.get('PYGEODESY_FOR_DOCS', None)  # for epydoc ...

# @module_property[_RO?] <https://GitHub.com/jtushman/proxy_tools/>
isLazy = None  # see @var isLazy above


class LazyImportError(ImportError):
    '''Raised if C{lazy import} is not supported, disabled or failed some other way.
    '''
    def __init__(self, *name_value, **txt):
        from pygeodesy.errors import _error_init
        _error_init(ImportError, self, name_value, **txt)


class _NamedEnum_RO(dict):
    '''(INTERNAL) C{Read_Only} enum-like C{dict} sub-class.
    '''
#   _name = NN  # also first kwd, __init__(_name=...)

    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:
            t = '%s %s' % (_dot_(self._name, attr), _doesn_t_exist_)  # PYCHOK _name
            raise AttributeError(t)

    def __setattr__(self, attr, value):
        t = 'Read_Only %s = %r' % (_dot_(self._name, attr), value)  # PYCHOK _name
        raise TypeError(t)

    def enums(self):
        for k, v in dict.items(self):
            if not k.startswith(_UNDERSCORE_):  # skip _name
                yield k, v


_ALL_INIT = 'pygeodesy_abspath', 'version'

# __all__ value for most modules, accessible as _ALL_LAZY.<module>
_ALL_LAZY = _NamedEnum_RO(_name='_ALL_LAZY',
                      azimuthal=('AzimuthalError', 'Equidistant', 'EquidistantKarney', 'Gnomonic', 'GnomonicKarney',
                                 'LambertEqualArea', 'Orthographic', 'Stereographic',
                                 'equidistant', 'gnomonic'),
                          bases=(),  # module and for backward compatibility only
                         basics=('EPS', 'EPS1', 'EPS1_2', 'EPS_2', 'INF', 'MANTIS', 'MAX', 'MIN',  # constants
                                 'NAN', 'NEG0', 'PI', 'PI2', 'PI_2', 'PI_4', 'R_M',
                                 'clips', 'halfs2',
                                 'isfinite', 'isinf', 'isint', 'isnan', 'isneg0', 'isscalar', 'issequence', 'isstr', 'issubclassof',
                                 'len2', 'map1', 'map2', 'property_doc_', 'property_RO'),
                          clipy=('ClipError',
                                 'clipCS3', 'clipSH', 'clipSH3'),
                            css=('CassiniSoldner', 'Css', 'CSSError', 'toCss'),
                          datum=('R_M', 'R_MA', 'R_MB', 'R_KM', 'R_NM', 'R_SM', 'R_FM', 'R_VM',
                                 'Datum',  'Ellipsoid',  'Transform',
                                 'Datums', 'Ellipsoids', 'Transforms'),
                     deprecated=('OK',  # DEPRECATED contants
                                 'HeightIDW', 'HeightIDW2', 'HeightIDW3', 'RefFrameError',  # DEPRECATED classes
                                 'anStr', 'areaof', 'bounds', 'clipDMS', 'clipStr', 'decodeEPSG2', 'encodeEPSG',  # most of the DEPRECATED functions, ...
                                 'equirectangular3', 'enStr2', 'false2f', 'falsed2f', 'fStr', 'fStrzs', 'hypot3',  # ... except ellipsoidal, spherical flavors
                                 'inStr', 'isenclosedby', 'nearestOn3', 'nearestOn4', 'parseUTM', 'perimeterof', 'polygon',
                                 'scalar', 'simplify2', 'toUtm', 'unStr', 'utmZoneBand2'),
                            dms=('F_D',   'F_DM',   'F_DMS',   'F_DEG',   'F_MIN',   'F_SEC',   'F__E',   'F__F',   'F__G',   'F_RAD',
                                 'F_D_',  'F_DM_',  'F_DMS_',  'F_DEG_',  'F_MIN_',  'F_SEC_',  'F__E_',  'F__F_',  'F__G_',  'F_RAD_',
                                 'F_D__', 'F_DM__', 'F_DMS__', 'F_DEG__', 'F_MIN__', 'F_SEC__', 'F__E__', 'F__F__', 'F__G__', 'F_RAD__',
                                 'S_DEG', 'S_MIN', 'S_SEC', 'S_RAD', 'S_SEP', 'ParseError',
                                 'bearingDMS', 'clipDegrees', 'clipRadians', 'compassDMS', 'compassPoint',
                                 'degDMS', 'latDMS', 'latlonDMS', 'lonDMS', 'normDMS',
                                 'parseDDDMMSS', 'parseDMS', 'parseDMS2', 'parse3llh', 'parseRad', 'precision', 'toDMS'),
                           ecef=('EcefCartesian', 'EcefError', 'EcefKarney', 'EcefMatrix', 'EcefSudano', 'EcefVeness', 'EcefYou'),
                     elevations=('elevation2', 'geoidHeight2'),
              ellipsoidalKarney=(),  # module only
             ellipsoidalNvector=(),  # module only
            ellipsoidalVincenty=('VincentyError',),  # nothing else
                       elliptic=('Elliptic', 'EllipticError'),
                           epsg=('Epsg', 'EPSGError'),
                         errors=('CrossError', 'IntersectionError', 'LenError', 'LimitError', 'PointsError',
                                 'RangeError', 'SciPyError', 'SciPyWarning',
                                 'crosserrors', 'exception_chaining', 'limiterrors', 'rangerrors'),
                            etm=('Etm', 'ETMError', 'ExactTransverseMercator',
                                 'parseETM5', 'toEtm8'),
                          fmath=('Fdot', 'Fhorner', 'Fpolynomial', 'Fsum',
                                 'cbrt', 'cbrt2',
                                 'favg', 'fdot', 'fdot3', 'fmean', 'fhorner', 'fidw', 'fpolynomial',
                                 'fpowers', 'fprod', 'frange', 'freduce', 'fsum', 'fsum_',
                                 'hypot', 'hypot_', 'hypot1', 'hypot2', 'sqrt3'),
                          formy=('antipode', 'antipode_', 'bearing', 'bearing_',
                                 'compassAngle', 'cosineForsytheAndoyerLambert', 'cosineForsytheAndoyerLambert_',
                                 'cosineAndoyerLambert', 'cosineAndoyerLambert_', 'cosineLaw', 'cosineLaw_',
                                 'euclidean', 'euclidean_', 'equirectangular', 'equirectangular_',
                                 'flatLocal', 'flatLocal_', 'flatPolar', 'flatPolar_',
                                 'haversine', 'haversine_', 'heightOf', 'horizon', 'hubeny', 'hubeny_',
                                 'intersections2', 'isantipode', 'isantipode_',
                                 'latlon2n_xyz', 'n_xyz2latlon', 'n_xyz2philam',
                                 'philam2n_xyz', 'points2', 'thomas', 'thomas_', 'vincentys', 'vincentys_'),
                        frechet=('Frechet', 'FrechetDegrees', 'FrechetError', 'FrechetRadians',
                                 'FrechetCosineAndoyerLambert', 'FrechetCosineForsytheAndoyerLambert',
                                 'FrechetCosineLaw', 'FrechetDistanceTo', 'FrechetEquirectangular',
                                 'FrechetEuclidean', 'FrechetFlatLocal', 'FrechetFlatPolar', 'FrechetHaversine',
                                 'FrechetHubeny', 'FrechetKarney', 'FrechetThomas', 'FrechetVincentys',
                                 'fractional', 'frechet_'),
                           gars=('Garef', 'GARSError'),
                        geohash=('Geohash', 'GeohashError'),
                         geoids=('GeoidError', 'GeoidG2012B', 'GeoidKarney', 'GeoidPGM', 'egmGeoidHeights', 'PGMError'),
                      hausdorff=('Hausdorff', 'HausdorffDegrees', 'HausdorffError', 'HausdorffRadians',
                                 'HausdorffCosineAndoyerLambert', 'HausdorffCosineForsytheAndoyerLambert',
                                 'HausdorffCosineLaw', 'HausdorffDistanceTo', 'HausdorffEquirectangular',
                                 'HausdorffEuclidean', 'HausdorffFlatLocal', 'HausdorffFlatPolar', 'HausdorffHaversine',
                                 'HausdorffHubeny', 'HausdorffKarney', 'HausdorffThomas', 'HausdorffVincentys',
                                 'hausdorff_', 'randomrangenerator'),
                        heights=('HeightError',
                                 'HeightIDWcosineAndoyerLambert', 'HeightIDWcosineForsytheAndoyerLambert',
                                 'HeightIDWcosineLaw', 'HeightIDWdistanceTo', 'HeightIDWequirectangular',
                                 'HeightIDWeuclidean', 'HeightIDWflatLocal', 'HeightIDWflatPolar', 'HeightIDWhaversine',
                                 'HeightIDWhubeny', 'HeightIDWkarney', 'HeightIDWthomas', 'HeightIDWvincentys',
                                 'HeightCubic', 'HeightLinear', 'HeightLSQBiSpline', 'HeightSmoothBiSpline'),
                        interns=('NN',),
                         karney=(),  # module only
                         lazily=('LazyImportError', 'isLazy'),
                            lcc=('Conic', 'Conics', 'Lcc', 'LCCError', 'toLcc'),
                           mgrs=('Mgrs', 'MGRSError', 'parseMGRS', 'toMgrs'),
                          named=('callername', 'classname', 'classnaming', 'modulename', 'nameof', 'notImplemented', 'notOverloaded'),
                        nvector=(),  # module and for backward compatibility only
                           osgr=('Osgr', 'OSGRError', 'parseOSGR', 'toOsgr'),
                         points=('LatLon_', 'LatLon2psxy', 'Numpy2LatLon', 'Tuple2LatLon',
                                 _areaOf_, 'boundsOf', 'centroidOf',
                                 _isclockwise_, 'isconvex', 'isconvex_', 'isenclosedBy', _ispolar_,
                                 'nearestOn5', _perimeterOf_),
               sphericalNvector=(),  # module only
          sphericalTrigonometry=(),  # module only
                       simplify=('simplify1', 'simplifyRDP', 'simplifyRDPm', 'simplifyRW', 'simplifyVW', 'simplifyVWm'),
                        streprs=('anstr', 'attrs', 'enstr2', 'fstr', 'fstrzs', 'hstr', 'instr', 'pairs', 'reprs', 'strs', 'unstr'),
                            trf=('RefFrame', 'RefFrames', 'TRFError', 'date2epoch', 'epoch2date'),
                          units=('Band', 'Bearing', 'Bearing_', 'Degrees', 'Distance', 'Easting',
                                 'Feet', 'Float', 'Float_', 'Height', 'Int', 'Int_',
                                 'Lam', 'Lam_', 'Lat', 'Lon', 'Meter', 'Northing', 'Number_',
                                 'Phi', 'Phi_', 'Precision_', 'Radians',
                                 'Radius', 'Radius_', 'Scalar', 'Scalar_', 'Str', 'UnitError', 'Zone'),
                            ups=('Ups', 'UPSError', 'parseUPS5', 'toUps8', 'upsZoneBand5'),
                          utily=('acos1', 'asin1', 'atan2d',
                                 'degrees', 'degrees90', 'degrees180', 'degrees360', 'degrees2m',
                                 'ft2m',
                                 'isNumpy2', 'isPoints2', 'isTuple2', 'iterNumpy2', 'iterNumpy2over',
                                 'm2degrees', 'm2ft', 'm2km', 'm2NM', 'm2SM',
                                 'radians', 'radiansPI', 'radiansPI2', 'radiansPI_2',
                                 'sincos2', 'sincos2d', 'splice', 'tan_2', 'tanPI_2_2',
                                 'unroll180', 'unrollPI',
                                 'wrap90', 'wrap180', 'wrap360', 'wrapPI_2','wrapPI', 'wrapPI2'),
                            utm=('Utm', 'UTMError', 'parseUTM5', 'toUtm8', 'utmZoneBand5'),
                         utmups=('UtmUps', 'UTMUPSError', 'parseUTMUPS5', 'toUtmUps8',
                                 'utmupsValidate', 'utmupsValidateOK', 'utmupsZoneBand5'),
                       vector3d=('Vector3d', 'VectorError'),
                    webmercator=('Wm', 'WebMercatorError', 'parseWM', 'toWm'),
                           wgrs=('Georef', 'WGRSError'))

# DEPRECATED __all__ names overloading those in _ALL_LAZY.deprecated where
# the new name is fully backward compatible in signature and return value
_ALL_OVERRIDING = _NamedEnum_RO(_name='_ALL_OVERRIDING',  # all DEPRECATED
                               basics=('clips as clipStr',),
                                fmath=('hypot_ as hypot3',),
                                formy=('points2 as polygon',),
                              heights=('HeightIDWequirectangular as HeightIDW2', 'HeightIDWeuclidean as HeightIDW',
                                       'HeightIDWhaversine as HeightIDW3'),
                               points=('areaOf as areaof',
                                       'isenclosedBy as isenclosedby', 'perimeterOf as perimeterof'),
                             simplify=('simplifyRW as simplify2',),
                              streprs=('anstr as anStr', 'enstr2 as enStr2', 'fstr as fStr', 'fstrzs as fStrzs',
                                       'instr as inStr', 'unstr as unStr'))

__all__ = _ALL_LAZY.lazily
__version__ = '20.08.04'


def _ALL_OTHER(*objs):
    '''(INTERNAL) List local objects for __all__.
    '''
    from pygeodesy import interns  # PYCHOK import

    def _dun(o):
        n = _dunder_name(o).rsplit(_DOT_, 1)[-1]
        u = _UNDERSCORE_ + n + _UNDERSCORE_
        return getattr(interns, u, n)

    return tuple(map(_dun, objs))


if _FOR_DOCS:
    _ALL_DOCS = _ALL_OTHER
    # (INTERNAL) Only export B{C{objs.__name__}} when making the
    # docs to force C{epydoc} to include certain classes, methods,
    # functions and other names in the documentation.  Using the
    # C{epydoc --private ...} command line option tends to include
    # too much internal documentation.
else:
    def _ALL_DOCS(*unused):
        return ()


def _all_imports(**more):
    '''(INTERNAL) Build C{dict} of all lazy imports.
    '''
    # imports naming conventions stored below - [<key>] = <from>:
    #  import <module>                        - [<module>] = <module>
    #  from <module> import <attr>            - [<attr>] = <module>
    #  from pygeodesy import <attr>           - [<attr>] = <attr>
    #  from <module> import <attr> as <name>  - [<name>] = <module>.<attr>
    imports = {}
    for _all_ in (_ALL_LAZY, _ALL_OVERRIDING, more):
        for mod, attrs in _all_.items():
            if isinstance(attrs, tuple) and not mod.startswith(_UNDERSCORE_):
                if mod not in imports:
                    imports[mod] = mod
                elif imports[mod] != mod:
                    t = _item_sq('imports', 'mod'), imports[mod], mod
                    raise AssertionError('%s: %r, not %r' % t)
                for attr in attrs:
                    attr, _, _as_ = attr.partition(' as ')
                    if _as_:
                        imports[_as_] = mod + _DOT_ + attr
                    else:
                        imports[attr] = mod
    return imports


def _all_missing2(_all_):
    '''(INTERNAL) Get deltas between pygeodesy.__all__ and lazily._all_imports.
    '''
    _alzy = _all_imports(**_NamedEnum_RO((a, ()) for a in _ALL_INIT))
    return ((_dot_('lazily', _all_imports.__name__), _COMMA_SPACE_.join(a for a in _all_ if a not in _alzy)),
            (_dot_('pygeodesy', '__all__'),          _COMMA_SPACE_.join(a for a in _alzy if a not in _all_)))


def _lazy_import2(_package_):  # MCCABE 16
    '''Check for and set up lazy importing.

       @arg _package_: The name of the package (C{str}) performing
                       the imports, to help facilitate resolving
                       relative imports, usually C{__package__}.

       @return: 2-Tuple C{(package, getattr)} of the importing package
                for easy reference within itself and the callable to
                be set to `__getattr__`.

       @raise LazyImportError: Lazy import not supported, an import
                               failed or a module name or attribute
                               name is invalid or does not exist.

       @note: This is the original function U{modutil.lazy_import
              <https://GitHub.com/brettcannon/modutil/blob/master/modutil.py>}
              modified to handle the C{__all__} and C{__dir__} attributes
              and call C{importlib.import_module(<module>.<name>, ...)}
              without causing a C{ModuleNotFoundError}.

       @see: The original U{modutil<https://PyPi.org/project/modutil>} and
             U{PEP 562<https://www.Python.org/dev/peps/pep-0562>}.
    '''
    import sys

    if sys.version_info[:2] < (3, 7):  # not supported
        t = 'no ' + _dot_(_package_, _lazy_import2.__name__)
        raise LazyImportError(t, txt='Python ' + sys.version.split()[0])

    import_module, package, parent = _lazy_init3(_package_)

    if isLazy > 2:  # trim import path names
        import os  # PYCHOK re-import
        cwdir = os.getcwd()
        cwdir = cwdir[:-len(os.path.basename(cwdir))]
        del os
    else:  # no import path names
        cwdir = NN

    import_ = _dot_(_package_, NN)  # namespace
    imports = _all_imports()

    def __getattr__(name):  # __getattr__ only for Python 3.7+
        # only called once for each undefined pygeodesy attribute
        if name in imports:
            # importlib.import_module() implicitly sets sub-modules
            # on this module as appropriate for direct imports (see
            # note in the _lazy_import.__doc__ above).
            mod, _, attr = imports[name].partition(_DOT_)
            if mod not in imports:
                raise LazyImportError('no module', txt=_dot_(parent, mod))
            imported = import_module(import_ + mod, parent)  # XXX '.' + mod
            if imported.__package__ not in (parent, '__main__', ''):
                raise LazyImportError(_dot_(mod, '__package__'), imported.__package__)
            # import the module or module attribute
            if attr:
                imported = getattr(imported, attr, _Missing)
            elif name != mod:
                imported = getattr(imported, name, _Missing)
            if imported is _Missing:
                raise LazyImportError('no attribute', txt=_dot_(mod, attr or name))

        elif name in ('__all__',):  # XXX '__dir__', '__members__'?
            imported = _ALL_INIT + tuple(imports.keys())
            mod = NN
        else:
            raise LazyImportError('no module or attribute', txt=_dot_(parent, name))

        setattr(package, name, imported)
        if isLazy > 1:
            z = NN
            if mod and mod != name:
                z = ' from .%s' % (mod,)
            if isLazy > 2:
                # sys._getframe(1) ... 'importlib._bootstrap' line 1032,
                # may throw a ValueError('call stack not deep enough')
                try:
                    f = sys._getframe(2)  # get importing file and line
                    n = f.f_code.co_filename
                    if cwdir and n.startswith(cwdir):
                        n = n[len(cwdir):]
                    z = '%s by %s line %d' % (z, n, f.f_lineno)
                except ValueError:  # PYCHOK no cover
                    pass
            print('# lazily imported %s%s' % (_dot_(parent, name), z))

        return imported  # __getattr__

    return package, __getattr__  # _lazy_import2


def _lazy_init3(_package_):
    '''(INTERNAL) Try to initialize lazy import.

       @arg _package_: The name of the package (C{str}) performing
                       the imports, to help facilitate resolving
                       relative imports, usually C{__package__}.

       @return: 3-Tuple C{(import_module, package, parent)} of module
                C{importlib.import_module}, the importing C{package}
                for easy reference within itself and the package name,
                aka the C{parent}.

       @raise LazyImportError: Lazy import not supported, an import
                               failed or a module name or attribute
                               name is invalid or does not exist.

       @note: Global C{isLazy} is set accordingly.
    '''
    global isLazy

    z = _environ.get('PYGEODESY_LAZY_IMPORT', None)
    if z is None:  # PYGEODESY_LAZY_IMPORT not set
        isLazy = 1  # ... but on by default on 3.7
    else:
        z = z.strip()  # like PYTHONVERBOSE et.al.
        isLazy = int(z) if z.isdigit() else (1 if z else 0)
    if isLazy < 1:  # not enabled
        raise LazyImportError('disabled', txt='%s=%r' % ('PYGEODESY_LAZY_IMPORT', z))
    if _environ.get('PYTHONVERBOSE', None):
        isLazy += 1

    try:  # to initialize
        from importlib import import_module

        package = import_module(_package_)
        parent = package.__spec__.parent  # __spec__ only in Python 3.7+
        if parent != _package_:  # assertion
            raise AttributeError('parent %r vs %r' % (parent, _package_))
    except (AttributeError, ImportError) as x:  # PYCHOK no cover
        isLazy = False  # failed
        raise LazyImportError(_lazy_init3.__name__, _package_, txt=str(x))

    return import_module, package, parent

# **) MIT License
#
# Copyright (C) 2018-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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
