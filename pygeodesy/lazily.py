
# -*- coding: utf-8 -*-

u'''Lazily import C{pygeodesy} modules and attributes, based on
U{lazy_import<https://modutil.ReadTheDocs.io/en/latest/#lazy_import>}
from I{Brett Cannon}'s U{modutil<https://PyPI.org/project/modutil>}.

C{Lazy import} is I{supported only for }U{Python 3.7+
<https://Snarky.CA/lazy-importing-in-python-3-7>} and is I{enabled by
default in }U{PyGeodesy 18.11.10+<https://PyPI.org/project/PyGeodesy>}
I{ and later}.

To I{enable} C{lazy import}, set C{env} variable C{PYGEODESY_LAZY_IMPORT}
to C{1}, C{2}, C{3} or higher prior to C{import pygeodesy}.  To I{disable}
C{lazy import}, set C{env} variable C{PYGEODESY_LAZY_IMPORT} to C{0} or
an empty string.  Use C{2} or higher to print a message for each lazily
imported module and attribute, similar to C{env} variable C{PYTHONVERBOSE}
showing imports.  Using C{3} or higher also shows the importing file name
and line number.

@note: C{Lazy import} applies only to top-level modules of C{pygeodesy}.
A C{lazy import} of a top-level module inherently loads all sub-modules
imported by that top-level module.
'''
from pygeodesy.interns import MISSING, NN, __all__ as _interns_a_l_l_, \
                             _areaOf_, _attribute_, _COMMASPACE_, \
                             _doesn_t_exist_, _DOT_, _dunder_name, \
                             _enabled_, _EQUALSPACED_, _immutable_, \
                             _isclockwise_, _ispolar_, _module_, \
                             _NL_, _no_, _not_, _or_, _perimeterOf_, \
                             _Python_, _pygeodesy_abspath_, _sep_, \
                             _SPACE_, _UNDER_, _version_

from os import environ as _env  # in .geodsolve
from os.path import basename as _basename
import sys as _sys  # in .props
try:
    from importlib import import_module
except ImportError:  # no import_module in Python 2.6-

    def import_module(name, package=None):  # PYCHOK
        raise LazyImportError(name=name, package=package,
                              txt=_no_(import_module.__name__))

_a_l_l_            = '__all__'
_FOR_DOCS          = _env.get('PYGEODESY_FOR_DOCS', NN)  # for epydoc ...
_imports_          = 'imports'
_p_a_c_k_a_g_e_    = '__package__'
_PYGEODESY_LAZY_IMPORT_  = 'PYGEODESY_LAZY_IMPORT'
_PYTHON_X_DEV      =  getattr(_sys, '_xoptions', {}).get('dev',  # Python 3.2
                     _env.get('PYTHONDEVMODE', NN))  # PYCHOK exported
_sub_packages      = 'deprecated' , 'geodesicx'
_sys_version_info2 = _sys.version_info[:2]  # in .fmath, .geodsolve

# @module_property[_RO?] <https://GitHub.com/jtushman/proxy_tools/>
isLazy = None  # see @var isLazy


class LazyImportError(ImportError):
    '''Raised if C{lazy import} is not supported, disabled or failed some other way.
    '''
    def __init__(self, *name_value, **txt):
        from pygeodesy.errors import _error_init
        _error_init(ImportError, self, name_value, **txt)


class _Dict(dict):
    '''(INTERNAL) Imports C{dict}.
    '''
    def add(self, key, value, *values):
        '''Add C{[key] = value}, typically C{[attr] = mod}.

           @raise AssertionError: The B{C{key}} already
                                  exists with different
                                  B{C{value}}.
        '''
        if key in self:
            val = self[key]  # duplicate OK
            if val != value and val not in values:  # PYCHOK no cover
                from pygeodesy.streprs import Fmt as _Fmt
                t = _Fmt.SQUARE(_imports_, key), val, value
                raise AssertionError('%s: %r, not %r' % t)
        else:
            self[key] = value


class _NamedEnum_RO(dict):
    '''(INTERNAL) C{Read_Only} enum-like C{dict} sub-class.
    '''
#   _name = NN  # also first kwd, __init__(_name=...)

    def _DOT_(self, attr):
        return _DOT_(self._name, attr)  # PYCHOK _name

    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:  # PYCHOK no cover
            from pygeodesy.errors import _AttributeError
            raise _AttributeError(self._DOT_(attr), txt=_doesn_t_exist_)

    def __setattr__(self, attr, value):  # PYCHOK no cover
        from pygeodesy.errors import _TypeError
        t = _EQUALSPACED_(self._DOT_(attr), value)
        raise _TypeError(t, txt=_immutable_)

    def enums(self):
        for k, v in dict.items(self):
            if not k.startswith(_UNDER_):  # skip _name
                yield k, v


_ALL_INIT = _pygeodesy_abspath_, _version_

# __all__ value for most modules, accessible as _ALL_LAZY.<module>
_ALL_LAZY = _NamedEnum_RO(_name='_ALL_LAZY',
                         albers=('AlbersEqualArea', 'AlbersEqualArea2', 'AlbersEqualArea4',
                                 'AlbersEqualAreaCylindrical', 'AlbersEqualAreaNorth', 'AlbersEqualAreaSouth',
                                 'AlbersError', 'Albers7Tuple'),
                      azimuthal=('AzimuthalError', 'Azimuthal7Tuple',
                                 'Equidistant', 'EquidistantExact', 'EquidistantGeodSolve', 'EquidistantKarney',
                                 'Gnomonic', 'GnomonicExact', 'GnomonicGeodSolve', 'GnomonicKarney',
                                 'LambertEqualArea', 'Orthographic', 'Stereographic',
                                 'equidistant', 'gnomonic'),
                         basics=('clips', 'copysign0', 'copytype', 'halfs2',
                                 'isbool', 'isclass', 'isfinite', 'isidentifier', 'isinf', 'isint', 'iskeyword',
                                 'isnan', 'isneg0', 'isodd', 'isscalar', 'issequence', 'isstr', 'issubclassof',
                                 'len2', 'map1', 'map2', 'neg', 'neg_',
                                 'signOf', 'splice', 'ub2str', 'unsign0'),
                          clipy=('ClipError',
                                 'ClipCS4Tuple', 'ClipLB6Tuple', 'ClipSH3Tuple',
                                 'clipCS4', 'clipLB6', 'clipSH', 'clipSH3'),
                            css=('CassiniSoldner', 'Css', 'CSSError', 'toCss',
                                 'EasNorAziRk4Tuple', 'LatLonAziRk4Tuple'),
                         datums=('Datum', 'Datums', 'Transform', 'Transforms'),
                     deprecated=('OK',  # DEPRECATED constants
                                 'bases', 'datum', 'nvector',  # DEPRECATED modules
                                 'ClipCS3Tuple', 'EcefCartesian', 'HeightIDW', 'HeightIDW2', 'HeightIDW3', 'RefFrameError', 'UtmUps4Tuple',  # DEPRECATED classes
                                 'anStr', 'areaof', 'bounds', 'clipCS3', 'clipDMS', 'clipStr', 'decodeEPSG2', 'encodeEPSG',  # most of the DEPRECATED functions, ...
                                 'equirectangular3', 'enStr2', 'false2f', 'falsed2f', 'fStr', 'fStrzs',  # ... except ellipsoidal, spherical flavors
                                 'hypot3', 'inStr', 'isenclosedby', 'joined', 'joined_',
                                 'nearestOn3', 'nearestOn4', 'parseUTM', 'perimeterof', 'polygon',
                                 'scalar', 'simplify2', 'toUtm', 'unStr', 'utmZoneBand2'),
                            dms=('F_D',   'F_DM',   'F_DMS',   'F_DEG',   'F_MIN',   'F_SEC',   'F__E',   'F__F',   'F__G',   'F_RAD',
                                 'F_D_',  'F_DM_',  'F_DMS_',  'F_DEG_',  'F_MIN_',  'F_SEC_',  'F__E_',  'F__F_',  'F__G_',  'F_RAD_',
                                 'F_D__', 'F_DM__', 'F_DMS__', 'F_DEG__', 'F_MIN__', 'F_SEC__', 'F__E__', 'F__F__', 'F__G__', 'F_RAD__',
                                 'S_DEG', 'S_MIN', 'S_SEC', 'S_RAD', 'S_SEP', 'ParseError',
                                 'bearingDMS', 'clipDegrees', 'clipRadians', 'compassDMS', 'compassPoint',
                                 'degDMS', 'latDMS', 'latlonDMS', 'lonDMS', 'normDMS',
                                 'parseDDDMMSS', 'parseDMS', 'parseDMS2', 'parse3llh', 'parseRad', 'precision', 'toDMS'),
                           ecef=('EcefError', 'EcefFarrell21', 'EcefFarrell22', 'EcefKarney', 'EcefMatrix', 'EcefSudano', 'Ecef9Tuple', 'EcefVeness', 'EcefYou'),
                     elevations=('elevation2', 'geoidHeight2',
                                 'Elevation2Tuple', 'GeoidHeight2Tuple'),
               ellipsoidalExact=(),  # module only
           ellipsoidalGeodSolve=(),  # module only
              ellipsoidalKarney=(),  # module only
             ellipsoidalNvector=('Ned3Tuple',),  # nothing else
            ellipsoidalVincenty=('VincentyError',),  # nothing else
                     ellipsoids=('R_M', 'R_MA', 'R_MB', 'R_KM', 'R_NM', 'R_SM', 'R_FM', 'R_GM', 'R_VM',
                                 'a_f2Tuple', 'Circle4Tuple', 'Curvature2Tuple',
                                 'Ellipsoid', 'Ellipsoid2', 'Ellipsoids',
                                 'a_b2e', 'a_b2e2', 'a_b2e22', 'a_b2e32', 'a_b2f', 'a_b2f_', 'a_b2f2', 'a_b2n',
                                 'a_f2b', 'a_f_2b', 'b_f2a', 'b_f_2a',
                                 'f2e2', 'f2e22', 'f2e32', 'f_2f', 'f2f_', 'f2f2', 'f2n', 'n2e2', 'n2f', 'n2f_'),
                       elliptic=('Elliptic', 'EllipticError', 'Elliptic3Tuple'),
                           epsg=('Epsg', 'EPSGError'),
                         errors=('CrossError', 'IntersectionError', 'NumPyError', 'LenError', 'LimitError', 'PointsError',
                                 'RangeError', 'SciPyError', 'SciPyWarning', 'TRFError', 'UnitError', 'VectorError',
                                 'crosserrors', 'exception_chaining', 'limiterrors', 'rangerrors'),
                            etm=('Etm', 'ETMError', 'ExactTransverseMercator',
                                 'EasNorExact4Tuple', 'LatLonExact4Tuple',
                                 'parseETM5', 'toEtm8'),
                          fmath=('Fdot', 'Fhorner', 'Fpolynomial', 'Fsum',
                                 'cbrt', 'cbrt2', 'euclid', 'euclid_',
                                 'facos1', 'fasin1', 'fatan', 'fatan1', 'fatan2', 'favg',
                                 'fdot', 'fdot3', 'fmean', 'fmean_', 'fhorner', 'fidw', 'fpolynomial',
                                 'fpowers', 'fprod', 'frange', 'freduce', 'fsum', 'fsum_',
                                 'hypot', 'hypot_', 'hypot1', 'hypot2', 'hypot2_',
                                 'norm2', 'norm_', 'sqrt0', 'sqrt3'),
                          formy=('antipode', 'antipode_', 'bearing', 'bearing_',
                                 'compassAngle', 'cosineForsytheAndoyerLambert', 'cosineForsytheAndoyerLambert_',
                                 'cosineAndoyerLambert', 'cosineAndoyerLambert_', 'cosineLaw', 'cosineLaw_',
                                 'equirectangular', 'equirectangular_', 'euclidean', 'euclidean_',
                                 'excessAbc', 'excessGirard', 'excessLHuilier',
                                 'excessKarney', 'excessKarney_', 'excessQuad', 'excessQuad_',
                                 'flatLocal', 'flatLocal_', 'flatPolar', 'flatPolar_',
                                 'haversine', 'haversine_', 'heightOf', 'horizon', 'hubeny', 'hubeny_',
                                 'intersections2', 'isantipode', 'isantipode_',
                                 'latlon2n_xyz', 'n_xyz2latlon', 'n_xyz2philam',
                                 'philam2n_xyz', 'radical2', 'thomas', 'thomas_', 'vincentys', 'vincentys_',
                                 'Radical2Tuple'),
                        frechet=('Frechet', 'FrechetDegrees', 'FrechetError', 'FrechetRadians',
                                 'FrechetCosineAndoyerLambert', 'FrechetCosineForsytheAndoyerLambert',
                                 'FrechetCosineLaw', 'FrechetDistanceTo', 'FrechetEquirectangular',
                                 'FrechetEuclidean', 'FrechetExact', 'FrechetFlatLocal', 'FrechetFlatPolar',
                                 'FrechetHaversine', 'FrechetHubeny', 'FrechetKarney', 'FrechetThomas',
                                 'FrechetVincentys', 'Frechet6Tuple',
                                 'frechet_'),
                           gars=('Garef', 'GARSError'),
                      geodesicx=('gx', 'gxarea', 'gxline',  # modules
                                 'Caps', 'GeodesicAreaExact', 'GeodesicExact', 'GeodesicLineExact', 'PolygonArea'),
                      geodsolve=('GeodesicSolve', 'GeodesicLineSolve'),
                        geohash=('Geohash', 'GeohashError', 'Neighbors8Dict', 'Resolutions2Tuple'),
                         geoids=('GeoidError', 'GeoidG2012B', 'GeoidKarney', 'GeoidPGM', 'egmGeoidHeights',
                                 'PGMError', 'GeoidHeight5Tuple'),
                      hausdorff=('Hausdorff', 'HausdorffDegrees', 'HausdorffError', 'HausdorffRadians',
                                 'HausdorffCosineAndoyerLambert', 'HausdorffCosineForsytheAndoyerLambert',
                                 'HausdorffCosineLaw', 'HausdorffDistanceTo', 'HausdorffEquirectangular',
                                 'HausdorffEuclidean', 'HausdorffExact', 'HausdorffFlatLocal', 'HausdorffFlatPolar',
                                 'HausdorffHaversine', 'HausdorffHubeny', 'HausdorffKarney', 'HausdorffThomas',
                                 'HausdorffVincentys', 'Hausdorff6Tuple',
                                 'hausdorff_', 'randomrangenerator'),
                        heights=('HeightError',
                                 'HeightIDWcosineAndoyerLambert', 'HeightIDWcosineForsytheAndoyerLambert',
                                 'HeightIDWcosineLaw', 'HeightIDWdistanceTo', 'HeightIDWequirectangular',
                                 'HeightIDWeuclidean', 'HeightIDWflatLocal', 'HeightIDWflatPolar', 'HeightIDWhaversine',
                                 'HeightIDWhubeny', 'HeightIDWkarney', 'HeightIDWthomas', 'HeightIDWvincentys',
                                 'HeightCubic', 'HeightLinear', 'HeightLSQBiSpline', 'HeightSmoothBiSpline'),
                        interns=_interns_a_l_l_,
                          iters=('LatLon2PsxyIter', 'PointsIter', 'points2',
                                 'isNumpy2', 'isPoints2', 'isTuple2', 'iterNumpy2', 'iterNumpy2over'),
                         karney=('Direct9Tuple', 'GDict', 'GeodesicError', 'GeodSolve12Tuple', 'Inverse10Tuple'),
                         lazily=('LazyImportError', 'isLazy', 'print_', 'printf'),
                            lcc=('Conic', 'Conics', 'Lcc', 'LCCError', 'toLcc'),
                            ltp=('Frustum', 'LocalCartesian', 'LocalError', 'Ltp'),
                      ltpTuples=('Aer', 'Aer4Tuple', 'Enu', 'Enu4Tuple', 'Footprint5Tuple', 'Local9Tuple',
                                 'Ned', 'Ned4Tuple', 'XyzLocal', 'Xyz4Tuple'),
                           mgrs=('Mgrs', 'MGRSError', 'parseMGRS', 'toMgrs', 'Mgrs4Tuple', 'Mgrs6Tuple'),
                          named=('callername', 'classname', 'classnaming', 'modulename', 'nameof', 'notImplemented', 'notOverloaded'),
                    namedTuples=('Bearing2Tuple', 'Bounds2Tuple', 'Bounds4Tuple',
                                 'Destination2Tuple', 'Destination3Tuple',
                                 'Distance2Tuple', 'Distance3Tuple', 'Distance4Tuple',
                                 'EasNor2Tuple', 'EasNor3Tuple', 'Intersection3Tuple',
                                 'LatLon2Tuple', 'LatLon3Tuple', 'LatLon4Tuple',
                                 'LatLonDatum3Tuple', 'LatLonDatum5Tuple',
                                 'LatLonPrec3Tuple', 'LatLonPrec5Tuple',
                                 'NearestOn3Tuple',
                                 'PhiLam2Tuple', 'PhiLam3Tuple', 'PhiLam4Tuple', 'Point3Tuple', 'Points2Tuple',
                                 'Triangle7Tuple', 'Triangle8Tuple', 'Trilaterate5Tuple',
                                 'UtmUps2Tuple', 'UtmUps5Tuple', 'UtmUps8Tuple', 'UtmUpsLatLon5Tuple',
                                 'Vector2Tuple', 'Vector3Tuple', 'Vector4Tuple'),
                           osgr=('Osgr', 'OSGRError', 'parseOSGR', 'toOsgr'),
                         points=('LatLon_', 'LatLon2psxy', 'NearestOn5Tuple', 'Numpy2LatLon', 'Shape2Tuple', 'Tuple2LatLon',
                                 _areaOf_, 'boundsOf', 'centroidOf', 'fractional',
                                 _isclockwise_, 'isconvex', 'isconvex_', 'isenclosedBy', _ispolar_,
                                 'luneOf', 'nearestOn5', _perimeterOf_, 'quadOf'),
                          props=('Property', 'Property_RO', 'property_RO', 'property_doc_',
                                 'deprecated_class', 'deprecated_function', 'deprecated_method',
                                 'deprecated_Property_RO', 'deprecated_property_RO', 'DeprecationWarnings'),
               sphericalNvector=(),  # module only
          sphericalTrigonometry=(),  # module only
                       simplify=('simplify1', 'simplifyRDP', 'simplifyRDPm', 'simplifyRW', 'simplifyVW', 'simplifyVWm'),
                        streprs=('anstr', 'attrs', 'enstr2', 'fstr', 'fstrzs', 'hstr', 'instr', 'pairs', 'reprs', 'strs', 'unstr'),
                            trf=('RefFrame', 'RefFrames', 'Transform7Tuple',
                                 'date2epoch', 'epoch2date', 'trfXform'),
                          units=('Band', 'Bearing', 'Bearing_', 'Bool',
                                 'Degrees', 'Degrees_', 'Degrees2', 'Distance', 'Distance_', 'Easting', 'Epoch',
                                 'Feet', 'FIx', 'Float', 'Float_', 'Height', 'Int', 'Int_',
                                 'Lam', 'Lam_', 'Lat', 'Lat_', 'Lon', 'Lon_',
                                 'Meter', 'Meter_', 'Meter2', 'Meter3', 'Northing', 'Number_',
                                 'Phi', 'Phi_', 'Precision_', 'Radians', 'Radians_', 'Radians2',
                                 'Radius', 'Radius_', 'Scalar', 'Scalar_', 'Str', 'Zone'),
                            ups=('Ups', 'UPSError', 'parseUPS5', 'toUps8', 'upsZoneBand5'),
                          utily=('acos1', 'acre2ha', 'acre2m2', 'asin1', 'atand', 'atan2b', 'atan2d',
                                 'chain2m', 'circle4',
                                 'degrees', 'degrees90', 'degrees180', 'degrees360', 'degrees2grades', 'degrees2m',
                                 'fathom2m', 'ft2m', 'furlong2m',
                                 'grades', 'grades400', 'grades2degrees', 'grades2radians',
                                 'm2degrees', 'm2ft', 'm2km', 'm2NM', 'm2radians', 'm2SM', 'm2yard',
                                 'radians', 'radiansPI', 'radiansPI2', 'radiansPI_2', 'radians2m',
                                 'sincos2', 'sincos2d', 'tan_2', 'tanPI_2_2',
                                 'unroll180', 'unrollPI',
                                 'wrap90', 'wrap180', 'wrap360', 'wrapPI_2','wrapPI', 'wrapPI2',
                                 'yard2m'),
                            utm=('Utm', 'UTMError', 'parseUTM5', 'toUtm8', 'utmZoneBand5'),
                         utmups=('UtmUps', 'UTMUPSError', 'parseUTMUPS5', 'toUtmUps8',
                                 'utmupsValidate', 'utmupsValidateOK', 'utmupsZoneBand5'),
                       vector3d=('Vector3d', 'intersection3d3', 'iscolinearWith', 'parse3d',
                                 'trilaterate2d2', 'trilaterate3d2'),
                    webmercator=('Wm', 'WebMercatorError', 'parseWM', 'toWm', 'EasNorRadius3Tuple'),
                           wgrs=('Georef', 'WGRSError'))

# DEPRECATED __all__ names overloading those in _ALL_LAZY.deprecated where
# the new name is fully backward compatible in signature and return value
_ALL_OVERRIDDEN = _NamedEnum_RO(_name='_ALL_OVERRIDING',  # all DEPRECATED
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
__version__ = '21.06.30'


def _ALL_OTHER(*objs):
    '''(INTERNAL) List local objects for __all__.
    '''
    from pygeodesy import interns  # PYCHOK import

    def _dun(o):
        n = _dunder_name(o).rsplit(_DOT_, 1)[-1]
        u = NN(_UNDER_, n, _UNDER_)
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
    imports = _Dict()
    imports_add = imports.add

    for ALL in (_ALL_LAZY, _ALL_OVERRIDDEN, more):
        for mod, attrs in ALL.items():
            if isinstance(attrs, tuple) and not mod.startswith(_UNDER_):
                imports_add(mod, mod)
                for attr in attrs:
                    attr, _, as_attr = attr.partition(' as ')
                    if as_attr:
                        imports_add(as_attr, _DOT_(mod, attr), *_sub_packages)
                    else:
                        imports_add(attr, mod)
    return imports


def _all_missing2(_all_):
    '''(INTERNAL) Get diffs between pygeodesy.__all__ and lazily._all_imports.
    '''
    _alzy = _all_imports(**_NamedEnum_RO((a, ()) for a in _ALL_INIT))
    return ((_DOT_('lazily', _all_imports.__name__), _COMMASPACE_.join(a for a in _all_ if a not in _alzy)),
            (_DOT_('pygeodesy', _a_l_l_),            _COMMASPACE_.join(a for a in _alzy if a not in _all_)))


def _caller3(up):  # in .named
    '''(INTERNAL) Get 3-tuple C{(caller name, file name, line number)}
       for the caller B{C{up}} stack frames in the Python call stack.
    '''
    # sys._getframe(1) ... 'importlib._bootstrap' line 1032,
    # may throw a ValueError('call stack not deep enough')
    f = _sys._getframe(up + 1)
    return (f.f_code.co_name,  # caller name
           _basename(f.f_code.co_filename),  # file name
            f.f_lineno)  # line number


def _lazy_import2(_pygeodesy_):  # MCCABE 15
    '''Check for and set up C{lazy import}.

       @arg _pygeodesy_: The name of the package (C{str}) performing
                         the imports, to help facilitate resolving
                         relative imports, usually C{__package__}.

       @return: 2-Tuple C{(package, getattr)} of the importing package
                for easy reference within itself and the callable to
                be set to `__getattr__`.

       @raise LazyImportError: Lazy import not supported or not enabled,
                               an import failed or the package name or
                               module name or attribute name is invalid
                               or does not exist.

       @note: This is the original function U{modutil.lazy_import
              <https://GitHub.com/brettcannon/modutil/blob/master/modutil.py>}
              modified to handle the C{__all__} and C{__dir__} attributes
              and call C{importlib.import_module(<module>.<name>, ...)}
              without causing a C{ModuleNotFoundError}.

       @see: The original U{modutil<https://PyPi.org/project/modutil>},
             U{PEP 562<https://www.Python.org/dev/peps/pep-0562>} and the
             U{new way<https://Snarky.CA/lazy-importing-in-python-3-7/>}.
    '''
    if _sys_version_info2 < (3, 7):  # not supported before 3.7
        t = _no_(_DOT_(_pygeodesy_, _lazy_import2.__name__))
        raise LazyImportError(t, txt=_Python_(_sys))

    package, parent = _lazy_init2(_pygeodesy_)

    packages = (parent, '__main__', NN) + tuple(
               _DOT_(parent, s) for s in _sub_packages)
    imports  = _all_imports()

    def __getattr__(name):  # __getattr__ only for Python 3.7+
        # only called once for each undefined pygeodesy attribute
        if name in imports:
            # importlib.import_module() implicitly sets sub-modules
            # on this module as appropriate for direct imports (see
            # note in the _lazy_import.__doc__ above).
            mod, _, attr = imports[name].partition(_DOT_)
            if mod not in imports:
                raise LazyImportError(_no_(_module_), txt=_DOT_(parent, mod))
            imported = import_module(_DOT_(_pygeodesy_, mod), parent)
            pkg = getattr(imported, _p_a_c_k_a_g_e_, None)
            if pkg not in packages:  # invalid package
                raise LazyImportError(_DOT_(mod, _p_a_c_k_a_g_e_), pkg)
            # import the module or module attribute
            if attr:
                imported = getattr(imported, attr, MISSING)
            elif name != mod:
                imported = getattr(imported, name, MISSING)
            if imported is MISSING:
                raise LazyImportError(_no_(_attribute_),
                                      txt=_DOT_(mod, attr or name))

        elif name in (_a_l_l_,):  # XXX '_d_i_r_', '_m_e_m_b_e_r_s_'?
            imported = _ALL_INIT + tuple(imports.keys())
            mod = NN
        else:
            raise LazyImportError(_no_(_module_, _or_, _attribute_),
                                  txt=_DOT_(parent, name))

        setattr(package, name, imported)
        if isLazy > 1:
            z = NN
            if mod and mod != name:
                z = ' from .%s' % (mod,)
            if isLazy > 2:
                try:  # see C{_caller3}
                    _, f, s = _caller3(2)
                    z = '%s by %s line %d' % (z, f, s)
                except ValueError:  # PYCHOK no cover
                    pass
            printf('# lazily imported %s%s', _DOT_(parent, name), z)

        return imported  # __getattr__

    return package, __getattr__  # _lazy_import2


def _lazy_init2(_pygeodesy_):
    '''(INTERNAL) Try to initialize lazy import.

       @arg _pygeodesy_: The name of the package (C{str}) performing
                         the imports, to help facilitate resolving
                         relative imports, usually C{__package__}.

       @return: 3-Tuple C{(import_module, package, parent)} of module
                C{importlib.import_module}, the importing C{package}
                for easy reference within itself, always C{pygeodesy}
                and the package name, aka the C{parent}, always
                C{'pygeodesy'}.

       @raise LazyImportError: Lazy import not supported or not enabled,
                               an import failed or the package name is
                               invalid or does not exist.

       @note: Global C{isLazy} is set accordingly.
    '''
    global isLazy

    z = _env.get(_PYGEODESY_LAZY_IMPORT_, None)
    if z is None:  # _PYGEODESY_LAZY_IMPORT_ not set
        isLazy = 1  # ... but on by default on 3.7
    else:
        z = z.strip()  # like PYTHONVERBOSE et.al.
        isLazy = int(z) if z.isdigit() else (1 if z else 0)
    if isLazy < 1:  # not enabled
        raise LazyImportError(_PYGEODESY_LAZY_IMPORT_, repr(z), txt=_not_(_enabled_))
    if _env.get('PYTHONVERBOSE', None):  # PYCHOK no cover
        isLazy += 1

    try:  # to initialize in Python 3+
        package = import_module(_pygeodesy_)
        parent = package.__spec__.parent  # __spec__ only in Python 3.7+
        if parent != _pygeodesy_:  # assert
            t = _COMMASPACE_(parent, _not_(_pygeodesy_))
            raise AttributeError(_EQUALSPACED_('parent', t))

    except (AttributeError, ImportError) as x:  # PYCHOK no cover
        isLazy = False  # failed
        raise LazyImportError(_lazy_init2.__name__, _pygeodesy_, txt=str(x))

    return package, parent


def print_(*args, **nl_nt_prefix_end_file_flush_sep):
    '''Python 3-style C{print} function.

       @arg args: Values to be converted to C{str} and
                  concatenated (C{any} types).
       @kwarg nl=0: Number of leading blank lines (C{int}).
       @kwarg nt=0: Number of additional , trailing blank lines (C{int}).
       @kwarg prefix=NN: To be inserted before the formatted text (C{str}).

       @note: Python 3+ keyword arguments C{end}, C{file} and C{flush}
              are silently ignored.
    '''
    sep = nl_nt_prefix_end_file_flush_sep.get(_sep_, _SPACE_)
    txt = sep.join(map(str, args))
    printf(txt, **nl_nt_prefix_end_file_flush_sep)


def printf(fmt, *args, **nl_nt_prefix_end_file_flush_sep):
    '''C-style C{printf} function.

       @arg fmt: C-style formating text (C{str}).
       @arg args: Values to be formatted (C{any} types).
       @kwarg nl=0: Number of leading blank lines (C{int}).
       @kwarg nt=0: Number of additional , trailing blank lines (C{int}).
       @kwarg prefix=NN: To be inserted before the formatted text (C{str}).

       @note: Python 3+ keyword arguments C{end}, C{file}, C{flush}
              and C{sep} are silently ignored.
    '''
    def _kwds(nl=0, nt=0, prefix=NN, **kwds):  # XXX end?
        nl = (_NL_ * nl) if nl else NN
        nt = (_NL_ * nt) if nt else NN
        return nl, nt, prefix, kwds

    nl, nt, prefix, _ = _kwds(**nl_nt_prefix_end_file_flush_sep)
    if args:
        fmt %= args
#   elif kwds:
#       fmt %= kwds
    print(NN(nl, prefix, fmt, nt))


# **) MIT License
#
# Copyright (C) 2018-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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
