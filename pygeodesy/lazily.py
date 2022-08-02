
# -*- coding: utf-8 -*-

u'''Lazily import C{pygeodesy} modules and attributes, based on
U{lazy_import<https://modutil.ReadTheDocs.io/en/latest/#lazy_import>}
from I{Brett Cannon}'s U{modutil<https://PyPI.org/project/modutil>}.

C{Lazy import} is I{supported only for }U{Python 3.7+
<https://Snarky.Ca/lazy-importing-in-python-3-7>} and is I{enabled by
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
                             _areaOf_, _attribute_, _by_, _COLONSPACE_, \
                             _COMMASPACE_, _doesn_t_exist_, _DOT_, \
                             _dunder_name, _enabled_, _EQUALSPACED_, \
                             _from_, _immutable_, _isclockwise_, \
                             _ispolar_, _module_, _NL_, _no_, _not_, _or_, \
                             _perimeterOf_, _Python_, _pygeodesy_abspath_, \
                             _SPACE_, _UNDER_, _version_

from os import getenv as _getenv  # in .errors, .geodsolve, .props, .units
from os.path import basename as _basename
import sys as _sys  # in .props

_a_l_l_                 = '__all__'
_FOR_DOCS               = _getenv('PYGEODESY_FOR_DOCS', NN)  # for epydoc ...
_from_DOT__             = _SPACE_(NN, _from_, _DOT_)
_imports_               = 'imports'
_lazily_                = 'lazily'
_lazily_imported__      = _SPACE_('#', _lazily_, 'imported', NN)
_line_                  = 'line'
_p_a_c_k_a_g_e_         = '__package__'
_pygeodesy_             = 'pygeodesy'
_PYGEODESY_LAZY_IMPORT_ = 'PYGEODESY_LAZY_IMPORT'
_PYTHON_X_DEV           =  getattr(_sys, '_xoptions', {}).get('dev',  # Python 3.2+
                          _getenv('PYTHONDEVMODE', NN))  # PYCHOK exported
_sub_packages           = 'deprecated', 'geodesicx'
_sys_version_info2      = _sys.version_info[:2]  # in .basics, .fmath, ...

# @module_property[_RO?] <https://GitHub.com/jtushman/proxy_tools/>
isLazy = None  # see @var isLazy in .__init__

try:
    from importlib import import_module
except ImportError:  # Python 2.6-

    def import_module(name, package=None):
        raise LazyImportError(name=name, package=package,
                              txt=_no_(import_module.__name__))


class LazyImportError(ImportError):
    '''Raised if C{lazy import} is not supported, disabled or failed some other way.
    '''
    def __init__(self, *name_value, **txt):
        _ALL_MODS.errors._error_init(ImportError, self, name_value, **txt)


class _Dict(dict):
    '''(INTERNAL) Imports C{dict}.
    '''
    def add(self, key, value, *values):
        '''Add C{[key] = value}, typically C{[attr] = mod}.

           @raise AssertionError: The B{C{key}} already exists
                                  with a different B{C{value}}.
        '''
        if key in self:
            val = self[key]  # duplicate OK
            if val != value and val not in values:  # PYCHOK no cover
                k = _ALL_MODS.streprs.Fmt.SQUARE(_imports_, key)
                t = _COLONSPACE_(k,       repr(val))
                t = _COMMASPACE_(t, _not_(repr(value)))
                raise AssertionError(t)
        else:
            self[key] = value


class _NamedEnum_RO(dict):
    '''(INTERNAL) C{Read_Only} enum-like C{dict} sub-class.
    '''
#   _name = NN  # also first kwd, __init__(_name=...)

    def _DOT_(self, attr):  # PYCHOK no cover
        return _DOT_(self._name, attr)  # PYCHOK _name

    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:
            t = self._DOT_(attr)
            raise AttributeError(_COLONSPACE_(t, _doesn_t_exist_))

    def __setattr__(self, attr, value):  # PYCHOK no cover
        t = _EQUALSPACED_(self._DOT_(attr), value)
        raise TypeError(_COLONSPACE_(t, _immutable_))

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
                                 'isbool', 'isclass', 'isclose', 'iscomplex',
                                 'isfinite', 'isidentifier', 'isinf', 'isint', 'isint0', 'iskeyword',
                                 'isnan', 'isnear0', 'isnear1', 'isneg0', 'isninf', 'isnon0', 'isodd',
                                 'isscalar', 'issequence', 'isstr', 'issubclassof', 'istuplist',
                                 'len2', 'map1', 'map2', 'neg', 'neg_', 'remainder',
                                 'signBit', 'signOf', 'splice', 'str2ub', 'ub2str', 'unsigned0'),
                   cartesianBase=(),  # module only
                           clipy=('ClipError',
                                 'ClipCS4Tuple', 'ClipLB6Tuple', 'ClipSH3Tuple',
                                 'clipCS4', 'clipLB6', 'clipSH', 'clipSH3'),
                            css=('CassiniSoldner', 'Css', 'CSSError', 'toCss',
                                 'EasNorAziRk4Tuple', 'EasNorAziRkEqu6Tuple', 'LatLonAziRk4Tuple'),
                         datums=('Datum', 'Datums', 'Transform', 'Transforms'),
                     deprecated=('EPS1_2', 'OK',  # DEPRECATED constants
                                 'bases', 'datum', 'nvector',  # DEPRECATED modules
                                 'ClipCS3Tuple', 'EcefCartesian', 'EasNorExact4Tuple', 'HeightIDW', 'HeightIDW2', 'HeightIDW3',  # DEPRECATED classes
                                 'LatLonExact4Tuple', 'Ned3Tuple', 'RefFrameError', 'Rhumb7Tuple', 'Transform7Tuple', 'UtmUps4Tuple',
                                 'anStr', 'areaof', 'bounds', 'clipCS3', 'clipDMS', 'clipStr', 'collins',   # most of the DEPRECATED functions, ...
                                 'decodeEPSG2', 'encodeEPSG', 'equirectangular3', 'enStr2',   # ... except ellipsoidal, spherical flavors
                                 'false2f', 'falsed2f', 'fStr', 'fStrzs', 'hypot3', 'inStr', 'isenclosedby', 'joined', 'joined_',
                                 'nearestOn3', 'nearestOn4', 'parseUTM', 'perimeterof', 'polygon',
                                 'scalar', 'simplify2', 'tienstra', 'toUtm', 'unsign0', 'unStr', 'utmZoneBand2'),
                            dms=('F_D',   'F_DM',   'F_DMS',   'F_DEG',   'F_MIN',   'F_SEC',   'F_D60',   'F__E',   'F__F',   'F__G',   'F_RAD',
                                 'F_D_',  'F_DM_',  'F_DMS_',  'F_DEG_',  'F_MIN_',  'F_SEC_',  'F_D60_',  'F__E_',  'F__F_',  'F__G_',  'F_RAD_',
                                 'F_D__', 'F_DM__', 'F_DMS__', 'F_DEG__', 'F_MIN__', 'F_SEC__', 'F_D60__', 'F__E__', 'F__F__', 'F__G__', 'F_RAD__',
                                 'S_DEG', 'S_MIN', 'S_SEC', 'S_DMS', 'S_RAD', 'S_SEP',
                                 'bearingDMS', 'clipDegrees', 'clipRadians', 'compassDMS', 'compassPoint',
                                 'degDMS', 'latDMS', 'latlonDMS', 'latlonDMS_', 'lonDMS', 'normDMS',
                                 'parseDDDMMSS', 'parseDMS', 'parseDMS2', 'parse3llh', 'parseRad', 'precision', 'toDMS'),
                           ecef=('EcefError', 'EcefFarrell21', 'EcefFarrell22', 'EcefKarney', 'EcefMatrix',
                                 'EcefSudano', 'Ecef9Tuple', 'EcefVeness', 'EcefYou'),
                     elevations=('Elevation2Tuple', 'GeoidHeight2Tuple',
                                 'elevation2', 'geoidHeight2'),
                ellipsoidalBase=(),  # module only
              ellipsoidalBaseDI=(),  # module only
               ellipsoidalExact=(),  # module only
           ellipsoidalGeodSolve=(),  # module only
              ellipsoidalKarney=(),  # module only
             ellipsoidalNvector=(),  # module only
            ellipsoidalVincenty=('VincentyError',),  # nothing else
                     ellipsoids=('R_M', 'R_MA', 'R_MB', 'R_KM', 'R_NM', 'R_SM', 'R_FM', 'R_GM', 'R_VM',
                                 'a_f2Tuple', 'Circle4Tuple', 'Curvature2Tuple',
                                 'Ellipsoid', 'Ellipsoid2', 'Ellipsoids',
                                 'a_b2e', 'a_b2e2', 'a_b2e22', 'a_b2e32', 'a_b2f', 'a_b2f_', 'a_b2f2', 'a_b2n',
                                 'a_f2b', 'a_f_2b', 'b_f2a', 'b_f_2a',
                                 'f2e2', 'f2e22', 'f2e32', 'f_2f', 'f2f_', 'f2f2', 'f2n', 'n2e2', 'n2f', 'n2f_'),
                       elliptic=('Elliptic', 'EllipticError', 'Elliptic3Tuple'),
                           epsg=('Epsg', 'EPSGError'),
                         errors=('CrossError', 'IntersectionError', 'NumPyError', 'LenError', 'LimitError',
                                 'MGRSError', 'ParseError', 'PointsError', 'RangeError', 'SciPyError', 'SciPyWarning',
                                 'TRFError', 'TriangleError', 'UnitError', 'VectorError',
                                 'crosserrors', 'exception_chaining', 'limiterrors', 'rangerrors'),
                            etm=('Etm', 'ETMError', 'ExactTransverseMercator',
                                 'parseETM5', 'toEtm8'),
                          fmath=('Fdot', 'Fhorner', 'Fpolynomial',
                                 'cbrt', 'cbrt2', 'euclid', 'euclid_',
                                 'facos1', 'fasin1', 'fatan', 'fatan1', 'fatan2', 'favg',
                                 'fdot', 'fdot3', 'fmean', 'fmean_', 'fhorner', 'fidw', 'fpolynomial',
                                 'fpowers', 'fprod', 'frange', 'freduce', 'fremainder',
                                 'hypot', 'hypot_', 'hypot1', 'hypot2', 'hypot2_',
                                 'norm2', 'norm_', 'sqrt0', 'sqrt3', 'sqrt_a'),
                          formy=('Radical2Tuple',
                                 'antipode', 'antipode_', 'bearing', 'bearing_',
                                 'compassAngle', 'cosineForsytheAndoyerLambert', 'cosineForsytheAndoyerLambert_',
                                 'cosineAndoyerLambert', 'cosineAndoyerLambert_', 'cosineLaw', 'cosineLaw_',
                                 'equirectangular', 'equirectangular_', 'euclidean', 'euclidean_',
                                 'excessAbc', 'excessGirard', 'excessLHuilier',
                                 'excessKarney', 'excessKarney_', 'excessQuad', 'excessQuad_',
                                 'flatLocal', 'flatLocal_', 'flatPolar', 'flatPolar_',
                                 'hartzell', 'haversine', 'haversine_', 'heightOf', 'horizon', 'hubeny', 'hubeny_',
                                 'intersections2', 'isantipode', 'isantipode_',
                                 'latlon2n_xyz', 'n_xyz2latlon', 'n_xyz2philam',
                                 'opposing', 'opposing_', 'philam2n_xyz',
                                 'radical2', 'thomas', 'thomas_', 'vincentys', 'vincentys_'),
                        frechet=('Frechet', 'FrechetDegrees', 'FrechetError', 'FrechetRadians',
                                 'FrechetCosineAndoyerLambert', 'FrechetCosineForsytheAndoyerLambert',
                                 'FrechetCosineLaw', 'FrechetDistanceTo', 'FrechetEquirectangular',
                                 'FrechetEuclidean', 'FrechetExact', 'FrechetFlatLocal', 'FrechetFlatPolar',
                                 'FrechetHaversine', 'FrechetHubeny', 'FrechetKarney', 'FrechetThomas',
                                 'FrechetVincentys', 'Frechet6Tuple',
                                 'frechet_'),
                         fstats=('Fcook', 'Flinear', 'Fwelford'),
                          fsums=('Fsum', 'Fsum2Tuple', 'ResidualError',
                                 'fsum', 'fsum_', 'fsum1', 'fsum1_',),
                           gars=('Garef', 'GARSError'),
                      geodesicx=('gx', 'gxarea', 'gxline',  # modules
                                 'GeodesicAreaExact', 'GeodesicExact', 'GeodesicLineExact', 'PolygonArea'),
                      geodsolve=('GeodesicSolve', 'GeodesicLineSolve', 'GeodSolve12Tuple'),
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
                         karney=('Area3Tuple', 'Caps', 'Direct9Tuple', 'GDict', 'GeodesicError', 'Inverse10Tuple'),
                            ktm=('KTMError', 'KTransverseMercator'),
                     latlonBase=(),  # module only
                         lazily=('LazyImportError', 'isLazy', 'print_', 'printf'),
                            lcc=('Conic', 'Conics', 'Lcc', 'LCCError', 'toLcc'),
                            ltp=('Attitude', 'AttitudeError', 'Frustum', 'LocalCartesian', 'LocalError', 'Ltp', 'tyr3d'),
                      ltpTuples=('Aer', 'Aer4Tuple', 'Attitude4Tuple', 'Enu', 'Enu4Tuple', 'Footprint5Tuple',
                                 'Local9Tuple', 'Ned', 'Ned4Tuple', 'XyzLocal', 'Xyz4Tuple'),
                           mgrs=('Mgrs', 'parseMGRS', 'toMgrs', 'Mgrs4Tuple', 'Mgrs6Tuple'),
                          named=('callername', 'classname', 'classnaming', 'modulename',
                                 'nameof', 'notImplemented', 'notOverloaded'),
                    namedTuples=('Bearing2Tuple', 'Bounds2Tuple', 'Bounds4Tuple',
                                 'Destination2Tuple', 'Destination3Tuple',
                                 'Distance2Tuple', 'Distance3Tuple', 'Distance4Tuple',
                                 'EasNor2Tuple', 'EasNor3Tuple', 'Forward4Tuple', 'Intersection3Tuple',
                                 'LatLon2Tuple', 'LatLon3Tuple', 'LatLon4Tuple',
                                 'LatLonDatum3Tuple', 'LatLonDatum5Tuple',
                                 'LatLonPrec3Tuple', 'LatLonPrec5Tuple',
                                 'NearestOn2Tuple', 'NearestOn3Tuple', 'NearestOn4Tuple',
                                 'NearestOn5Tuple', 'NearestOn6Tuple', 'NearestOn8Tuple',
                                 'PhiLam2Tuple', 'PhiLam3Tuple', 'PhiLam4Tuple', 'Point3Tuple', 'Points2Tuple',
                                 'Reverse4Tuple', 'Triangle7Tuple', 'Triangle8Tuple', 'Trilaterate5Tuple',
                                 'UtmUps2Tuple', 'UtmUps5Tuple', 'UtmUps8Tuple', 'UtmUpsLatLon5Tuple',
                                 'Vector2Tuple', 'Vector3Tuple', 'Vector4Tuple'),
                    nvectorBase=(),  # module only
                           osgr=('Osgr', 'OSGRError', 'parseOSGR', 'toOsgr'),
                         points=('LatLon_', 'LatLon2psxy', 'Numpy2LatLon', 'Shape2Tuple', 'Tuple2LatLon',
                                 _areaOf_, 'boundsOf', 'centroidOf', 'fractional',
                                 _isclockwise_, 'isconvex', 'isconvex_', 'isenclosedBy', _ispolar_,
                                 'luneOf', 'nearestOn5', _perimeterOf_, 'quadOf'),
                          props=('Property', 'Property_RO', 'property_RO', 'property_doc_',
                                 'deprecated_class', 'deprecated_function', 'deprecated_method',
                                 'deprecated_Property_RO', 'deprecated_property_RO', 'DeprecationWarnings'),
                     resections=('Collins5Tuple', 'ResectionError', 'Survey3Tuple', 'Tienstra7Tuple',
                                 'TriAngle4Tuple', 'TriSide2Tuple', 'TriSide4Tuple',
                                 'cassini', 'collins5', 'pierlot', 'tienstra7',
                                 'snellius3', 'wildberger3',
                                 'triAngle', 'triAngle4', 'triSide', 'triSide2', 'triSide4'),
                     rhumbsolve=('RhumbSolve', 'RhumbLineSolve', 'RhumbSolve7Tuple'),
                         rhumbx=('Rhumb', 'RhumbError', 'RhumbLine', 'RhumbOrder2Tuple', 'Rhumb8Tuple'),
                  sphericalBase=(),  # module only
               sphericalNvector=(),  # module only
          sphericalTrigonometry=(),  # module only
                       simplify=('simplify1', 'simplifyRDP', 'simplifyRDPm', 'simplifyRW', 'simplifyVW', 'simplifyVWm'),
                      solveBase=(),  # module only
                        streprs=('anstr', 'attrs', 'enstr2', 'fstr', 'fstrzs', 'hstr', 'instr', 'pairs', 'reprs', 'strs', 'unstr'),
                            trf=('Helmert7Tuple', 'RefFrame', 'RefFrames',
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
                                 'chain2m', 'circle4', 'cot', 'cot_', 'cotd', 'cotd_',
                                 'degrees', 'degrees90', 'degrees180', 'degrees360', 'degrees2grades', 'degrees2m',
#                                                                                    'degrees2grades as degrees2gons',
                                 'fathom2m', 'ft2m', 'furlong2m',
                                 'grades', 'grades400', 'grades2degrees', 'grades2radians',
#                                'grades as gons', 'grades400 as gons400', 'grades2degrees as gons2degrees', 'grades2radians as gons2radians',
                                 'm2chain', 'm2degrees', 'm2fathom', 'm2ft', 'm2furlong',
                                 'm2km', 'm2NM', 'm2radians', 'm2SM', 'm2toise', 'm2yard',
                                 'radians', 'radiansPI', 'radiansPI2', 'radiansPI_2', 'radians2m',
                                 'sincos2', 'sincos2_', 'sincos2d', 'sincos2d_',
                                 'tand', 'tand_', 'tan_2', 'tanPI_2_2', 'toise2m',
                                 'unroll180', 'unrollPI',
                                 'wrap90', 'wrap180', 'wrap360', 'wrapPI_2','wrapPI', 'wrapPI2',
                                 'yard2m'),
                            utm=('Utm', 'UTMError', 'parseUTM5', 'toUtm8', 'utmZoneBand5'),
                         utmups=('UtmUps', 'UTMUPSError', 'parseUTMUPS5', 'toUtmUps8',
                                 'utmupsValidate', 'utmupsValidateOK', 'utmupsZoneBand5'),
                     utmupsBase=(),  # module only
                       vector2d=('Circin6Tuple', 'Circum3Tuple', 'Circum4Tuple', 'Meeus2Tuple', 'Radii11Tuple', 'Soddy4Tuple',
                                 'circin6', 'circum3', 'circum4_', 'meeus2', 'radii11', 'soddy4'),
                       vector3d=('Vector3d', 'intersection3d3', 'iscolinearWith', 'nearestOn', 'nearestOn6', 'parse3d',
                                 'trilaterate2d2', 'trilaterate3d2'),
                   vector3dBase=(),  # module only
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


class _ALL_MODS(object):
    '''(INTERNAL) Memoize import of any L{pygeodesy} module.
    '''
    def _DOT_(self, name):  # PYCHOK no cover
        return _DOT_(self.__class__.__name__, name)

    def __getattr__(self, name):
        '''Get a C{pygeodesy} module or attribute by B{C{name}}.

           @arg name: Unqualified module or attribute name (C{str}).
        '''
        try:
            n = _DOT_(_pygeodesy_, name)
            return self.getmodule(n)
        except ImportError as x:
            raise LazyImportError(str(x), txt=_doesn_t_exist_)

    def __setattr__(self, name, value):  # PYCHOK no cover
        t = _EQUALSPACED_(self._DOT_(name), repr(value))
        raise AttributeError(_COLONSPACE_(t, _immutable_))

    def getattr(self, module_name, name, *dflt):
        '''Get an attribute of a C{pygeodesy} module.

           @arg module_name: Un- or qualified module name (C{str}).
           @arg name: Attribute name (C{str}).

           @return: The C{pygeodesy} module's attribute.

           @raise AttributeError: No attribute with that B{C{name}}.

           @raise ImportError: Importing B{C{module_name}} failed.
        '''
        n = module_name
        if n.split(_DOT_, 1)[0] != _pygeodesy_:
            n = _DOT_(_pygeodesy_, n)
        m = self.getmodule(n)
        return m if name in (None, NN) else (
               getattr(m, name, *dflt) if dflt else getattr(m, name))

    def getmodule(self, name):
        '''Get a C{pygeodesy} module.

           @arg name: Qualified module name (C{str}).

           @return: The C{pygeodesy} module.

           @raise ImportError: Importing module B{C{name}} failed.
        '''
        try:
            return _sys.modules[name]
        except KeyError:
            return import_module(name, _pygeodesy_)

#   def _imported(self, name, module):  # in _lazy_import2 below
#       try:
#           if name not in _sys.modules:
#               _sys.modules[name] = module
#       except (AttributeError, KeyError):
#           pass
#       return module

    def items(self):  # no module named 'items'
        '''Yield the modules imported so far.
        '''
        for n, m in _sys.modules.items():
            yield n, m

_ALL_MODS = _ALL_MODS()  # PYCHOK singleton

__all__ = _ALL_LAZY.lazily
__version__ = '22.08.02'


def _ALL_OTHER(*objs):
    '''(INTERNAL) Get class and function B{C{objs}} for __all__.
    '''
    _interns = _ALL_MODS.interns  # from pygeodesy import interns

    def _dun(o):
        n = _dunder_name(o).rsplit(_DOT_, 1)[-1]
        i =  NN(_UNDER_, n, _UNDER_)  # intern'd
        return getattr(_interns, i, n)

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
    _add = imports.add

    for ALL in (_ALL_LAZY, _ALL_OVERRIDDEN, more):
        for mod, attrs in ALL.items():
            if isinstance(attrs, tuple) and not mod.startswith(_UNDER_):
                _add(mod, mod)
                for attr in attrs:
                    attr, _, as_attr = attr.partition(' as ')
                    if as_attr:
                        _add(as_attr, _DOT_(mod, attr), *_sub_packages)
                    else:
                        _add(attr, mod)
    return imports


def _all_missing2(_all_):
    '''(INTERNAL) Get diffs between pygeodesy.__all__ and lazily._all_imports.
    '''
    def _diff(one, two):
        return _COMMASPACE_.join(a for a in one if a not in two)

    _alzy = _all_imports(**_NamedEnum_RO((a, ()) for a in _ALL_INIT))
    return ((_DOT_(_lazily_, _all_imports.__name__), _diff(_all_, _alzy)),
            (_DOT_(_pygeodesy_, _a_l_l_),            _diff(_alzy, _all_)))


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


def _lazy_import2(package_name):  # MCCABE 14
    '''Check for and set up C{lazy import}.

       @arg package_name: The name of the package (C{str}) performing
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
             U{new way<https://Snarky.Ca/lazy-importing-in-python-3-7/>}.
    '''
    if package_name != _pygeodesy_ or _sys_version_info2 < (3, 7):  # not supported before 3.7
        t = _no_(_DOT_(package_name, _lazy_import2.__name__))  # PYCHOK no cover
        raise LazyImportError(t, txt=_Python_(_sys.version))

    package, parent = _lazy_init2(_pygeodesy_)

    packages = (parent, '__main__', NN) + tuple(
               _DOT_(parent, s) for s in _sub_packages)
    imports  = _all_imports()

    def __getattr__(name):  # __getattr__ only for Python 3.7+
        # only called once for each undefined pygeodesy attribute
        if name in imports:
            # importlib.import_module() implicitly sets sub-modules
            # on this module as appropriate for direct imports (see
            # note in the _lazy_import2.__doc__ above).
            mod, _, attr = imports[name].partition(_DOT_)
            if mod not in imports:
                raise LazyImportError(_no_(_module_), txt=_DOT_(parent, mod))
            imported = import_module(_DOT_(_pygeodesy_, mod), parent)
            pkg = getattr(imported, _p_a_c_k_a_g_e_, None)
            if pkg not in packages:  # invalid package
                raise LazyImportError(_DOT_(mod, _p_a_c_k_a_g_e_), pkg)
            # _ALL_MODS._imported(mod, imported)
            # import the module or module attribute
            if attr:
                imported = getattr(imported, attr, MISSING)
            elif name != mod:
                imported = getattr(imported, name, MISSING)
            if imported is MISSING:  # PYCHOK no cover
                t = _DOT_(mod, attr or name)
                raise LazyImportError(_no_(_attribute_), txt=t)

        elif name in (_a_l_l_,):  # XXX '_d_i_r_', '_m_e_m_b_e_r_s_'?
            imported = _ALL_INIT + tuple(imports.keys())
            mod = NN
        else:  # PYCHOK no cover
            t = _no_(_module_, _or_, _attribute_)
            raise LazyImportError(t, txt=_DOT_(parent, name))

        setattr(package, name, imported)
        if isLazy > 1:
            t = NN(_lazily_imported__, _DOT_(parent, name))
            if mod and mod != name:
                t = NN(t, _from_DOT__, mod)
            if isLazy > 2:
                try:  # see C{_caller3}
                    _, f, s = _caller3(2)
                    t = _SPACE_(t, _by_, f, _line_, s)
                except ValueError:
                    pass
            print(t)  # XXX printf

        return imported  # __getattr__

    return package, __getattr__  # _lazy_import2


def _lazy_init2(package_name):
    '''(INTERNAL) Try to initialize lazy import.

       @arg package_name: The name of the package (C{str}) performing
                          the imports, to help facilitate resolving
                          relative imports, usually C{__package__}.

       @return: 2-Tuple C{(package, parent)} of the importing C{package}
                for easy reference within itself and its name aka the
                C{parent}, same as B{C{package_name}}.

       @raise LazyImportError: Lazy import not supported or not enabled,
                               an import failed or the package name is
                               invalid or does not exist.

       @note: Global C{isLazy} is set accordingly.
    '''
    global isLazy

    z = _getenv(_PYGEODESY_LAZY_IMPORT_, None)
    if z is None:  # _PYGEODESY_LAZY_IMPORT_ not set
        isLazy = 1  # ... but only by default on 3.7
    else:
        z = z.strip()  # like PYTHONVERBOSE et.al.
        isLazy = int(z) if z.isdigit() else (1 if z else 0)
    if isLazy < 1:  # not enabled
        raise LazyImportError(_PYGEODESY_LAZY_IMPORT_, repr(z), txt=_not_(_enabled_))
    if _getenv('PYTHONVERBOSE', None):  # PYCHOK no cover
        isLazy += 1

    try:  # to initialize in Python 3+
        package = import_module(package_name)
        parent = package.__spec__.parent  # __spec__ only in Python 3.7+
        if parent != package_name:  # assert
            t = _COMMASPACE_(parent, _not_(package_name))  # PYCHOK no cover
            raise AttributeError(_EQUALSPACED_('parent', t))

    except (AttributeError, ImportError) as x:
        isLazy = False  # failed
        raise LazyImportError(_lazy_init2.__name__, package_name, txt=str(x))

    return package, parent


def _pairs(*args, **kwds):  # in .errors, .ktm
    # from pygeodesy.streprs import pairs
    return _ALL_MODS.streprs.pairs(*args, **kwds)


def print_(*args, **nl_nt_prefix_end_file_flush_sep):  # PYCHOK no cover
    '''Python 3+ C{print}-like formatting and printing.

       @arg args: Arguments to be converted to C{str} and joined by B{C{sep}}
                  (C{any} type, all positional).
       @kwarg nl_nt_prefix_end_file_flush_sep: Keyword arguments C{B{nl}=0}
                 for the number of leading blank lines (C{int}), C{B{nt}=0}
                 the number of trailing blank lines (C{int}), C{B{prefix}=NN}
                 to be inserted before the formatted text (C{str}) and Python
                 3+ C{print} keyword arguments C{B{end}}, C{B{sep}}, C{B{file}}
                 and C{B{flush}}.

       @return: Number of bytes written.
    '''
    return printf(NN, *args, **nl_nt_prefix_end_file_flush_sep)


def printf(fmt, *args, **nl_nt_prefix_end_file_flush_sep_kwds):
    '''C{Printf-style} and Python 3+ C{print}-like formatting and printing.

       @arg fmt: U{Printf-style<https://Docs.Python.org/3/library/stdtypes.html#
                 printf-style-string-formatting>} format specification (C{str}).
       @arg args: Arguments to be formatted (C{any} types, all positional).
       @kwarg nl_nt_prefix_end_file_flush_sep_kwds: Keyword arguments C{B{nl}=0}
                 for the number of leading blank lines (C{int}), C{B{nt}=0} the
                 number of trailing blank lines (C{int}), C{B{prefix}=NN} to
                 be inserted before the formatted text (C{str}) and Python 3+
                 C{print} keyword arguments C{B{end}}, C{B{sep}}, C{B{file}} and
                 C{B{flush}}.  Any remaining C{B{kwds}} are U{printf-style
                 <https://Docs.Python.org/3/library/stdtypes.html#printf-style-string-formatting>}
                 keyword arguments to be formatted, I{iff no B{C{args}} are present}.

       @return: Number of bytes written.
    '''
    b, e, s, f, fl, p, kwds = _xprint6(**nl_nt_prefix_end_file_flush_sep_kwds)
    try:
        if args:
            t = (fmt % args) if fmt else s.join(map(str, args))
        elif kwds:  # PYCHOK no cover
            t = (fmt % kwds) if fmt else s.join(_pairs(kwds, prec=p))
        else:  # PYCHOK no cover
            t =  fmt
    except Exception as x:
        _E, t = _ALL_MODS.errors._xError2(x)
        unstr = _ALL_MODS.streprs.unstr
        raise _E(t, txt=unstr(printf.__name__, fmt, *args, **
                        nl_nt_prefix_end_file_flush_sep_kwds))
    n = f.write(NN(b, t, e))
    if fl:
        f.flush()
    return n


def _xprint6(nl=0, nt=0, prec=6, prefix=NN, sep=_SPACE_, file=_sys.stdout,
                                            end=_NL_, flush=False, **kwds):
    '''(INTERNAL) Unravel the C{printf} and remaining keyword arguments.
    '''
    if nl > 0:
        prefix = NN(_NL_ * nl, prefix)
    if nt > 0:
        end = NN(end, _NL_ * nt)
    return prefix, end, sep, file, flush, prec, kwds


if __name__ == '__main__':

    from timeit import timeit

    def t1():
        from pygeodesy.trf import RefFrame
        return RefFrame

    def t2():
        return _ALL_MODS.trf.RefFrame

    t1(); t2()  # PYCHOK prime each

    t1 = timeit(t1, number=1000000)
    t2 = timeit(t2, number=1000000)
    v = _Python_(_sys.version)
    printf('%.8f import vs %.8f _ALL_MODS: %.3fX, %s', t1, t2, t2 / t1, v)

# % python3.11 -W ignore -m pygeodesy.lazily
# 0.31956017 import vs 0.52837638 _ALL_MODS: 1.653X, Python 3.11.0b5

# % python3.10 -W ignore -m pygeodesy.lazily
# 0.31828208 import vs 0.58981700 _ALL_MODS: 1.853X, Python 3.10.5

# % python2 -m pygeodesy.lazily
# 1.19996715 import vs 1.39310884 _ALL_MODS: 1.161X, Python 2.7.18

# **) MIT License
#
# Copyright (C) 2018-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
