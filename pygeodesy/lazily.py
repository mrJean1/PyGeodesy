
# -*- coding: utf-8 -*-

u'''Lazily import C{pygeodesy} modules and attributes, based on
U{lazy_import<https://modutil.ReadTheDocs.io/en/latest/#lazy_import>}
from I{Brett Cannon}'s U{modutil<https://PyPI.org/project/modutil>}.

C{Lazy import} is I{supported only for }U{Python 3.7+
<https://Snarky.Ca/lazy-importing-in-python-3-7>} and is I{enabled by
default} in U{PyGeodesy 18.11.10<https://PyPI.org/project/PyGeodesy>}
I{and newer}.

To I{enable} C{lazy import}, set C{env} variable C{PYGEODESY_LAZY_IMPORT}
to C{1}, C{2}, C{3} or higher prior to C{import pygeodesy}.  To I{disable}
C{lazy import}, set C{env} variable C{PYGEODESY_LAZY_IMPORT} to C{0} or
an empty string.  Use C{2} or higher to print a message for each lazily
imported module and attribute, similar to C{env} variable C{PYTHONVERBOSE}
showing imports.  Using C{3} or higher also shows the importing file name
and line number.

@note: C{Lazy import} applies only to top-level modules of C{pygeodesy}.
       The C{lazy import} of a top-level module invariably loads all
       sub-modules imported by that top-level module.

@note: C{Lazy import} raises a L{LazyAttributeError} or L{LazyImportError}
       depending on the cause of the error and such errors can occur late,
       after all initial imports.
'''

# from pygeodesy.errors import _xError2  # _ALL_MODS
from pygeodesy.interns import MISSING, NN, __all__ as _interns__all__, _areaOf_, \
                             _attribute_, _by_, _COLONSPACE_, _COMMASPACE_, \
                             _doesn_t_exist_, _DOT_, _enabled_, _EQUALSPACED_, \
                             _from_, _HASH_, _immutable_, _isclockwise_, _ispolar_, \
                             _NL_, _no_, _NorthPole_, _not_, _or_, _pygeodesy_, \
                             _line_, _module_, _pygeodesy_abspath_, _Python_, _QUOTE1_, \
                             _QUOTE2_, _SouthPole_, _SPACE_, _sub_packages, _UNDER_, \
                             _version_, _dunder_nameof, _headof, _tailof  # _DEPRECATED_
from pygeodesy.interns import _intern  # PYCHOK used!
# from pygeodesy.streprs import Fmt, pairs, unstr  # _ALL_MODS
from pygeodesy import _isfrozen  # handle as w/o lazy import

from os import getenv as _getenv  # in .errors, .geodsolve, .props, .units
from os.path import basename as _basename
import sys as _sys  # in .basics._sizeof
try:
    from importlib import import_module
except ImportError:  # Python 2.6-
    def import_module(name, *package):
        t = _ALL_MODS.streprs.unstr(import_module, name, *package)
        raise LazyImportError(t, txt=_doesn_t_exist_)

_a_l_l_                 = '__all__'  # .__main__
__as__                  = ' as '
_FOR_DOCS               = _getenv('PYGEODESY_FOR_DOCS', NN)  # for epydoc ...
_from_DOT__             = _SPACE_(NN, _from_, _DOT_)
_i0                     = ()  # PYCHOK empty tuple
_init__all__            = _FOR_DOCS or _getenv('PYGEODESY_INIT__ALL__', _a_l_l_) == _a_l_l_  # PYCHOK expoted
_lazily_                = 'lazily'
_lazily_imported__      = _SPACE_(_HASH_, _lazily_, 'imported', NN)
_p_a_c_k_a_g_e_         = '__package__'
_PYGEODESY_GEOCONVERT_  = 'PYGEODESY_GEOCONVERT'  # PYCHOK .mgrs, test.bases
_PYGEODESY_GEODSOLVE_   = 'PYGEODESY_GEODSOLVE'   # PYCHOK .geodsolve, test.bases
_PYGEODESY_LAZY_IMPORT_ = 'PYGEODESY_LAZY_IMPORT'
_PYGEODESY_RHUMBSOLVE_  = 'PYGEODESY_RHUMBSOLVE'  # PYCHOK .rhumb.solve, test.bases
_PYTHON_X_DEV           =  getattr(_sys, '_xoptions', {}).get('dev',  # Python 3.2+
                          _getenv('PYTHONDEVMODE', NN))  # PYCHOK exported
_sys_version_info2      = _sys.version_info[:2]  # in .basics, .fmath, ...
_unlazy = _unLazy0      = _isfrozen or _sys_version_info2 < (3, 7)  # PYCHOK mod.__getattr__ 3.7+
_WARNINGS_X_DEV         = _getenv('PYGEODESY_WARNINGS', NN) and (
                          _PYTHON_X_DEV or bool(_sys.warnoptions))  # PYCHOK .props
# @module_property[_RO?] <https://GitHub.com/jtushman/proxy_tools/>
isLazy                  = None  # see @var isLazy in .__init__


class LazyAttributeError(AttributeError):
    '''Raised if a C{lazily imported} attribute is missing or invalid.
    '''
    def __init__(self, *name_value, **txt):
        _ALL_MODS.errors._error_init(AttributeError, self, name_value, **txt)


class LazyImportError(ImportError):
    '''Raised if C{lazy import} is not supported, disabled or failed some other way.
    '''
    def __init__(self, *name_value, **txt):
        _ALL_MODS.errors._error_init(ImportError, self, name_value, **txt)


class _Dict(dict):
    '''(INTERNAL) Imports C{dict}.
    '''
    _name = NN

    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:
            return dict.__getattr__(self, attr)

#   def __setattr__(self, attr, value):
#       if attr in self:
#           self[attr] = value
#       else:
#           dict.__setattr__(self, attr, value)

    def add(self, name, mod_, *subs):
        '''Add a C{[name] = mod_} item.

           @raise AssertionError: The B{C{name}} already exists
                                  with a different B{C{mod_}}.
        '''
        if name in self:
            sub = self[name]  # duplicate OK
            if sub != mod_ and sub not in subs:  # PYCHOK no cover
                t = _DOT_(self._name, name)
                t = _COLONSPACE_(t,       repr(sub))
                t = _COMMASPACE_(t, _not_(repr(mod_)))
                raise AssertionError(t)
        else:
            self[name] = mod_

    def _NAME(self, which):
        self._name = _intern(which.__name__.upper())


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
            raise LazyAttributeError(t, txt=_doesn_t_exist_)

    def __setattr__(self, attr, value):  # PYCHOK no cover
        t = _EQUALSPACED_(self._DOT_(attr), value)
        raise LazyAttributeError(t, txt=_immutable_)

    def enums(self):
        # Yield all C{(mod_, tuple)} pairs
        for m, t in dict.items(self):
            n = m.replace(_UNDER_, _DOT_)
            if n != m:
                if m.startswith(_UNDER_):
                    t = None  # skip _name= ...
                else:
                    u = m.rstrip(_UNDER_)
                    if u != m:
                        u = len(u)
                        n = n[:u] + m[u:]
            if isinstance(t, tuple):
                yield n, t

    def fill_D(self, _D, which):
        # Fill C{_Dict _D}.
        _D._NAME(which)
        _a = _D.add
        for m, t in self.enums():
            _a(m, _DOT_(m, NN, NN))  # import module
            for a in t:
                a, _, as_ = a.partition(__as__)
                if as_:  # import attr as attr_
                    _a(as_, _DOT_(m, a, NN), *_sub_packages)
                else:
                    _a(a, m)
        return _D


def _i(*names):
    '''(INTERNAL) Intern all C{names}.
    '''
    return tuple(map(_intern, names)) if names else _i0


def _ALL_ATTRS(*attrs):
    '''(INTERNAL) Unravel all exported module attributes.
    '''
    t = ()
    for attr in attrs:
        t += tuple(map(_attrof, attr))
    return t


_ALL_INIT = _i(_pygeodesy_abspath_, _version_)

# __all__ value for most modules, accessible as _ALL_LAZY.<module>
_ALL_LAZY = _NamedEnum_RO(_name='_ALL_LAZY',
                         albers=_i('AlbersEqualArea', 'AlbersEqualArea2', 'AlbersEqualArea4',
                                    'AlbersEqualAreaCylindrical', 'AlbersEqualAreaNorth', 'AlbersEqualAreaSouth',
                                    'AlbersError', 'Albers7Tuple'),
                       auxilats=_i(),  # module only
                      azimuthal=_i('AzimuthalError', 'Azimuthal7Tuple',
                                   'Equidistant', 'EquidistantExact', 'EquidistantGeodSolve', 'EquidistantKarney',
                                   'Gnomonic', 'GnomonicExact', 'GnomonicGeodSolve', 'GnomonicKarney',
                                   'LambertEqualArea', 'Orthographic', 'Stereographic',
                                   'equidistant', 'gnomonic'),
                         basics=_i('clips', 'copysign0', 'copytype', 'halfs2',
                                   'int1s', 'isbool', 'isCartesian', 'isclass', 'iscomplex', 'isDEPRECATED', 'isfloat',
                                   'isidentifier', 'isinstanceof', 'isint', 'iskeyword', 'isLatLon', 'islistuple',
                                   'isNvector', 'isodd', 'isscalar', 'issequence', 'isstr', 'issubclassof',
                                   'len2', 'map1', 'map2', 'neg', 'neg_',
                                   'signBit', 'signOf', 'splice', 'str2ub', 'ub2str', 'unsigned0'),
                       booleans=_i('BooleanFHP', 'BooleanGH', 'LatLonFHP', 'LatLonGH',
                                   'isBoolean'),
                  cartesianBase=_i('RadiusThetaPhi3Tuple', 'rtp2xyz', 'rtp2xyz_', 'xyz2rtp', 'xyz2rtp_'),
                          clipy=_i('ClipCS4Tuple', 'ClipFHP4Tuple', 'ClipGH4Tuple', 'ClipLB6Tuple', 'ClipSH3Tuple',
                                   'clipCS4', 'clipFHP4', 'clipGH4', 'clipLB6', 'clipSH', 'clipSH3'),
                            css=_i('CassiniSoldner', 'Css', 'CSSError', 'toCss',
                                   'EasNorAziRk4Tuple', 'EasNorAziRkEqu6Tuple', 'LatLonAziRk4Tuple'),
                      constants=_i('DIG', 'EPS', 'EPS0', 'EPS02', 'EPS1', 'EPS2', 'EPS4', 'EPS_2',
                                   'INF', 'INT0', 'MANT_DIG', 'MAX', 'MAX_EXP', 'MIN', 'MIN_EXP', 'NAN', 'NEG0', 'NINF',
                                   'PI', 'PI2', 'PI_2', 'PI3', 'PI_3', 'PI3_2', 'PI4', 'PI_4',
                                   'R_FM', 'R_GM', 'R_KM', 'R_M', 'R_MA', 'R_MB', 'R_NM', 'R_QM', 'R_SM', 'R_VM',
                                   'float_', 'float0_', 'isclose', 'isfinite', 'isinf', 'isint0',
                                   'isnan', 'isnear0', 'isnear1', 'isnear90', 'isneg0', 'isninf', 'isnon0',
                                   'remainder'),
                         datums=_i('Datum', 'Datums', 'Transform', 'Transforms'),
#                    deprecated=_i(),  # module only
                            dms=_i('F_D',   'F_DM',   'F_DMS',   'F_DEG',   'F_MIN',   'F_SEC',   'F_D60',   'F__E',   'F__F',   'F__G',   'F_RAD',
                                   'F_D_',  'F_DM_',  'F_DMS_',  'F_DEG_',  'F_MIN_',  'F_SEC_',  'F_D60_',  'F__E_',  'F__F_',  'F__G_',  'F_RAD_',
                                   'F_D__', 'F_DM__', 'F_DMS__', 'F_DEG__', 'F_MIN__', 'F_SEC__', 'F_D60__', 'F__E__', 'F__F__', 'F__G__', 'F_RAD__',
                                   'S_DEG', 'S_MIN', 'S_SEC', 'S_DMS', 'S_RAD', 'S_SEP',
                                   'bearingDMS', 'clipDegrees', 'clipRadians', 'compassDMS', 'compassPoint',
                                   'degDMS', 'latDMS', 'latlonDMS', 'latlonDMS_', 'lonDMS', 'normDMS',
                                   'parseDDDMMSS', 'parseDMS', 'parseDMS2', 'parse3llh', 'parseRad', 'precision', 'toDMS'),
                           ecef=_i('EcefError', 'EcefFarrell21', 'EcefFarrell22', 'EcefKarney', 'EcefMatrix',
                                   'EcefSudano', 'Ecef9Tuple', 'EcefVeness', 'EcefYou'),
                     elevations=_i('Elevation2Tuple', 'GeoidHeight2Tuple',
                                   'elevation2', 'geoidHeight2'),
                ellipsoidalBase=_i(),  # module only
              ellipsoidalBaseDI=_i(),  # module only
               ellipsoidalExact=_i(),  # module only
           ellipsoidalGeodSolve=_i(),  # module only
              ellipsoidalKarney=_i(),  # module only
             ellipsoidalNvector=_i(),  # module only
            ellipsoidalVincenty=_i('VincentyError',),  # nothing else
                     ellipsoids=_i('a_f2Tuple', 'Circle4Tuple', 'Curvature2Tuple',
                                   'Ellipsoid', 'Ellipsoid2', 'Ellipsoids',
                                   'a_b2e', 'a_b2e2', 'a_b2e22', 'a_b2e32', 'a_b2f', 'a_b2f_', 'a_b2f2', 'a_b2n',
                                   'a_f2b', 'a_f_2b', 'b_f2a', 'b_f_2a',
                                   'e2f', 'e22f',
                                   'f2e2', 'f2e22', 'f2e32', 'f_2f', 'f2f_', 'f2f2', 'f2n', 'n2e2', 'n2f', 'n2f_'),
                       elliptic=_i('Elliptic', 'EllipticError', 'Elliptic3Tuple'),
                           epsg=_i('Epsg', 'EPSGError'),
                         errors=_i('AuxError', 'ClipError', 'CrossError', 'GeodesicError', 'IntersectionError',
                                   'NumPyError', 'LenError', 'LimitError', 'MGRSError',
                                   'ParseError', 'PointsError', 'RangeError', 'RhumbError',
                                   'SciPyError', 'SciPyWarning', 'TRFError', 'TriangleError', 'UnitError', 'VectorError',
                                   'crosserrors', 'exception_chaining', 'isError', 'itemsorted',
                                   'limiterrors', 'rangerrors'),
                            etm=_i('Etm', 'ETMError', 'ExactTransverseMercator',
                                   'parseETM5', 'toEtm8'),
                          fmath=_i('Fdot', 'Fhorner', 'Fhypot', 'Fpolynomial', 'Fpowers', 'Fn_rt', 'Fcbrt', 'Fsqrt',
                                   'bqrt', 'cbrt', 'cbrt2', 'euclid', 'euclid_',
                                   'facos1', 'fasin1', 'fatan', 'fatan1', 'fatan2', 'favg',
                                   'fdot', 'fdot3', 'fmean', 'fmean_', 'fhorner', 'fidw', 'fpolynomial',
                                   'fpowers', 'fprod', 'frange', 'freduce', 'fremainder',
                                   'hypot', 'hypot_', 'hypot1', 'hypot2', 'hypot2_',
                                   'norm2', 'norm_', 'sqrt0', 'sqrt3', 'sqrt_a', 'zcrt', 'zqrt'),
                          formy=_i('Radical2Tuple',
                                   'antipode', 'antipode_', 'bearing', 'bearing_',
                                   'compassAngle', 'cosineForsytheAndoyerLambert', 'cosineForsytheAndoyerLambert_',
                                   'cosineAndoyerLambert', 'cosineAndoyerLambert_', 'cosineLaw', 'cosineLaw_',
                                   'equirectangular', 'equirectangular_', 'euclidean', 'euclidean_',
                                   'excessAbc_', 'excessCagnoli_', 'excessGirard_', 'excessLHuilier_',
                                   'excessKarney', 'excessKarney_', 'excessQuad', 'excessQuad_',
                                   'flatLocal', 'flatLocal_', 'flatPolar', 'flatPolar_',
                                   'hartzell', 'haversine', 'haversine_', 'heightOf', 'heightOrthometric', 'horizon', 'hubeny', 'hubeny_',
                                   'intersection2', 'intersections2', 'isantipode', 'isantipode_', 'isnormal', 'isnormal_',
                                   'latlon2n_xyz', 'normal', 'normal_', 'n_xyz2latlon', 'n_xyz2philam',
                                   'opposing', 'opposing_', 'philam2n_xyz', 'radical2',
                                   'thomas', 'thomas_', 'vincentys', 'vincentys_'),
                        frechet=_i('Frechet', 'FrechetDegrees', 'FrechetError', 'FrechetRadians',
                                   'FrechetCosineAndoyerLambert', 'FrechetCosineForsytheAndoyerLambert',
                                   'FrechetCosineLaw', 'FrechetDistanceTo', 'FrechetEquirectangular',
                                   'FrechetEuclidean', 'FrechetExact', 'FrechetFlatLocal', 'FrechetFlatPolar',
                                   'FrechetHaversine', 'FrechetHubeny', 'FrechetKarney', 'FrechetThomas',
                                   'FrechetVincentys', 'Frechet6Tuple',
                                   'frechet_'),
                         fstats=_i('Fcook', 'Flinear', 'Fwelford'),
                          fsums=_i('Fsum', 'Fsum2Tuple', 'ResidualError',
                                   'fsum', 'fsum_', 'fsumf_', 'fsum1', 'fsum1_', 'fsum1f_'),
                           gars=_i('Garef', 'GARSError'),
                      geodesicw=_i('Geodesic', 'GeodesicLine', 'Geodesic_WGS84'),
                      geodesicx=_i('gx', 'gxarea', 'gxbases', 'gxline',  # modules, see _sub_packages
                                   'GeodesicAreaExact', 'GeodesicExact', 'GeodesicLineExact', 'PolygonArea'),
                      geodsolve=_i('GeodesicSolve', 'GeodesicLineSolve', 'GeodSolve12Tuple'),
                        geohash=_i('Geohash', 'GeohashError', 'Neighbors8Dict', 'Resolutions2Tuple'),
                         geoids=_i('GeoidError', 'GeoidG2012B', 'GeoidKarney', 'GeoidPGM', 'egmGeoidHeights',
                                   'PGMError', 'GeoidHeight5Tuple'),
                      hausdorff=_i('Hausdorff', 'HausdorffDegrees', 'HausdorffError', 'HausdorffRadians',
                                   'HausdorffCosineAndoyerLambert', 'HausdorffCosineForsytheAndoyerLambert',
                                   'HausdorffCosineLaw', 'HausdorffDistanceTo', 'HausdorffEquirectangular',
                                   'HausdorffEuclidean', 'HausdorffExact', 'HausdorffFlatLocal', 'HausdorffFlatPolar',
                                   'HausdorffHaversine', 'HausdorffHubeny', 'HausdorffKarney', 'HausdorffThomas',
                                   'HausdorffVincentys', 'Hausdorff6Tuple',
                                   'hausdorff_', 'randomrangenerator'),
                        heights=_i('HeightCubic', 'HeightError',
                                   'HeightIDWcosineAndoyerLambert', 'HeightIDWcosineForsytheAndoyerLambert',
                                   'HeightIDWcosineLaw', 'HeightIDWdistanceTo', 'HeightIDWequirectangular',
                                   'HeightIDWeuclidean', 'HeightIDWexact', 'HeightIDWflatLocal', 'HeightIDWflatPolar',
                                   'HeightIDWhaversine', 'HeightIDWhubeny', 'HeightIDWkarney', 'HeightIDWthomas',
                                   'HeightIDWvincentys', 'HeightLinear', 'HeightLSQBiSpline', 'HeightSmoothBiSpline'),
                        interns=_interns__all__,
                          iters=_i('LatLon2PsxyIter', 'PointsIter', 'points2',
                                   'isNumpy2', 'isPoints2', 'isTuple2', 'iterNumpy2', 'iterNumpy2over'),
                         karney=_i('Area3Tuple', 'Caps', 'Direct9Tuple', 'GDict', 'Inverse10Tuple', 'Rhumb8Tuple'),
                            ktm=_i('KTMError', 'KTransverseMercator'),
                     latlonBase=_i(),  # module only
                         lazily=_i('LazyAttributeError', 'LazyImportError', 'isLazy', 'print_', 'printf'),
                            lcc=_i('Conic', 'Conics', 'Lcc', 'LCCError', 'toLcc'),
                            ltp=_i('Attitude', 'AttitudeError', 'ChLV', 'ChLVa', 'ChLVe', 'Frustum',
                                   'LocalCartesian', 'LocalError', 'Ltp', 'tyr3d'),
                      ltpTuples=_i('Aer', 'Aer4Tuple', 'Attitude4Tuple',
                                   'ChLVEN2Tuple', 'ChLV9Tuple', 'ChLVYX2Tuple', 'ChLVyx2Tuple',
                                   'Enu', 'Enu4Tuple', 'Footprint5Tuple', 'Local9Tuple', 'Los',
                                   'Ned', 'Ned4Tuple', 'Uvw', 'Uvw3Tuple', 'XyzLocal', 'Xyz4Tuple'),
                           mgrs=_i('Mgrs', 'parseMGRS', 'toMgrs', 'Mgrs4Tuple', 'Mgrs6Tuple'),
                          named=_i('ADict',
                                   'callername', 'classname', 'classnaming', 'modulename',
                                   'nameof', 'notImplemented', 'notOverloaded'),
                    namedTuples=_i('Bearing2Tuple', 'Bounds2Tuple', 'Bounds4Tuple',
                                   'Destination2Tuple', 'Destination3Tuple',
                                   'Distance2Tuple', 'Distance3Tuple', 'Distance4Tuple',
                                   'EasNor2Tuple', 'EasNor3Tuple', 'Forward4Tuple', 'Intersection3Tuple',
                                   'LatLon2Tuple', 'LatLon3Tuple', 'LatLon4Tuple',
                                   'LatLonDatum3Tuple', 'LatLonDatum5Tuple',
                                   'LatLonPrec3Tuple', 'LatLonPrec5Tuple',
                                   'NearestOn2Tuple', 'NearestOn3Tuple', 'NearestOn6Tuple', 'NearestOn8Tuple',
                                   'PhiLam2Tuple', 'PhiLam3Tuple', 'PhiLam4Tuple', 'Point3Tuple', 'Points2Tuple',
                                   'Reverse4Tuple', 'Triangle7Tuple', 'Triangle8Tuple', 'Trilaterate5Tuple',
                                   'UtmUps2Tuple', 'UtmUps5Tuple', 'UtmUps8Tuple', 'UtmUpsLatLon5Tuple',
                                   'Vector2Tuple', 'Vector3Tuple', 'Vector4Tuple'),
                    nvectorBase=_i(_NorthPole_, _SouthPole_),
                           osgr=_i('Osgr', 'OSGRError', 'parseOSGR', 'toOsgr'),
                         points=_i('LatLon_', 'LatLon2psxy', 'Numpy2LatLon', 'Shape2Tuple', 'Tuple2LatLon',
                                   _areaOf_, 'boundsOf', 'centroidOf', 'fractional',
                                   _isclockwise_, 'isconvex', 'isconvex_', 'isenclosedBy', _ispolar_,
                                   'luneOf', 'nearestOn5', 'perimeterOf', 'quadOf'),
                          props=_i('Property', 'Property_RO', 'property_RO', 'property_doc_',
                                   'deprecated_class', 'deprecated_function', 'deprecated_method',
                                   'deprecated_Property_RO', 'deprecated_property_RO', 'DeprecationWarnings'),
                     resections=_i('Collins5Tuple', 'ResectionError', 'Survey3Tuple', 'Tienstra7Tuple',
                                   'TriAngle5Tuple', 'TriSide2Tuple', 'TriSide4Tuple',
                                   'cassini', 'collins5', 'pierlot', 'pierlotx', 'tienstra7',
                                   'snellius3', 'wildberger3',
                                   'triAngle', 'triAngle5', 'triArea', 'triSide', 'triSide2', 'triSide4'),
                          rhumb=_i(),  # module only
                     rhumb_aux_=_i('RhumbAux', 'RhumbLineAux'),
                      rhumb_ekx=_i('Rhumb', 'RhumbLine'),
                    rhumb_solve=_i('RhumbSolve', 'RhumbLineSolve', 'RhumbSolve7Tuple'),
                  sphericalBase=_i(),  # module only
               sphericalNvector=_i(),  # module only
          sphericalTrigonometry=_i(),  # module only
                       simplify=_i('simplify1', 'simplifyRDP', 'simplifyRDPm', 'simplifyRW', 'simplifyVW', 'simplifyVWm'),
                      solveBase=_i(),  # module only
                        streprs=_i('anstr', 'attrs', 'enstr2', 'fstr', 'fstrzs', 'hstr', 'instr',
                                   'lrstrip', 'pairs', 'reprs', 'strs', 'unstr'),
                            trf=_i('RefFrame', 'RefFrames', 'TransformXform', 'TRFXform', 'TRFXform7Tuple',
                                   'date2epoch', 'epoch2date', 'trfTransform0', 'trfXform'),
                      triaxials=_i('BetaOmega2Tuple', 'BetaOmega3Tuple', 'Jacobi2Tuple',
                                   'JacobiConformal', 'JacobiConformalSpherical',
                                   'Triaxial', 'Triaxial_', 'TriaxialError', 'Triaxials', 'hartzell4'),
                          units=_i('Band', 'Bearing', 'Bearing_', 'Bool',
                                   'Degrees', 'Degrees_', 'Degrees2', 'Distance', 'Distance_', 'Easting', 'Epoch',
                                   'Feet', 'FIx', 'Float_', 'Height', 'Height_', 'HeightX', 'Int_',
                                   'Lam', 'Lam_', 'Lat', 'Lat_', 'Lon', 'Lon_',
                                   'Meter', 'Meter_', 'Meter2', 'Meter3', 'Northing', 'Number_',
                                   'Phi', 'Phi_', 'Precision_', 'Radians', 'Radians_', 'Radians2',
                                   'Radius_', 'Scalar', 'Scalar_', 'Zone'),
                      unitsBase=_i('Float', 'Int', 'Radius', 'Str'),
                            ups=_i('Ups', 'UPSError', 'parseUPS5', 'toUps8', 'upsZoneBand5'),
                          utily=_i('acos1', 'acre2ha', 'acre2m2', 'asin1', 'atan1', 'atan1d', 'atan2b', 'atan2d',
                                   'chain2m', 'circle4', 'cot', 'cot_', 'cotd', 'cotd_',
                                   'degrees', 'degrees90', 'degrees180', 'degrees360', 'degrees2grades', 'degrees2m',
#                                                                                      'degrees2grades as degrees2gons',
                                   'fathom2m', 'ft2m', 'furlong2m',
                                   'grades', 'grades400', 'grades2degrees', 'grades2radians',
#                                  'grades as gons', 'grades400 as gons400', 'grades2degrees as gons2degrees', 'grades2radians as gons2radians',
                                   'km2m', 'm2chain', 'm2degrees', 'm2fathom', 'm2ft', 'm2furlong',
                                   'm2km', 'm2NM', 'm2radians', 'm2SM', 'm2toise', 'm2yard', 'NM2m',
                                   'radians', 'radiansPI', 'radiansPI2', 'radiansPI_2', 'radians2m',
                                   'sincos2', 'SinCos2', 'sincos2_', 'sincos2d', 'sincos2d_', 'sincostan3', 'SM2m',
                                   'tand', 'tand_', 'tan_2', 'tanPI_2_2', 'toise2m', 'truncate',
                                   'unroll180', 'unrollPI',
                                   'wrap90', 'wrap180', 'wrap360', 'wrapPI_2', 'wrapPI', 'wrapPI2', 'wrap_normal',
                                   'yard2m'),
                            utm=_i('Utm', 'UTMError', 'parseUTM5', 'toUtm8', 'utmZoneBand5'),
                         utmups=_i('UtmUps', 'UTMUPSError', 'parseUTMUPS5', 'toUtmUps8',
                                   'utmupsValidate', 'utmupsValidateOK', 'utmupsZoneBand5'),
                     utmupsBase=_i(),  # module only
                       vector2d=_i('Circin6Tuple', 'Circum3Tuple', 'Circum4Tuple', 'Meeus2Tuple', 'Radii11Tuple', 'Soddy4Tuple',
                                   'circin6', 'circum3', 'circum4_', 'meeus2', 'radii11', 'soddy4'),
                       vector3d=_i('Vector3d', 'intersection3d3', 'iscolinearWith', 'nearestOn', 'nearestOn6', 'parse3d',
                                   'trilaterate2d2', 'trilaterate3d2'),
                   vector3dBase=_i(),  # module only
                    webmercator=_i('Wm', 'WebMercatorError', 'parseWM', 'toWm', 'EasNorRadius3Tuple'),
                           wgrs=_i('Georef', 'WGRSError'),)

_ALL_DEPRECATED = _NamedEnum_RO(_name='_ALL_DEPRECATED',
                           deprecated=_i('bases', 'datum', 'nvector',  # DEPRECATED modules and ...
                                         'rhumbaux', 'rhumbBase', 'rhumbsolve', 'rhumbx'),  # ... names
                     deprecated_bases=_i('LatLonHeightBase', 'points2'),
                   deprecated_classes=_i('ClipCS3Tuple', 'EasNorExact4Tuple', 'EcefCartesian',
                                         'HeightIDW', 'HeightIDW2', 'HeightIDW3', 'Helmert7Tuple',
                                         'LatLonExact4Tuple', 'NearestOn4Tuple', 'Ned3Tuple',
                                         'RefFrameError', 'Rhumb7Tuple', 'RhumbOrder2Tuple',
                                         'Transform7Tuple', 'TriAngle4Tuple', 'UtmUps4Tuple'),
                 deprecated_consterns=_i('EPS1_2', 'MANTIS', 'OK'),
                     deprecated_datum=_i('Curvature2Tuple', 'Datum',  'Ellipsoid',  'Transform',  # assert
                                                            'Datums', 'Ellipsoids', 'Transforms',
                                         'R_M', 'R_MA', 'R_MB', 'R_KM', 'R_NM', 'R_SM', 'R_FM', 'R_VM'),
                 deprecated_functions=_i('anStr', 'areaof', 'atand', 'bounds',  # most of the DEPRECATED functions, except ...
                                         'clipCS3', 'clipDMS', 'clipStr', 'collins', 'copysign',  # ...  ellipsoidal, spherical flavors
                                         'decodeEPSG2', 'encodeEPSG', 'enStr2', 'equirectangular3',
                                         'excessAbc', 'excessGirard', 'excessLHuilier',
                                         'false2f', 'falsed2f', 'float0', 'fStr', 'fStrzs', 'hypot3',
                                         'inStr', 'isenclosedby', 'istuplist',
                                         'joined', 'joined_', 'nearestOn3', 'nearestOn4',
                                         'parseUTM', 'perimeterof', 'polygon', 'scalar', 'simplify2',
                                         'tienstra', 'toUtm', 'trfTransforms', 'triAngle4',
                                         'unsign0', 'unStr', 'utmZoneBand2'),
                   deprecated_nvector=_i('LatLonNvectorBase', 'Nvector', 'sumOf', 'NorthPole', 'SouthPole'),)


class _ALL_MODS(object):
    '''(INTERNAL) Memoize import of any L{pygeodesy} module.
    '''
    def _DOT_(self, name):  # PYCHOK no cover
        return _DOT_(self.__class__.__name__, name)

    def __getattr__(self, name):
        '''Get a C{pygeodesy} module or attribute by B{C{name}}.

           @arg name: Qualified module or attribute name (C{str}).

           @raise ImportError: Importing module B{C{name}} failed.

           @raise AttributeError: No attribute named B{C{name}}.
        '''
        m = self.getmodule(name)
        return m if _tailof(m.__name__) == name else \
               getattr(m, _tailof(name))

    def __setattr__(self, attr, value):  # PYCHOK no cover
        t = _EQUALSPACED_(self._DOT_(attr), repr(value))
        raise AttributeError(_COLONSPACE_(t, _immutable_))

    def getattr(self, mod, *attr_dflt):  # , parent=_pygeodesy_
        '''Get an attribute of/or a C{pygeodesy} module.

           @arg mod: Qualified module name (C{str}).
           @arg attr_dflt: Optional attribute name (C{str}) and
                           optional default value (any C{type}).

           @return: The C{pygeodesy} module's attribute value.

           @raise ImportError: Importing module B{C{mod}} failed.

           @raise AttributeError: No attribute named B{C{attr}}.
        '''
        v = self.getmodule(mod)
        if attr_dflt:
            v = getattr(v, *attr_dflt)
        return v

    def getmodule(self, name, parent=_pygeodesy_):
        '''Get a C{pygeodesy} module.

           @arg name: Qualified module name (C{str}).

           @return: The C{pygeodesy} module.

           @raise ImportError: Importing module B{C{name}} failed.
        '''
        if _headof(name) != parent:
            name = _DOT_(parent, name)
        try:
            return _sys.modules[name]
        except KeyError:
            return import_module(name, parent)

    def items(self):  # no module named 'items'
        '''Yield the modules imported so far.
        '''
        for n, m in _sys.modules.items():
            yield n, m

    @property  # property_RO
    def name(self):
        return self.__class__.__name__

_ALL_MODS = _ALL_MODS()  # PYCHOK singleton

__all__ = _ALL_LAZY.lazily
__version__ = '24.02.02'


def _ALL_OTHER(*objs):
    '''(INTERNAL) Get class and function B{C{objs}} for __all__.
    '''
    _interns = _ALL_MODS.interns  # from pygeodesy import interns

    def _interned(o):  # intern'd base name
        n = _tailof(_dunder_nameof(o))
        i =  NN(_UNDER_, n, _UNDER_)  # intern'd
        return getattr(_interns, i, n)

    return tuple(map(_interned, objs))  # map2


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


def _all_deprecates():
    '''(INTERNAL) Build C{dict} of all deprecated imports.
    '''
    _D = _ALL_DEPRECATES
    if not _D:
        _ALL_DEPRECATED.fill_D(_D, _all_deprecates)  # see _all_imports()
    return _D

_ALL_DEPRECATES = _Dict()  # PYCHOK _ALL_DEPRECATED.imports()


def _all_imports():
    '''(INTERNAL) Build C{dict} of all lazy imports.
    '''
    # imports naming conventions stored below - [<key>] = <from>:
    #  import <module>                        - [<module>] = <module>
    #  from <module> import <attr>            - [<attr>] = <module>
    #  from pygeodesy import <attr>           - [<attr>] = <attr>
    #  from <module> import <attr> as <name>  - [<name>] = <module>.<attr>.
    _D = _ALL_IMPORTS
    if not _D:
        _ALL_LAZY.fill_D(_D, _all_imports)  # see _all_deprecates()
    return _D

_ALL_IMPORTS = _Dict()  # PYCHOK _ALL_LAZY.imports()


def _all_missing2(_all_):
    '''(INTERNAL) Get diffs between pygeodesy.__all__ and lazily._all_imports.
    '''
    def _diff(one, two):
        return tuple(sorted(a for a in one if a not in two))

    _alzy = _Dict((a, a) for a in _ALL_INIT)
    _alzy.update(_all_imports())  # without _all_backups!
    return ((_DOT_(_lazily_, _all_imports.__name__), _diff(_all_, _alzy)),
            (_DOT_(_pygeodesy_, _a_l_l_),     _diff(_alzy.keys(), _all_)))


def _attrof(attr_as):  # .testDeprecated
    a_, _, as_ = attr_as.partition(__as__)
    return as_ or a_.rstrip(_DOT_)


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


def _lazy_attr(unused):  # PYCHOK overwritten in _lazy_import
    pass


# def _lazy_attributes(_name_):
#     '''(INTERNAL) Return a function to C{B{_name_}.__getattr__(attr)}
#        on lazily imported modules and sub-modules.
#     '''
#     if _unlazy:
#         raise AssertionError(_COMMASPACE_(_name_, _not_(_DEPRECATED_)))
#
#     def _getattr(attr, *dflt):
#         try:  # a module name
#             return _ALL_MODS.getmodule(attr)
#         except (AttributeError, ImportError):
#             return _ALL_MODS.getattr(_name_, attr, *dflt)
#
#     return _getattr


def _lazy_import2(pack):  # MCCABE 14
    '''Check for and set up C{lazy import}.

       @arg pack: The name of the package (C{str}) performing the imports,
                  to help resolving relative imports, usually C{__package__}.

       @return: 2-Tuple C{(package, getattr)} of the importing package for
                easy reference within itself and the callable to be set to
                C{package.__getattr__}.

       @raise LazyAttributeError: The package, module or attribute name is
                                  invalid or does not exist.

       @raise LazyImportError: Lazy import not supported or not enabled or
                               an import failed.

       @note: This is I{Brett Cannon}'s function U{modutil.lazy_import
              <https://GitHub.com/brettcannon/modutil/blob/master/modutil.py>}
              modified to handle the C{__all__} and C{__dir__} attributes and
              call C{importlib.import_module(<module>.<name>, ...)} without
              causing a C{ModuleNotFoundError}.

       @see: The original U{modutil<https://PyPI.org/project/modutil>},
             U{PEP 562<https://www.Python.org/dev/peps/pep-0562>} and the
             U{new way<https://Snarky.Ca/lazy-importing-in-python-3-7/>}.
    '''
    if pack != _pygeodesy_ or _unlazy:  # new in 3.7
        t = _no_(_DOT_(pack, _lazy_import2.__name__))  # PYCHOK no cover
        raise LazyImportError(t, txt=_Python_(_sys.version))

    package, parent = _lazy_init2(pack)  # _pygeodesy_

    subpacks   = set((parent, '__main__', NN) + tuple(
                     _DOT_(parent, s) for s in _sub_packages))
    imports    = _all_imports()
    deprecates = _all_deprecates()

    def __getattr__(name):  # __getattr__ only for Python 3.7+
        # only called once for each undefined pygeodesy attribute
        mod = imports.get(name, NN) or deprecates.get(name, NN)
        if mod:
            # importlib.import_module() implicitly sets sub-modules
            # on this module as appropriate for direct imports (see
            # note in the _lazy_import2.__doc__ above).
            if mod.endswith(_DOT_):  # import mod[.attr] as name
                mod, _, attr = mod[:-1].rpartition(_DOT_)
            else:  # from mod import name
                attr = name
            try:
                t = _DOT_(pack, mod)
                imported = import_module(t, parent)
            except ImportError:
                # <https://GitHub.com/mrJean1/PyGeodesy/issues/76>
                raise LazyImportError(_no_(_module_), txt=t)
            t = getattr(imported, _p_a_c_k_a_g_e_, None)
            if t not in subpacks:  # invalid module package
                raise LazyImportError(_DOT_(mod, _p_a_c_k_a_g_e_), t)
            if attr:  # get the attribute
                imported = getattr(imported, attr, MISSING)
                if imported is MISSING:  # PYCHOK no cover
                    t = _DOT_(mod, attr)
                    # <https://GitHub.com/mrJean1/PyGeodesy/issues/76>
                    raise LazyAttributeError(_no_(_attribute_), txt=t)

        elif name in (_a_l_l_,):  # XXX '_d_i_r_', '_m_e_m_b_e_r_s_'?
            imported = _ALL_INIT + tuple(imports.keys())
        else:  # PYCHOK no cover
            t = _no_(_module_, _or_, _attribute_)
            # <https://GitHub.com/mrJean1/PyGeodesy/issues/76>
            raise LazyAttributeError(t, txt=_DOT_(parent, name))

        setattr(package, name, imported)
        if isLazy > 1:
            t = NN(_lazily_imported__, _DOT_(parent, name))
            if mod and _tailof(mod) != name:
                t = NN(t, _from_DOT__, mod)
            if isLazy > 2:
                try:  # see C{_caller3}
                    _, f, s = _caller3(2)
                    t = _SPACE_(t, _by_, f, _line_, s)
                except ValueError:
                    pass
            printf(t)  # XXX print

        return imported  # __getattr__

    global _lazy_attr
    _lazy_attr = __getattr__

    return package, __getattr__  # _lazy_import2


# def _lazy_import_all(_name_):
#     '''(INTERNAL) Return a function mimicking C{from B{__name__} import *},
#        of all items, see .deprecated.__init__
#     '''
#     if _unlazy:
#         raise AssertionError(_COMMASPACE_(_name_, _not_(_DEPRECATED_)))
#
#     _getattr = _lazy_attributes(_name_)  # _name_.__getattr__
#     _import_start = _lazy_import_star(_name_, ALL_=_ALL_IMPORTS)
#
#     def _import_all(attr, *dflt):
#         return _import_star(_name_) if attr == _a_l_l_ else \
#                _getattr(attr, *dflt)
#
#     return _import_all


def _lazy_import_as(_name_):
    '''(INTERNAL) Return a function to C{import B{__name__}.mod as mod}
       I{of modules only}, see .deprecated, .rhumb or get an attribute
       lazily exported by C{__name__}.
    '''
    if _unlazy:
        return None

    def _import_as(mod):
        try:
            return _ALL_MODS.getmodule(_DOT_(_name_, mod))
        except ImportError:
            return _lazy_attr(mod)

    return _import_as


# def _lazy_import_star(_name_, ALL_=_ALL_DEPRECATES):
#     '''(INTERNAL) Return a function to mimick C{from B{__name__} import *},
#        of all DEPRECATED items, see .deprecated, .testDeprecated
#     '''
#     if _unlazy:
#         raise AssertionError(_COMMASPACE_(_name_, _not_(_DEPRECATED_)))
#
#     def _import_star(_into_):
#         '''Do C{from B{__name__} import *} inside module C{B{__into__}}.
#         '''
#         d  =  dict()
#         nm = _tailof(_name_)
#         _g = _ALL_MODS.getattr  # pygeodesy.__getattr__
#         _h = _headof
#         for a, m in ALL_.items():
#             if _h(m) == nm:
#                 try:
#                     d[a] = _g(m, a)
#                 except (AttributeError, ImportError):
#                     pass
#         _sys.modules[_into_].__dict__.update(d)
#         return d.keys()  # imported names
#
#     return _import_star


# def _lazy_subs(_name_, force=_FOR_DOCS, over=False):
#     '''(INTERNAL) Return the names of a package's sub-packages and
#        update the package's C{__dict__} accordingly.
#     '''
#     sm = dict()
#     if force and _name_ != '__main__':
#         nm = _tailof(_name_)
#         _a = _ALL_MODS.getattr
#         _m = _ALL_MODS.getmodule
#         d  = _a(_name_, '__dict__', {})
#         for n in _a(_name_, _a_l_l_, ()):
#             try:  # n is a class name, get its mod name
#                 m = _a(_name_, n).__module__
#                 n, s = m.split(_DOT_)[-2:]
#                 if n == nm and s not in sm:
#                     # like  import m as s
#                     m = _m(m)
#                     sm[s] = m if over else d.get(s, m)
#             except (AttributeError, ImportError, ValueError) as x:
#                 pass
#         d.update(sm)
#
#     return _ALL_OTHER(*sm.values())


def _lazy_init2(pack):
    '''(INTERNAL) Initialize lazy import and set globals C{isLazy} and C{_unLazy0}.

       @arg pack: The name of the package (C{str}) performing the imports,
                  to help resolving relative imports, usually C{__package__}.

       @return: 2-Tuple C{(package, parent)} with the importing C{package}
                for easy reference within itself and its name aka the
                C{parent}, same as B{C{pack}}.

       @raise LazyImportError: Lazy import not supported or not enabled,
                               an import failed or the package name is
                               invalid or does not exist.

       @note: Global C{isLazy} is set accordingly.
    '''
    global isLazy, _unLazy0

    z = _getenv(_PYGEODESY_LAZY_IMPORT_, None)
    if z is None:  # _PYGEODESY_LAZY_IMPORT_ not set
        isLazy = 1  # ... but only by default on 3.7
    else:
        z = z.strip()  # like PYTHONVERBOSE et.al.
        isLazy = int(z) if z.isdigit() else (1 if z else 0)

    _unLazy0 = _unlazy or not isLazy  # pre-3.7 or w/o lazy import

    if isLazy < 1:  # not enabled
        raise LazyImportError(_PYGEODESY_LAZY_IMPORT_, repr(z), txt=_not_(_enabled_))
    if _getenv('PYTHONVERBOSE', None):  # PYCHOK no cover
        isLazy += 1

    try:  # to initialize in Python 3+
        package = import_module(pack)
        parent = package.__spec__.parent  # __spec__ only in Python 3.7+
        if parent != pack:  # assert
            t = _COMMASPACE_(parent, _not_(pack))  # PYCHOK no cover
            raise AttributeError(_EQUALSPACED_('parent', t))

    except (AttributeError, ImportError) as x:
        isLazy = False  # failed
        raise LazyImportError(_lazy_init2.__name__, pack, cause=x)

    return package, parent


def _pairs(*args, **kwds):  # in .ktm
    # from pygeodesy.streprs import pairs
    return _ALL_MODS.streprs.pairs(*args, **kwds)


def print_(*args, **nl_nt_prefix_end_file_flush_sep):  # PYCHOK no cover
    '''Python 3+ C{print}-like formatting and printing.

       @arg args: Arguments to be converted to C{str} and joined by B{C{sep}}
                  (any C{type}, all positional).
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
       @arg args: Arguments to be formatted (any C{type}, all positional).
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
    b, e, s, f, fl, p, kwds = _xprint7(**nl_nt_prefix_end_file_flush_sep_kwds)
    try:
        if args:
            t = (fmt % args) if fmt else s.join(map(str, args))
        elif kwds:  # PYCHOK no cover
            t = (fmt % kwds) if fmt else s.join(_pairs(kwds, prec=p))
        else:  # PYCHOK no cover
            t =  fmt
    except Exception as x:
        _E, s = _ALL_MODS.errors._xError2(x)
        unstr = _ALL_MODS.streprs.unstr
        t = unstr(printf, fmt, *args, **nl_nt_prefix_end_file_flush_sep_kwds)
        raise _E(s, txt=t, cause=x)
    try:
        n = f.write(NN(b, t, e))
    except UnicodeEncodeError:  # XXX only Windows
        t = t.replace('\u2032', _QUOTE1_).replace('\u2033', _QUOTE2_)
        n = f.write(NN(b, t, e))
    if fl:  # PYCHOK no cover
        f.flush()
    return n


def _xprint7(nl=0, nt=0, prec=6, prefix=NN, sep=_SPACE_, file=_sys.stdout,
                                            end=_NL_, flush=False, **kwds):
    '''(INTERNAL) Unravel the C{printf} and remaining keyword arguments.
    '''
    if nl > 0:
        prefix = NN(_NL_ * nl, prefix)
    if nt > 0:
        end = NN(end, _NL_ * nt)
    return prefix, end, sep, file, flush, prec, kwds


# del _i, _i0, _intern

if __name__ == '__main__':

    from timeit import timeit

    def t1():
        from pygeodesy.trf import RefFrame
        return RefFrame

    def t2():
        return _ALL_MODS.trf.RefFrame

    assert t1() is t2()  # prime each

    t1 = timeit(t1, number=1000000)
    t2 = timeit(t2, number=1000000)
    v = _Python_(_sys.version)
    printf('%.8f import vs %.8f _ALL_MODS: %.3fX, %s', t1, t2, t2 / t1, v)
    del t1, t2, v

# python3.12 -m pygeodesy.lazily
# 0.13352763 import vs 0.70804508 _ALL_MODS: 5.303X, Python 3.12.0

# % python3.11 -W ignore -m pygeodesy.lazily
# 0.37998008 import vs 0.79537812 _ALL_MODS: 2.093X, Python 3.11.5

# % python3.10 -W ignore -m pygeodesy.lazily
# 0.39046367 import vs 0.90492925 _ALL_MODS: 2.318X, Python 3.10.8

# % python2 -m pygeodesy.lazily
# 1.17563510 import vs 2.02626395 _ALL_MODS: 1.724X, Python 2.7.18

# **) MIT License
#
# Copyright (C) 2018-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
