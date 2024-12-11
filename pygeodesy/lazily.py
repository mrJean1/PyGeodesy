
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

from pygeodesy import internals as _internals, interns as _interns, \
                     _isfrozen  # DON'T _lazy_import2
# from pygeodesy.errors import _error_init, _xkwds_item2  # _ALL_MODS
from pygeodesy.internals import _caller3, _DUNDER_nameof, _getPYGEODESY, _headof, \
                                _is_DUNDER_main, printf, _tailof, _versions
from pygeodesy.interns import NN, _attribute_, _by_, _COLONSPACE_, _COMMASPACE_, \
                             _doesn_t_exist_, _DOT_, _DUNDER_all_, _EQUALSPACED_, \
                             _from_, _HASH_, _immutable_, _line_, _module_, _no_, \
                             _not_, _or_, _pygeodesy_abspath_, _pygeodesy_,  _sys, \
                             _SUB_PACKAGES, _UNDER_, _version_, _intern  # function
try:
    from importlib import import_module
except ImportError as x:  # Python 2.6-
    raise ImportError(_COLONSPACE_(x, _versions()))
# import sys as _sys  # from .interns

_a0                = ()  # PYCHOK empty tuple
_asSPACED_         = ' as '
_FOR_DOCS          = _getPYGEODESY('FOR_DOCS')  # for epydoc ...
_init__all__       = _FOR_DOCS or _getPYGEODESY('_init__all__', _DUNDER_all_) == _DUNDER_all_  # PYCHOK exported
_lazily_           = 'lazily'
_PYTHON_X_DEV      =  getattr(_sys.flags, 'dev_mode', False)  # PYCHOK Python 3.2+
_unlazy = _unLazy0 = _isfrozen or _internals._MODS.sys_version_info2 < (3, 7)  # PYCHOK mod.__getattr__ 3.7+
_WARNINGS_X_DEV    = _getPYGEODESY('WARNINGS') and (_PYTHON_X_DEV or bool(_sys.warnoptions))  # PYCHOK .props

# @module_property[_RO?] <https://GitHub.com/jtushman/proxy_tools/> <https://discuss.Python.org/t/47379>
isLazy = None  # see @var isLazy in .__init__


class LazyAttributeError(AttributeError):
    '''Raised if a C{lazily imported} attribute is missing or invalid.
    '''
    def __init__(self, *args, **kwds):
        _ALL_MODS.errors._error_init(AttributeError, self, args, **kwds)


class LazyImportError(ImportError):
    '''Raised if C{lazy import} is not supported, disabled or failed some other way.
    '''
    def __init__(self, *args, **kwds):
        _ALL_MODS.errors._error_init(ImportError, self, args, **kwds)


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
        if name in self:  # PYCHOK no cover
            sub = self[name]  # duplicate OK
            if sub != mod_ and sub not in subs:
                t = _DOT_(self._name, name)
                t = _COLONSPACE_(t,       repr(sub))
                t = _COMMASPACE_(t, _not_(repr(mod_)))
                raise AssertionError(t)
        else:
            self[name] = mod_

    def _NAME(self, which):
        self._name = _intern(_DUNDER_nameof(which).upper())


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
        t = _EQUALSPACED_(self._DOT_(attr), repr(value))
        raise LazyAttributeError(_immutable_, txt=t)

    def enums(self):
        # Yield all C{(mod_, tuple)} pairs
        for m, t in dict.items(self):
            n = m.replace(_UNDER_, _DOT_)
            if n != m:
                if m.startswith(_UNDER_):
                    continue  # skip _name= ...
                u = m.rstrip(_UNDER_)
                if u != m:
                    u = len(u)
                    n = n[:u] + m[u:]
            yield n, t

    def fill_D(self, _D, which):
        # Fill C{_Dict _D}.
        _D._NAME(which)
        _a = _D.add
        for m, t in self.enums():
            _a(m, _DOT_(m, NN, NN))  # import module
            for a in t:
                a, _, as_ = a.partition(_asSPACED_)
                if as_:  # import attr as attr_
                    _a(as_, _DOT_(m, a, NN), *_SUB_PACKAGES)
                else:
                    _a(a, m)
        return _D


def _a(*names):
    '''(INTERNAL) Intern all C{names}.
    '''
    return tuple(map(_intern, names)) if names else _a0


def _ALL_ATTRS(*attrs):
    '''(INTERNAL) Unravel all exported module attributes.
    '''
    t = ()
    for attr in attrs:
        t += tuple(map(_getattras, attr))
    return t


_ALL_INIT = _a(_pygeodesy_abspath_, _version_)

# __all__ value for most modules, accessible as _ALL_LAZY.<module>
_ALL_LAZY = _NamedEnum_RO(_name='_ALL_LAZY',
                         albers=_a('AlbersEqualArea', 'AlbersEqualArea2', 'AlbersEqualArea4',
                                   'AlbersEqualAreaCylindrical', 'AlbersEqualAreaNorth', 'AlbersEqualAreaSouth',
                                   'AlbersError', 'Albers7Tuple'),
                       auxilats=_a(),  # module only
                      azimuthal=_a('AzimuthalError', 'Azimuthal7Tuple',
                                   'Equidistant', 'EquidistantExact', 'EquidistantGeodSolve', 'EquidistantKarney',
                                   'Gnomonic', 'GnomonicExact', 'GnomonicGeodSolve', 'GnomonicKarney',
                                   'LambertEqualArea', 'Orthographic', 'Stereographic',
                                   'equidistant', 'gnomonic'),
                         basics=_a('clips', 'copysign0', 'copytype', 'halfs2',
                                   'int1s', 'isbool', 'isCartesian', 'isclass', 'iscomplex', 'isDEPRECATED', 'isfloat',
                                   'isidentifier', 'isinstanceof', 'isint', 'isiterable', 'isiterablen', 'iskeyword',
                                   'isLatLon', 'islistuple', 'isNvector', 'isodd', 'isscalar', 'issequence', 'isstr',
                                   'issubclassof', 'itemsorted',
                                   'len2', 'map1', 'map2', 'neg', 'neg_',
                                   'signBit', 'signOf', 'splice', 'str2ub', 'ub2str', 'unsigned0'),
                       booleans=_a('BooleanFHP', 'BooleanGH', 'LatLonFHP', 'LatLonGH',
                                   'isBoolean'),
                  cartesianBase=_a('RadiusThetaPhi3Tuple', 'rtp2xyz', 'rtp2xyz_', 'xyz2rtp', 'xyz2rtp_'),
                          clipy=_a('ClipCS4Tuple', 'ClipFHP4Tuple', 'ClipGH4Tuple', 'ClipLB6Tuple', 'ClipSH3Tuple',
                                   'clipCS4', 'clipFHP4', 'clipGH4', 'clipLB6', 'clipSH', 'clipSH3'),
                            css=_a('CassiniSoldner', 'Css', 'CSSError', 'toCss',
                                   'EasNorAziRk4Tuple', 'EasNorAziRkEqu6Tuple', 'LatLonAziRk4Tuple'),
                      constants=_a('DIG', 'EPS', 'EPS0', 'EPS02', 'EPS1', 'EPS2', 'EPS4', 'EPS_2',
                                   'INF', 'INT0', 'MANT_DIG', 'MAX', 'MAX_EXP', 'MIN', 'MIN_EXP', 'NAN', 'NEG0', 'NINF',
                                   'PI', 'PI2', 'PI_2', 'PI3', 'PI_3', 'PI3_2', 'PI4', 'PI_4',
                                   'R_FM', 'R_GM', 'R_KM', 'R_M', 'R_MA', 'R_MB', 'R_NM', 'R_QM', 'R_SM', 'R_VM',
                                   'float_', 'float0_', 'isclose', 'isfinite', 'isinf', 'isint0',
                                   'isnan', 'isnear0', 'isnear1', 'isnear90', 'isneg0', 'isninf', 'isnon0',
                                   'remainder'),
                         datums=_a('Datum', 'Datums', 'Transform', 'Transforms'),
#                    deprecated=_a(),  # module only
                            dms=_a('F_D',   'F_DM',   'F_DMS',   'F_DEG',   'F_MIN',   'F_SEC',   'F_D60',   'F__E',   'F__F',   'F__G',   'F_RAD',
                                   'F_D_',  'F_DM_',  'F_DMS_',  'F_DEG_',  'F_MIN_',  'F_SEC_',  'F_D60_',  'F__E_',  'F__F_',  'F__G_',  'F_RAD_',
                                   'F_D__', 'F_DM__', 'F_DMS__', 'F_DEG__', 'F_MIN__', 'F_SEC__', 'F_D60__', 'F__E__', 'F__F__', 'F__G__', 'F_RAD__',
                                   'S_DEG', 'S_MIN', 'S_SEC', 'S_DMS', 'S_RAD', 'S_SEP',
                                   'bearingDMS', 'clipDegrees', 'clipRadians', 'compassDMS', 'compassPoint',
                                   'degDMS', 'latDMS', 'latlonDMS', 'latlonDMS_', 'lonDMS', 'normDMS',
                                   'parseDDDMMSS', 'parseDMS', 'parseDMS2', 'parse3llh', 'parseRad', 'precision', 'toDMS'),
                           ecef=_a('EcefError', 'EcefFarrell21', 'EcefFarrell22', 'EcefKarney', 'EcefMatrix',
                                   'EcefSudano', 'Ecef9Tuple', 'EcefVeness', 'EcefYou'),
                     elevations=_a('Elevation2Tuple', 'GeoidHeight2Tuple',
                                   'elevation2', 'geoidHeight2'),
                ellipsoidalBase=_a(),  # module only
              ellipsoidalBaseDI=_a(),  # module only
               ellipsoidalExact=_a(),  # module only
           ellipsoidalGeodSolve=_a(),  # module only
              ellipsoidalKarney=_a(),  # module only
             ellipsoidalNvector=_a(),  # module only
            ellipsoidalVincenty=_a('VincentyError',),  # nothing else
                     ellipsoids=_a('a_f2Tuple', 'Circle4Tuple', 'Curvature2Tuple',
                                   'Ellipsoid', 'Ellipsoid2', 'Ellipsoids',
                                   'a_b2e', 'a_b2e2', 'a_b2e22', 'a_b2e32', 'a_b2f', 'a_b2f_', 'a_b2f2', 'a_b2n',
                                   'a_f2b', 'a_f_2b', 'b_f2a', 'b_f_2a',
                                   'e2f', 'e22f',
                                   'f2e2', 'f2e22', 'f2e32', 'f_2f', 'f2f_', 'f2f2', 'f2n', 'n2e2', 'n2f', 'n2f_'),
                       elliptic=_a('Elliptic', 'EllipticError', 'Elliptic3Tuple'),
                           epsg=_a('Epsg', 'EPSGError'),
                         errors=_a('AuxError', 'ClipError', 'CrossError', 'GeodesicError', 'IntersectionError',
                                   'NumPyError', 'LenError', 'LimitError', 'MGRSError',
                                   'ParseError', 'PointsError', 'RangeError', 'RhumbError',
                                   'SciPyError', 'SciPyWarning', 'TRFError', 'TriangleError', 'UnitError', 'VectorError',
                                   'crosserrors', 'exception_chaining', 'isError', 'limiterrors', 'rangerrors'),
                            etm=_a('Etm', 'ETMError', 'ExactTransverseMercator',
                                   'parseETM5', 'toEtm8'),
                          fmath=_a('Fdot', 'Fhorner', 'Fhypot', 'Fpolynomial', 'Fpowers', 'Fcbrt', 'Froot', 'Fsqrt',
                                   'bqrt', 'cbrt', 'cbrt2', 'euclid', 'euclid_',
                                   'facos1', 'fasin1', 'fatan', 'fatan1', 'fatan2', 'favg',
                                   'fdot', 'fdot_', 'fdot3', 'fma', 'fmean', 'fmean_', 'fhorner', 'fidw', 'f2mul_',
                                   'fpolynomial', 'fpowers', 'fprod', 'frandoms', 'frange', 'freduce', 'fremainder',
                                   'hypot', 'hypot_', 'hypot1', 'hypot2', 'hypot2_',
                                   'norm2', 'norm_', 'sqrt0', 'sqrt3', 'sqrt_a', 'zcrt', 'zqrt'),
                          formy=_a('Radical2Tuple',
                                   'angle2chord', 'antipode', 'antipode_', 'bearing', 'bearing_', 'chord2angle',
                                   'compassAngle', 'cosineForsytheAndoyerLambert', 'cosineForsytheAndoyerLambert_',
                                   'cosineAndoyerLambert', 'cosineAndoyerLambert_', 'cosineLaw', 'cosineLaw_',
                                   'equirectangular', 'equirectangular4', 'euclidean', 'euclidean_',
                                   'excessAbc_', 'excessCagnoli_', 'excessGirard_', 'excessLHuilier_',
                                   'excessKarney', 'excessKarney_', 'excessQuad', 'excessQuad_',
                                   'flatLocal', 'flatLocal_', 'flatPolar', 'flatPolar_',
                                   'hartzell', 'haversine', 'haversine_', 'heightOf', 'heightOrthometric', 'horizon', 'hubeny', 'hubeny_',
                                   'intersection2', 'intersections2', 'isantipode', 'isantipode_', 'isnormal', 'isnormal_',
                                   'normal', 'normal_', 'opposing', 'opposing_', 'radical2',
                                   'thomas', 'thomas_', 'vincentys', 'vincentys_'),
                        frechet=_a('Frechet', 'FrechetDegrees', 'FrechetError', 'FrechetRadians',
                                   'FrechetCosineAndoyerLambert', 'FrechetCosineForsytheAndoyerLambert',
                                   'FrechetCosineLaw', 'FrechetDistanceTo', 'FrechetEquirectangular',
                                   'FrechetEuclidean', 'FrechetExact', 'FrechetFlatLocal', 'FrechetFlatPolar',
                                   'FrechetHaversine', 'FrechetHubeny', 'FrechetKarney', 'FrechetThomas',
                                   'FrechetVincentys', 'Frechet6Tuple',
                                   'frechet_'),
                         fstats=_a('Fcook', 'Flinear', 'Fwelford'),
                          fsums=_a('Fsum', 'DivMod2Tuple', 'Fsum2Tuple', 'ResidualError',
                                   'f2product', 'fsum', 'fsum_', 'fsumf_', 'fsum1', 'fsum1_', 'fsum1f_', 'nonfiniterrors'),
                           gars=_a('Garef', 'GARSError'),
                      geodesici=_a('Intersectool', 'Intersectool5Tuple', 'Intersect7Tuple',
                                   'Intersector',  'Intersector5Tuple',  'Middle5Tuple', 'XDict'),
                      geodesicw=_a('Geodesic', 'GeodesicLine', 'Geodesic_WGS84'),
                      geodesicx=_a('gx', 'gxarea', 'gxbases', 'gxline',  # modules
                                   'GeodesicAreaExact', 'GeodesicExact', 'GeodesicLineExact', 'PolygonArea'),
                      geodsolve=_a('GeodesicSolve', 'GeodesicLineSolve', 'GeodSolve12Tuple'),
                        geohash=_a('Geohash', 'Geohashed', 'GeohashError', 'Neighbors8Dict', 'Resolutions2Tuple', 'Sizes3Tuple'),
                         geoids=_a('GeoidError', 'GeoidG2012B', 'GeoidKarney', 'GeoidPGM', 'egmGeoidHeights',
                                   'PGMError', 'GeoidHeight5Tuple'),
                      hausdorff=_a('Hausdorff', 'HausdorffDegrees', 'HausdorffError', 'HausdorffRadians',
                                   'HausdorffCosineAndoyerLambert', 'HausdorffCosineForsytheAndoyerLambert',
                                   'HausdorffCosineLaw', 'HausdorffDistanceTo', 'HausdorffEquirectangular',
                                   'HausdorffEuclidean', 'HausdorffExact', 'HausdorffFlatLocal', 'HausdorffFlatPolar',
                                   'HausdorffHaversine', 'HausdorffHubeny', 'HausdorffKarney', 'HausdorffThomas',
                                   'HausdorffVincentys', 'Hausdorff6Tuple',
                                   'hausdorff_', 'randomrangenerator'),
                        heights=_a('HeightCubic', 'HeightError',
                                   'HeightIDWcosineAndoyerLambert', 'HeightIDWcosineForsytheAndoyerLambert',
                                   'HeightIDWcosineLaw', 'HeightIDWdistanceTo', 'HeightIDWequirectangular',
                                   'HeightIDWeuclidean', 'HeightIDWexact', 'HeightIDWflatLocal', 'HeightIDWflatPolar',
                                   'HeightIDWhaversine', 'HeightIDWhubeny', 'HeightIDWkarney', 'HeightIDWthomas',
                                   'HeightIDWvincentys', 'HeightLinear', 'HeightLSQBiSpline', 'HeightSmoothBiSpline'),
                      internals=_internals.__all__,
                        interns=_interns.__all__,
                          iters=_a('LatLon2PsxyIter', 'PointsIter', 'points2',
                                   'isNumpy2', 'isPoints2', 'isTuple2', 'iterNumpy2', 'iterNumpy2over'),
                         karney=_a('Area3Tuple', 'Caps', 'Direct9Tuple', 'GDict', 'Inverse10Tuple', 'Rhumb8Tuple'),
                            ktm=_a('KTMError', 'KTransverseMercator'),
                     latlonBase=_a('latlon2n_xyz', 'philam2n_xyz'),
                         lazily=_a('LazyAttributeError', 'LazyImportError', 'isLazy'),
                            lcc=_a('Conic', 'Conics', 'Lcc', 'LCCError', 'toLcc'),
                            ltp=_a('Attitude', 'AttitudeError', 'ChLV', 'ChLVa', 'ChLVe', 'Frustum',
                                   'LocalCartesian', 'LocalError', 'Ltp', 'tyr3d'),
                      ltpTuples=_a('Aer', 'Aer4Tuple', 'Attitude4Tuple',
                                   'ChLVEN2Tuple', 'ChLV9Tuple', 'ChLVYX2Tuple', 'ChLVyx2Tuple',
                                   'Enu', 'Enu4Tuple', 'Footprint5Tuple', 'Local9Tuple', 'Los',
                                   'Ned', 'Ned4Tuple', 'Uvw', 'Uvw3Tuple', 'XyzLocal', 'Xyz4Tuple'),
                           mgrs=_a('Mgrs', 'parseMGRS', 'toMgrs', 'Mgrs4Tuple', 'Mgrs6Tuple'),
                          named=_a('ADict',
                                   'callername', 'classname', 'classnaming', 'modulename',
                                   'nameof', 'notImplemented', 'notOverloaded'),
                    namedTuples=_a('Bearing2Tuple', 'Bounds2Tuple', 'Bounds4Tuple',
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
                    nvectorBase=_a('NorthPole', 'SouthPole', 'n_xyz2latlon', 'n_xyz2philam'),
                           osgr=_a('Osgr', 'OSGRError', 'parseOSGR', 'toOsgr'),
                         points=_a('LatLon_', 'LatLon2psxy', 'Numpy2LatLon', 'Shape2Tuple', 'Tuple2LatLon',
                                   'areaOf', 'boundsOf', 'centroidOf', 'fractional',
                                   'isclockwise', 'isconvex', 'isconvex_', 'isenclosedBy', 'ispolar',
                                   'luneOf', 'nearestOn5', 'perimeterOf', 'quadOf'),
                          props=_a('Property', 'Property_RO', 'property_doc_',
                                   'property_RO', 'property_ROnce', 'property_ROver',
                                   'deprecated_class', 'deprecated_function', 'deprecated_method',
                                   'deprecated_Property_RO', 'deprecated_property_RO', 'DeprecationWarnings'),
                     resections=_a('Collins5Tuple', 'ResectionError', 'Survey3Tuple', 'Tienstra7Tuple',
                                   'TriAngle5Tuple', 'TriSide2Tuple', 'TriSide4Tuple',
                                   'cassini', 'collins5', 'pierlot', 'pierlotx', 'tienstra7',
                                   'snellius3', 'wildberger3',
                                   'triAngle', 'triAngle5', 'triArea', 'triSide', 'triSide2', 'triSide4'),
                          rhumb=_a(),  # module only
                     rhumb_aux_=_a('RhumbAux', 'RhumbLineAux'),
                      rhumb_ekx=_a('Rhumb', 'RhumbLine'),
                    rhumb_solve=_a('RhumbSolve', 'RhumbLineSolve', 'RhumbSolve7Tuple'),
                  sphericalBase=_a(),  # module only
               sphericalNvector=_a(),  # module only
          sphericalTrigonometry=_a(),  # module only
                       simplify=_a('simplify1', 'simplifyRDP', 'simplifyRW', 'simplifyVW'),
                      solveBase=_a(),  # module only
                        streprs=_a('anstr', 'attrs', 'enstr2', 'fstr', 'fstrzs', 'hstr', 'instr',
                                   'lrstrip', 'pairs', 'reprs', 'strs', 'unstr'),
                            trf=_a('RefFrame', 'RefFrames', 'TransformXform', 'TRFXform', 'TRFXform7Tuple',
                                   'date2epoch', 'epoch2date', 'trfTransform0', 'trfTransforms', 'trfXform'),
                      triaxials=_a('BetaOmega2Tuple', 'BetaOmega3Tuple', 'Jacobi2Tuple',
                                   'JacobiConformal', 'JacobiConformalSpherical',
                                   'Triaxial', 'Triaxial_', 'TriaxialError', 'Triaxials', 'hartzell4'),
                          units=_a('Azimuth', 'Band', 'Bearing', 'Bearing_', 'Bool',
                                   'Degrees', 'Degrees_', 'Degrees2', 'Distance', 'Distance_', 'Easting', 'Epoch',
                                   'Feet', 'FIx', 'Float_', 'Height', 'Height_', 'HeightX', 'Int_',
                                   'Lam', 'Lamd', 'Lat', 'Lat_', 'Lon', 'Lon_',
                                   'Meter', 'Meter_', 'Meter2', 'Meter3', 'Northing', 'Number_',
                                   'Phi', 'Phid', 'Precision_', 'Radians', 'Radians_', 'Radians2',
                                   'Radius_', 'Scalar', 'Scalar_', 'Zone'),
                      unitsBase=_a('Float', 'Int', 'Radius', 'Str'),
                            ups=_a('Ups', 'UPSError', 'parseUPS5', 'toUps8', 'upsZoneBand5'),
                          utily=_a('acos1', 'acre2ha', 'acre2m2', 'asin1', 'atan1', 'atan1d', 'atan2b', 'atan2d',
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
                                   'tan', 'tan_', 'tand', 'tand_', 'tan_2', 'tanPI_2_2', 'toise2m', 'truncate',
                                   'unroll180', 'unrollPI',
                                   'wrap90', 'wrap180', 'wrap360', 'wrapPI_2', 'wrapPI', 'wrapPI2', 'wrap_normal',
                                   'yard2m'),
                            utm=_a('Utm', 'UTMError', 'parseUTM5', 'toUtm8', 'utmZoneBand5'),
                         utmups=_a('UtmUps', 'UTMUPSError', 'parseUTMUPS5', 'toUtmUps8',
                                   'utmupsValidate', 'utmupsValidateOK', 'utmupsZoneBand5'),
                     utmupsBase=_a(),  # module only
                       vector2d=_a('Circin6Tuple', 'Circum3Tuple', 'Circum4Tuple', 'Meeus2Tuple', 'Radii11Tuple', 'Soddy4Tuple', 'Triaxum5Tuple',
                                   'circin6', 'circum3', 'circum4', 'circum4_', 'meeus2', 'radii11', 'soddy4', 'triaxum5', 'trilaterate2d2'),
                       vector3d=_a('Vector3d', 'intersection3d3', 'iscolinearWith', 'nearestOn', 'nearestOn6', 'parse3d',
                                   'trilaterate3d2'),
                   vector3dBase=_a(),  # module only
                    webmercator=_a('Wm', 'WebMercatorError', 'parseWM', 'toWm', 'EasNorRadius3Tuple'),
                           wgrs=_a('Georef', 'WGRSError'),)

_ALL_DEPRECATED = _NamedEnum_RO(_name='_ALL_DEPRECATED',
                           deprecated=_a('bases', 'datum', 'nvector',  # DEPRECATED modules and ...
                                         'rhumbaux', 'rhumbBase', 'rhumbsolve', 'rhumbx'),  # ... names
                     deprecated_bases=_a('LatLonHeightBase', 'points2'),
                   deprecated_classes=_a('ClipCS3Tuple', 'EasNorExact4Tuple', 'EcefCartesian', 'Fn_rt',
                                         'HeightIDW', 'HeightIDW2', 'HeightIDW3', 'Helmert7Tuple',
                                         'Lam_', 'LatLonExact4Tuple', 'NearestOn4Tuple', 'Ned3Tuple',
                                         'Phi_', 'RefFrameError', 'Rhumb7Tuple', 'RhumbOrder2Tuple',
                                         'Transform7Tuple', 'TriAngle4Tuple', 'UtmUps4Tuple', 'XDist'),
                 deprecated_consterns=_a('EPS1_2', 'MANTIS', 'OK'),
                     deprecated_datum=_a('Curvature2Tuple', 'Datum',  'Ellipsoid',  'Transform',  # assert
                                                            'Datums', 'Ellipsoids', 'Transforms',
                                         'R_FM', 'R_KM', 'R_M', 'R_MA', 'R_MB', 'R_NM', 'R_SM', 'R_VM'),
                 deprecated_functions=_a('anStr', 'areaof', 'atand', 'bounds',  # most of the DEPRECATED functions, except ellipsoidal ...
                                         'clipCS3', 'clipDMS', 'clipStr', 'collins', 'copysign',  # ... and spherical flavors
                                         'decodeEPSG2', 'encodeEPSG', 'enStr2', 'equirectangular_', 'equirectangular3',
                                         'excessAbc', 'excessGirard', 'excessLHuilier',
                                         'false2f', 'falsed2f', 'float0', 'fStr', 'fStrzs', 'Fsum2product',
                                         'hypot3', 'inStr', 'isenclosedby', 'istuplist',
                                         'joined', 'joined_', 'nearestOn3', 'nearestOn4',
                                         'parseUTM', 'perimeterof', 'polygon',
                                         'scalar', 'simplify2', 'simplifyRDPm', 'simplifyVWm',
                                         'tienstra', 'toUtm', 'triAngle4',
                                         'unsign0', 'unStr', 'utmZoneBand2'),
                   deprecated_nvector=_a('LatLonNvectorBase', 'Nvector', 'sumOf', 'NorthPole', 'SouthPole'),)


class _ALL_MODS(_internals._MODS_Base):
    '''(INTERNAL) Memoized import of any L{pygeodesy} module.
    '''
    def __getattr__(self, name):
        '''Get a C{pygeodesy} module or attribute by B{C{name}}.

           @arg name: Un/qualified module or qualified attribute name (C{str}).

           @raise ImportError: Importing module B{C{name}} failed.

           @raise AttributeError: No attribute named B{C{name}}.
        '''
        try:
            v = _lazy_dict[name]  # package.__dict__
        except KeyError:
            v = _lazy_module(name)  # package.__getattr__
            if _tailof(_DUNDER_nameof(v)) != name:
                try:
                    v = getattr(v, _tailof(name))
                except AttributeError:
                    pass  # XXX LazyAttributeError?
        return v

    def getattr(self, name, *attr_dflt):  # , parent=_pygeodesy_
        '''Get an attribute of/or a C{pygeodesy} module.

           @arg name: Un/qualified module name (C{str}).
           @arg attr_dflt: Optional attribute name (C{str}) and
                           optional default value (any C{type}).

           @return: The C{pygeodesy} module's attribute value.

           @raise ImportError: Importing module B{C{name}} failed.

           @raise AttributeError: No attribute named B{C{attr}}.
        '''
        v = self.getmodule(name)
        if attr_dflt:
            v = getattr(v, *attr_dflt)
        return v

    def getmodule(self, name, parent=_pygeodesy_):
        '''Get a C{pygeodesy} module or the C{__main__}.

           @arg name: Un/qualified module name (C{str}).

           @return: The C{pygeodesy} module.

           @raise ImportError: Importing module B{C{name}} failed.
        '''
        if _headof(name) != parent and not _is_DUNDER_main(name):
            name = _DOT_(parent, name)
        try:
            return _sys.modules[name]
        except KeyError:
            return _getmodule(name, parent)

    def imported(self, name):
        '''Return module or package C{name} if already imported.
        '''
        return _sys.modules.get(name, None)

    def into(self, **mod_DUNDER_name):
        '''Lazily import module C{mod} into module C{DUNDER_name}
           and set C{DUNDER_name._mod} to module C{mod}, I{once}.
        '''
        class _Into(object):

            def __getattr__(unused, name):
                mod, dun = self.errors._xkwds_item2(mod_DUNDER_name)
                _mod = _UNDER_(NN, mod)
                d =  self.getmodule(dun)  # '__main__' OK
                i = _getattribute(d, _mod, dun)
                assert isinstance(i, _Into)
                m =  self.getmodule(mod)
                setattr(d, _mod, m)  # overwrite C{d._mod}
                return getattr(m, name)

        return _Into()

#   @_Property_RO
#   def _isBoolean(self):
#       '''(INTERNAL) Get function C(.booleans.isBoolean}, I{once}.
#       '''
#       return self.booleans.isBoolean

    def items(self):  # no module named 'items'
        '''Yield the modules imported so far.
        '''
        for n, m in _sys.modules.items():
            if _headof(n) == _pygeodesy_:
                yield n, m

_internals._MODS = _ALL_MODS = _ALL_MODS()  # PYCHOK singleton

__all__ = _ALL_LAZY.lazily
__version__ = '24.11.28'


def _ALL_OTHER(*objs):
    '''(INTERNAL) Get class and function B{C{objs}} for __all__.
    '''
    def _interned(o):  # intern'd base name
        n = _tailof(_DUNDER_nameof(o))
        i =  NN(_UNDER_, n, _UNDER_)  # intern'd
        return getattr(_interns, i, n)

    return tuple(map(_interned, objs))  # map2


if _FOR_DOCS:  # PYCHOK no cover
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
    return ((_DOT_(_lazily_, _all_imports.__name__),  _diff(_all_, _alzy)),
            (_DOT_(_pygeodesy_, _DUNDER_all_), _diff(_alzy.keys(), _all_)))


def _getattras(attr_as):  # test/testDeprecated
    '''(INTERNAL) Get the C{"as name"} or C{"name"} of a lazy entry.
    '''
    a_, _, as_ = attr_as.partition(_asSPACED_)
    return as_ or a_.rstrip(_DOT_)


def _getattribute(m, name, mod=_pygeodesy_):
    '''(INTERNAL) Get attr C{m.name}.
    '''
    try:
        return getattr(m, name)
    except AttributeError:
        name = _DOT_(mod, name)
        # <https://GitHub.com/mrJean1/PyGeodesy/issues/76>
        raise LazyAttributeError(_no_(_attribute_), txt=name)


def _getmodule(name, *parent):
    '''(INTERNAL) Wrapper for C{import_module}.
    '''
    try:
        return import_module(name, parent)
    except ImportError:
        # <https://GitHub.com/mrJean1/PyGeodesy/issues/76>
        raise LazyImportError(_no_(_module_), txt=name)


# def _lazy_attributes(DUNDER_name):
#     '''(INTERNAL) Return a function to C{B{__name__}.__getattr__(attr)}
#        on lazily imported modules and sub-modules.
#     '''
#     if _unlazy:
#         raise AssertionError(_COMMASPACE_(DUNDER_name, _not_(_DEPRECATED_)))
#
#     def _getattr(attr, *dflt):
#         try:  # a module name
#             return _ALL_MODS.getmodule(attr)
#         except (AttributeError, ImportError):
#             return _ALL_MODS.getattr(DUNDER_name, attr, *dflt)
#
#     return _getattr


_lazy_dict = {}  # PYCHOK overwritten by _lazy_import2


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
    _DOT_   = _interns._DOT_
    _SPACE_ = _interns._SPACE_

    if pack != _pygeodesy_ or _unlazy:  # Python 3.7+ # PYCHOK no cover
        t = _DOT_(pack, _DUNDER_nameof(_lazy_import2))
        raise LazyImportError(_no_(t), txt=_versions())

    package, parent   = _lazy_init2(pack)  # _pygeodesy_

    _DUNDER_package_  = '__package__'
    _lazily_imported_ = _SPACE_(_HASH_, _lazily_, 'imported', parent)

    sub_packages =  set((parent, NN) + tuple(
                   _DOT_(parent, s) for s in _SUB_PACKAGES))
    imports      = _all_imports()
    deprecates   = _all_deprecates()

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
            v = _getmodule(_DOT_(pack, mod), parent)
            t =  getattr(v, _DUNDER_package_, None)
            if t not in sub_packages:  # invalid module package
                raise LazyImportError(_DOT_(mod, _DUNDER_package_), t)
            if attr:  # get mod.attr
                v = _getattribute(v, attr, mod)

        elif name in (_DUNDER_all_,):  # XXX _DUNDER_dir_, _DUNDER_members_?
            v = _ALL_INIT + tuple(imports.keys())
        else:  # PYCHOK no cover
            t = _no_(_module_, _or_, _attribute_)
            # <https://GitHub.com/mrJean1/PyGeodesy/issues/76>
            raise LazyAttributeError(t, txt=_DOT_(parent, name))

        setattr(package, name, v)  # package.__dict__[name] = val
        if isLazy > 1:
            t = _DOT_(_lazily_imported_, name)
            if mod and _tailof(mod) != name:
                t = _SPACE_(t, _from_, _DOT_(NN, mod))
            if isLazy > 2:
                try:  # see C{_caller3}
                    _, f, s = _caller3(2)
                    t = _SPACE_(t, _by_, f, _line_, s)
                except ValueError:
                    pass
            printf(t)  # XXX print

        return v  # __getattr__

    global _lazy_dict, _lazy_module
    _lazy_dict   = package.__dict__
    _lazy_module = __getattr__

    return package, __getattr__  # _lazy_import2


# def _lazy_import_all(DUNDER_name):
#     '''(INTERNAL) Return a function mimicking C{from B{__name__} import *},
#        of all items, see .deprecated.__init__
#     '''
#     if _unlazy:
#         raise AssertionError(_COMMASPACE_(DUNDER_name, _not_(_DEPRECATED_)))
#
#     _getattr = _lazy_attributes(DUNDER_name)  # __name__.__getattr__
#     _import_start = _lazy_import_star(DUNDER_name, ALL_=_ALL_IMPORTS)
#
#     def _import_all(attr, *dflt):
#         return _import_star(DUNDER_name) if attr == _DUNDER_all_ else \
#                _getattr(attr, *dflt)
#
#     return _import_all


def _lazy_import_as(DUNDER_name):
    '''(INTERNAL) Return a function to C{import B{__name__}.mod as mod}
       I{of modules only}, see .deprecated, .rhumb or get an attribute
       lazily exported by C{__name__}.
    '''
    if _unlazy:
        return None

    def _import_as(mod):
        try:
            return _ALL_MODS.getmodule(_DOT_(DUNDER_name, mod))
        except ImportError:
            return _lazy_module(mod)

    return _import_as


# def _lazy_import_star(DUNDER_name, ALL_=_ALL_DEPRECATES):
#     '''(INTERNAL) Return a function to mimick C{from B{__name__} import *},
#        of all DEPRECATED items, see .deprecated, .testDeprecated
#     '''
#     if _unlazy:
#         raise AssertionError(_COMMASPACE_(DUNDER_name, _not_(_DEPRECATED_)))
#
#     def _import_star(_into_):
#         '''Do C{from B{__name__} import *} inside module C{B{__into__}}.
#         '''
#         d  =  dict()
#         nm = _tailof(DUNDER_name)
#         _g = _ALL_MODS.getattr  # pygeodesy.__getattr__
#         for a, m in ALL_.items():
#             if _headof(m) == nm:
#                 try:
#                     d[a] = _g(m, a)
#                 except (AttributeError, ImportError):
#                     pass
#         _sys.modules[_into_].__dict__.update(d)
#         return d.keys()  # imported names
#
#     return _import_star


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

    z = _getPYGEODESY('LAZY_IMPORT', None)
    if z is None:  # not set, but ...
        isLazy = 1  # ... default on 3.7+
    else:
        z = z.strip()  # like PYTHONVERBOSE et.al.
        isLazy = int(z) if z.isdigit() else (1 if z else 0)

    _unLazy0 = _unlazy or not isLazy  # pre-3.7 or w/o lazy import

    if isLazy < 1:  # invalid, not enabled
        e = _internals._PYGEODESY('LAZY_IMPORT')
        raise LazyImportError(e, repr(z), txt_not_='enabled')
    if _sys.flags.verbose:  # PYCHOK no cover
        isLazy += 1

    try:  # to initialize in Python 3+
        package = import_module(pack)
        parent = package.__spec__.parent  # __spec__ only in Python 3.7+
        if parent != pack:  # assert
            t = _COMMASPACE_(parent, _not_(pack))  # PYCHOK no cover
            raise AttributeError(_EQUALSPACED_('parent', t))

    except (AttributeError, ImportError) as x:
        isLazy = False  # failed
        z = _DUNDER_nameof(_lazy_init2)
        raise LazyImportError(z, pack, cause=x)

    return package, parent


def _lazy_module(name):  # overwritten by _lazy_import2
    '''(INTERNAL) Get or import a C{pygeodesy} module.
    '''
    try:  # most likely ... module has been imported
        m = _ALL_MODS.getmodule(name)
    except (AttributeError, ImportError) as x:
        raise LazyImportError(name, cause=x)
    _lazy_dict[name] = m  # cache
    return m


# def _lazy_subs(DUNDER_name, force=_FOR_DOCS, over=False):
#     '''(INTERNAL) Return the names of a package's sub-packages and
#        update the package's C{__dict__} accordingly.
#     '''
#     sm = dict()
#     if force and not _is_DUNDER_main(DUNDER_name):
#         nm = _tailof(DUNDER_name)
#         _a = _ALL_MODS.getattr
#         _m = _ALL_MODS.getmodule
#         d  = _a(DUNDER_name, _DUNDER_dict_, {})
#         for n in _a(DUNDER_name, _DUNDER_all_, ()):
#             try:  # n is a class name, get its mod name
#                 m = _a(DUNDER_name, n).__module__
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

# del _a, _a0

if _is_DUNDER_main(__name__):  # PYCHOK no cover

    def _main():
        from timeit import timeit

        def t1():
            from pygeodesy.trf import RefFrame
            return RefFrame

        def t2():
            return _ALL_MODS.trf.RefFrame

        assert t1() is t2()  # prime each

        t1 =  timeit(t1, number=1000000)
        t2 =  timeit(t2, number=1000000)
        A  = _DUNDER_nameof(_ALL_MODS.__class__)
        v  = _versions()
        printf('%.6f import vs %.6f %s: %.2fX, %s', t1, t2, A, (t1 / t2), v)

    _main()

# % python3.13 -W ignore -m pygeodesy.lazily
# 0.106602 import vs 0.078136 _ALL_MODS: 1.36X, pygeodesy 24.10.24 Python 3.13.0 64bit arm64 macOS 14.6.1

# % python3.12 -W ignore -m pygeodesy.lazily
# 0.138844 import vs 0.080458 _ALL_MODS: 1.73X, pygeodesy 24.10.24 Python 3.12.7 64bit arm64 macOS 14.6.1

# % python3.11 -W ignore -m pygeodesy.lazily
# 0.387520 import vs 0.254229 _ALL_MODS: 1.52X, pygeodesy 24.10.24 Python 3.11.5 64bit arm64 macOS 14.6.1

# % python3.10 -W ignore -m pygeodesy.lazily
# 0.371269 import vs 0.272897 _ALL_MODS: 1.36X, pygeodesy 24.10.24 Python 3.10.8 64bit arm64 macOS 14.6.1

# % python3.8 -W ignore -m pygeodesy.lazily
# 0.555572 import vs 0.370304 _ALL_MODS: 1.50X, pygeodesy 24.10.24 Python 3.8.10 64bit arm64_x86_64 macOS 10.16

# % python2 -m pygeodesy.lazily
# 1.160292 import vs 0.490279 _ALL_MODS: 2.37X, pygeodesy 24.10.24 Python 2.7.18 64bit arm64_x86_64 macOS 10.16

# **) MIT License
#
# Copyright (C) 2018-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
