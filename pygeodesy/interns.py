# -*- coding: utf-8 -*-

u'''Single C{str}ing constants, C{intern}'ed across C{pygeodesy}
modules and function L{pygeodesy.machine}.
'''
import sys as _sys

_COMMASPACE_  = ', '  # overriden below
_pf2List      = []    # cached _platform2 list
_Py3List      = []    # cached _pythonarchine list

_sub_packages = 'auxilats', 'deprecated', 'geodesicx'  # PYCHOK in .lazily,
# ... make._dist, MANIFEST, setup.setup, test.bases, test.testModules


class _Dash(str):
    '''(INTERNAL) Extended C{str} for prefix_DASH_.
    '''
    def __call__(self, *args):
        '''Join C{self} plus all B{C{args}} like C{'-'.join((self,) + B{args})}.
        '''
        return _DASH_(self, *args)  # re-callable


class _Int(int):
    '''(INTERNAL) Unique C{int}.
    '''
    pass


class Str_(str):
    '''Extended, I{callable} C{str} class, not nameable.

       @see: Nameable and callable class L{pygeodesy.Str}.
    '''
    def join_(self, *args):
        '''Join all positional B{C{args}} like C{self.join(B{args})}.

           @return: All B{C{args}} joined by this instance (L{Str_}).

           @note: An other L{Str_} instance is returned to make the
                  result re-callable.
        '''
        return Str_(str.join(self, map(str, args)))  # re-callable

    __call__ = join_

NN = Str_('')  # PYCHOK Nomen Nescio <https://Wiktionary.org/wiki/N.N.>


class _Prefix(Str_):
    '''(INTERNAL) Extended C{str} for prefix.
    '''
    def __call__(self, *args):
        '''Join C{self} plus all B{C{args}} like C{" ".join((self,) + B{args})}.
        '''
        return _SPACE_.join_(self, *args)  # re-callable


class _PyPy__(str):  # overwritten by singleton below
    '''(INTERNAL) Extended C{str} for C{"PyPy"} and version.
    '''
    def __call__(self, version=NN):
        '''Return C{"PyPy <version>"} or C{NN}.
        '''
        v = version or _sys.version
        if _PyPy__ in v:
            v = v.split(_PyPy__)[1].split(None, 1)[0]  # == _DOT_.join(_sys.pypy_version_info[:3])
            return NN(_PyPy__, v)
        else:
            return NN


class _Python_(str):  # overwritten by singleton below
    '''(INTERNAL) Extended C{str} for C{"Python"} and version.
    '''
    def __call__(self, version=NN):
        '''Return C{"Python <version>"}.
        '''
        v = version or _sys.version
        return _SPACE_(self, v.split(None, 1)[0])


class _Range(str):
    '''(INTERNAL) Extended C{str} for C{range} strings.
    '''
    def __call__(self, lo, hi, lopen=False, ropen=False,
                               prec=0, sep=_COMMASPACE_):
        '''Return the range as C{"(lo, hi)"}, C{"(lo, hi]"},
           C{"[lo, hi)"} or C{"[lo, hi]"}.
        '''
        from pygeodesy.streprs import Fmt
        r = NN(Fmt.f(lo, prec=prec), sep,
               Fmt.f(hi, prec=prec))
        f = (Fmt.PAREN if ropen else Fmt.LOPEN) if lopen else \
            (Fmt.ROPEN if ropen else Fmt.SQUARE)
        return f(r)


class _Slicer(str):
    '''(INTERNAL) String slicer C{.fromX} or C{.tillY}.
    '''
    def __getattr__(self, name):  # .fromX, .tillY
        n = len(name) - 4
        if n > 0:
            # assert len('till') == len(_from_) == 4
            if name.startswith('till'):
                i = self.find(name[4:])
                return self if i < 0 else _Slicer(self[:i + n])
            elif name.startswith(_from_):
                i = self.find(name[4:])
                return self if i < 0 else _Slicer(self[i:])
        return str.__getattr__(self, name)  # PYCHOK no cover


class MISSING(object):
    '''(INTERNAL) Singleton C{str}.
    '''
    def toRepr(self, **unused):
        return self.__class__.__name__

    __repr__ = __str__ = toStr = toRepr

MISSING          = MISSING()  # PYCHOK singleton
MISSING.__name__ = str(MISSING)

# __DUNDER__-style names would get mangled in classes
_0_                   = '0'                  # PYCHOK 'zero'
_0to9_                = '0123456789'         # PYCHOK OK
_1_                   = '1'                  # PYCHOK OK
_2_                   = '2'                  # PYCHOK OK
_3_                   = '3'                  # PYCHOK OK
_4_                   = '4'                  # PYCHOK OK
_a_                   = 'a'                  # PYCHOK OK
_A_                   = 'A'                  # PYCHOK OK
_a12_                 = 'a12'                # PYCHOK OK
_area_                = 'area'               # PYCHOK OK
_Airy1830_            = 'Airy1830'           # PYCHOK OK
_AiryModified_        = 'AiryModified'       # PYCHOK OK
_ambiguous_           = 'ambiguous'          # PYCHOK OK
_AMPERSAND_      = Str_('&')                 # PYCHOK OK
_an_                  = 'an'                 # PYCHOK OK
_and_                 = 'and'                # PYCHOK OK
# _AND_               = _AMPERSAND_          # PYCHOK OK
_angle_               = 'angle'              # PYCHOK OK
_antipodal_           = 'antipodal'          # PYCHOK OK
_areaOf_              = 'areaOf'             # PYCHOK OK
_arg_                 = 'arg'                # PYCHOK OK
_at_                  = 'at'                 # PYCHOK OK
_AT_             = Str_('@')                 # PYCHOK OK
_AtoZnoIO_    = _Slicer('ABCDEFGHJKLMNPQRSTUVWXYZ')  # PYCHOK in .gars, .mgrs and .wgrs
_attribute_           = 'attribute'          # PYCHOK OK
_azi1_                = 'azi1'               # PYCHOK OK
_azi12_               = 'azi12'              # PYCHOK OK
_azi2_                = 'azi2'               # PYCHOK OK
_azimuth_             = 'azimuth'            # PYCHOK OK
_b_                   = 'b'                  # PYCHOK OK
_B_                   = 'B'                  # PYCHOK OK
_BACKSLASH_      = Str_('\\')                # PYCHOK OK
_band_                = 'band'               # PYCHOK OK
_BANG_           = Str_('!')                 # PYCHOK OK
_BAR_            = Str_('|')                 # PYCHOK OK
_bearing_             = 'bearing'            # PYCHOK OK
_Bessel1841_          = 'Bessel1841'         # PYCHOK OK
_beta_                = 'beta'               # PYCHOK OK
_by_                  = 'by'                 # PYCHOK OK
_c_                   = 'c'                  # PYCHOK OK
_C_                   = 'C'                  # PYCHOK OK
_cartesian_           = 'cartesian'          # PYCHOK OK
_center_              = 'center'             # PYCHOK OK
# _CIRCUMFLEX_   = Str_('^')                 # PYCHOK OK
_Clarke1866_          = 'Clarke1866'         # PYCHOK OK
_Clarke1880IGN_       = 'Clarke1880IGN'      # PYCHOK OK
_clip_                = 'clip'               # PYCHOK OK
_clipid_              = 'clipid'             # PYCHOK OK
_coincident_          = 'coincident'         # PYCHOK OK
_colinear_            = 'colinear'           # PYCHOK OK
_COLON_          = Str_(':')                 # PYCHOK OK
_COLONSPACE_     = Str_(': ')                # PYCHOK OK
_COMMA_          = Str_(',')                 # PYCHOK OK
_COMMASPACE_     = Str_(_COMMASPACE_)        # PYCHOK OK
_composite_           = 'composite'          # PYCHOK OK
_concentric_          = 'concentric'         # PYCHOK OK
_convergence_ = _Prefix('convergence')       # PYCHOK OK
_conversion_          = 'conversion'         # PYCHOK OK
_convex_              = 'convex'             # PYCHOK OK
_cubic_               = 'cubic'              # PYCHOK OK
_d_                   = 'd'                  # PYCHOK OK
_D_                   = 'D'                  # PYCHOK OK
_DASH_           = Str_('-')                 # PYCHOK == _MINUS_
_datum_               = 'datum'              # PYCHOK OK
_decode3_             = 'decode3'            # PYCHOK OK
_deg_                 = 'deg'                # PYCHOK OK
_degrees_             = 'degrees'            # PYCHOK OK
_degrees2_            = 'degrees2'           # PYCHOK SQUARED
_delta_               = 'delta'              # PYCHOK OK
_DEPRECATED_          = 'DEPRECATED'         # PYCHOK OK
_DEQUALSPACED_   = Str_(' == ')              # PYCHOK OK
_distance_            = 'distance'           # PYCHOK OK
_distant_     = _Prefix('distant')           # PYCHOK OK
_doesn_t_exist_       = "doesn't exist"      # PYCHOK OK
_DOT_            = Str_('.')                 # PYCHOK OK
_down_                = 'down'               # PYCHOK OK
_e_                   = 'e'                  # PYCHOK OK
_E_                   = 'E'                  # PYCHOK OK
_earth_               = 'earth'              # PYCHOK OK
_east_                = 'east'               # PYCHOK OK
_easting_             = 'easting'            # PYCHOK OK
_ecef_                = 'ecef'               # PYCHOK OK
_edge_                = 'edge'               # PYCHOK OK
_elevation_           = 'elevation'          # PYCHOK OK
_ELLIPSIS_       = Str_('...')               # PYCHOK OK
_ELLIPSIS4_      = Str_('....')              # PYCHOK OK
# _ELLIPSISPACED_ = Str_(' ... ')            # PYCHOK <https://www.ThePunctuationGuide.com/ellipses.html>
_ellipsoid_           = 'ellipsoid'          # PYCHOK OK
_ellipsoidal_         = 'ellipsoidal'        # PYCHOK OK
_enabled_             = 'enabled'            # PYCHOK OK
_encode_              = 'encode'             # PYCHOK OK
_end_                 = 'end'                # PYCHOK OK
_epoch_               = 'epoch'              # PYCHOK OK
_eps_                 = 'eps'                # PYCHOK OK
_EQUAL_          = Str_('=')                 # PYCHOK OK
_EQUALSPACED_    = Str_(' = ')               # PYCHOK OK
_Error_               = 'Error'              # PYCHOK OK
_exceed_PI_radians_   = 'exceed PI radians'  # PYCHOK OK
_exceeds_     = _Prefix('exceeds')           # PYCHOK OK
# _EXCLAMATION_       = _BANG_               # PYCHOK OK
_exists_              = 'exists'             # PYCHOK OK
_f_                   = 'f'                  # PYCHOK OK
_F_                   = 'F'                  # PYCHOK OK
_feet_                = 'feet'               # PYCHOK OK
_few_                 = 'few'                # PYCHOK OK
_fi_                  = 'fi'                 # PYCHOK OK
_finite_              = 'finite'             # PYCHOK OK
_from_                = 'from'               # PYCHOK OK
_g_                   = 'g'                  # PYCHOK OK
_gamma_               = 'gamma'              # PYCHOK OK
_GRS80_               = 'GRS80'              # PYCHOK OK
_h_                   = 'h'                  # PYCHOK OK
_H_                   = 'H'                  # PYCHOK OK
_HASH_                = '#'                  # PYCHOK OK
_height_              = 'height'             # PYCHOK OK
_hemipole_            = 'hemipole'           # PYCHOK OK
_i_                   = 'i'                  # PYCHOK OK
_iadd_op_             = '+='                 # PYCHOK OK
_immutable_           = 'immutable'          # PYCHOK OK
_in_                  = 'in'                 # PYCHOK OK
_incompatible_        = 'incompatible'       # PYCHOK OK
_INF_                 = 'INF'                # PYCHOK OK
_infinite_            = 'infinite'           # PYCHOK _not_finite_
_initial_             = 'initial'            # PYCHOK OK
_inside_              = 'inside'             # PYCHOK OK
_insufficient_        = 'insufficient'       # PYCHOK OK
_intersection_        = 'intersection'       # PYCHOK OK
_Intl1924_            = 'Intl1924'           # PYCHOK OK
_invalid_             = 'invalid'            # PYCHOK OK
_invokation_          = 'invokation'         # PYCHOK OK
_isclockwise_         = 'isclockwise'        # PYCHOK OK
_ispolar_             = 'ispolar'            # PYCHOK OK
_j_                   = 'j'                  # PYCHOK OK
_k0_                  = 'k0'                 # PYCHOK OK
_kind_                = 'kind'               # PYCHOK OK
_knots_               = 'knots'              # PYCHOK OK
_Krassovski1940_      = 'Krassovski1940'     # PYCHOK OK
_Krassowsky1940_      = 'Krassowsky1940'     # PYCHOK OK
_LANGLE_              = '<'                  # PYCHOK OK
_lam_                 = 'lam'                # PYCHOK OK
_lat_                 = 'lat'                # PYCHOK OK
_lat0_                = 'lat0'               # PYCHOK OK
_lat1_                = 'lat1'               # PYCHOK OK
_lat2_                = 'lat2'               # PYCHOK OK
_latlon_              = 'latlon'             # PYCHOK OK
_LatLon_              = 'LatLon'             # PYCHOK OK
_LCURLY_              = '{'                  # PYCHOK LBRACE
_len_                 = 'len'                # PYCHOK OK
_limit_               = 'limit'              # PYCHOK OK
_line_                = 'line'               # PYCHOK OK
_linear_              = 'linear'             # PYCHOK OK
_LPAREN_              = '('                  # PYCHOK OK
_lon_                 = 'lon'                # PYCHOK OK
_lon0_                = 'lon0'               # PYCHOK OK
_lon1_                = 'lon1'               # PYCHOK OK
_lon2_                = 'lon2'               # PYCHOK OK
_low_                 = 'low'                # PYCHOK OK
_LSQUARE_             = '['                  # PYCHOK LBRACK
_ltp_                 = 'ltp'                # PYCHOK OK
_m_                   = 'm'                  # PYCHOK OK
_M_                   = 'M'                  # PYCHOK OK
_m12_                 = 'm12'                # PYCHOK OK
_M12_                 = 'M12'                # PYCHOK OK
_M21_                 = 'M21'                # PYCHOK OK
_MANT_DIG_            = 'MANT_DIG'           # PYCHOK OK
_MAX_                 = 'MAX'                # PYCHOK OK
_mean_                = 'mean'               # PYCHOK OK
_meanOf_              = 'meanOf'             # PYCHOK OK
_meridional_          = 'meridional'         # PYCHOK OK
_meter_               = 'meter'              # PYCHOK OK
_meter2_              = 'meter2'             # PYCHOK SQUARED
_MGRS_                = 'MGRS'               # PYCHOK OK
_MIN_                 = 'MIN'                # PYCHOK OK
_MINUS_               = _DASH_               # PYCHOK OK
_module_              = 'module'             # PYCHOK OK
_n_                   = 'n'                  # PYCHOK OK
_N_                   = 'N'                  # PYCHOK OK
_n_a_                 = 'n/a'                # PYCHOK OK
_N_A_                 = 'N/A'                # PYCHOK OK
_NAD27_               = 'NAD27'              # PYCHOK OK
_NAD83_               = 'NAD83'              # PYCHOK OK
_name_                = 'name'               # PYCHOK OK
_NAN_                 = 'NAN'                # PYCHOK OK
_near_          = _Dash('near')              # PYCHOK OK
_nearestOn2_          = 'nearestOn2'         # PYCHOK OK
_negative_            = 'negative'           # PYCHOK OK
_NL_             = Str_('\n')                # PYCHOK OK
_NLATvar_        = Str_(_NL_ + '@var ')      # PYCHOK OK
_NLHASH_         = Str_(_NL_ + '# ')         # PYCHOK OK
# _NLNL_              = _DNL_                # PYCHOK OK
_NN_                  = 'NN'                 # PYCHOK OK
_no_          = _Prefix('no')                # PYCHOK OK
_north_               = 'north'              # PYCHOK OK
_northing_            = 'northing'           # PYCHOK OK
_NorthPole_           = 'NorthPole'          # PYCHOK OK
_not_         = _Prefix('not')               # PYCHOK OK
_NOTEQUAL_            = _BANG_ + _EQUAL_     # PYCHOK OK
_not_finite_          = 'not finite'         # PYCHOK _not_(_finite_), _infinite_
_not_scalar_          = 'not scalar'         # PYCHOK _not_(_scalar_)
_NTF_                 = 'NTF'                # PYCHOK OK
_null_                = 'null'               # PYCHOK OK
_number_              = 'number'             # PYCHOK OK
_numpy_               = 'numpy'              # PYCHOK OK
_Nv00_                = 'Nv00'               # PYCHOK OK
_of_                  = 'of'                 # PYCHOK OK
_on_                  = 'on'                 # PYCHOK OK
_opposite_            = 'opposite'           # PYCHOK OK
_or_                  = 'or'                 # PYCHOK OK
_other_               = 'other'              # PYCHOK OK
_outside_             = 'outside'            # PYCHOK OK
_overlap_             = 'overlap'            # PYCHOK OK
_parallel_            = 'parallel'           # PYCHOK OK
_PERCENT_             = '%'                  # PYCHOK OK
_PERCENTDOTSTAR_      = '%.*'                # PYCHOK _DOT_(_PERCENT_, _STAR_)
_phi_                 = 'phi'                # PYCHOK OK
_PLUS_           = Str_('+')                 # PYCHOK OK
_PLUSMINUS_           = _PLUS_ + _MINUS_     # PYCHOK OK
_point_               = 'point'              # PYCHOK OK
_points_              = 'points'             # PYCHOK OK
_pole_                = 'pole'               # PYCHOK OK
_precision_           = 'precision'          # PYCHOK OK
_prime_vertical_      = 'prime_vertical'     # PYCHOK OK
_pygeodesy_           = 'pygeodesy'          # PYCHOK OK
_pygeodesy_abspath_   = 'pygeodesy_abspath'  # PYCHOK OK
_PyPy__       = _PyPy__('PyPy ')             # PYCHOK + _SPACE_
_Python_     = _Python_('Python')            # PYCHOK singleton
_python_              = 'python'             # PYCHOK OK
_QUOTE1_              = "'"                  # PYCHOK OK
_QUOTE2_              = '"'                  # PYCHOK OK
# _QUOTE3_            = "'''"                # PYCHOK OK
# _QUOTE6_            = '"""'                # PYCHOK OK
_R_                   = 'R'                  # PYCHOK OK
_radians_             = 'radians'            # PYCHOK OK
_radians2_            = 'radians2'           # PYCHOK SQUARED
_radius_              = 'radius'             # PYCHOK OK
_radius1_             = 'radius1'            # PYCHOK OK
_radius2_             = 'radius2'            # PYCHOK OK
_range_        = _Range('range')             # PYCHOK OK
_RANGLE_              = '>'                  # PYCHOK OK
_RCURLY_              = '}'                  # PYCHOK RBRACE
_reciprocal_          = 'reciprocal'         # PYCHOK OK
_reframe_             = 'reframe'            # PYCHOK OK
_resolution_          = 'resolution'         # PYCHOK OK
_rIn_                 = 'rIn'                # PYCHOK OK
_RPAREN_              = ')'                  # PYCHOK OK
_RSQUARE_             = ']'                  # PYCHOK RBRACK
_s_                   = 's'                  # PYCHOK OK
_S_                   = 'S'                  # PYCHOK OK
_s12_                 = 's12'                # PYCHOK OK
_S12_                 = 'S12'                # PYCHOK OK
_scalar_              = 'scalar'             # PYCHOK OK
_scale_               = 'scale'              # PYCHOK OK
_scale0_              = 'scale0'             # PYCHOK OK
_scipy_               = 'scipy'              # PYCHOK OK
_semi_circular_       = 'semi-circular'      # PYCHOK OK
_sep_                 = 'sep'                # PYCHOK OK
_singular_            = 'singular'           # PYCHOK OK
_SLASH_          = Str_('/')                 # PYCHOK OK
_small_               = 'small'              # PYCHOK OK
_Sphere_              = 'Sphere'             # PYCHOK OK
_spherical_           = 'spherical'          # PYCHOK OK
_SouthPole_           = 'SouthPole'          # PYCHOK OK
_SPACE_          = Str_(' ')                 # PYCHOK OK
_specified_           = 'specified'          # PYCHOK OK
_STAR_           = Str_('*')                 # PYCHOK OK
_start_               = 'start'              # PYCHOK OK
_std_                 = 'std'                # PYCHOK OK
_stdev_               = 'stdev'              # PYCHOK OK
_supported_           = 'supported'          # PYCHOK OK
_sx_                  = 'sx'                 # PYCHOK OK
_sy_                  = 'sy'                 # PYCHOK OK
_sz_                  = 'sz'                 # PYCHOK OK
_tbd_                 = 'tbd'                # PYCHOK OK
_TILDE_               = '~'                  # PYCHOK OK
_to_                  = 'to'                 # PYCHOK OK
_tolerance_   = _Prefix('tolerance')         # PYCHOK OK
_too_         = _Prefix('too')               # PYCHOK OK
_transform_           = 'transform'          # PYCHOK OK
_tx_                  = 'tx'                 # PYCHOK OK
_ty_                  = 'ty'                 # PYCHOK OK
_tz_                  = 'tz'                 # PYCHOK OK
_UNDER_          = Str_('_')                 # PYCHOK OK
_units_               = 'units'              # PYCHOK OK
_UNUSED_              = 'UNUSED'             # PYCHOK OK
_up_                  = 'up'                 # PYCHOK OK
_UPS_                 = 'UPS'                # PYCHOK OK
_utf_8_               = 'utf-8'              # PYCHOK OK
_UTM_                 = 'UTM'                # PYCHOK OK
_V_                   = 'V'                  # PYCHOK OK
_valid_               = 'valid'              # PYCHOK OK
_value_               = 'value'              # PYCHOK OK
_version_             = 'version'            # PYCHOK OK
_vs_                  = 'vs'                 # PYCHOK OK
_W_                   = 'W'                  # PYCHOK OK
_WGS72_               = 'WGS72'              # PYCHOK OK
_WGS84_               = 'WGS84'              # PYCHOK OK
_width_               = 'width'              # PYCHOK OK
_with_                = 'with'               # PYCHOK OK
_x_                   = 'x'                  # PYCHOK OK
_X_                   = 'X'                  # PYCHOK OK
_xyz_                 = 'xyz'                # PYCHOK OK
_y_                   = 'y'                  # PYCHOK OK
_Y_                   = 'Y'                  # PYCHOK OK
_z_                   = 'z'                  # PYCHOK OK
_Z_                   = 'Z'                  # PYCHOK OK
_zone_                = 'zone'               # PYCHOK OK

_EW_       = _E_  + _W_   # PYCHOK common cardinals
_NE_       = _N_  + _E_   # PYCHOK positive ones
_NS_       = _N_  + _S_   # PYCHOK OK
_NSEW_     = _NS_ + _EW_  # PYCHOK OK
_NW_       = _N_  + _W_   # PYCHOK OK
_SE_       = _S_  + _E_   # PYCHOK OK
_SW_       = _S_  + _W_   # PYCHOK negative ones
# _NESW_   = _NE_ + _SW_  # PYCHOK clockwise

_DDOT_     = Str_(_DOT_   * 2)  # PYCHOK OK
# _DEQUAL_ = Str_(_EQUAL_ * 2)  # PYCHOK OK
_DNL_      = Str_(_NL_    * 2)  # PYCHOK OK
# _DSLASH_ = Str_(_SLASH_ * 2)  # PYCHOK OK
# _DSTAR_  = Str_(_STAR_  * 2)  # PYCHOK OK
_DUNDER_   = Str_(_UNDER_ * 2)  # PYCHOK OK

_LR_PAIRS  = {_LANGLE_:  _RANGLE_,
              _LCURLY_:  _RCURLY_,
              _LPAREN_:  _RPAREN_,
              _LSQUARE_: _RSQUARE_}  # PYCHOK OK


def _dunder_nameof(inst, *dflt):
    '''(INTERNAL) Get the double_underscore __name__ attr.
    '''
    try:
        return inst.__name__
    except AttributeError:
        pass
    return dflt[0] if dflt else inst.__class__.__name__


def _enquote(strs, quote=_QUOTE2_):  # in .basics, .solveBase
    '''(INTERNAL) Enquote a string containing whitespace.
    '''
    if len(strs.split()) > 1:
        strs = NN(quote, strs, quote)
    return strs


def _is(a, b):  # PYCHOK no cover
    '''(INTERNAL) C{a is b}? in C{PyPy}
    '''
    return (a == b) if _isPyPy() else (a is b)


def _isPyPy():
    '''(INTERNAL) Is this C{PyPy}?
    '''
    # platform.python_implementation() == 'PyPy'
    return _pythonarchine()[0].startswith(_PyPy__)


def _load_lib(name):
    '''(INTERNAL) Load a C{dylib}, B{C{name}} must startwith('lib').
    '''
    # macOS 11+ (aka 10.16) no longer provides direct loading of
    # system libraries.  As a result, C{ctypes.util.find_library}
    # will not find any library, unless previously installed by a
    # low-level dlopen(name) call (with the library base C{name}).
    from ctypes import CDLL
    from ctypes.util import find_library

    ns = find_library(name), name
    if _sys.platform[:6] == 'darwin':  # and os.name == 'posix'
        from ctypes import _dlopen, DEFAULT_MODE
        from os.path import join
        ns += (_DOT_(name, 'dylib'),
               _DOT_(name, 'framework'), join(
               _DOT_(name, 'framework'), name))
    else:  # non-macOS
        DEFAULT_MODE = 0

        def _dlopen(*unused):
            return True

    for n in ns:
        try:
            if n and _dlopen(n, DEFAULT_MODE):  # pre-load handle
                lib = CDLL(n)  # == ctypes.cdll.LoadLibrary(n)
                if lib._name:  # has a qualified name
                    return lib
        except (AttributeError, OSError):
            pass

    return None  # raise OSError


def machine():
    '''Return standard C{platform.machine}, but distinguishing Intel I{native}
       from Intel I{emulation} on Apple Silicon (on macOS only).

       @return: Machine C{'arm64'} for Apple Silicon I{native}, C{'x86_64'}
                for Intel I{native}, C{"arm64_x86_64"} for Intel I{emulation},
                etc. (C{str} with C{comma}s replaced by C{underscore}s).
    '''
    return _platform2()[1]


def _platform2(sep=NN):
    '''(INTERNAL) Get platform architecture and machine as C{2-list} or C{str}.
    '''
    L = _pf2List
    if not L:
        import platform
        m = platform.machine()  # ARM64, arm64, x86_64, iPhone13,2, etc.
        m = m.replace(_COMMA_, _UNDER_)
        if m.lower() == 'x86_64':  # on Intel or Rosetta2 ...
            v = platform.mac_ver()  # ... and only on macOS ...
            if v and _version2(v[0]) > (10, 15):  # ... 11+ aka 10.16
                # <https://Developer.Apple.com/forums/thread/659846>
                if _sysctl_uint('sysctl.proc_translated') == 1:  # and \
#                  _sysctl_uint('hw.optional.arm64') == 1:  # PYCHOK indent
                    m = _UNDER_('arm64', m)  # Apple Si emulating Intel x86-64
        L[:] = [platform.architecture()[0],  # bits
                m]  # arm64, arm64_x86_64, x86_64, etc.
    return sep.join(L) if sep else L  # 2-list()


def _pythonarchine(sep=NN):  # in test/bases.py versions
    '''(INTERNAL) Get PyPy and Python versions and C{_platform2} as C{3- or 4-list} or C{str}.
    '''
    L = _Py3List
    if not L:
        v    = _sys.version
        L[:] = [_Python_(v)] + _platform2()
        pypy = _PyPy__(v)
        if pypy:
            L.insert(0, pypy)
    return sep.join(L) if sep else L  # 3- or 4-list


def _sysctl_uint(name):
    '''(INTERNAL) Get an unsigned int sysctl item by name, use on macOS ONLY!
    '''
    libc = _load_lib('libc')
    if libc:  # <https://StackOverflow.com/questions/759892/python-ctypes-and-sysctl>
        from ctypes import byref, c_char_p, c_size_t, c_uint, sizeof  # get_errno
        n = name if str is bytes else bytes(name, _utf_8_)  # PYCHOK isPython2 = str is bytes
        u = c_uint(0)
        z = c_size_t(sizeof(u))
        r = libc.sysctlbyname(c_char_p(n), byref(u), byref(z), None, c_size_t(0))
    else:  # could find or load 'libc'
        r = -2
    return int(r if r else u.value)  # -1 ENOENT error, -2 no libc


def _under(name):  # PYCHOK in .datums, .auxilats, .ups, .utm, .utmupsBase, ...
    '''(INTERNAL) Prefix C{name} with I{underscore}.
    '''
    return name if name.startswith(_UNDER_) else NN(_UNDER_, name)


def _usage(file_py, *args):  # in .etm
    '''(INTERNAL) Build "usage: python -m ..." cmd line for module B{C{file_py}}.
    '''
    import os
    m = os.path.dirname(file_py).replace(os.getcwd(), _ELLIPSIS_) \
                                .replace(os.sep, _DOT_).strip()
    b, x = os.path.splitext(os.path.basename(file_py))
    if x == '.py' and b != '__main__':
        m = _DOT_(m or _pygeodesy_, b)
    p = NN(_python_, _sys.version_info[0])
    return NN('usage', _SPACE_(_COLON_, p, '-m', _enquote(m), *args))


def _version2(version, n=2):
    '''(INTERNAL) Split C{B{version} str} into a C{1-, 2- or 3-tuple} of C{int}s.
    '''
    def _ints(vs):
        for v in vs:
            if v:
                try:
                    yield int(v.strip())
                except (TypeError, ValueError):
                    pass

    t = tuple(_ints(version.split(_DOT_, 2)))
    if len(t) < n:
        t += (0,) * n
    return tuple(t[:n])


__all__ = (_NN_,  # not MISSING!
            Str_.__name__,  # classes
            machine.__name__)  # in .lazily
__version__ = '23.09.13'

if __name__ == '__main__':

    from pygeodesy import itemsorted, printf

    t = b = 0
    for n, v in itemsorted(locals(), asorted=False, reverse=True):
        if n.endswith(_UNDER_) and n.startswith(_UNDER_) and \
                               not n.startswith(_DUNDER_):
            t += 1
            b += len(v)
            m  = n[1:-1]
            if m != v and m.replace(_UNDER_, _SPACE_) != v:
                printf('%4d: %s = %r', t, n, v)
    n = len(locals())
    printf('%4d (%d) names, %s chars total, %.2f chars avg', t, n, b, float(b) / t, nl=1)

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
