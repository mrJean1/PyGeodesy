# -*- coding: utf-8 -*-

u'''Single C{str}ing constants, C{intern}'ed across C{pygeodesy}
modules and function L{pygeodesy.machine}.
'''
_COMMASPACE_ = ', '   # overriden below
_pl2List     =  None  # cached _platform2 lists
_Py3List     =  None  # cached _pythonarchine lists


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


class _Prefix(Str_):
    '''(INTERNAL) Extended C{str} for prefix.
    '''
    def __call__(self, *args):
        '''Join C{self} plus all B{C{args}} like C{" ".join((self,) + B{args})}.
        '''
        return _SPACE_.join_(self, *args)  # re-callable


class _Python_(str):  # overwritten by singleton below
    '''(INTERNAL) Extended C{str} for C{"Python"} and version.
    '''
    def __call__(self, version=None):
        '''Return C{"Python <version>"}.
        '''
        if not version:
            from sys import version
        return _SPACE_(self, version.split(None, 1)[0])


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
        if lopen:
            r = Fmt.PAREN(r) if ropen else Fmt.LOPEN(r)
        else:
            r = Fmt.ROPEN(r) if ropen else Fmt.SQUARE(r)
        return r


class _Slicer(str):
    '''(INTERNAL) String slicer C{.fromX} or C{.tillY}.
    '''
    def __getattr__(self, name):  # .fromX, .tillY
        if name.startswith(_till_):
            i = self.find(name[len(_till_):])
            if 0 < i < len(self):
                return _Slicer(self[:i + 1])
        elif name.startswith(_from_):
            i = self.find(name[len(_from_):])
            if 0 < (i + 1) < len(self):
                return _Slicer(self[i:])
        else:
            return getattr(str, name)
        return self


class MISSING(object):
    '''(INTERNAL) Singleton C{str}.
    '''
    def toRepr(self, **unused):
        return self.__class__.__name__

    __repr__ = toRepr
    __str__  = toRepr
    toStr    = toRepr

MISSING          = MISSING()  # PYCHOK singleton
MISSING.__name__ = str(MISSING)

NN = Str_('')  # Nomen Nescio <https://Wiktionary.org/wiki/N.N.>

# __DUNDER__-style names would get mangled in classes
_0_                   = '0'                  # PYCHOK 'zero'
_0to9_                = '0123456789'         # PYCHOK expected
_1_                   = '1'                  # PYCHOK expected
_2_                   = '2'                  # PYCHOK expected
_3_                   = '3'                  # PYCHOK expected
_4_                   = '4'                  # PYCHOK expected
_a_                   = 'a'                  # PYCHOK expected
_A_                   = 'A'                  # PYCHOK expected
_a12_                 = 'a12'                # PYCHOK expected
_area_                = 'area'               # PYCHOK expected
_Airy1830_            = 'Airy1830'           # PYCHOK expected
_AiryModified_        = 'AiryModified'       # PYCHOK expected
_ambiguous_           = 'ambiguous'          # PYCHOK expected
_AMPERSAND_      = Str_('&')                 # PYCHOK expected
_an_                  = 'an'                 # PYCHOK expected
_and_                 = 'and'                # PYCHOK expected
# _AND_               = _AMPERSAND_          # PYCHOK expected
_angle_               = 'angle'              # PYCHOK expected
_antipodal_           = 'antipodal'          # PYCHOK expected
_areaOf_              = 'areaOf'             # PYCHOK expected
_arg_                 = 'arg'                # PYCHOK expected
_at_                  = 'at'                 # PYCHOK expected
_AT_             = Str_('@')                 # PYCHOK expected
_AtoZnoIO_    = _Slicer('ABCDEFGHJKLMNPQRSTUVWXYZ')  # PYCHOK in C{gars}, C{mgrs} and C{wgrs}
_attribute_           = 'attribute'          # PYCHOK expected
_azi1_                = 'azi1'               # PYCHOK expected
_azi12_               = 'azi12'              # PYCHOK expected
_azi2_                = 'azi2'               # PYCHOK expected
_azimuth_             = 'azimuth'            # PYCHOK expected
_b_                   = 'b'                  # PYCHOK expected
_B_                   = 'B'                  # PYCHOK expected
_BACKSLASH_      = Str_('\\')                # PYCHOK expected
_band_                = 'band'               # PYCHOK expected
_BANG_           = Str_('!')                 # PYCHOK expected
_BAR_            = Str_('|')                 # PYCHOK expected
_bearing_             = 'bearing'            # PYCHOK expected
_Bessel1841_          = 'Bessel1841'         # PYCHOK expected
_by_                  = 'by'                 # PYCHOK expected
_c_                   = 'c'                  # PYCHOK expected
_C_                   = 'C'                  # PYCHOK expected
_cartesian_           = 'cartesian'          # PYCHOK expected
_center_              = 'center'             # PYCHOK expected
# _CIRCUMFLEX_        = '^'                  # PYCHOK expected
_Clarke1866_          = 'Clarke1866'         # PYCHOK expected
_Clarke1880IGN_       = 'Clarke1880IGN'      # PYCHOK expected
_coincident_          = 'coincident'         # PYCHOK expected
_colinear_            = 'colinear'           # PYCHOK expected
_COLON_          = Str_(':')                 # PYCHOK expected
_COLONSPACE_     = Str_(': ')                # PYCHOK expected
_COMMA_          = Str_(',')                 # PYCHOK expected
_COMMASPACE_     = Str_(_COMMASPACE_)        # PYCHOK expected
_concentric_          = 'concentric'         # PYCHOK expected
_convergence_ = _Prefix('convergence')       # PYCHOK expected
_conversion_          = 'conversion'         # PYCHOK expected
_convex_              = 'convex'             # PYCHOK expected
_cubic_               = 'cubic'              # PYCHOK expected
_d_                   = 'd'                  # PYCHOK expected
_D_                   = 'D'                  # PYCHOK expected
_DASH_           = Str_('-')                 # PYCHOK == _MINUS_
_datum_               = 'datum'              # PYCHOK expected
_decode3_             = 'decode3'            # PYCHOK expected
_deg_                 = 'deg'                # PYCHOK expected
_degrees_             = 'degrees'            # PYCHOK expected
_degrees2_            = 'degrees2'           # PYCHOK SQUARED
_DEPRECATED_          = 'DEPRECATED'         # PYCHOK expected
_DEQUALSPACED_   = Str_(' == ')              # PYCHOK expected
_distance_            = 'distance'           # PYCHOK expected
_distanceTo_          = 'distanceTo'         # PYCHOK expected
_distant_     = _Prefix('distant')           # PYCHOK expected
_doesn_t_exist_       = "doesn't exist"      # PYCHOK expected
_DOT_            = Str_('.')                 # PYCHOK expected
_down_                = 'down'               # PYCHOK expected
_e_                   = 'e'                  # PYCHOK expected
_E_                   = 'E'                  # PYCHOK expected
_east_                = 'east'               # PYCHOK expected
_easting_             = 'easting'            # PYCHOK expected
_ecef_                = 'ecef'               # PYCHOK expected
_edge_                = 'edge'               # PYCHOK expected
_elevation_           = 'elevation'          # PYCHOK expected
_ELLIPSIS_       = Str_('...')               # PYCHOK expected
_ELLIPSIS4_      = Str_('....')              # PYCHOK expected
# _ELLIPSISPACED_ = Str_(' ... ')            # PYCHOK <https://www.ThePunctuationGuide.com/ellipses.html>
_ellipsoid_           = 'ellipsoid'          # PYCHOK expected
_ellipsoidal_         = 'ellipsoidal'        # PYCHOK expected
_enabled_             = 'enabled'            # PYCHOK expected
_encode_              = 'encode'             # PYCHOK expected
_end_                 = 'end'                # PYCHOK expected
_epoch_               = 'epoch'              # PYCHOK expected
_EQUAL_          = Str_('=')                 # PYCHOK expected
_EQUALSPACED_    = Str_(' = ')               # PYCHOK expected
_exceed_PI_radians_   = 'exceed PI radians'  # PYCHOK expected
_exceeds_     = _Prefix('exceeds')           # PYCHOK expected
_exists_              = 'exists'             # PYCHOK expected
_f_                   = 'f'                  # PYCHOK expected
_F_                   = 'F'                  # PYCHOK expected
_feet_                = 'feet'               # PYCHOK expected
_few_                 = 'few'                # PYCHOK expected
_fi_                  = 'fi'                 # PYCHOK expected
_finite_              = 'finite'             # PYCHOK expected
_from_                = 'from'               # PYCHOK expected
_g_                   = 'g'                  # PYCHOK expected
_gamma_               = 'gamma'              # PYCHOK expected
_GRS80_               = 'GRS80'              # PYCHOK expected
_h_                   = 'h'                  # PYCHOK expected
_H_                   = 'H'                  # PYCHOK expected
_height_              = 'height'             # PYCHOK expected
_hemipole_            = 'hemipole'           # PYCHOK expected
_i_                   = 'i'                  # PYCHOK expected
_I_                   = 'I'                  # PYCHOK expected
_iadd_                = '+='                 # PYCHOK expected
_immutable_           = 'immutable'          # PYCHOK expected
_in_                  = 'in'                 # PYCHOK expected
_incompatible_        = 'incompatible'       # PYCHOK expected
_INF_                 = 'INF'                # PYCHOK expected
_infinite_            = 'infinite'           # PYCHOK _not_finite_
_initial_             = 'initial'            # PYCHOK expected
_inside_              = 'inside'             # PYCHOK expected
_intersection_        = 'intersection'       # PYCHOK expected
_Intl1924_            = 'Intl1924'           # PYCHOK expected
_invalid_             = 'invalid'            # PYCHOK expected
_isclockwise_         = 'isclockwise'        # PYCHOK expected
_ispolar_             = 'ispolar'            # PYCHOK expected
_j_                   = 'j'                  # PYCHOK expected
_k0_                  = 'k0'                 # PYCHOK expected
_kind_                = 'kind'               # PYCHOK expected
_knots_               = 'knots'              # PYCHOK expected
_Krassovski1940_      = 'Krassovski1940'     # PYCHOK expected
_Krassowsky1940_      = 'Krassowsky1940'     # PYCHOK expected
_LANGLE_              = '<'                  # PYCHOK expected
_lam_                 = 'lam'                # PYCHOK expected
_lat_                 = 'lat'                # PYCHOK expected
_lat0_                = 'lat0'               # PYCHOK expected
_lat1_                = 'lat1'               # PYCHOK expected
_lat2_                = 'lat2'               # PYCHOK expected
_latlon_              = 'latlon'             # PYCHOK expected
_LatLon_              = 'LatLon'             # PYCHOK expected
_LCURLY_              = '{'                  # PYCHOK LBRACE
_len_                 = 'len'                # PYCHOK expected
_linear_              = 'linear'             # PYCHOK expected
#_LPAREN_             = '('                  # PYCHOK expected
_lon_                 = 'lon'                # PYCHOK expected
_lon0_                = 'lon0'               # PYCHOK expected
_lon1_                = 'lon1'               # PYCHOK expected
_lon2_                = 'lon2'               # PYCHOK expected
#_LSQUARE_            = '['                  # PYCHOK LBRACK
_ltp_                 = 'ltp'                # PYCHOK expected
_m_                   = 'm'                  # PYCHOK expected
_M_                   = 'M'                  # PYCHOK expected
_m12_                 = 'm12'                # PYCHOK expected
_M12_                 = 'M12'                # PYCHOK expected
_M21_                 = 'M21'                # PYCHOK expected
_MANT_DIG_            = 'MANT_DIG'           # PYCHOK expected
_MAX_                 = 'MAX'                # PYCHOK expected
_mean_                = 'mean'               # PYCHOK expected
_meanOf_              = 'meanOf'             # PYCHOK expected
_meridional_          = 'meridional'         # PYCHOK expected
_meter_               = 'meter'              # PYCHOK expected
_meter2_              = 'meter2'             # PYCHOK SQUARED
_MGRS_                = 'MGRS'               # PYCHOK expected
_MIN_                 = 'MIN'                # PYCHOK expected
_MINUS_               = _DASH_               # PYCHOK expected
_module_              = 'module'             # PYCHOK expected
_n_                   = 'n'                  # PYCHOK expected
_N_                   = 'N'                  # PYCHOK expected
_n_a_                 = 'n/a'                # PYCHOK expected
_N_A_                 = 'N/A'                # PYCHOK expected
_NAD27_               = 'NAD27'              # PYCHOK expected
_NAD83_               = 'NAD83'              # PYCHOK expected
_name_                = 'name'               # PYCHOK expected
_NAN_                 = 'NAN'                # PYCHOK expected
_near_          = _Dash('near')              # PYCHOK expected
_nearestOn2_          = 'nearestOn2'         # PYCHOK expected
_negative_            = 'negative'           # PYCHOK expected
_NL_             = Str_('\n')                # PYCHOK expected
_NLATvar_        = Str_(_NL_ + '@var ')      # PYCHOK expected
_NLHASH_         = Str_(_NL_ + '# ')         # PYCHOK expected
# _NLNL_              = _DNL_                # PYCHOK expected
_NN_                  = 'NN'                 # PYCHOK expected
_no_          = _Prefix('no')                # PYCHOK expected
_north_               = 'north'              # PYCHOK expected
_northing_            = 'northing'           # PYCHOK expected
_NorthPole_           = 'NorthPole'          # PYCHOK expected
_not_         = _Prefix('not')               # PYCHOK expected
_not_finite_          = 'not finite'         # PYCHOK _not_(_finite_), _infinite_
_not_scalar_          = 'not scalar'         # PYCHOK _not_(_scalar_)
_NTF_                 = 'NTF'                # PYCHOK expected
_null_                = 'null'               # PYCHOK expected
_number_              = 'number'             # PYCHOK expected
_numpy_               = 'numpy'              # PYCHOK expected
_Nv00_                = 'Nv00'               # PYCHOK expected
_O_                   = 'O'                  # PYCHOK letter "Oh", not zero
_on_                  = 'on'                 # PYCHOK expected
_opposite_            = 'opposite'           # PYCHOK expected
_or_                  = 'or'                 # PYCHOK expected
_other_               = 'other'              # PYCHOK expected
_outside_             = 'outside'            # PYCHOK expected
_overlap_             = 'overlap'            # PYCHOK expected
_PERCENT_             = '%'                  # PYCHOK expected
_PERCENTDOTSTAR_      = '%.*'                # PYCHOK _DOT_(_PERCENT_, _STAR_)
_perimeterOf_         = 'perimeterOf'        # PYCHOK expected
_phi_                 = 'phi'                # PYCHOK expected
_PLUS_           = Str_('+')                 # PYCHOK expected
_PLUSMINUS_           = _PLUS_ + _MINUS_     # PYCHOK expected
_point_               = 'point'              # PYCHOK expected
_points_              = 'points'             # PYCHOK expected
_pole_                = 'pole'               # PYCHOK expected
_precision_           = 'precision'          # PYCHOK expected
_prime_vertical_      = 'prime_vertical'     # PYCHOK expected
_pygeodesy_abspath_   = 'pygeodesy_abspath'  # PYCHOK expected
_PyPy__               = 'PyPy '              # PYCHOK + _SPACE_
_Python_     = _Python_('Python')            # PYCHOK singleton
_python_              = 'python'             # PYCHOK expected
_QUOTE1_              = "'"                  # PYCHOK expected
_QUOTE2_              = '"'                  # PYCHOK expected
# _QUOTE3_            = "'''"                # PYCHOK expected
# _QUOTE6_            = '"""'                # PYCHOK expected
_R_                   = 'R'                  # PYCHOK expected
_radians_             = 'radians'            # PYCHOK expected
_radians2_            = 'radians2'           # PYCHOK SQUARED
_radius_              = 'radius'             # PYCHOK expected
_radius1_             = 'radius1'            # PYCHOK expected
_radius2_             = 'radius2'            # PYCHOK expected
_range_        = _Range('range')             # PYCHOK expected
_RANGLE_              = '>'                  # PYCHOK expected
_RCURLY_              = '}'                  # PYCHOK RBRACE
_reciprocal_          = 'reciprocal'         # PYCHOK expected
_reframe_             = 'reframe'            # PYCHOK expected
_resolution_          = 'resolution'         # PYCHOK expected
_rIn_                 = 'rIn'                # PYCHOK expected
#_RPAREN_             = ')'                  # PYCHOK expected
#_RSQUARE_            = ']'                  # PYCHOK RBRACK
_s_                   = 's'                  # PYCHOK expected
_S_                   = 'S'                  # PYCHOK expected
_s12_                 = 's12'                # PYCHOK expected
_S12_                 = 'S12'                # PYCHOK expected
_scalar_              = 'scalar'             # PYCHOK expected
_scale_               = 'scale'              # PYCHOK expected
_scipy_               = 'scipy'              # PYCHOK expected
_semi_circular_       = 'semi-circular'      # PYCHOK expected
_sep_                 = 'sep'                # PYCHOK expected
_sets_                = 'sets'               # PYCHOK expected
_singular_            = 'singular'           # PYCHOK expected
_SLASH_          = Str_('/')                 # PYCHOK expected
_small_               = 'small'              # PYCHOK expected
_Sphere_              = 'Sphere'             # PYCHOK expected
_spherical_           = 'spherical'          # PYCHOK expected
_SouthPole_           = 'SouthPole'          # PYCHOK expected
_SPACE_          = Str_(' ')                 # PYCHOK expected
_specified_           = 'specified'          # PYCHOK expected
_STAR_           = Str_('*')                 # PYCHOK expected
_start_               = 'start'              # PYCHOK expected
_std_                 = 'std'                # PYCHOK expected
_stdev_               = 'stdev'              # PYCHOK expected
_supported_           = 'supported'          # PYCHOK expected
_sx_                  = 'sx'                 # PYCHOK expected
_sy_                  = 'sy'                 # PYCHOK expected
_sz_                  = 'sz'                 # PYCHOK expected
_tbd_                 = 'tbd'                # PYCHOK expected
_TILDE_               = '~'                  # PYCHOK expected
_till_                = 'till'               # PYCHOK expected
_to_                  = 'to'                 # PYCHOK expected
_tolerance_   = _Prefix('tolerance')         # PYCHOK expected
_too_         = _Prefix('too')               # PYCHOK expected
_transform_           = 'transform'          # PYCHOK expected
_tx_                  = 'tx'                 # PYCHOK expected
_ty_                  = 'ty'                 # PYCHOK expected
_tz_                  = 'tz'                 # PYCHOK expected
_UNDER_          = Str_('_')                 # PYCHOK expected
_units_               = 'units'              # PYCHOK expected
_up_                  = 'up'                 # PYCHOK expected
_UPS_                 = 'UPS'                # PYCHOK expected
_utf_8_               = 'utf-8'              # PYCHOK expected
_UTM_                 = 'UTM'                # PYCHOK expected
_V_                   = 'V'                  # PYCHOK expected
_valid_               = 'valid'              # PYCHOK expected
_value_               = 'value'              # PYCHOK expected
_version_             = 'version'            # PYCHOK expected
_vs_                  = 'vs'                 # PYCHOK expected
_W_                   = 'W'                  # PYCHOK expected
_WGS72_               = 'WGS72'              # PYCHOK expected
_WGS84_               = 'WGS84'              # PYCHOK expected
_width_               = 'width'              # PYCHOK expected
_x_                   = 'x'                  # PYCHOK expected
_X_                   = 'X'                  # PYCHOK expected
_xyz_                 = 'xyz'                # PYCHOK expected
_y_                   = 'y'                  # PYCHOK expected
_Y_                   = 'Y'                  # PYCHOK expected
_z_                   = 'z'                  # PYCHOK expected
_Z_                   = 'Z'                  # PYCHOK expected
_zone_                = 'zone'               # PYCHOK expected

_EW_                  = _E_  + _W_           # PYCHOK common cardinals
_NE_                  = _N_  + _E_           # PYCHOK positive ones
_NS_                  = _N_  + _S_           # PYCHOK expected
_NSEW_                = _NS_ + _EW_          # PYCHOK expected
_NW_                  = _N_  + _W_           # PYCHOK expected
_SE_                  = _S_  + _E_           # PYCHOK expected
_SW_                  = _S_  + _W_           # PYCHOK negative ones
# _NESW_              = _NE_ + _SW_          # PYCHOK clockwise

_DDOT_           = Str_(_DOT_   * 2)         # PYCHOK expected
# _DEQUAL_       = Str_(_EQUAL_ * 2)         # PYCHOK expected
_DNL_            = Str_(_NL_    * 2)         # PYCHOK expected
# _DSLASH_       = Str_(_SLASH_ * 2)         # PYCHOK expected
# _DSTAR_        = Str_(_STAR_  * 2)         # PYCHOK expected
_DUNDER_         = Str_(_UNDER_ * 2)         # PYCHOK expected


def _dunder_name(inst, *dflt):
    '''(INTERNAL) Get the double_underscore __name__ attr.
    '''
    try:
        return inst.__name__
    except AttributeError:
        pass
    return dflt[0] if dflt else inst.__class__.__name__


def _enquote(strs, quote=_QUOTE2_):  # in .basics
    '''(INTERNAL) Enquote a string containing whitespace.
    '''
    if len(strs.split()) > 1:
        strs = NN(quote, strs, quote)
    return strs


def _load_lib(name):
    '''(INTERNAL) Load a C{dylib}, B{C{name}} must startwith('lib').
    '''
    # macOS 11+ (aka 10.16) no longer provides direct loading of
    # system libraries, instead it installs the library after a
    # low-level dlopen(name) call where name is the library base
    # name.  As a result, ctypes.util.find_library can not find
    # any library not previously dlopen'ed.
    from ctypes import CDLL
    from ctypes.util import find_library
    from sys import platform

    ns = find_library(name), name
    if platform[:6] == 'darwin':  # and os.name == 'posix'
        from ctypes import _dlopen
        from os.path import join
        ns += (_DOT_(name, 'dylib'),
               _DOT_(name, 'framework'), join(
               _DOT_(name, 'framework'), name))
    else:  # non-macOS
        def _dlopen(unused):
            return True

    for n in ns:
        try:
            if n and _dlopen(n):  # handle
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
    global _pl2List
    if _pl2List is None:
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
        _pl2List = [platform.architecture()[0],  # bits
                    m]  # arm64, arm64_x86_64, x86_64, etc.
    return sep.join(_pl2List) if sep else _pl2List  # 2-list()


def _pythonarchine(sep=NN):  # in test/base.py versions
    '''(INTERNAL) Get PyPy and Python versions and C{_platform2} as C{3-} or C{4-list} or C{str}.
    '''
    global _Py3List
    if _Py3List is None:
        from sys import version  # XXX shadows?
        _Py3List = [_Python_(version)] + _platform2()
        if _PyPy__ in version:  # see test/base.py
            v = version.split(_PyPy__)[1].split(None, 1)[0]
            _Py3List.insert(0, NN(_PyPy__, v))
    return sep.join(_Py3List) if sep else _Py3List  # 3- or 4-list


def _splituple(strs, *sep_splits):  # in .basics, .mgrs
    '''(INTERNAL) Split a C{comma}- or C{whitespace}-separated
       string into a C{tuple} of stripped strings.
    '''
    t = (strs.split(*sep_splits) if sep_splits else
         strs.replace(_COMMA_, _SPACE_).split()) if strs else ()
    return tuple(s.strip() for s in t if s)


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


def _under_name(name):  # in .datums
    '''(INTERNAL) Prefix C{name} with I{underscore}.
    '''
    return name if name.startswith(_UNDER_) else NN(_UNDER_, name)


def _usage(file_py, *args):  # in .etm
    '''(INTERNAL) Build "usage: python -m ..." cmd line for module B{C{file_py}}.
    '''
    import os, sys  # PYCHOK imports
    p = NN(_python_, sys.version_info[0])
    m = os.path.dirname(file_py).replace(os.getcwd(), _ELLIPSIS_) \
                                .replace(os.sep, _DOT_).strip()
    b, x = os.path.splitext(os.path.basename(file_py))
    if x == '.py' and b != '__main__':
        m = _DOT_(m, b)
    return NN('usage', _SPACE_(_COLON_, p, '-m', _enquote(m), *args))


def _version2(version, n=2):
    '''(INTERNAL) Split C{B{version} str} into 1-, 2- or 3-tuple of C{int}s.
    '''
    def _int(*vs):
        for v in vs:
            try:
                yield int(v)
            except (TypeError, ValueError):
                pass

    t = tuple(_int(*_splituple(version, _DOT_, 2))) + (0,) * n
    return t[:n]


__all__ = (_NN_,  # not MISSING!
            Str_.__name__,  # classes
            machine.__name__)  # in .lazily
__version__ = '22.09.16'

if __name__ == '__main__':

    from pygeodesy.errors import itemsorted
    from pygeodesy.lazily import printf

    t = b = 0
    for n, v in itemsorted(locals(), asorted=False, reverse=True):
        if n.startswith(_UNDER_) and n.endswith(_UNDER_) and \
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
# Copyright (C) 2016-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
