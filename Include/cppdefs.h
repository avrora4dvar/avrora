/* Option selectors. Uncomment if used without Build
#undef  NONLINEAR
#define TANGENT
#undef  ADJOINT
#undef  MPI
*/

/*
List of available CPP options, AVRORA
*/

#undef  UV_ADV
#undef  UV_COR
#undef  UV_VIS2
#undef  TS_DIF2
#undef  TS_FIXED
#undef  NONLIN_EOS
#undef  CURVGRID

#undef  DJ_GRADPS /* Shchepetkins' splines, NONLINEAR only */
#undef  DJ_STD2   /* Standard 2nd order density Jacobian */

#undef  NS_PERIODIC

#undef NORTHERN_WALL
#undef SOUTHERN_WALL
#undef WESTERN_WALL
#undef EASTERN_WALL

#undef NORTH_M2FLATHER
#undef SOUTH_M2FLATHER
#undef EAST_M2FLATHER
#undef WEST_M2FLATHER


#undef  TRIDIAG_ROMS
#undef  TRIDIAG_AVRORA
#undef  STEP3D_OLD
#undef  EXCLUDE
#undef  UV_COR_MASK

#undef ANA_FWD_MIX /* if defined, const fwd Akv and Akt */

/* CASES */

#if defined NONLINEAR

# define  UV_ADV
# define  UV_COR
# define CURVGRID
# define UV_VIS2
# define TS_DIF2
# undef  TS_FIXED
# define NONLIN_EOS

# define DJ_GRADPS /* Shchepetkins' splines, NONLINEAR only */

# define NS_PERIODIC

# define TRIDIAG_AVRORA

#endif /* NONLINEAR */

# if defined TANGENT || defined ADJOINT
#define UV_ADV
#define CURVGRID
#define UV_COR
#define UV_VIS2
#define TS_DIF2
#undef  TS_FIXED
#undef  NONLIN_EOS
#undef  ANA_FWD_MIX /* if defined, const fwd Akv and Akt */

#undef  DJ_GRADPS /* Shchepetkins' splines, NONLINEAR only */
#define DJ_STD2   /* Standard 2nd order density Jacobian */

#undef  NS_PERIODIC

#undef  NORTHERN_WALL
#undef  SOUTHERN_WALL
#define EASTERN_WALL
#undef  WESTERN_WALL

#define NORTH_M2FLATHER
#define SOUTH_M2FLATHER
#undef  EAST_M2FLATHER
#define WEST_M2FLATHER

#undef  STEP3D_NEW
#undef  TRIDIAG_ROMS
#undef  TRIDIAG_AVRORA
#undef  STEP3D_OLD
#undef  EXCLUDE
#undef  UV_COR_MASK

# endif /* TANGENT || ADJOINT */
