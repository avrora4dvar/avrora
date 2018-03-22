#define INSTANT_SSH
#undef  INSTANT_SST
#define AVG

#ifdef INSTANT_SSH
# undef  DAILY_SSH
#else
# define DAILY_SSH
#endif

#ifdef INSTANT_SST
# undef  DAILY_SST
#else
# define DAILY_SST
#endif

#ifdef AVG
# undef DAILY_SST
#endif
