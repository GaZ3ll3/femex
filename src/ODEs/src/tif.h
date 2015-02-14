#ifndef tifh
#define tifh

#include <stdint.h>
#include "mex.h"

/* PREFIX: tif   (type interface) */

/* Fortran Integers:   */
#ifdef MX_COMPAT_32
typedef  int32_t Fint;
#define FINT_MAX INT32_MAX
#else
typedef  int64_t Fint;
#define FINT_MAX INT64_MAX
#endif


#ifdef FORTRANNOUNDER
#ifdef FORTRANUPP
#ifdef FORTRANLEADUNDER
/* _FORTRANNAME */
#define TIFTEST_ _TIFTEST
#else
/* FORTRANNAME */
#define TIFTEST_ TIFTEST
#endif
#else
#ifdef FORTRANLEADUNDER
/* _fortranname */
#define TIFTEST_ _tiftest
#else
/* fortranname */
#define TIFTEST_ tiftest
#endif
#endif
#else
#ifdef FORTRANUPP
#ifdef FORTRANLEADUNDER
/* _FORTRANNAME_ */
#define TIFTEST_ _TIFTEST_
#else
/* FORTRANNAME_ */
#endif
#else
#ifdef FORTRANLEADUNDER
/* _fortranname_ */
#define TIFTEST_ _tiftest_
#else
/* fortranname_ */
#define TIFTEST_ tiftest_
#endif
#endif
#endif

extern void TIFTEST_ (mwSize *n);

void tif_checkFint (void);

/* size_t -> mwSize */
mwSize tif_st2mws (size_t s);

/* mwSize -> Fint
 * because mwSize can be UNSIGNED and Fint is SIGNED
 * not every mwSize value fits in an Fint */
Fint tif_mws2Fint (mwSize s);

/* does this value (of type size_t) fit in a mwSize? */
char tif_fits_st_in_mws (size_t s);

/* print size of types and other environment stuff */
void tif_printInfos(void);

#endif
