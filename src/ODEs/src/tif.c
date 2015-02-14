#include "tif.h"
#include <stdio.h>

char tif_fits_st_in_mws (size_t s) {
  return (s>MWSIZE_MAX?0:1);
}

mwSize tif_st2mws (size_t s) {
  static char msg[255];
  if (s>MWSIZE_MAX) {
    snprintf(msg,255,
      "Found size_t value %" FMT_SIZE_T "u too large for mwSize with "
      "MWSIZE_MAX=%" FMT_SIZE_T "u.",s,(size_t)MWSIZE_MAX);
    mexErrMsgTxt(msg);
    return 0;
  }
  else {
    return (mwSize)s;
  }
}

Fint tif_mws2Fint (mwSize s) {
  static char msg[255];
  if (s>FINT_MAX) {
    snprintf(msg,255,
      "Found mwSize value %li too large for Fortran's (signed) integers, "
      "because FINT_MAX=%" FMT_SIZE_T "u.", (long)s,(size_t)FINT_MAX);
    mexErrMsgTxt(msg);
    return 0;
  }
  else {
    return (Fint)s;
  }
}

void tif_printInfos(void) {
  mexPrintf("\n");
  mexPrintf("  sizeof(int)      = %" FMT_SIZE_T "u\n",sizeof(int));
  mexPrintf("  sizeof(long)     = %" FMT_SIZE_T "u\n",sizeof(long));
  mexPrintf("  sizeof(size_t)   = %" FMT_SIZE_T "u\n",sizeof(size_t));
  mexPrintf("  sizeof(mwSize)   = %" FMT_SIZE_T "u\n",sizeof(mwSize));
  mexPrintf("  sizeof(Fint)     = %" FMT_SIZE_T "u\n",sizeof(Fint));
  mexPrintf("  MWSIZE_MAX       = 0x%" FMT_SIZE_T "x = %" FMT_SIZE_T "u\n",
    (size_t)MWSIZE_MAX,(size_t)MWSIZE_MAX);
  mexPrintf("  FINT_MAX         = 0x%" FMT_SIZE_T "x = %" FMT_SIZE_T "u\n",
    (size_t)FINT_MAX,(size_t)FINT_MAX);
  mexPrintf("  MX_COMPAT_32     = "
#ifdef MX_COMPAT_32
  "defined    => compatibleArrayDims"
#else
  "undefined  => largeArrayDims"
#endif
  "\n");
  mexPrintf("  FORTRANNOUNDER   = "
#ifdef FORTRANNOUNDER
  "defined"
#else
  "undefined"
#endif
  "\n");
  mexPrintf("  FORTRANUPP       = "
#ifdef FORTRANUPP
  "defined"
#else
  "undefined"
#endif
  "\n");
  mexPrintf("  FORTRANLEADUNDER = "
#ifdef FORTRANLEADUNDER
  "defined"
#else
  "undefined"
#endif
  "\n");
  mexPrintf("  FMT_SIZE_T       = %s\n",FMT_SIZE_T);
  mexPrintf("\n");
}

void tif_checkFint (void) {
  uint8_t bytearr[16];
  size_t i;

  for (i=0; i<16; i++) { bytearr[i]=0xff; } /* 16 bytes with all bits set */

  TIFTEST_((mwSize*)(&bytearr));            /* let Fortran return 0 integer */

  /* test, how many bytes are zero */
  i=0;
  while (i<16 && 0x00==bytearr[i]) { i++; }

  if (i!=sizeof(mwSize)) {
    tif_printInfos();
    mexPrintf("There is a problem concerning the size of Fortran-Integers.\n");
#ifdef MX_COMPAT_32
    mexPrintf("MX_COMPAT_32 was DEFINED.\n");
#else
    mexPrintf("MX_COMPAT_32 was NOT DEFINED.\n");
#endif
    mexPrintf("Hence I assumed a Fortran-integer "
              "has a size of %" FMT_SIZE_T "u bytes.\n",sizeof(mwSize));
    mexPrintf("But a test-call indicates that a Fortran-integer "
              "has a size of %" FMT_SIZE_T "u bytes.\n",i);
    mexPrintf("I'm confused.\n");
    mexErrMsgTxt("Inconsistent MX_COMPAT_32 setting and Fortran-integer size");
  }
}

