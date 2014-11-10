#ifndef _MUMPS_MATCHING_H
#define _MUMPS_MATCHING_H

/**************** include nameshsl.h if necessary *********/
#include "long_integer.h"
#include "../../../ilupack/include/long_integer.h"
#include "../../../ilupack/include/namesmumps_matching.h"

void dmumps_match(integer *, integer *, integer *, 
		  integer *, integer *, doubleprecision *, 
		  doubleprecision *, doubleprecision *);
void smumps_match(integer *, integer *, integer *, 
		  integer *, integer *, real *, 
		  real *, real *);
void cmumps_match(integer *, integer *, integer *, 
		  integer *, integer *, complex *, 
		  real *, real *);
void zmumps_match(integer *, integer *, integer *, 
		  integer *, integer *, doublecomplex *, 
		  doubleprecision *, doubleprecision *);

#endif /* _MUMPS_MATCHING_H */
