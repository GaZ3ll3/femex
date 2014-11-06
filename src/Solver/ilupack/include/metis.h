#ifndef _METIS_H_
#define _METIS_H_
/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: metis.h,v 1.1 1998/11/27 17:59:21 karypis Exp $
 */


#include <stdio.h>
#ifdef __STDC__
#include <stdlib.h>
#else
#include <malloc.h>
#endif
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#include "../../../ilupack/include/long_integer.h"

#include "../../../ilupack/include/metis_defs.h"
#include "../../../ilupack/include/metis_struct.h"
#include "../../../ilupack/include/metis_macros.h"
#include "../../../ilupack/include/metis_rename.h"
#include "../../../ilupack/include/metis_proto.h"

#endif
