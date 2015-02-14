/*
 * DOP853MEX-Interface by C. Ludwig
 * Version:  $Id: dopri5Mex.c 1241 2014-11-04 19:29:12Z luchr $ */
#define DOPRI5MexVersion "04. November 2014"
/*
 * English:
 *   Questions, requests, problems, notes, 
 *   remarks and bugs to
 *  ludwig@ma.tum.de
 * 
 * German:
 *   Fragen, Wünsche, Probleme, Anregungen, 
 *   Anmerkungen, Bemerkungen und Bugs an
 *  ludwig@ma.tum.de
 */
#include <math.h>
#include "mex.h"
#include "string.h"
#include "options.h"
#include "dopri5Mex.h"

#include "dopriMex.c"

