#ifndef dopri5Mexh
#define dopri5Mexh

#include "dopriMex.h"

extern void DOPRI5_ (Fint *n, 
  DOPRIRightSide fcn, double *tStart, double *x, double *tEnd,
	double  *rtol, double *atol, Fint *itol,
	DOPRISolout solout, Fint *iout,
	double *work, Fint *lwork, Fint *iwork, Fint *liwork,
	double *rpar, Fint *ipar, Fint *idid);
	
extern double CONTD5_ (Fint *i,
  double *x, double *con,Fint *icomp,Fint *nd);

#endif
