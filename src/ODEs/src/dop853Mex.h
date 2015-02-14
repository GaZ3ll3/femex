#ifndef dop853Mexh
#define dop853Mexh

#include "dopriMex.h"

extern void DOP853_ (Fint *n, 
  DOPRIRightSide fcn, double *tStart, double *x, double *tEnd,
	double  *rtol, double *atol, Fint *itol,
	DOPRISolout solout, Fint *iout,
	double *work, Fint *lwork, Fint *iwork, Fint *liwork,
	double *rpar, Fint *ipar, Fint *idid);
	
extern double CONTD8_ (Fint *i,
  double *x,double *con,Fint *icomp,Fint *nd);

#endif
