#ifndef dopriMexh
#define dopriMexh

#include "tif.h"

/* Optionsettings */
#define OPT_WARNMISS "OptWarnMiss"
#define OPT_WARNTYPE "OptWarnType"
#define OPT_WARNSIZE "OptWarnSize"

/* H and Tols */
#define OPT_INITIALSS "InitialStep"
#define OPT_RTOL "RelTol"
#define OPT_ATOL "AbsTol"

/* Output */
#define OPT_OUTPUTFUNCTION "OutputFcn"
#define OPT_OUTPUTCALLMODE "OutputCallMode"

/* StepSizeSelection */
#define OPT_RHO "rho"
#define OPT_SSMINSEL "StepSizeMinSelection"
#define OPT_SSMAXSEL "StepSizeMaxSelection"
#define OPT_SSBETA "StepSizeBeta"
#define OPT_MAXSS "MaxStep"
#define OPT_MAXSTEPS "MaxNumberOfSteps"

/* Rest */
#define OPT_IGPIDO "IncludeGridPointsInDenseOutput"
#define OPT_EPS "eps"
#define OPT_STEST "StiffTestAfterStep"
#define OPT_FUNCCALLMETHOD "FuncCallMethod"

struct ListElement
{ /* Verkettete Liste mit (t,x)-Paaren, also Vektoren der Länge d+1 */
  double* values;
  struct ListElement *next;
};
typedef struct ListElement SListElement;
typedef SListElement* PListElement;

struct ParameterGlobal 
{ /* globale Variablen usw. */
  mwSize d;               /* Dimension des Systems NICHT in size_t! */
  mwSize tLength;         /* Länge des t-Vektors */  
  double* tPointer;       /* Zeiger auf alle Zeitwerte */
  double direction;       /* sign(tEnd-tStart) */  
  int funcCallMethod;     /* Methode, um Matlab-Funktionen aufzurufen */ 
};
typedef struct ParameterGlobal SParameterGlobal;

struct ParameterOptions
{ /* Parameter für Optionen */    
  const mxArray *opt;  /* Options */
  char optCreated;     /* Flag, ob Options selbst erzeugt wurde */
};
typedef struct ParameterOptions SParameterOptions;

struct ParameterDOPRI
{ /* Parameter für DOPRI */
  double tStart;    /* Startzeitpunkt */
  double tEnd;      /* Endzeitpunkt */
  double *xStart;   /* Startwert */
  double *RTOL;     /* releative Toleranz */
  double *ATOL;     /* absolute Toleranz */
  Fint ITOL;        /* Switch für RTOL und ATOL */
  Fint IOUT;        /* Switch für SOLOUT */
  char denseFlag;   /* Flag, ob dense output */
  double *WORK;     /* Double-Arbeits-Array */
  Fint LWORK;       /* Länge von WORK */
  Fint *IWORK;      /* Integer-Arbeits-Array */
  Fint LIWORK;      /* Länge von IWORK */
  double *RPAR;     /* Zusatz double-array */
  Fint *IPAR;       /* Zusatz int-array */
  Fint IDID;        /* Status */
};
typedef struct ParameterDOPRI SParameterDOPRI;

struct ParameterRightSide
{ /* Parameter für rechte Seite f */
  char *rightSideFcn;            /* Funktionsname für rechte Seite, 
                                    falls als String angegeben */
  const mxArray *rightSideFcnH;  /* Funktionshandle oder inline-function für rightSide */
  mxArray *tArg;                 /* Zum Aufruf von rightSideFcn: t */
  mxArray *xArg;                 /* Zum Aufruf von rightSideFcn: x */
};
typedef struct ParameterRightSide SParameterRightSide;

struct ParameterOutput
{ /* Parameter zum Speichern der DOPRI-Ausgabe */
  SListElement txList;        /* Start der txListe */
  PListElement lastTXElement; /* letzter Eintrag in txListe */
  mwSize numberOfElements;    /* Anzahl der Einträge in txListe */
  mwSize tPos;                /* Position im t-Vektor */
  int includeGrid;            /* Flag, dense Ausgabe mit Gridpoints */
  char* outputFcn;            /* Outputfunction, falls als String angegeben */
  const mxArray *outputFcnH;  /* Outputfunction: handle oder inline-function */
  int outputCallMode;         /* Modus für outputFcn-Aufruf */
  mxArray *tArg;              /* Zum Aufruf von outputFcn: t */
  mxArray *xArg;              /* Zum Aufruf von outputFcn: x */  
  mxArray *emptyArg;          /* Zum Aufruf von outputFcn: empty array [] */
  mxArray *toldArg ;          /* Zum Aufruf von outputFcn: tOld */
};
typedef struct ParameterOutput SParameterOutput;

struct DopriDense
{ /* Argumente zum Aufruf von CONTD */
  double *con;
  Fint *icomp;
  Fint *nd;
};
typedef struct DopriDense SDopriDense;

typedef void (*DOPRIRightSide)(Fint *n, double *t,
  double *x, double *f, double *rpar, Fint *ipar);

void DOPRIRightSideFunc (Fint *n, double *t,
  double *x, double *f, double *rpar, Fint *ipar);

typedef void (*DOPRISolout)(Fint *nr, double *told,
  double *t, double *x, Fint *n,
	double *con, Fint *icomp, Fint *nd,
	double *rpar, Fint *ipar, Fint *irtrn);

void DOPRISoloutFunc (Fint *nr, double *told,
  double *t, double *x, Fint *n,
	double *con, Fint *icomp, Fint *nd,
	double *rpar, Fint *ipar, Fint *irtrn);

#ifdef FORTRANNOUNDER
/* Fotran functions without underscore */
#ifdef FORTRANUPP
/* Fotran functions without underscore  & UPPERCASE letters */
#ifdef FORTRANLEADUNDER
/* Fortran functions without underscore & UPPERCASE letters & leading underscore */
#define DOP853_ _DOP853
#define CONTD8_ _CONTD8
#define DOPRI5_ _DOPRI5
#define CONTD5_ _CONTD5
#else
#define DOP853_ DOP853
#define CONTD8_ CONTD8
#define DOPRI5_ DOPRI5
#define CONTD5_ CONTD5
#endif
#else
/* Fotran functions without underscore  & lowercase letters */
#ifdef FORTRANLEADUNDER
/* Fortran functions without underscore & lowercase letters & leading underscore */
#define DOP853_ _dop853
#define CONTD8_ _contd8
#define DOPRI5_ _dopri5
#define CONTD5_ _contd5
#else
/* Fortran functions without underscore & lowercase letters & without leading underscore */
#define DOP853_ dop853
#define CONTD8_ contd8
#define DOPRI5_ dopri5
#define CONTD5_ contd5
#endif
#endif
#else
/* Fortran functions with underscore */
#ifdef FORTRANUPP
/* Fortran functions with underscore & UPPERCASE letters */
#ifdef FORTRANLEADUNDER
/* Fortran functions with underscore & UPPERCASE letters & leading underscore */
#define DOP853_ _DOP853_
#define CONTD8_ _CONTD8_
#define DOPRI5_ _DOPRI5_
#define CONTD5_ _CONTD5_
#else
/* Fortran functions with underscore & UPPERCASE letters */
#endif
#else
/* Fortran functions with underscore & lowercase letters */
#ifdef FORTRANLEADUNDER
/* Fortran functions with underscore & lowercase letters  & leading underscore*/
#define DOP853_ _dop853_
#define CONTD8_ _contd8_
#define DOPRI5_ _dopri5_
#define CONTD5_ _contd5_
#else
/* Fortran functions with underscore & lowercase letters  & without leading underscore*/
#define DOP853_ dop853_
#define CONTD8_ contd8_
#define DOPRI5_ dopri5_
#define CONTD5_ contd5_
#endif
#endif
#endif


#endif
