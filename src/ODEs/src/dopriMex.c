/* Common code for dopri5Mex and dop853Mex */

static SOptionSettings optSet;
static SParameterOptions paramOpt;
static SParameterGlobal paramGlobal;
static SParameterDOPRI paramDOPRI;
static SParameterRightSide paramRightSide;
static SParameterOutput paramOutput;

static SDopriDense dopriDense;
static char isInnerCall=0;

static char ismxArrayString (const mxArray *arr) {
  if (arr==NULL) return 0;
  if (mxIsChar(arr)) return 1;
  return 0;
}

static char ismxArrayFunction (const mxArray *arr) {
  if (arr==NULL) return 0;
  return (mxGetClassID(arr)==mxFUNCTION_CLASS)?1:0;
}

static char ismxArrayInline (const mxArray *arr) {
  if (arr==NULL) return 0;
  return (strcmp(mxGetClassName(arr),"inline")==0)?1:0;
}

static void clearList (PListElement current) {
  PListElement next;
  
  while (current!=NULL) {
    next=current->next;
    if (current->values!=NULL)
      {mxFree(current->values);current->values=NULL;}
    mxFree(current);
    current=next;
  }
}

static void initVars (void) {
  /* Option settings */
  optSet.warnMiss=0;optSet.warnType=1;optSet.warnSize=1;
  
  /* Parameters for Options*/
  paramOpt.opt=NULL;paramOpt.optCreated=0;
  
  /* global parameters */
  paramGlobal.tPointer=NULL;
  isInnerCall=0;
  
  /* parameters for DOPRI */
  paramDOPRI.xStart=NULL;
  paramDOPRI.RTOL=NULL;paramDOPRI.ATOL=NULL;
  paramDOPRI.WORK=NULL;paramDOPRI.IWORK=NULL;
  paramDOPRI.RPAR=NULL;paramDOPRI.IPAR=NULL;
  
  /* parameters for RightSide */
  paramRightSide.rightSideFcn=NULL;
  paramRightSide.rightSideFcnH=NULL;
  paramRightSide.tArg=NULL;paramRightSide.xArg=NULL;
  
  /* parameters for Output */
  paramOutput.txList.values=NULL;paramOutput.txList.next=NULL;
  paramOutput.lastTXElement=&paramOutput.txList;
  paramOutput.numberOfElements=0;
  paramOutput.outputFcn=NULL;
  paramOutput.outputFcnH=NULL;
  paramOutput.tArg=NULL;paramOutput.xArg=NULL;paramOutput.emptyArg=NULL;
  paramOutput.toldArg=NULL;
}

static void doneVars (void) {
  /* Parameters for Options*/
  if ((paramOpt.optCreated) && (paramOpt.opt!=NULL))
    {mxDestroyArray((mxArray*)paramOpt.opt);paramOpt.opt=NULL;}
    
  /* global parameters */
  /* It's not allowed to free paramGlobal.tPointer because
     it belogns to the caller */
  /* paramGlobal.tPointer darf nicht entsorgt werden, denn
     er gehört dem Aufrufer */
     
  /* parameters for DOPRI */
  if (paramDOPRI.xStart!=NULL)
    {mxFree(paramDOPRI.xStart);paramDOPRI.xStart=NULL;}
  if (paramDOPRI.RTOL!=NULL)
    {mxFree(paramDOPRI.RTOL);paramDOPRI.RTOL=NULL;}
  if (paramDOPRI.ATOL!=NULL)
    {mxFree(paramDOPRI.ATOL);paramDOPRI.ATOL=NULL;}
  if (paramDOPRI.WORK!=NULL)
    {mxFree(paramDOPRI.WORK);paramDOPRI.WORK=NULL;}
  if (paramDOPRI.IWORK!=NULL)
    {mxFree(paramDOPRI.IWORK);paramDOPRI.IWORK=NULL;}
  if (paramDOPRI.RPAR!=NULL)
    {mxFree(paramDOPRI.RPAR);paramDOPRI.RPAR=NULL;}
  if (paramDOPRI.IPAR!=NULL)
    {mxFree(paramDOPRI.IPAR);paramDOPRI.IPAR=NULL;}
    
  /* parameters for RightSide */
  if (paramRightSide.rightSideFcn!=NULL)
    {mxFree(paramRightSide.rightSideFcn);paramRightSide.rightSideFcn=NULL;}
  paramRightSide.rightSideFcnH=NULL; /* gehört Aufrufer => nicht freigeben */
  if (paramRightSide.tArg!=NULL)
    {mxDestroyArray(paramRightSide.tArg);paramRightSide.tArg=NULL;}
  if (paramRightSide.xArg!=NULL)
    {mxDestroyArray(paramRightSide.xArg);paramRightSide.xArg=NULL;}
    
  /* parameters for Output */
  clearList(paramOutput.txList.next);paramOutput.txList.next=NULL;
  paramOutput.lastTXElement=&paramOutput.txList;
  if (paramOutput.outputFcn!=NULL)
    {mxFree(paramOutput.outputFcn);paramOutput.outputFcn=NULL;}
  paramOutput.outputFcnH=NULL; /* gehört Aufrufer => nicht freigeben */
  if (paramOutput.tArg!=NULL)
    {mxDestroyArray(paramOutput.tArg);paramOutput.tArg=NULL;}
  if (paramOutput.xArg!=NULL)
    {mxDestroyArray(paramOutput.xArg);paramOutput.xArg=NULL;}
  if (paramOutput.emptyArg!=NULL) 
    {mxDestroyArray(paramOutput.emptyArg);paramOutput.emptyArg=NULL;}
  if (paramOutput.toldArg!=NULL) 
    {mxDestroyArray(paramOutput.toldArg);paramOutput.toldArg=NULL;}
}

static void stopMexFunctionImpl (int errNo,
  size_t i1, size_t i2, size_t i3, int i4, double d1, char doneFlag) {
  const char *msg;

  (void)i4;

  isInnerCall=0;
#ifdef DOPRI5MexVersion
  mexPrintf("Error (%i) [Version: %s]\n",errNo,DOPRI5MexVersion);
#endif
#ifdef DOP853MexVersion
  mexPrintf("Error (%i) [Version: %s]\n",errNo,DOP853MexVersion);
#endif
  switch (errNo) {
    #include "errors.c"
    default: msg="Unknown errornumber (Unbekannte Fehlernummer)";break;
  }
  
  if (doneFlag) {doneVars();}
  mexErrMsgTxt(msg);
}

static void stopMexFunction (int errNo,
  size_t i1, size_t i2, size_t i3, int i4, double d1) {

  stopMexFunctionImpl(errNo,i1,i2,i3,i4,d1,1);
}


static void addTXtoList (double t, double *x) {
  double* dpointer;
  PListElement target;
  
  target=mxMalloc( sizeof(SListElement) );
  target->next=NULL;paramOutput.lastTXElement->next=target;
  paramOutput.lastTXElement=target;paramOutput.numberOfElements++;
  target->values=mxMalloc((size_t)(paramGlobal.d+1)*sizeof(double));
  dpointer=target->values;
  *dpointer=t;dpointer++;
  memcpy(dpointer,x,(size_t)paramGlobal.d*sizeof(double));  
}

static void checkNumberOfArgs (int nlhs, int nrhs) {
  if ((nlhs<2) || (nlhs>4)) stopMexFunction(1,(size_t)nlhs,0,0,0,0);    
    
  if ((nrhs!=3) && (nrhs!=4)) stopMexFunction(2,(size_t)nrhs,0,0,0,0);
}

static void processArgs (int nrhs, const mxArray* prhs[]) {
  size_t len;
  mwSize buflen;
  double *dpointer;
  
  /* 1st arg: right side */
  if (ismxArrayString(prhs[0])) {
    if ((mxGetNumberOfDimensions(prhs[0])!=2) || (mxGetM(prhs[0])!=1))
      stopMexFunction(4,(size_t)mxGetNumberOfDimensions(prhs[0]),mxGetM(prhs[0]),0,0,0);

    paramRightSide.rightSideFcn=mxArrayToString(prhs[0]);
    paramRightSide.rightSideFcnH=prhs[0];
  } else 
  if ( (ismxArrayFunction(prhs[0])) || (ismxArrayInline(prhs[0])) ) {
    if (((size_t)mxGetNumberOfDimensions(prhs[0])!=2) || (mxGetM(prhs[0])!=1) ||
        (mxGetN(prhs[0])!=1)) 
      stopMexFunction(15,(size_t)mxGetNumberOfDimensions(prhs[0]),
                      mxGetM(prhs[0]),mxGetN(prhs[0]),0,0);
    paramRightSide.rightSideFcn=NULL; /* kein String */ 
    paramRightSide.rightSideFcnH=prhs[0];
  } else {
    stopMexFunction(3,0,0,0,0,0);
  }
  
  /* 2nd arg: row-vector containing time-values */
  if (!mxIsDouble(prhs[1])) stopMexFunction(5,0,0,0,0,0);
  if (((size_t)mxGetNumberOfDimensions(prhs[1])!=2) || (mxIsSparse(prhs[1])) || 
      (mxGetM(prhs[1])!=1))
    stopMexFunction(6,(size_t)mxGetNumberOfDimensions(prhs[1]),mxGetM(prhs[1]),0,0,0);
  len=mxGetN(prhs[1]);
  if (!tif_fits_st_in_mws(len)) 
    stopMexFunction(20,len,(size_t)MWSIZE_MAX,0,0,0);
  paramGlobal.tLength=tif_st2mws(len);
  if (paramGlobal.tLength<2) stopMexFunction(7,0,0,0,0,0);
  if (paramGlobal.tLength>2) {
    paramDOPRI.denseFlag=1;paramDOPRI.IOUT=2;
  } else {
    paramDOPRI.denseFlag=0;paramDOPRI.IOUT=1;
  }
  dpointer=mxGetPr(prhs[1]);paramGlobal.tPointer=dpointer;
  paramDOPRI.tStart=dpointer[0];
  paramDOPRI.tEnd=dpointer[paramGlobal.tLength-1];
  if (paramDOPRI.tStart==paramDOPRI.tEnd) stopMexFunction(12,0,0,0,0,0);
  paramGlobal.direction=(paramDOPRI.tEnd-paramDOPRI.tStart)>0?1.0:-1.0;
  for (buflen=1; buflen<paramGlobal.tLength; buflen++,dpointer++)
    if (paramGlobal.direction*(dpointer[0]-dpointer[1])>0)
      stopMexFunction(14,paramGlobal.direction>0.0?1:0,0,0,0,0);        
  
  /* 3rd arg: start vector */
  if (!mxIsDouble(prhs[2])) stopMexFunction(8,0,0,0,0,0);
  if (((size_t)mxGetNumberOfDimensions(prhs[2])!=2) || (mxGetN(prhs[2])!=1) ||
      (mxIsSparse(prhs[2])))
    stopMexFunction(9,(size_t)mxGetNumberOfDimensions(prhs[2]),mxGetN(prhs[2]),0,0,0);
  len=mxGetM(prhs[2]);
  if (!tif_fits_st_in_mws(len)) 
    stopMexFunction(19,len,(size_t)MWSIZE_MAX,0,0,0);
  paramGlobal.d=tif_st2mws(len);
  if (paramGlobal.d<1) stopMexFunction(10,0,0,0,0,0);
  dpointer=mxGetPr(prhs[2]);
  paramDOPRI.xStart=mxMalloc((size_t)paramGlobal.d*sizeof(double));
  memcpy(paramDOPRI.xStart,dpointer,(size_t)paramGlobal.d*sizeof(double));
  /* little remark: A COPY of the startvector is made, because
     it will be passed to Fortran code. To be sure, that
     it will not be changed there, a copy will be passed.
     x0 belongs to the caller and it is not allowed 
     (due to Matlab-Contract) to change it. */
  /* kleine Anmerkung: der Startvektor wird KOPIERT, da er an den 
     Fortran code weitergegeben wird. Um ganz sicher zu gehen,
     dass er nicht verändert wird, wird eine Kopie übergeben.
     x0 gehört ja dem Aufrufer und darf nach Matlab-Konvention
     nicht verändert werden. */
  
  /* 4th arg: struct with options */
  if (nrhs==4) {
    if (!mxIsStruct(prhs[3])) stopMexFunction(11,0,0,0,0,0);
    paramOpt.opt=prhs[3];paramOpt.optCreated=0;
  } else {
    paramOpt.opt=mxCreateStructMatrix(1,1,0,NULL);
    paramOpt.optCreated=1;
  }  
}

static void extractOptionSettings (void) {
  optSet.warnMiss=opt_getIntFromOpt(paramOpt.opt,&optSet,OPT_WARNMISS,0);
  optSet.warnType=opt_getIntFromOpt(paramOpt.opt,&optSet,OPT_WARNTYPE,1);
  optSet.warnSize=opt_getIntFromOpt(paramOpt.opt,&optSet,OPT_WARNSIZE,1);
}

static void extractTOLs (void) {
  size_t len1,len2;
  mwSize m1,n1,m2,n2;
  char res1,res2;
  char takeScalar;
  
  res1=opt_getSizeOfOptField(paramOpt.opt,OPT_RTOL,&len1,&len2);
  if (!tif_fits_st_in_mws(len1)) stopMexFunction(21,len1,0,0,0,0);
  if (!tif_fits_st_in_mws(len2)) stopMexFunction(21,len2,0,0,0,0);
  m1=tif_st2mws(len1); n1=tif_st2mws(len2);

  res2=opt_getSizeOfOptField(paramOpt.opt,OPT_ATOL,&len1,&len2);
  if (!tif_fits_st_in_mws(len1)) stopMexFunction(22,len1,0,0,0,0);
  if (!tif_fits_st_in_mws(len2)) stopMexFunction(22,len2,0,0,0,0);
  m2=tif_st2mws(len1); n2=tif_st2mws(len2);
  
  takeScalar=0;
  if (((res1!=0) || (res2!=0)) ||
      ((res1==0) && (m1==1)) ||
      ((res2==0) && (m2==1))) {
    takeScalar=1; 
  } else {
    if ((n1!=1) || (n2!=1)) takeScalar=1;
    if (m1!=m2) takeScalar=1;
    if (m1!=paramGlobal.d) takeScalar=1;
  }

  if (takeScalar) {
    paramDOPRI.RTOL=mxMalloc(1*sizeof(double));
    paramDOPRI.ATOL=mxMalloc(1*sizeof(double));
    *paramDOPRI.RTOL=opt_getDoubleFromOpt(paramOpt.opt,&optSet,OPT_RTOL,1e-3);
    *paramDOPRI.ATOL=opt_getDoubleFromOpt(paramOpt.opt,&optSet,OPT_ATOL,1e-6);
    paramDOPRI.ITOL=0;
  } else {
    paramDOPRI.RTOL=mxMalloc((size_t)paramGlobal.d*sizeof(double));
    paramDOPRI.ATOL=mxMalloc((size_t)paramGlobal.d*sizeof(double));
    opt_getDoubleVectorFromOpt(paramOpt.opt,&optSet,OPT_RTOL,
      (size_t)paramGlobal.d,1,paramDOPRI.RTOL);
    opt_getDoubleVectorFromOpt(paramOpt.opt,&optSet,OPT_ATOL,
      (size_t)paramGlobal.d,1,paramDOPRI.ATOL);
    paramDOPRI.ITOL=1;
  }
}

static void extractOutput (void) {
  mxArray *arr;

  if (paramDOPRI.denseFlag) {
    paramOutput.includeGrid=opt_getIntFromOpt(paramOpt.opt,&optSet,
      OPT_IGPIDO,0);
  }

  paramOutput.outputFcnH=NULL;paramOutput.outputFcn=NULL;
  arr=mxGetField(paramOpt.opt,0,OPT_OUTPUTFUNCTION);
  if ((arr!=NULL) && (!mxIsEmpty(arr))) {
    if ( (ismxArrayFunction(arr)) || (ismxArrayInline(arr)) ) {
      if (((size_t)mxGetNumberOfDimensions(arr)!=2) || (mxGetM(arr)!=1) ||
          (mxGetN(arr)!=1) || mxIsSparse(arr)) 
        stopMexFunction(401,(size_t)mxGetNumberOfDimensions(arr),
                        mxGetM(arr),mxGetN(arr),0,0);
      paramOutput.outputFcnH=arr;
    } else {
      paramOutput.outputFcn=opt_getStringFromOpt(paramOpt.opt,&optSet,
        OPT_OUTPUTFUNCTION,NULL);
      if (paramOutput.outputFcn!=NULL) {
        if (strlen(paramOutput.outputFcn)==0) {
          mxFree(paramOutput.outputFcn);paramOutput.outputFcn=NULL;
        } else {
          paramOutput.outputFcnH=arr;
        }
      }
    }
  }
  
  paramOutput.outputCallMode=1;
  if ((paramOutput.outputFcnH!=NULL) && (paramDOPRI.denseFlag)) {
    paramOutput.outputCallMode=opt_getIntFromOpt(paramOpt.opt,&optSet,
      OPT_OUTPUTCALLMODE,1);
    if ((paramOutput.outputCallMode<1) || (paramOutput.outputCallMode>3))
      paramOutput.outputCallMode=1;
  }
}

static void extractGlobalOptions (void) {
  paramGlobal.funcCallMethod=opt_getIntFromOpt(paramOpt.opt,&optSet,
    OPT_FUNCCALLMETHOD,1);
  switch (paramGlobal.funcCallMethod) {
    case 0: /* use mexCallMATLAB direct */
      if (paramRightSide.rightSideFcn==NULL) stopMexFunction(17,0,0,0,0,0);
      if ((paramOutput.outputFcnH!=NULL) && (paramOutput.outputFcn==NULL))
        stopMexFunction(18,0,0,0,0,0);
      break;
    case 1: /* use mexCallMATLAB to call feval */
      break;
    default:
      stopMexFunction(16,0,0,0,0,0);
      break;
  }
}

static void extractOptionsPart1 (void) {
  extractOptionSettings();
  extractTOLs();
  extractOutput();  
  extractGlobalOptions();
}

static void extractIWORKOpt (void) {
  int i;
  
  i=opt_getIntFromOpt(paramOpt.opt,&optSet,OPT_MAXSTEPS,100000);
  if (i<=0) stopMexFunction(101,0,0,0,i,0);
  paramDOPRI.IWORK[1-1]=i;
  
  /* it's the only set of coefficients */
  /* gibt momentan nur einen Koeff-Satz */
  paramDOPRI.IWORK[2-1]=1; /* gibt momentan nur einen Koeff-Satz */

  /* DOPRI: please don't try to write */
  /* DOPRI soll gar nicht erst versuchen zu schreiben: */
  paramDOPRI.IWORK[3-1]=-1; 
    
  i=opt_getIntFromOpt(paramOpt.opt,&optSet,OPT_STEST,1000);
  if (i==0) stopMexFunction(102,0,0,0,0,0);
  paramDOPRI.IWORK[4-1]=i;
  
  /* paramDOPRI.IWORK[5-1] was prepared in prepareWorkArrays */
  /* paramDOPRI.IWORK[5-1] wurde schon in prepareWorkArrays erledigt */  
}

static void extractWORKOpt (void) {
  double d;
  
  d=opt_getDoubleFromOpt(paramOpt.opt,&optSet,OPT_EPS,2.3e-16);
  if (!((d>=1e-35) && (d<1.0))) stopMexFunction(103,0,0,0,0,d);
  paramDOPRI.WORK[1-1]=d;
  
  d=opt_getDoubleFromOpt(paramOpt.opt,&optSet,OPT_RHO,0.9);
  if (!((d>1e-4) && (d<1.0))) stopMexFunction(104,0,0,0,0,d);
  paramDOPRI.WORK[2-1]=d;
  
#ifdef DOPRI5MexVersion
  d=opt_getDoubleFromOpt(paramOpt.opt,&optSet,OPT_SSMINSEL,0.2);
#endif
#ifdef DOP853MexVersion
  d=opt_getDoubleFromOpt(paramOpt.opt,&optSet,OPT_SSMINSEL,0.333);
#endif
  if (!(d>0.0)) stopMexFunction(105,0,0,0,0,d);
  paramDOPRI.WORK[3-1]=d;
  
#ifdef DOPRI5MexVersion
  d=opt_getDoubleFromOpt(paramOpt.opt,&optSet,OPT_SSMAXSEL,10.0);
#endif
#ifdef DOP853MexVersion
  d=opt_getDoubleFromOpt(paramOpt.opt,&optSet,OPT_SSMAXSEL,6.0);
#endif
  if (!(d>0.0)) stopMexFunction(106,0,0,0,0,d);
  paramDOPRI.WORK[4-1]=d;
  
#ifdef DOPRI5MexVersion
  d=opt_getDoubleFromOpt(paramOpt.opt,&optSet,OPT_SSBETA,0.04);
#endif
#ifdef DOP853MexVersion
  d=opt_getDoubleFromOpt(paramOpt.opt,&optSet,OPT_SSBETA,0.0);
#endif
  if (!(d<=0.2)) stopMexFunction(107,0,0,0,0,d);
  paramDOPRI.WORK[5-1]=d;
  
  d=opt_getDoubleFromOpt(paramOpt.opt,&optSet,OPT_MAXSS,
    paramDOPRI.tEnd-paramDOPRI.tStart);
  if (!(d!=0.0)) stopMexFunction(108,0,0,0,0,d);
  paramDOPRI.WORK[6-1]=d;
  
  d=opt_getDoubleFromOpt(paramOpt.opt,&optSet,OPT_INITIALSS,0.0);
  paramDOPRI.WORK[7-1]=d;
}


static void extractOptionsPart2 (void) {
  extractIWORKOpt();
  extractWORKOpt();
}

static void prepareHelpMxArrays (void) {
  paramRightSide.tArg=mxCreateDoubleMatrix(1,1,mxREAL);
  paramRightSide.xArg=mxCreateDoubleMatrix(paramGlobal.d,1,mxREAL);
  
  if ((paramOutput.outputFcnH!=NULL) || (paramDOPRI.denseFlag)) {
    paramOutput.tArg=mxCreateDoubleMatrix(1,1,mxREAL);
    paramOutput.xArg=mxCreateDoubleMatrix(paramGlobal.d,1,mxREAL);
    paramOutput.emptyArg=mxCreateDoubleMatrix(1,0,mxREAL);
    paramOutput.toldArg=mxCreateDoubleMatrix(1,1,mxREAL);
  }
}

static void prepareWorkArrays (void) {
  Fint i;
  mwSize nrdense;
  
  if (paramDOPRI.denseFlag) nrdense=paramGlobal.d; else nrdense=0;

#ifdef DOPRI5MexVersion
  paramDOPRI.LWORK=tif_mws2Fint(8*paramGlobal.d + 5*nrdense + 21);
#endif
#ifdef DOP853MexVersion
  paramDOPRI.LWORK=tif_mws2Fint(11*paramGlobal.d + 8*nrdense + 21);
#endif  
  paramDOPRI.WORK=mxMalloc((size_t)paramDOPRI.LWORK*sizeof(double));
  /* Std-Values */
  for (i=1-1; i<=20-1; i++) paramDOPRI.WORK[i]=0.0;
  
  paramDOPRI.LIWORK=tif_mws2Fint(nrdense+21);
  paramDOPRI.IWORK=mxMalloc((size_t)paramDOPRI.LIWORK*sizeof(Fint));  
  
  /* Std-Values */
  for (i=1-1; i<=20-1; i++) paramDOPRI.IWORK[i]=0;
  
  paramDOPRI.IWORK[5-1]=tif_mws2Fint(nrdense);
}

void DOPRIRightSideFunc (Fint *n, double *t,
  double *x, double *f, double *rpar, Fint *ipar) {

  size_t d;
  mxArray *rhs[3];
  mxArray *lhs[1];

  (void)rpar;
  (void)ipar;
  
  d=(size_t)(*n);lhs[0]=NULL;

  /* Call by value */
  /* Use always the SAME Matlab tArg and xArg */
  /* IMPORTANT NOTICE (if the right side is also a MEX-File)
     the right side must not "garble" the passed mxArrays:
     the size must not be changed and the memory must not
     be freed. 
     Hence: the right side has to take care, that the
     passed mxArrays have the same size and enough memory
     when returning. The values in the memory(-block)
     may be overwritten.
     If the right side is an m-file MATLAB obeys this
     restriction automatically. */
  /* WICHTIGE Anmerkung (falls rechte Seite auch MEX-File ist)
     die rechte Seite darf die übergebenen mxArrays nicht
     "verstümmeln": die Größe darf nicht verändert werden
     und auch der Speicherplatz darf nicht freigegeben werden.
     Also: die rechte Seite muss sicherstellen, dass die
     übergebenen mxArrays am Ende wieder diesselbe Größe
     und ausreichend Speicher haben. Der Speicher selbst
     darf natürlich zu rechenzwecken überschrieben werden.
     Bei m-Files achtet MATLAB automatisch darauf.
  */
  *mxGetPr(paramRightSide.tArg)= *t;
  memcpy(mxGetPr(paramRightSide.xArg),x,d*sizeof(double));
  switch (paramGlobal.funcCallMethod) {
    case 0: 
      rhs[0]=paramRightSide.tArg; 
      rhs[1]=paramRightSide.xArg;
        
      /* Call User's right side */
      mexCallMATLAB(1,lhs,2,rhs,paramRightSide.rightSideFcn);
      break;
    case 1:
      rhs[0]=(mxArray*)paramRightSide.rightSideFcnH;
      rhs[1]=paramRightSide.tArg;
      rhs[2]=paramRightSide.xArg;

      /* Call User's right side */
      mexCallMATLAB(1,lhs,3,rhs,"feval");
      break;
    default: stopMexFunction(1002,0,0,0,0,0);break;
  }
  
  /* check return values */
  if (lhs[0]==NULL) stopMexFunction(301,0,0,0,0,0);
  if ((!mxIsDouble(lhs[0])) || (mxGetNumberOfDimensions(lhs[0])!=2))
    stopMexFunction(302,0,0,0,0,0);

  if (!(((mxGetM(lhs[0])==d) && (mxGetN(lhs[0])==1)) ||
        ((mxGetM(lhs[0])==1) && (mxGetN(lhs[0])==d))))
    stopMexFunction(303,mxGetM(lhs[0]),mxGetN(lhs[0]),0,0,0);
    
  /* copy back */
  memcpy(f,mxGetPr(lhs[0]),d*sizeof(double));
  
  /* free memory */
  mxDestroyArray(lhs[0]);  
}

static int callOutputFcn (int reason, double tOld, double t, double *x, int doCopy) {
  int doPoint,erg,offset;
  double *dpointer;
  mxArray* rhs[5];
  mxArray* lhs[1];
    
  if (paramOutput.outputFcnH==NULL) return 0;
  doPoint=0;lhs[0]=NULL;erg=0;
  switch (reason)
  {
    case 1: /* Init-case */
      offset=0;
      switch (paramGlobal.funcCallMethod) {
        case 0: break;
        case 1:
          rhs[offset++]=(mxArray*)paramOutput.outputFcnH;
          break;
        default: stopMexFunction(1002,0,0,0,0,0);break;
      }
      rhs[offset]=mxCreateDoubleMatrix(1,2,mxREAL);
      dpointer=mxGetPr(rhs[offset]);
      *dpointer=paramDOPRI.tStart;dpointer++;
      *dpointer=paramDOPRI.tEnd;
      offset++;
      
      rhs[offset]=mxCreateDoubleMatrix(paramGlobal.d,1,mxREAL);
      memcpy(mxGetPr(rhs[offset]),paramDOPRI.xStart,(size_t)paramGlobal.d*sizeof(double));
      offset++;

      rhs[offset++]=mxCreateString("init");

      switch (paramGlobal.funcCallMethod) {
        case 0: 
          mexCallMATLAB(0,lhs,offset,rhs,paramOutput.outputFcn);
          while (--offset>=0) {mxDestroyArray(rhs[offset]);}
          break;
        case 1:
          mexCallMATLAB(0,lhs,offset,rhs,"feval");
          while (--offset>=1) {mxDestroyArray(rhs[offset]);}
          break;
        default: stopMexFunction(1002,0,0,0,0,0);break;
      }
      break;
    case 2: /* Done-case */
      offset=0;
      switch (paramGlobal.funcCallMethod) {
        case 0: break;
        case 1:
          rhs[offset++]=(mxArray*)paramOutput.outputFcnH;
          break;
        default: stopMexFunction(1002,0,0,0,0,0);break;
      }
      rhs[offset++]=mxCreateDoubleMatrix(0,0,mxREAL);
      rhs[offset++]=mxCreateDoubleMatrix(0,0,mxREAL);
      rhs[offset++]=mxCreateString("done");
      
      switch (paramGlobal.funcCallMethod) {
        case 0: 
          mexCallMATLAB(0,lhs,offset,rhs,paramOutput.outputFcn);
          while (--offset>=0) {mxDestroyArray(rhs[offset]);}
          break;
        case 1:
          mexCallMATLAB(0,lhs,offset,rhs,"feval");
          while (--offset>=1) {mxDestroyArray(rhs[offset]);}
          break;
        default: stopMexFunction(1002,0,0,0,0,0);break;
      }
      break;
    case 3: /* new Grid-node */
      if ((paramDOPRI.denseFlag) && (paramOutput.outputCallMode==1)) break;
      doPoint=1;
      break;
    case 4: /* new dense-node */
      if ((paramDOPRI.denseFlag) && (paramOutput.outputCallMode==2)) break;
      doPoint=1;
      break;
    case 5: /* new dense and Grid-node */
      doPoint=1;
      break;
    default: stopMexFunction(1001,0,0,0,0,0);
  }
  
  if (doPoint) {
    *mxGetPr(paramOutput.tArg)=t;
    *mxGetPr(paramOutput.toldArg)=tOld;
    if (doCopy)
      memcpy(mxGetPr(paramOutput.xArg),x,(size_t)paramGlobal.d*sizeof(double));

    offset=0;
    switch (paramGlobal.funcCallMethod) {
      case 0: break;
      case 1:
        rhs[offset++]=(mxArray*)paramOutput.outputFcnH;
        break;
      default: stopMexFunction(1002,0,0,0,0,0);break;
    }
    rhs[offset++]=paramOutput.tArg;rhs[offset++]=paramOutput.xArg;
    rhs[offset++]=paramOutput.emptyArg;
    rhs[offset++]=paramOutput.toldArg;

    switch (paramGlobal.funcCallMethod) {
      case 0:
        mexCallMATLAB(1,lhs,offset,rhs,paramOutput.outputFcn);
        break;
      case 1:
        mexCallMATLAB(1,lhs,offset,rhs,"feval");
        break;
      default: stopMexFunction(1002,0,0,0,0,0);break;
    }

    if (lhs[0]==NULL) {
      erg=0; 
    } else {
      erg=(int)mxGetScalar(lhs[0]);
      mxDestroyArray(lhs[0]);
    }
  }
  return erg;
}

void DOPRISoloutFunc (Fint *nr, double *told,
  double *t, double *x, Fint *n,
	double *con, Fint *icomp, Fint *nd,
	double *rpar, Fint *ipar, Fint *irtrn) {
  Fint erg,i,ilast;
  double *dpointer;
  double tdense;
  char alreadySaved;
  
  (void)(n); (void)rpar; (void)ipar;
  if (*nr==1) paramOutput.tPos=0;
  
  erg=0;
  if (paramDOPRI.denseFlag) {
    dopriDense.con=con;dopriDense.icomp=icomp;dopriDense.nd=nd;
    isInnerCall=1;
    alreadySaved=0;
    while (1) {
      if (paramOutput.tPos>=paramGlobal.tLength) break;
      tdense=paramGlobal.tPointer[paramOutput.tPos];
      if (tdense==*t) {
        addTXtoList(*t,x);paramOutput.tPos++;
        alreadySaved=1;
        erg=callOutputFcn(5,*told,*t,x,1);if ((erg!=0) && (erg!=1)) erg=0;
        if (erg==1) {erg=-1;break;} /* Stop, NOW */
        continue;
      }
      if (paramGlobal.direction*(tdense-(*t))>0) break;
      
      dpointer=mxGetPr(paramOutput.xArg);
      
      ilast=tif_mws2Fint(paramGlobal.d);
      for (i=1; i<=ilast; i++,dpointer++) {
        *dpointer=
#ifdef DOPRI5MexVersion
                  CONTD5_
#endif
#ifdef DOP853MexVersion
                  CONTD8_ 
#endif
        (&i,&tdense,con,icomp,nd);
      }
        
      dpointer=mxGetPr(paramOutput.xArg);
      addTXtoList(tdense,dpointer);
      erg=callOutputFcn(4,*told,tdense,dpointer,0);if ((erg!=0) && (erg!=1)) erg=0;
      
      if (erg==1) {erg=-1;break;} /* Stop, NOW */
      
      paramOutput.tPos++;
    }
    if (!alreadySaved) {
      if (paramOutput.includeGrid) addTXtoList(*t,x);
      erg=callOutputFcn(3,*told,*t,x,1);if ((erg!=0) && (erg!=1)) erg=0;
      if (erg==1) {erg=-1;} /* Stop, NOW */
    }
    isInnerCall=0;
    dopriDense.con=NULL;dopriDense.icomp=NULL;dopriDense.nd=NULL;
  } else {
    addTXtoList(*t,x);
    erg=callOutputFcn(3,*told,*t,x,1);if ((erg!=0) && (erg!=1)) erg=0;
    if (erg==1) {erg=-1;} /* Stop, NOW */
  }
  
  *irtrn=erg;  
}

static void createTXArrays (mxArray* plhs[]) {
  mwSize count,i,n,d;
  PListElement current;
  double *tPointer, *xPointer, *vPointer;
  
  d=paramGlobal.d;
  n=paramOutput.numberOfElements;
  plhs[0]=mxCreateDoubleMatrix(n,1,mxREAL);
  plhs[1]=mxCreateDoubleMatrix(n,d,mxREAL);
  
  tPointer=mxGetPr(plhs[0]);
  xPointer=mxGetPr(plhs[1]);
  
  current=paramOutput.txList.next;
  count=0;
  while (current!=NULL) {
    vPointer=current->values;
    *tPointer= *vPointer; tPointer++; vPointer++;
    for (i=0; i<d; i++,vPointer++) xPointer[count+i*n]= *vPointer;
    
    current=current->next;count++;
  }
}

static void createStatVector (mxArray* plhs[]) {
  double* dpointer;

  plhs[2]=mxCreateDoubleMatrix(1,5,mxREAL);
  dpointer=mxGetPr(plhs[2]);
  
  *dpointer=(double)paramDOPRI.IDID;dpointer++;
  *dpointer=(double)paramDOPRI.IWORK[17-1];dpointer++;
  *dpointer=(double)paramDOPRI.IWORK[18-1];dpointer++;
  *dpointer=(double)paramDOPRI.IWORK[19-1];dpointer++;
  *dpointer=(double)paramDOPRI.IWORK[20-1];
}

static void createHPred (mxArray* plhs[]) {  
  plhs[3]=mxCreateDoubleMatrix(1,1,mxREAL);
  
  *mxGetPr(plhs[3]) = paramDOPRI.WORK[7-1];
}

static void createOutput (int nlhs, mxArray* plhs[]) {
  createTXArrays(plhs);
  
  if (nlhs>=3) createStatVector(plhs);
  if (nlhs>=4) createHPred(plhs);
}

static void doInnerCall (int nlhs, mxArray *plhs[],
                         int nrhs, const mxArray *prhs[]) {
  /* So jetzt wird es trickreich: Wenn wir jetzt
       stopMexFunction(...);
     aufrufen, so ist das für Matlab tödlich!
     Begründung:
     Schauen wir uns jetzt einmal den Stack in diesem Moment an:
       dopriMex (InnerCall)
       Matlab m-file
       dopridMex (NormalCall)
       Matlab (bsp. Command-Line)
     Wenn ich jetzt stopMexFunction mache, dann wird ALLES
     von dopriMex aufgeräumt und dann mexErrMsgTxt aufgerufen.
     Matlab will jetzt gepfelgt den kompletten Stack von oben
     nach unten abarbeiten und den Speicher aufräumen. Aber 
     da wären wir ja schon schneller gewesen, so dass Matlab
     auf dem Stack ungültige Speicherhandles findet ... und sich
     damit verabschiedet!! 
     Deshalb beenden wir uns jetzt, OHNE unseren Speicher 
     aufzuräumen. Damit verlassen wir uns darauf, dass Matlab
     unseren ganzen Sch... aufräumt! Dies sollte problemlos funktionieren,
     aber ich habe mich noch NIE auf sowas verlassen, bis jetzt ... */

  Fint i,ilast;
  double t;
  double *dpointer;
  const mxArray *tArg;
  mxArray *xArg;

  if ((nlhs!=1) || (nrhs!=1)) stopMexFunctionImpl(501,0,0,0,0,0.0,0);

  tArg=prhs[0];
  if ((tArg==NULL) || (mxIsEmpty(tArg))) stopMexFunctionImpl(502,0,0,0,0,0.0,0);
  if (!mxIsDouble(tArg) || (mxIsSparse(tArg)) || 
     (mxGetNumberOfDimensions(tArg)!=2) || 
     (mxGetM(tArg)!=1) || (mxGetN(tArg)!=1)) 
    stopMexFunctionImpl(503,0,0,0,0,0.0,0);
  t= *mxGetPr(tArg);

  xArg=mxCreateDoubleMatrix(1,paramGlobal.d,mxREAL);
  dpointer=mxGetPr(xArg);
  ilast=tif_mws2Fint(paramGlobal.d);
  for (i=1; i<=ilast; i++,dpointer++) {
    *dpointer=
#ifdef DOPRI5MexVersion
                  CONTD5_
#endif
#ifdef DOP853MexVersion
                  CONTD8_
#endif
      (&i,&t,dopriDense.con,dopriDense.icomp,dopriDense.nd);
  }

  plhs[0]=xArg;
}

void mexFunction (int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[]) {
  Fint d;

  /* check for inner call during Solout routine */
  if (isInnerCall && (nrhs==1)) {doInnerCall(nlhs,plhs,nrhs,prhs);return;}

  tif_checkFint();

  initVars();
  
  checkNumberOfArgs(nlhs,nrhs);
  processArgs(nrhs,prhs);
  
  extractOptionsPart1();
  prepareWorkArrays();
  extractOptionsPart2();
  prepareHelpMxArrays();
  
  d=tif_mws2Fint(paramGlobal.d);
  callOutputFcn(1,0.0,0.0,NULL,0);
#ifdef DOPRI5MexVersion
  DOPRI5_ 
#endif
#ifdef DOP853MexVersion
  DOP853_ 
#endif
  (&d,&DOPRIRightSideFunc,&paramDOPRI.tStart,
    paramDOPRI.xStart,&paramDOPRI.tEnd,
    paramDOPRI.RTOL,paramDOPRI.ATOL,&paramDOPRI.ITOL,
    &DOPRISoloutFunc,&paramDOPRI.IOUT,
    paramDOPRI.WORK,&paramDOPRI.LWORK,
    paramDOPRI.IWORK,&paramDOPRI.LIWORK,
    paramDOPRI.RPAR,paramDOPRI.IPAR,&paramDOPRI.IDID);
  callOutputFcn(2,0.0,0.0,NULL,0);
  
  createOutput(nlhs,plhs);
  
  doneVars();
}
