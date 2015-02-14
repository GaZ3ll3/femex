#include <string.h>
#include <stdint.h>
#include "options.h"

char opt_getSizeOfOptField (const mxArray *opt, const char *name, size_t *m, size_t *n)
{
  mxArray *field;
  
  field=mxGetField(opt,(mwIndex)0,name);
  if (field==NULL) return (char)1;
  
  if (mxGetNumberOfDimensions(field)!=(size_t)2) return (char)2;
  
  if (mxIsEmpty(field)) return (char)3;
  
  *m=mxGetM(field);*n=mxGetN(field);
  return 0;
}

double opt_getDoubleFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, double defaultValue)
{
  mxArray *field;
  
  field=mxGetField(opt,(mwIndex)0,name);
  if ((field==NULL) ||  (mxIsEmpty(field)))
  {
    if (optSet->warnMiss) 
      mexPrintf("Option '%s' fehlt, nehme default %e\n",name,defaultValue);    
    return defaultValue;
  }
  if (!mxIsDouble(field))
  {
    if (optSet->warnType)
      mexPrintf("Option '%s' hat falschen Typ, nehme default %e\n",name,defaultValue);
    return defaultValue;
  }
  if ((mxGetNumberOfDimensions(field)!=(mwSize)2) ||
      (mxGetM(field)!=(size_t)1) || (mxGetN(field)!=(size_t)1))
  {
    if (optSet->warnSize)
      mexPrintf("Option '%s' hat falsche Größe, nehme default %e\n",name,defaultValue);
    return defaultValue;
  }
  return mxGetScalar(field);
}

char opt_getDoubleVectorFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, size_t m, size_t n, double *dpointer)
{
  size_t l;
  mxArray *field;
  
  mxAssert(m==1 || n==1,"getDoubleVectorFromOpt: m!=1 && n!=1");

  field=mxGetField(opt,(mwIndex)0,name);
  if ((field==NULL) ||  (mxIsEmpty(field)))
  {
    if (optSet->warnMiss) mexPrintf("Option '%s' fehlt\n",name);
    return (char)1;
  }
  if (!mxIsDouble(field))
  {
    if (optSet->warnType) mexPrintf("Option '%s' hat falschen Typ\n",name);
    return (char)2;
  }
  if ((mxGetM(field)!=(size_t)m) || (mxGetN(field)!=(size_t)n))
  {
    if (optSet->warnSize) mexPrintf("Option '%s' hat falsche Größe\n");
    return (char)3;
  }
  if (m>n) l=m; else l=n;
  
  memcpy(dpointer,mxGetPr(field),l*sizeof(double));
  
  return (char)0;
}

int opt_getIntFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, int defaultValue)
{
  mxArray *field;
  mxClassID classID;
  void *data;
  
  field=mxGetField(opt,(mwIndex)0,name);
  if ((field==NULL) ||  (mxIsEmpty(field)))
  {
    if (optSet->warnMiss) 
      mexPrintf("Option '%s' fehlt, nehme default %i\n",name,defaultValue);    
    return defaultValue;
  }
  
  classID=mxGetClassID(field);
  switch (classID)
  {
    case mxDOUBLE_CLASS: case mxINT8_CLASS:   case mxUINT8_CLASS:
    case mxINT16_CLASS:  case mxUINT16_CLASS: case mxINT32_CLASS:
    case mxUINT32_CLASS: case mxINT64_CLASS:  case mxUINT64_CLASS:
    break;
    default:
    if (optSet->warnType)
      mexPrintf("Option '%s' hat falschen Typ, nehme default %i\n",
        name,defaultValue);
    return defaultValue;
  }
    
  if ((mxGetNumberOfDimensions(field)!=(mwSize)2) ||
     (mxGetM(field)!=(size_t)1) || (mxGetN(field)!=(size_t)1))
  {
    if (optSet->warnSize)
      mexPrintf("Option '%s' hat falsche Größe, nehme default %i\n",
        name,defaultValue);
    return defaultValue;
  } 
  
  data=mxGetData(field);
  switch (classID)
  {
    case mxDOUBLE_CLASS: return (int)*((double*)  data);break;
    case mxINT8_CLASS:   return (int)*((int8_t*)  data);break;
    case mxUINT8_CLASS:  return (int)*((uint8_t*) data);break;
    case mxINT16_CLASS:  return (int)*((int16_t*) data);break;
    case mxUINT16_CLASS: return (int)*((uint16_t*)data);break;
    case mxINT32_CLASS:  return (int)*((int32_t*) data);break;
    case mxUINT32_CLASS: return (int)*((uint32_t*)data);break;
    case mxINT64_CLASS:  return (int)*((int64_t*) data);break;
    case mxUINT64_CLASS: return (int)*((uint64_t*)data);break;
    default:
      return defaultValue;
  }

  mxAssert(0,"getIntFromOpt: should not happen");
  return 0;
}

char* opt_getStringFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, char* defaultValue)
{
  static const char* noString = "<NULL>";
  
  mxArray *field;
  char takeDefault;
  size_t buflen;
  char* erg;

  takeDefault=(char)0;
  field=mxGetField(opt,(mwIndex)0,name);
  if ((field==NULL) ||  (mxIsEmpty(field)))
  {
    if (optSet->warnMiss) 
      mexPrintf("Option '%s' fehlt, nehme default '%s'\n",name,
        defaultValue==NULL?noString:defaultValue);
    takeDefault=(char)1;
  } 
  
  if ((!takeDefault) && (!mxIsChar(field)))
  {
    if (optSet->warnType)
      mexPrintf("Option '%s' hat falschen Typ, nehme default %s\n",
        name,defaultValue==NULL?noString:defaultValue);
    takeDefault=(char)1;
  }
  
  if (!takeDefault)
  {
    if ((mxGetNumberOfDimensions(field)!=(mwSize)2) || (mxGetM(field)!=(size_t)1))
    {
      if (optSet->warnSize)
        mexPrintf("Option '%s' hat falsche Größe, nehme default %s\n",
          name,defaultValue==NULL?noString:defaultValue);
      takeDefault=(char)1;
    }
  } 
  
  if (takeDefault)
  {
    if (defaultValue==NULL) return NULL;
    buflen=strlen(defaultValue)+(size_t)1;
    erg=mxMalloc(buflen);
    strcpy(erg,defaultValue);
  } 
  else
  {
    erg=mxArrayToString(field);
  }
  return erg;
}
