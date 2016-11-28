#ifndef optionsh
#define optionsh

#include "mex.h"

/* PREFIX: opt */

struct OptionSettings 
{ /* mögliche Einstellungen (Warnstufe, usw.) bei Optionen */
  int warnMiss;  /* Warnung, wenn Optionseintrag fehlt */
  int warnType;  /* Warnung, wenn Optionseintrag falschen Typ */
  int warnSize;  /* Warnung, wenn Optionseintrag falsche Größe */
};
typedef struct OptionSettings SOptionSettings;
typedef SOptionSettings* POptionSettings;

char opt_getSizeOfOptField (const mxArray *opt, const char *name, size_t *m, size_t *n);

double opt_getDoubleFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, double defaultValue);
  
char opt_getDoubleVectorFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, size_t m, size_t n, double *dpointer);
  
int opt_getIntFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, int defaultValue);

char* opt_getStringFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, char* defaultValue);


#endif
