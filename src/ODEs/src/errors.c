/* description of errors, included in dopri5Mex.c and dop853Mex */
/* Fehlerbeschreibungen, werden in dopri5Mex.c und dop853Mex.c included */

/* general Errors */
/* allgemeine Fehler */
case 1:
  tif_printInfos();
  mexPrintf("English: \n");
  mexPrintf("Only 2,3 or 4 output arguments are possible, but %" FMT_SIZE_T "u were requested.\n",i1);
  mexPrintf("German: \n");
  mexPrintf("Es werden 2,3 oder 4 Ausgabeargumente unterstützt. Verlangt wurden\n");
  mexPrintf("aber %" FMT_SIZE_T "u Ausgabeargumente.\n",i1);
  msg="Invalid number of output arguments (Ungültige Anzahl von Ausgabeargumenten)";break;
case 2:
  mexPrintf("English: \n");
  mexPrintf("Only 3 or 4 input arguments are possible, but %" FMT_SIZE_T "u were passed.\n",i1);
  mexPrintf("German: \n");
  mexPrintf("Es werden 3 oder 4 Eingabeargumente unterstützt. Es wurden aber\n");
  mexPrintf("%" FMT_SIZE_T "u Argumente übergeben.\n",i1);
  msg="Invalid number of input arguments (Ungültige Anzahl von Eingabeargumenten)";break;
case 3:
  mexPrintf("English: \n");
  mexPrintf("1st argument has to be a string (function name), a function handle\n");
  mexPrintf("or an inline function.\n");
  mexPrintf("German: \n");
  mexPrintf("1. Argument muss ein String (Funktionsname), ein Funktions-Handle\n");
  mexPrintf("oder eine inline-Funktion sein.\n");
  msg="1. Arg: Function (String,handle,inline) expected";break;
case 4:
  mexPrintf("English: \n");
  mexPrintf("1st argument has characters, but needs to be ONE string. ");
  if (i1!=2) 
    mexPrintf("Not a 0-, 1- or >=3-dimensional thing.");
  else
    if (i2!=1) mexPrintf("Not a matrix containing strings.");
  mexPrintf("\n");
  mexPrintf("German: \n");
  mexPrintf("1. Argument enthält Zeichen, muss aber EIN String sein. ");
  if (i1!=2) 
    mexPrintf("Kein 0-,1- oder >=3-dimensionales Gebilde."); 
  else
    if (i2!=1) mexPrintf("Keine Matrix aus Strings.");
  mexPrintf("\n");
  msg="1. Arg: only ONE string expected (nur EIN String erwartet)";break;
case 5:
  mexPrintf("English: \n");
  mexPrintf("2nd argument has to be a double row vector.\n");
  mexPrintf("German: \n");
  mexPrintf("2. Argument muss ein double-Zeilenvektor sein.\n");
  msg="2. Arg: double row vector expected (double-Vektor erwartet)";break;
case 6:
  mexPrintf("English: \n");
  mexPrintf("2nd argument has to be a double row vector. ");
  if (i1!=2) 
    mexPrintf("Not a 0-, 1- or >=3-dimensional thing.");
  else
    if (i2!=1) mexPrintf("Not a double matrix");
  mexPrintf("\n");
  mexPrintf("German: \n");
  mexPrintf("2. Argument muss ein double-Vektor sein. ");
  if (i1!=2)
    mexPrintf("Kein 0-,1- oder >=3-dimensionales Gebilde."); 
  else
    if (i2!=1) mexPrintf("Keine double-Matrix.");
  mexPrintf("\n");
  msg="2. Arg: (1,k) double vector expected ((1,k) double-Vektor erwartet)";break;
case 7:
  mexPrintf("English: \n");
  mexPrintf("2nd argument has to be a double vector with at least 2 components.\n");
  mexPrintf("The length of the vector found is %" FMT_SIZE_T "u.\n",paramGlobal.tLength);
  mexPrintf("German: \n");
  mexPrintf("2. Argument muss double-Vektor mit mindestens der\n");
  mexPrintf("Länge 2 sein. Der übergebene Vektor hat Länge %" FMT_SIZE_T "u.\n",
    paramGlobal.tLength);
  msg="2. Arg: (1,k) double vector (k>=2) expected ((1,k) double-Vektor (k>=2) erwartet)";break;
case 8:
  mexPrintf("English: \n");
  mexPrintf("3rd argument has to be a double column vector.\n");
  mexPrintf("German: \n");
  mexPrintf("3. Argument muss ein double-Spaltenvektor sein.\n");
  msg="3. Arg: double vector expected (double-Vektor erwartet)";break;
case 9:
  mexPrintf("English: \n");
  mexPrintf("3rd argument has to be a double column vector.\n");
  if (i1!=2) 
    mexPrintf("Not a 0-, 1- or >=3-dimensional thing.");
  else
    if (i2!=1) mexPrintf("Not a double matrix");
  mexPrintf("\n");
  mexPrintf("German: \n");
  mexPrintf("3. Argument muss ein double-Spaltenvektor sein. ");
  if (i1!=2)
    mexPrintf("Kein 0-,1- oder >=3-dimensionales Gebilde."); 
  else
    if (i2!=1) mexPrintf("Keine double-Matrix.");
  mexPrintf("\n");
  msg="3. Arg: (d,1) double vector expected ((d,1) double-Vektor erwartet)";break;
case 10:
  mexPrintf("English: \n");
  mexPrintf("3rd argument has to be a double column vector with at least one component.\n");
  mexPrintf("The length of the vector found was %" FMT_SIZE_T "u.\n",paramGlobal.d);
  mexPrintf("German: \n");
  mexPrintf("3. Argument muss ein double-Spaltenvektor mit mindestens\n");
  mexPrintf("der Länge 1 sein. Der übergebene Vektor hat Länge %" FMT_SIZE_T "u.\n",paramGlobal.d);
  msg="3. Arg: (d,1) double vector (d>=1) expected ((d,1) double-Vektor (d>=1) erwartet)";break;    
case 11:
  mexPrintf("English: \n");
  mexPrintf("4th argument has to be a struct.\n");
  mexPrintf("German: \n");
  mexPrintf("4. Argument muss eine struct sein.\n");
  msg="4. Arg: struct expected (struct erwartet)";break;
case 12:
  mexPrintf("English: \n");
  mexPrintf("3rd arg: start- and end-time are the same.\n");
  mexPrintf("German: \n");
  mexPrintf("3. Argument: Start- und Endzeitpunkt dürfen nicht\n");
  mexPrintf("übereinstimmen!\n");
  msg="3. Arg: tStart==tEnd";break;  
case 14:
  mexPrintf("English: \n");
  mexPrintf("2nd argument has to be ordered.\n");
  if (i1>0)
    mexPrintf("Because of %f=tStart<tEnd=%f, the 2nd argument has to be ascending.\n",
      paramDOPRI.tStart,paramDOPRI.tEnd);
  else
    mexPrintf("Because of %f=tStart>tEnd=%f, the 2nd argument has to be descending.\n",
      paramDOPRI.tStart,paramDOPRI.tEnd);
  mexPrintf("German: \n");
  mexPrintf("2. Argument muss sortiert sein.\n");
  if (i1>0)
    mexPrintf("Da %f=tStart<tEnd=%f ist, muss das zweite Argument steigend\n",
     paramDOPRI.tStart,paramDOPRI.tEnd); 
  else    
    mexPrintf("Da %f=tStart>tEnd=%f ist, muss das zweite Argument fallend\n",
     paramDOPRI.tStart,paramDOPRI.tEnd); 
  mexPrintf("sortiert sein.\n");
  msg="2. Arg: has to be ordered (muss sortiert sein)";break;
case 15:
  mexPrintf("English: \n");
  mexPrintf("1st argument contained inline-funcs or function handles, but I need only ONE. ");
  if (i1!=2) 
    mexPrintf("Not a 0-, 1- or >=3-dimensional thing.");
  else
    if ((i2!=1) || (i3!=1)) 
      mexPrintf("Not a matrix containing inline-funcs or function handles.");
  mexPrintf("\n");
  mexPrintf("German: \n");
  mexPrintf("1. Argument enthält zwar Inline-Funktionen oder Funktions-Handles, aber");
  mexPrintf("benötigt wird nur EIN Handle oder EINE Inline-Funktion. ");
  if (i1!=2) 
    mexPrintf("Kein 0-,1- oder >=3-dimensionales Gebilde."); 
  else
    if ((i2!=1) || (i3!=1)) 
      mexPrintf("Keine Matrix mit Inline-Funktionen oder Funktions-Handles .");
  mexPrintf("\n");
  msg="1. Arg: only ONE function expected (nur EINE Funktion erwartet)";break;
case 16:
  mexPrintf("English: \n");
  mexPrintf("Concering Option '%s':\n",OPT_FUNCCALLMETHOD);
  mexPrintf("Requirement: '%s'==0 or '%s'==1\n",OPT_FUNCCALLMETHOD,OPT_FUNCCALLMETHOD);
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_FUNCCALLMETHOD);
  mexPrintf("Es muss gelten: '%s'==0 or '%s'==1\n",OPT_FUNCCALLMETHOD,OPT_FUNCCALLMETHOD);
  msg="Invalid call method (ungültige Methode für Funktionsaufruf)";break;
case 17:
  mexPrintf("English: \n");
  mexPrintf("Concering Option '%s':\n",OPT_FUNCCALLMETHOD);
  mexPrintf("'%s'==0 was chosen. Hence all Matlab-functions must be given as string.\n",
    OPT_FUNCCALLMETHOD);
  mexPrintf("But the rightSide was not a string.\n");
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_FUNCCALLMETHOD);
  mexPrintf("Es wurde '%s'==0 gewählt. Also müssen alle Matlab-Funktionen als String angegeben werden.\n",
    OPT_FUNCCALLMETHOD);
  mexPrintf("Aber die rechte Seite war kein String mit Funktionsname.\n");
  msg="CallMethod 0 => all Matlab-Funcs must be given as Strings (Bei Aufrufmethode 0: müssen Funktionsnamen übergeben werden)";break;
case 18:
  mexPrintf("English: \n");
  mexPrintf("Concering Option '%s':\n",OPT_FUNCCALLMETHOD);
  mexPrintf("'%s'==0 was chosen. Hence the all Matlab-functions must be given as string.\n",
    OPT_FUNCCALLMETHOD);
  mexPrintf("But the Output Function was not a string.\n");
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_FUNCCALLMETHOD);
  mexPrintf("Es wurde '%s'==0 gewählt. Also müssen alle Matlab-Funktionen als String angegeben werden.\n",
    OPT_FUNCCALLMETHOD);
  mexPrintf("Aber die Output Funktion war kein String mit Funktionsname.\n");
  msg="CallMethod 0 => all Matlab-Funcs must be given as Strings (Bei Aufrufmethode 0: müssen Funktionsnamen übergeben werden)";break;
case 19:
  mexPrintf("English: \n");
  mexPrintf("3rd argument is too large (length is %" FMT_SIZE_T "u). I was compiled\n",
    i1);
  mexPrintf("with a maximal vector size of %" FMT_SIZE_T "u. Try recompiling\n",i2);
  mexPrintf("with largeArrayDims turned on.\n");
  mexPrintf("German: \n");
  mexPrintf("3. Argument ist zu groß (Länge ist % " FMT_SIZE_T "u). Ich wurde für \n",
    i1);
  mexPrintf("Vektoren mit maximal %" FMT_SIZE_T "u Einträge kompiliert. Versuchen\n",i2);
  mexPrintf("Sie ein erneutes Kompilieren mit aktivierten largeArrayDims.\n");
  msg="Dimension too large for mwSize (Dimension zu groß für mwSize)";break;
case 20:
  mexPrintf("English: \n");
  mexPrintf("2nd argument is too large (length is %" FMT_SIZE_T "u). I was compiled\n",
    i1);
  mexPrintf("with a maximal vector size of %" FMT_SIZE_T "u. Try recompiling\n",i2);
  mexPrintf("with largeArrayDims turned on.\n");
  mexPrintf("German: \n");
  mexPrintf("2. Argument ist zu groß (Länge ist % " FMT_SIZE_T "u). Ich wurde für \n",
    i1);
  mexPrintf("Vektoren mit maximal %" FMT_SIZE_T "u Einträge kompiliert. Versuchen\n",i2);
  mexPrintf("Sie ein erneutes Kompilieren mit aktivierten largeArrayDims.\n");
  msg="Dimension too large for mwSize (Dimension zu groß für mwSize)";break;
case 21:
  mexPrintf("English: \n");
  mexPrintf("Concering Option '%s':\n",OPT_RTOL);
  mexPrintf("'%s' was too large (length is %" FMT_SIZE_T "u) for mwSize.\n",
    OPT_RTOL,i1);
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_RTOL);
  mexPrintf("'%s'==0 war zu groß (Länge ist %" FMT_SIZE_T "u) für mwSize.\n",
    OPT_RTOL,i1);
  msg="Dimension too large for mwSize (Dimension zu groß für mwSize)";break;
case 22:
  mexPrintf("English: \n");
  mexPrintf("Concering Option '%s':\n",OPT_ATOL);
  mexPrintf("'%s' was too large (length is %" FMT_SIZE_T "u) for mwSize.\n",
    OPT_RTOL,i1);
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_ATOL);
  mexPrintf("'%s'==0 war zu groß (Länge ist %" FMT_SIZE_T "u) für mwSize.\n",
    OPT_RTOL,i1);
  msg="Dimension too large for mwSize (Dimension zu groß für mwSize)";break;

/* Errors concering WORK and IWORK-options */
/* Fehler bei WORK und IWORK-Options */
case 101:
  mexPrintf("English: \n");
  mexPrintf("Concerning Option '%s':\n",OPT_MAXSTEPS);
  mexPrintf("Requirement: 0<'%s'.\n",OPT_MAXSTEPS);
  mexPrintf("But I found '%s'=%i.\n",OPT_MAXSTEPS,i4);
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_MAXSTEPS);
  mexPrintf("Es muss gelten: 0<'%s'.\n",OPT_MAXSTEPS);
  mexPrintf("Gefunden wurde aber '%s'=%i.\n",OPT_MAXSTEPS,i4);
  msg="invalid maximal number of steps (ungültige maximale Schrittanzahl)";break;
case 102:
  mexPrintf("English: \n");
  mexPrintf("Concerning Option '%s':\n",OPT_STEST);
  mexPrintf("Requirement: '%s'!=0.\n",OPT_STEST);
  mexPrintf("But I found '%s'=%i.\n",OPT_STEST,i4);
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_STEST);
  mexPrintf("Es muss gelten: '%s'!=0.\n",OPT_STEST);
  mexPrintf("Gefunden wurde aber '%s'=%i.\n",OPT_STEST,i4);
  msg="invalid value for steffness test (ungültige Angabe für Stiffness-Test)";break;
case 103:
  mexPrintf("English: \n");
  mexPrintf("Concerning Option '%s':\n",OPT_EPS);
  mexPrintf("Requirement: 1e-35<'%s'<1.0.\n",OPT_EPS);
  mexPrintf("But I found '%s'=%e.\n",OPT_EPS,d1);
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_EPS);
  mexPrintf("Es muss gelten: 1e-35<'%s'<1.0.\n",OPT_EPS);
  mexPrintf("Gefunden wurde aber '%s'=%e.\n",OPT_EPS,d1);
  msg="invalid precision (ungültige Maschinengenauigkeit)";break;
case 104:
  mexPrintf("English: \n");
  mexPrintf("Concerning Option '%s':\n",OPT_RHO);
  mexPrintf("Requirement: 1e-4<'%s'<1.\n",OPT_RHO);
  mexPrintf("But I found '%s'=%e.\n",OPT_RHO,d1);
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_RHO);
  mexPrintf("Es muss gelten: 1e-4<'%s'<1.\n",OPT_RHO);
  mexPrintf("Gefunden wurde aber '%s'=%e.\n",OPT_RHO,d1);
  msg="invalid value for safety factor in step size prediction (ungültiger Sicherheitsfaktor bei der Schrittweitensteuerung)";break;
case 105:
  mexPrintf("English: \n");
  mexPrintf("Concerning Option '%s':\n",OPT_SSMINSEL);
  mexPrintf("Requirement: 0<'%s'.\n",OPT_SSMINSEL);
  mexPrintf("But I found '%s'=%e.\n",OPT_SSMINSEL,d1);
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_SSMINSEL);
  mexPrintf("Es muss gelten: 0<'%s'.\n",OPT_SSMINSEL);
  mexPrintf("Gefunden wurde aber '%s'=%e.\n",OPT_SSMINSEL,d1);
  msg="invalid lower bound for step size prediction (ungültige untere Grenze bei der Schrittweitensteuereung)";break;
case 106:
  mexPrintf("English: \n");
  mexPrintf("Concerning Option '%s':\n",OPT_SSMAXSEL);
  mexPrintf("Requirement: 0<'%s'.\n",OPT_SSMAXSEL);
  mexPrintf("But I found '%s'=%e.\n",OPT_SSMAXSEL,d1);
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_SSMAXSEL);
  mexPrintf("Es muss gelten: 0<'%s'.\n",OPT_SSMAXSEL);
  mexPrintf("Gefunden wurde aber '%s'=%e.\n",OPT_SSMAXSEL,d1);
  msg="invalid upper bound for step size prediciton (ungültige obere Grenze bei der Schrittweitensteuereung)";break;
case 107:
  mexPrintf("English: \n");
  mexPrintf("Concering Option '%s':\n",OPT_SSBETA);
  mexPrintf("Requirement: '%s'<=0.2.\n",OPT_SSBETA);
  mexPrintf("But I found '%s'=%e.\n",OPT_SSBETA,d1);
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_SSBETA);
  mexPrintf("Es muss gelten: '%s'<=0.2.\n",OPT_SSBETA);
  mexPrintf("Gefunden wurde aber '%s'=%e.\n",OPT_SSBETA,d1);
  msg="invalid SC-stabilizing factor (ungültiger SC-Stabilisierungsfaktor)";break;
case 108:
  mexPrintf("English: \n");
  mexPrintf("Concerning Option '%s':\n",OPT_MAXSS);
  mexPrintf("Requirement: '%s'!=0.\n",OPT_MAXSS);
  mexPrintf("But I found '%s'=%e.\n",OPT_MAXSS,d1);
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_MAXSS);
  mexPrintf("Es muss gelten: '%s'!=0.\n",OPT_MAXSS);
  mexPrintf("Gefunden wurde aber '%s'=%e.\n",OPT_MAXSS,d1);
  msg="invalid maximal number of steps (ungültige maximale Schrittweite)";break;

/* Errors concering rightSide */
/* Fehler im Zusammenhang mit der rechten Seite */
case 301:
  mexPrintf("English: \n");
  mexPrintf("Problem calling the right side.\n");
  mexPrintf("I have not got return values.\n");
  mexPrintf("German: \n");
  mexPrintf("Problem beim Aufruf der rechten Seite.\n");
  mexPrintf("Habe keine Rückgabe erhalten.\n");
  msg="Right side without return values (Rechte Seite ohne Rückgabe)";break;
case 302:
  mexPrintf("English: \n");
  mexPrintf("Problem calling the right side.\n");
  mexPrintf("Return value was not a double vector.\n");
  mexPrintf("German: \n");
  mexPrintf("Problem beim Aufruf der rechten Seite.\n");
  mexPrintf("Die Rückgabe war kein double-Vektor.\n");
  msg="Return value was not a double vector (Rückgabe der rechten Seite kein double-Vektor)";break;
case 303:
  mexPrintf("English: \n");
  mexPrintf("Problem calling the right side.\n");
  mexPrintf("The return value was not a (d,1) or (1,d) double Vector.\n");
  mexPrintf("A (%" FMT_SIZE_T "u,%" FMT_SIZE_T "u) double matrix was returned.\n",i1,i2);
  mexPrintf("German: \n");
  mexPrintf("Problem beim Aufruf der rechten Seite.\n");
  mexPrintf("Die Rückgabe war kein (d,1) oder (1,d) double-Vektor.\n");
  mexPrintf("Es wurde eine (%" FMT_SIZE_T "u,%" FMT_SIZE_T "u) double-Matrix zurückgegeben.\n",i1,i2);
  msg="Return value of right side was not a (d,1) or (1,d) double vector (Rückgabe der rechten Seite kein (d,1) oder (1,d) double-Vektor)";break;

/* Errors concering output options */
/* Fehler im Zusammenhang mit Output Options */
case 401:
  mexPrintf("English: \n");
  mexPrintf("Concering option '%s':\n",OPT_OUTPUTFUNCTION);
  mexPrintf("Option contained inline-funcs or function handles, but I need only ONE. ");
  if (i1!=2) 
    mexPrintf("Not a 0-, 1- or >=3-dimensional thing.");
  else
    if ((i2!=1) || (i3!=1)) 
      mexPrintf("Not a matrix containing inline-funcs or function handles.");
  mexPrintf("\n");
  mexPrintf("German: \n");
  mexPrintf("Zur Option '%s':\n",OPT_OUTPUTFUNCTION);
  mexPrintf("Option enthält zwar Inline-Funktionen oder Funktions-Handles, aber");
  mexPrintf("benötigt wird nur EIN Handle oder EINE Inline-Funktion. ");
  if (i1!=2) 
    mexPrintf("Kein 0-,1- oder >=3-dimensionales Gebilde."); 
  else
    if ((i2!=1) || (i3!=1)) 
      mexPrintf("Keine Matrix mit Inline-Funktionen oder Funktions-Handles .");
  mexPrintf("\n");
  msg="only ONE function expected (nur EINE Funktion erwartet)";break;

/* Errors in inner Call during SolOut routine */
/* Fehler während des inner Calls zum Zeitpunkt der SolOut Routine */
case 501:
  mexPrintf("English:\n");
  mexPrintf("Only one input and output argument is allowed\n");
  mexPrintf("in a call during the output function.\n");
  mexPrintf("\n");
  mexPrintf("German:\n");
  mexPrintf("Bei Aufrufen innerhalb der Output Function wird nur\n");
  mexPrintf("ein Eingabe- und ein Ausgabeargument unterstützt.\n");
  msg="wrong number of input and/or output parameters";break;
case 502:
  mexPrintf("English:\n");
  mexPrintf("1st input argument in the call during the ouput function is empty.\n");
  mexPrintf("\n");
  mexPrintf("German:\n");
  mexPrintf("Das erste Eingabeargument während des Aufrufs aus der Output Function ist leer.\n");
  msg="1st arg was empty at call inside output function";break;
case 503:
  mexPrintf("English:\n");
  mexPrintf("1st input argument in the call during the ouput function\n");
  mexPrintf("must be a scalar.\n");
  mexPrintf("\n");
  mexPrintf("German:\n");
  mexPrintf("Das erste Eingabeargument während des Aufrufs aus der\n");
  mexPrintf("Output Function muss ein Skalar sein.\n");
  msg="1st arg at call inside output function must be a scalar";break;

/* Errors that should never occur */
/* Fehler, die niemals auftreten sollten */
case 1001:
  msg="internal error: callOutputFcn: unknown reason";
  break;
case 1002:
  msg="internal error: unknown funcCallMethod";
  break;

