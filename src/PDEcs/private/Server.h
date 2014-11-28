/*
 * Server.h
 *
 *  Created on: Nov 22, 2014
 *      Author: lurker
 */

#ifndef SRC_PDECS_PRIVATE_SERVER_H_
#define SRC_PDECS_PRIVATE_SERVER_H_


#include <cstdlib>
#include <vector>
#include <cmath>

#include <iostream>
#include <iterator>
#include <string.h>

#include <mexplus.h>
#include <pprint.h>

#include "utils.h"

using namespace std;
using namespace mexplus;


#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

namespace MEX {

class Server {
public:
	Server();
	virtual ~Server();

	/*
	 * public methods
	 */

	// default inner product
	Real_t Inner_Prod(Real_t* u, Real_t* v, size_t n);
	Real_t Inner_Prod_omp(Real_t* u, Real_t* v, size_t n);


	// calculate gradient from inner product

	// Real_t* to int32_t* have to cast.
	void Inner_Prod(Real_t* u, Real_t* v, size_t n, Real_t* pI, Real_t* pJ, Real_t* pV);
	// this requires a cast of int32(MATLAB_PTR) operation in Matlab's environment.
	void Inner_Prod(Real_t* u, Real_t* v, size_t n, int32_t* pI, int32_t* pJ, Real_t* pV);




};

} /* namespace MEX */

#endif /* SRC_PDECS_PRIVATE_SERVER_H_ */
