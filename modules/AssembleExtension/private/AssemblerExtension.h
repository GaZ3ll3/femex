/*
 * AssemblerExtension.h
 *
 *  Created on: Feb 5, 2015
 *      Author: lurker
 */

#ifndef MODULES_ASSEMBLEEXTENSION_ASSEMBLEREXTENSION_H_
#define MODULES_ASSEMBLEEXTENSION_ASSEMBLEREXTENSION_H_

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

namespace Extension {

class AssemblerExtension {

public:
	AssemblerExtension();
	virtual ~AssemblerExtension();

	/*
	 * integral on int \phi_i_x \phi_j
	 */
	void AssembleGradXFunc(Real_t* &pI, Real_t* &pJ, Real_t* &pV,
			MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr Ref, MatlabPtr RefX,MatlabPtr RefY,
			MatlabPtr Weights, MatlabPtr Fcn);
	/*
	 * integral on int \phi_i_y \phi_j
	 */
	void AssembleGradYFunc(Real_t* &pI, Real_t* &pJ, Real_t* &pV,
			MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr Ref, MatlabPtr RefX, MatlabPtr RefY,
			MatlabPtr Weights, MatlabPtr Fcn);

};

}

#endif /* MODULES_ASSEMBLEEXTENSION_ASSEMBLEREXTENSION_H_ */
