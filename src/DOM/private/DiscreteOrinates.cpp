/*
 * DiscreteOrinates.cpp
 *
 *  Created on: Feb 21, 2015
 *      Author: lurker
 */

#include "DiscreteOrinates.h"

namespace Core {

DiscreteOrinates::DiscreteOrinates(MatlabPtr _nAngle) noexcept {
	nAngle = _nAngle;
	initAngle = 0.;

}

DiscreteOrinates::DiscreteOrinates(MatlabPtr _nAngle, MatlabPtr _initAngle) noexcept {
	nAngle = _nAngle;
	initAngle = _initAngle;
}

DiscreteOrinates::~DiscreteOrinates() {
#ifdef DEBUG
	mexPrintf("Discrete Orinates Method detached.\n");
#endif
}

void DiscreteOrinates::RayInt(Real_t*& output, MatlabPtr nodes, MatlabPtr elems,
		MatlabPtr neighbors, MatlabPtr weights, MatlabPtr Fcn){

	auto numberofnodes        = mxGetN(nodes);
	auto numberofelems        = mxGetN(elems);
	auto numberofnodesperelem = mxGetM(elems);


	auto pnodes               = mxGetPr(nodes);
	auto pelems               = (int32_t*)mxGetPr(elems);
	auto pneighbors           = (int32_t*)mxGetPr(neighbors);
	auto pweights             = mxGetPr(weights);
	auto interp               = mxGetPr(Fcn);


	mwSize vertex;
	Real_t theta;

	for (int32_t i = 0; i < nAngle; i++) {
		theta = initAngle + 2 * M_PIl/ nAngle;

		// reallocate the empty set.

		for (int32_t j = 0; j < numberofelems; j++) {
//			vertex_1 = pelems[numberofnodesperelem*i] - 1;
//			vertex_2 = pelems[numberofnodesperelem*i + 1] - 1;
//			vertex_3 = pelems[numberofnodesperelem*i + 2] - 1;

			/*
			 * calculate ray integral for each vertex.
			 *
			 * skip if already visited.
			 */
			for (int32_t k = 0; k < 3; k++) {
				vertex = pelems[numberofnodesperelem * i + k] - 1;

				// insert the vertex into the set.




			}
		}
	}
}

} /* namespace Core */
