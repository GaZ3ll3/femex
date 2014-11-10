/*
 * Solver.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: lurker
 */

#include "Solver.h"

namespace MEX {

Solver::Solver() {
	// TODO Auto-generated constructor stub

}

Solver::~Solver() {
	// TODO Auto-generated destructor stub
}

} /* namespace MEX */

using namespace MEX;

template class mexplus::Session<Solver>;


namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 0);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Solver>::create(new Solver()));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<Solver>::destroy(input.get(0));
}

}

MEX_DISPATCH
