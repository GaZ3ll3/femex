/*
 * Solver.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: lurker
 */

#include "Solver.h"

namespace MEX {

Solver::Solver() {

}

Solver::~Solver() {
#ifdef DEBUG
	std::cout << "Solver detached.\n" << std::endl;
#endif
}

void Solver::Reference(MatlabPtr& DX, MatlabPtr& DY, MatlabPtr Points){

	auto _numberofpoints = mxGetN(Points);

	auto Vander = mxCreateNumericMatrix(_numberofpoints,_numberofpoints,mxDOUBLE_CLASS, mxREAL);
	auto VanderX = mxCreateNumericMatrix(_numberofpoints, _numberofpoints, mxDOUBLE_CLASS, mxREAL);
	auto VanderY = mxCreateNumericMatrix(_numberofpoints, _numberofpoints, mxDOUBLE_CLASS, mxREAL);


	auto Vander_ptr  = mxGetPr(Vander);
	auto VanderX_ptr = mxGetPr(VanderX);
	auto VanderY_ptr = mxGetPr(VanderY);

	auto nodes_ptr   = mxGetPr(Points);

	int deg = round((sqrt(8*_numberofpoints + 1) - 3)/2);

	if ((deg + 1) * (deg + 2)/2 != _numberofpoints){
		mexErrMsgTxt("Invalid length of input nodes\n");
	}

	for (size_t col = 0; col < _numberofpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++){
				*Vander_ptr++ = pow(nodes_ptr[2*col], i - j)*pow(nodes_ptr[2*col + 1], j);
			}
		}
	}

	for (size_t col = 0; col < _numberofpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++) {
				if (j == i){
					*VanderX_ptr++ = 0.;
				}
				else{
					*VanderX_ptr++ = (i - j) * pow(nodes_ptr[2*col], i - j - 1)*pow(nodes_ptr[2*col + 1], j);
				}
			}
		}
	}



	for (size_t col = 0; col < _numberofpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++) {
				if (j == 0){
					*VanderY_ptr++ = 0.;
				}
				else{
					*VanderY_ptr++ = j* pow(nodes_ptr[2*col], i - j)*pow(nodes_ptr[2*col + 1], j - 1);
				}
			}
		}
	}


	mxArray* RHS_x[] = {Vander, VanderX};
	mexCallMATLAB(1, &DX, 2, RHS_x, "mldivide");

	mxArray* RHS_y[] = {Vander, VanderY};
	mexCallMATLAB(1, &DY, 2, RHS_y, "mldivide");

	mxDestroyArray(Vander);
	mxDestroyArray(VanderX);
	mxDestroyArray(VanderY);

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

MEX_DEFINE(reference) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs,  prhs, 2);
	OutputArguments output(nlhs, plhs, 2);

	Solver* solver = Session<Solver>::get(input.get(0));

	size_t numberofpoints = mxGetN(prhs[1]);

	plhs[0] = mxCreateNumericMatrix(numberofpoints, numberofpoints, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(numberofpoints, numberofpoints, mxDOUBLE_CLASS, mxREAL);

	solver->Reference(plhs[0], plhs[1], CAST(prhs[1]));
}

}

MEX_DISPATCH
