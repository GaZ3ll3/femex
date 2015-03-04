/*
 * Adjoint.cpp
 *
 *  Created on: Feb 28, 2015
 *      Author: lurker
 */

#include "Adjoint.h"

namespace Core {

Adjoint::Adjoint() noexcept{

	parameter.resize(0);
	solution.resize(0);
	adjoint_solution.resize(0);
	gradient.resize(0);

	J = 0.;
}

Adjoint::Adjoint(std::size_t param_size, std::size_t solution_size) noexcept {

	parameter.resize(param_size);
	solution.resize(solution_size);
	adjoint_solution.resize(solution_size);
	gradient.resize(param_size);

	J = 0.;
}

Adjoint::~Adjoint() {

	parameter.clear();
	solution.clear();
	adjoint_solution.clear();
	gradient.clear();

#ifdef DEBUG
	mexPrintf("Adjoint Scheme is detached.\n");
#endif
}

void Adjoint::Update_Gradient(){

}

void Adjoint::Update_Solution(MatlabPtr _solution){

}
void Adjoint::Update_Param(MatlabPtr _parameter){

}


void Adjoint::Update_J(){

}

void Adjoint::Update_Adj(){

}

} /* namespace Core */


using namespace Core;
template class mexplus::Session<Adjoint>;

namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Adjoint>::create(new Adjoint()));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<Adjoint>::destroy(input.get(0));
}
}
