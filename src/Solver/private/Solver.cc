/*
 * Solver.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: lurker
 */

#include "Solver.h"

namespace MEX {

Solver::Solver() {
	// default solver settings
	option._type = Solver_type::direct;
	option._restol = 1e-12;
	option._multilevel = false;
	option._nrestart = 30;
	option._maxIter = 500;
	option._smoother = Smoother::gsb;
	option._amg = AMG::ilu;
}


void Solver::set_amg(int32_t _amg){
	option._amg  = static_cast<AMG>(_amg);
}

void Solver::set_maxIter(int32_t _maxIter){
	option._maxIter = _maxIter;
}

void Solver::set_nrestart(int32_t _nrestart) {
	option._nrestart = _nrestart;
}

void Solver::set_restol(Real_t _restol) {
	option._restol = _restol;
}

void Solver::set_smoother(int32_t _smoother){
	option._smoother = static_cast<Smoother>(_smoother);
}

void Solver::set_type(int32_t _type) {
	option._type = static_cast<Solver_type>(_type);
}


int32_t Solver::get_amg() {
	return static_cast<int32_t>(option._amg);
}
int32_t Solver::get_maxIter() {
	return option._maxIter;
}

int32_t Solver::get_nrestart(){
	return option._nrestart;
}

Real_t Solver::get_restol() {
	return option._restol;
}

int32_t Solver::get_smoother(){
	return static_cast<int32_t>(option._smoother);
}

int32_t Solver::get_type() {
	return static_cast<int32_t>(option._type);
}


Solver::~Solver() {
#ifdef DEBUG
	std::cout << "Solver detached.\n" << std::endl;
#endif
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
