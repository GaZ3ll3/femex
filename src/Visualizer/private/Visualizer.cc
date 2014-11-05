/*
 * Visualizer.cpp
 *
 *  Created on: Nov 4, 2014
 *      Author: lurker
 */

#include "Visualizer.h"

namespace MEX {

Visualizer::Visualizer() {

}

Visualizer::~Visualizer() {
}

} /* namespace MEX */


using namespace MEX;

template class mexplus::Session<Visualizer>;


namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 0);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Visualizer>::create(new Visualizer()));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<Visualizer>::destroy(input.get(0));
}

}

MEX_DISPATCH
