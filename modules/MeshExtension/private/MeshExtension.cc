/*
 * MeshExtension.cpp
 *
 *  Created on: Feb 7, 2015
 *      Author: lurker
 */

#include "MeshExtension.h"

namespace Extension {

MeshExtension::MeshExtension() {

}

MeshExtension::~MeshExtension() {
#ifdef DEBUG
	mexPrintf("MeshExtension detached...\n");
#endif
}

} /* namespace Extension */

using namespace Extension;

template class mexplus::Session<MeshExtension>;

namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 0);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<MeshExtension>::create(new MeshExtension()));
}

MEX_DEFINE(delete)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<MeshExtension>::destroy(input.get(0));
}
}

MEX_DISPATCH
