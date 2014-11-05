/*
 * Boundary.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: lurker
 */

#include "Boundary.h"

namespace MEX {

Boundary::Boundary(MatlabPtr _edges) {

	b_edges.resize(mxGetM(_edges)*mxGetN(_edges));
	memcpy(&b_edges[0], mxGetPr(_edges), mxGetM(_edges)*mxGetN(_edges)*sizeof(int32_t));
	for (auto it = b_edges.begin(); it != b_edges.end(); it++) {
		b_edge_set.insert(*it);
	}
}

Boundary::~Boundary() {

	b_edges.clear();
	b_edge_set.clear();
#ifdef DEBUG
	mexPrintf("Boundary detached\n");
#endif
}

} /* namespace MEX */



using namespace MEX;

template class mexplus::Session<DirichletBC>;


namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<DirichletBC>::create(new DirichletBC(
			const_cast<MatlabPtr>(input.get(0))
						)));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<DirichletBC>::destroy(input.get(0));
}

MEX_DEFINE(report) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	auto boundary = Session<DirichletBC>::get(input.get(0));
	mexPrintf("Dirichlet Boundary Nodes Number: %d \n", boundary->b_edge_set.size());

}

// can be implemented by set_difference from <algorithm>
MEX_DEFINE(dofs) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 1);
	auto boundary = Session<DirichletBC>::get(input.get(0));
	auto _range   = input.get<int32_t>(1);

	plhs[0] = mxCreateNumericMatrix(_range - boundary->b_edge_set.size(), 1, mxINT32_CLASS, mxREAL);

	auto plhs_ptr = Matlab_Cast<int32_t>(plhs[0]);
	for (int32_t i = 1; i <= _range; ++i) {
		if (boundary->b_edge_set.find(i) != boundary->b_edge_set.end()) {
			continue;
		}
		else {
			*plhs_ptr++ = i;
		}
	}



}
} // namespace

MEX_DISPATCH
