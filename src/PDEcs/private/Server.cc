/*
 * Server.cpp
 *
 *  Created on: Nov 22, 2014
 *      Author: lurker
 */

#include "Server.h"

namespace MEX {

Server::Server() {
#ifdef DEBUG
	mexPrintf("PDEcs-Server attached\n");
#endif

}

Server::~Server() {
#ifdef DEBUG
	mexPrintf("PDEcs-Server detached\n");
#endif
}

// not optimized
Real_t Server::Inner_Prod(Real_t* u, Real_t* v, size_t n){



	if (n < 12) {
		Real_t sum  = u[0]*v[0];
		for (size_t i = 1; i < n; i++) {
			sum += u[i]  * v[i];
		}
		return sum;
	}

	else{
		/*
		 * manually unroll the loop,
		 *
		 * should have 2x performance
		 */

		auto m = n / 8; auto r = n % 8;

		register size_t j;
		Real_t worker_0 , worker_1, worker_2, worker_3, worker_4, worker_5, worker_6 ,worker_7;

		for (auto i = 0; i < m; i++) {
			j = 8 * i;
			worker_0 += u[j] * v[j];
			worker_1 += u[j + 1] * v[j + 1];
			worker_2 += u[j + 2] * v[j + 2];
			worker_3 += u[j + 3] * v[j + 3];
			worker_4 += u[j + 4] * v[j + 4];
			worker_5 += u[j + 5] * v[j + 5];
			worker_6 += u[j + 6] * v[j + 6];
			worker_7 += u[j + 7] * v[j + 7];
		}

		auto sum = worker_0 + worker_1 + worker_2 + worker_3 + worker_4 + worker_5 + worker_6 + worker_7;

		if (r > 0) {
			for (auto j = 0; j < r; j++) {
				sum += u[4 * m + j] * v [4 * m + j];
			}
		}
		return sum;
	}

}

} /* namespace MEX */


using namespace MEX;

template class mexplus::Session<Server>;


namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 0);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Server>::create(new Server()));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<Server>::destroy(input.get(0));
}

MEX_DEFINE(prod) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 1);

	auto server = Session<Server>::get(input.get(0));


	auto u_ptr = Matlab_Cast<Real_t>(CAST(prhs[1]));

	auto v_ptr = Matlab_Cast<Real_t>(CAST(prhs[2]));

	auto n     = mxGetM(prhs[1]) > mxGetM(prhs[2])? mxGetM(prhs[2]):mxGetM(prhs[1]);


	auto ret = server->Inner_Prod(u_ptr, v_ptr, n);
	plhs[0] = mxCreateDoubleScalar(ret);

}

} /* namespace */


MEX_DISPATCH
