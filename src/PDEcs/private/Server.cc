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

// not optimized, serial environment
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

		/*
		 *
		 * OpenMP only helps when there is
		 * not too much overheads.
		 *
		 */
		auto m = n / 8; auto r = n % 8;

		size_t j = 0;
		Real_t worker_0 = 0., worker_1= 0., worker_2 = 0., worker_3 = 0.,\
				worker_4 = 0. , worker_5 = 0., worker_6 = 0.,worker_7 = 0., w_sum = 0.;


			for (auto i = 0; i < m; i++) {
				worker_0 += u[j] * v[j];
				worker_1 += u[j + 1] * v[j + 1];
				worker_2 += u[j + 2] * v[j + 2];
				worker_3 += u[j + 3] * v[j + 3];
				worker_4 += u[j + 4] * v[j + 4];
				worker_5 += u[j + 5] * v[j + 5];
				worker_6 += u[j + 6] * v[j + 6];
				worker_7 += u[j + 7] * v[j + 7];

				j += 8;
			}

			w_sum = worker_0 + worker_1 + worker_2 + worker_3 + worker_4 + worker_5 + worker_6 + worker_7;

		if (r > 0) {
			for (auto k = 0; k < r; k++) {
				w_sum += u[j +k] * v [j + k];
			}
		}
		return w_sum;
	}

}

/*
 *  If there is no unroll-loops, it is still not fast enough
 *
 *  For smaller inner product, there is no much difference.
 */
Real_t Server::Inner_Prod_omp(Real_t* u, Real_t* v, size_t n){


	size_t chunk, thread_num, nthreads;
	size_t _start, _end, i;
	Real_t result;

#pragma omp parallel default(none) \
	shared(u,v,nthreads, n) private(chunk,i,thread_num,_start,_end)\
	reduction(+:result)
    {

		thread_num = omp_get_thread_num();
		nthreads = omp_get_num_threads();

		chunk = n / nthreads;
		_start =  thread_num * chunk;
		_end = (thread_num + 1) * chunk;

		if (thread_num == nthreads - 1) {
			_end = n;
		}
		for(i = _start; i < _end; i++) {
			result += u[i] * v[i];
		}
    }
    return result;
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

MEX_DEFINE(sprod) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 1);

	auto server = Session<Server>::get(input.get(0));
	auto u_ptr = Matlab_Cast<Real_t>(CAST(prhs[1]));
	auto v_ptr = Matlab_Cast<Real_t>(CAST(prhs[2]));
	auto n     = mxGetM(prhs[1]) > mxGetM(prhs[2])? mxGetM(prhs[2]):mxGetM(prhs[1]);
	auto ret = server->Inner_Prod(u_ptr, v_ptr, n);

	plhs[0] = mxCreateDoubleScalar(ret);

}

MEX_DEFINE(pprod) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 1);

	auto server = Session<Server>::get(input.get(0));
	auto u_ptr = Matlab_Cast<Real_t>(CAST(prhs[1]));
	auto v_ptr = Matlab_Cast<Real_t>(CAST(prhs[2]));
	auto n     = mxGetM(prhs[1]) > mxGetM(prhs[2])? mxGetM(prhs[2]):mxGetM(prhs[1]);
	auto ret = server->Inner_Prod_omp(u_ptr, v_ptr, n);

	plhs[0] = mxCreateDoubleScalar(ret);

}

} /* namespace */


MEX_DISPATCH
