/*
 * Treecode_.cpp
 *
 *  Created on: Oct 5, 2015
 *      Author: lurker
 */

#include "treecode.h"


#include <mexplus.h>
#include "utils.h"

using namespace std;
using namespace mexplus;


template class mexplus::Session<treecode>;


namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 4);
	OutputArguments output(nlhs, plhs, 1);

	auto px = mxGetPr(prhs[0]);
	auto py = mxGetPr(prhs[1]);
	auto pl = mxGetPr(prhs[2]);
	auto pv = mxGetPr(prhs[3]);

	output.set(0, Session<treecode>::create(new treecode(
			*px, *py, * pl , *pv
			)));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	auto root = Session<treecode>::get(input.get(0));

	Session<treecode>::destroy(input.get(0));
}

MEX_DEFINE(setAttribute) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])  {
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 0);

	auto attributes = mxGetPr(prhs[1]);

	auto root = Session<treecode>::get(input.get(0));
	for (auto it : root->root->points) {
		it->attribute = *attributes;
		attributes++;
	}

}

MEX_DEFINE(split) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);

	auto root = Session<treecode>::get(input.get(0));
	root->root->populate();
}

MEX_DEFINE(buildmatrix) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 1);

	auto root = Session<treecode>::get(input.get(0));
	auto num = root->size * root->size;
	auto theta_ptr = mxGetPr(prhs[1]);
	if (num == 0) {
		std::cout << "empty cell, stopped" << std::endl;
		return;
	}
	plhs[0] = mxCreateNumericMatrix(num, num,mxDOUBLE_CLASS, mxREAL);
	auto Radptr  = mxGetPr(plhs[0]);

//	omp_set_num_threads(omp_get_num_procs());
	omp_set_num_threads(2048);

	int i;

#pragma omp parallel for shared(theta_ptr, root, Radptr, num) private(i)
	for (i = 0; i < num; i++) {
		traversal(root, *theta_ptr, root->root->points[i], root->root, num, Radptr);
	}
}

MEX_DEFINE(apply) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 1);

	auto root = Session<treecode>::get(input.get(0));
	auto num = root->size * root->size;

	auto theta_ptr = mxGetPr(prhs[1]);
	auto rhs_ptr = mxGetPr(prhs[2]);


	if (num == 0) {
		std::cout << "empty cell, stopped" << std::endl;
		return;
	}

	/*
	 * already zero vector
	 */
	plhs[0] = mxCreateNumericMatrix(num, 1,mxDOUBLE_CLASS, mxREAL);
	auto Radptr  = mxGetPr(plhs[0]);

//	omp_set_num_threads(omp_get_num_procs());
	omp_set_num_threads(2048);

	int i;

#pragma omp parallel for shared(theta_ptr, root, Radptr, num, rhs_ptr) private(i)
	for (i = 0; i < num; i++) {
		traversal(root, *theta_ptr, root->root->points[i], root->root, num, Radptr, rhs_ptr);
	}
}

MEX_DEFINE(fast_apply) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 1);

	auto root = Session<treecode>::get(input.get(0));
	auto num = root->size * root->size;

	auto theta_ptr = mxGetPr(prhs[1]);
	auto rhs_ptr = mxGetPr(prhs[2]);


	if (num == 0) {
		std::cout << "empty cell, stopped" << std::endl;
		return;
	}

	/*
	 * already zero vector
	 */
	plhs[0] = mxCreateNumericMatrix(num, 1,mxDOUBLE_CLASS, mxREAL);
	auto Radptr  = mxGetPr(plhs[0]);

//	omp_set_num_threads(omp_get_num_procs());
	omp_set_num_threads(2048);

	int i;

#pragma omp parallel for shared(theta_ptr, root, Radptr, num, rhs_ptr) private(i)
	for (i = 0; i < num; i++) {
		fast_traversal(root, *theta_ptr, root->root->points[i], root->root, num, Radptr, rhs_ptr);
	}
}

MEX_DEFINE(preprocess) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 0);

	auto root = Session<treecode>::get(input.get(0));
	auto rhs_ptr = mxGetPr(prhs[1]);

	traversal_down(root->root, rhs_ptr);
	traversal_up(root->root);


}

MEX_DEFINE(buildop)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
// not implemented, since the source class only permits int32_t as size.
}


}

MEX_DISPATCH
