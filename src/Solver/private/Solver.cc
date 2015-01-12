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

void Solver::Gradient(MatlabPtr& GradX, MatlabPtr& GradY, MatlabPtr Solution, MatlabPtr Nodes,
		MatlabPtr Elems, MatlabPtr DX, MatlabPtr DY){

	auto numberofelem        = mxGetN(Elems);
	auto numberofnodeperelem = mxGetM(Elems);

	auto gradX_ptr = mxGetPr(GradX);
	auto gradY_ptr = mxGetPr(GradY);

	auto elem_ptr = (int32_t*)mxGetPr(Elems);
	auto node_ptr = mxGetPr(Nodes);

	auto sln_ptr  = mxGetPr(Solution);

	auto DX_ptr   = mxGetPr(DX);
	auto DY_ptr   = mxGetPr(DY);


	Real_t temp_ptr[numberofnodeperelem];
	Real_t gx_ptr[numberofnodeperelem];
	Real_t gy_ptr[numberofnodeperelem];

	size_t vertex_1, vertex_2, vertex_3;
	Real_t det;
	Real_t Jacobian[2][2];

	for (size_t i = 0; i < numberofelem; i++) {

		for (size_t j = 0; j < numberofnodeperelem; j++) {
			temp_ptr[j] = sln_ptr[elem_ptr[i * numberofnodeperelem + j] - 1];
		}


		for (size_t j = 0; j < numberofnodeperelem; j++) {
			gx_ptr[j] = 0.;
			gy_ptr[j] = 0.;
			/*
			 * since k is relatively small, thus the summation won't be incorrect too much.
			 */
			for (size_t k = 0; k < numberofnodeperelem; k++) {
				gx_ptr[j] += temp_ptr[k] * DX_ptr[j * numberofnodeperelem + k];
				gy_ptr[j] += temp_ptr[k] * DY_ptr[j * numberofnodeperelem + k];
			}
		}
		vertex_1 = elem_ptr[numberofnodeperelem*i] - 1;
		vertex_2 = elem_ptr[numberofnodeperelem*i + 1] - 1;
		vertex_3 = elem_ptr[numberofnodeperelem*i + 2] - 1;

		Jacobian[0][0] = node_ptr[2*vertex_3 + 1] - node_ptr[2*vertex_1 + 1];
		Jacobian[1][1] = node_ptr[2*vertex_2    ] - node_ptr[2*vertex_1    ];
		Jacobian[0][1] = node_ptr[2*vertex_1 + 1] - node_ptr[2*vertex_2 + 1];
		Jacobian[1][0] = node_ptr[2*vertex_1    ] - node_ptr[2*vertex_3    ];

		// Orientation corrected.
		det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

		for (size_t j = 0; j < numberofnodeperelem; j++) {
			gradX_ptr[i * numberofnodeperelem + j] =
					(Jacobian[0][0] * gx_ptr[j] + Jacobian[0][1] * gy_ptr[j])/det;
			gradY_ptr[i * numberofnodeperelem + j] =
					(Jacobian[1][0] * gx_ptr[j] + Jacobian[1][1] * gy_ptr[j])/det;
		}
	}
}

void Solver::Neumann(MatlabPtr& Neumann, MatlabPtr GradX, MatlabPtr GradY, MatlabPtr Nodes, MatlabPtr Elems,
		MatlabPtr Boundaries, MatlabPtr BoundaryId){

	auto numberofelem        = mxGetN(Elems);
	auto numberofnodeperelem = mxGetM(Elems);

	auto gradX_ptr = mxGetPr(GradX);
	auto gradY_ptr = mxGetPr(GradY);

	auto elem_ptr = (int32_t*)mxGetPr(Elems);
	auto node_ptr = mxGetPr(Nodes);

	auto neumann_ptr = mxGetPr(Neumann);

	auto boundary_ptr = (int32_t*)mxGetPr(Boundaries);
	auto id_ptr       = (int32_t*)mxGetPr(BoundaryId);


	auto numberofboundary = mxGetN(Boundaries);
	auto numberofnodeperboundary = mxGetM(Boundaries);



	/*
	 * Each boundary-segment contains 'degree + 1' end points,
	 * which produces 'degree + 1' Neumann data for each boundary.
	 *
	 * Neumann is a matrix with (number of boundary) x (degree + 1).
	 *
	 * each segment has same exterior normal vector, so it won't cost
	 * too much on calculating the dot product.
	 *
	 * complexity O(number of boundary) x O(degree + 1).
	 */

	// normal vector, uninitialized.
	Real_t normal[2], length;
	mwSize first, second, b_id;
	mwSize tri_u, tri_v, tri_w;
	mwSize edge[numberofnodeperboundary];

	for (size_t i = 0; i < numberofboundary; i++) {

		/*
		 * pick first two nodes.
		 */
		first  = boundary_ptr[i * numberofnodeperboundary] - 1;
		second = boundary_ptr[i * numberofnodeperboundary + 1] - 1;

		b_id  = id_ptr[i] - 1;

		tri_u = elem_ptr[b_id * numberofnodeperelem] - 1;
		tri_v = elem_ptr[b_id * numberofnodeperelem + 1] - 1;
		tri_w = elem_ptr[b_id * numberofnodeperelem + 2] - 1;

		/*
		 * match the edge
		 */
		// u - v
		if (tri_u == first && tri_v == second) {
			edge[0] = 1; edge[1] = 2;
			for (size_t k = 2; k < numberofnodeperboundary; k++) {
				edge[k] = 2 + k;
			}
		}
		// v - w
		else if (tri_v == first && tri_w == second) {
			edge[0] = 2; edge[1] = 3;
			for (size_t k = 2; k < numberofnodeperboundary; k++) {
				edge[k] = numberofnodeperboundary + k;
			}
		}
		// w - u
		else if (tri_w == first && tri_u == second){
			edge[0] = 3; edge[1] = 1;
			for (size_t k = 2; k < numberofnodeperboundary; k++) {
				edge[k] = 2 * numberofnodeperboundary - 2 + k ;
			}
		}
		else {
			// Error!
			mexErrMsgTxt("Solver::Neumann::The boundary does not match the element.\n");
		}


		/*
		 * normal direction
		 */
		normal[0] = node_ptr[2 * second + 1] - node_ptr[2 * first + 1];
		normal[1] = node_ptr[2 * first] - node_ptr[2 * second];

		length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1]);

		/*
		 * normalize
		 */
		normal[0] /= length;
		normal[1] /= length;

		for (size_t j = 0; j < numberofnodeperboundary; j++){
			neumann_ptr[i * numberofnodeperboundary + j] =
					gradX_ptr[b_id * numberofnodeperelem + edge[j] - 1] * normal[0] +
					gradY_ptr[b_id * numberofnodeperelem + edge[j] - 1] * normal[1];
		}
	}
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


MEX_DEFINE(grad) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs,  prhs, 6);
	OutputArguments output(nlhs, plhs, 2);

	Solver* solver = Session<Solver>::get(input.get(0));

	size_t numberofelem        = mxGetN(prhs[3]);
	size_t numberofnodeperelem = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix(numberofnodeperelem, numberofelem, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(numberofnodeperelem, numberofelem, mxDOUBLE_CLASS, mxREAL);

	solver->Gradient(plhs[0], plhs[1], CAST(prhs[1]),
			CAST(prhs[2]),CAST(prhs[3]),
			CAST(prhs[4]),CAST(prhs[5]));
}


MEX_DEFINE(neumann) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs,  prhs, 7);
	OutputArguments output(nlhs, plhs, 1);

	Solver* solver = Session<Solver>::get(input.get(0));

	size_t numberofboundary        = mxGetN(prhs[5]);
	size_t numberofnodeperboundary = mxGetM(prhs[5]);

	plhs[0] = mxCreateNumericMatrix(numberofnodeperboundary, numberofboundary, mxDOUBLE_CLASS, mxREAL);

	solver->Neumann(plhs[0], CAST(prhs[1]),
			CAST(prhs[2]),CAST(prhs[3]),
			CAST(prhs[4]),CAST(prhs[5]),
			CAST(prhs[6]));
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
