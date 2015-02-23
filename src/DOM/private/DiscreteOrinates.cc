/*
 * DiscreteOrinates.cpp
 *
 *  Created on: Feb 21, 2015
 *      Author: lurker
 */

#include "DiscreteOrinates.h"

namespace Core {

DiscreteOrinates::DiscreteOrinates(MatlabPtr _nAngle) noexcept {
	nAngle = static_cast<int32_t>(*mxGetPr(_nAngle));
	initAngle = 0.;

}

DiscreteOrinates::DiscreteOrinates(MatlabPtr _nAngle, MatlabPtr _initAngle) noexcept {
	nAngle = static_cast<int32_t>(*mxGetPr(_nAngle));
	initAngle = static_cast<int32_t>(*mxGetPr(_initAngle));
}

DiscreteOrinates::~DiscreteOrinates() {
#ifdef DEBUG
	mexPrintf("Discrete Orinates Method detached.\n");
#endif
}

void DiscreteOrinates::RayInt(Real_t*& output, MatlabPtr nodes, MatlabPtr elems,
		MatlabPtr neighbors, MatlabPtr edges, MatlabPtr weights, MatlabPtr Fcn){


	auto numberofnodes        = mxGetN(nodes);
	auto numberofelems        = mxGetN(elems);
	auto numberofnodesperelem = mxGetM(elems);
	auto numberofedges        = mxGetN(edges);
	auto numberofnodesperedge = mxGetM(edges);


	auto pnodes               = mxGetPr(nodes);
	auto pelems               = (int32_t*)mxGetPr(elems);
	auto pneighbors           = (int32_t*)mxGetPr(neighbors);
	auto pweights             = mxGetPr(weights);
	auto interp               = mxGetPr(Fcn);
	auto pedges               = (int32_t*)mxGetPr(edges);

	mwSize vertex, e_vl, e_vr;;
	Real_t theta , x1, x2, y1, y2,  a, b;
	std::unordered_set<int32_t> visited;


	Real_t ret, t, eta;
	int32_t index;

#ifdef ADJACENT
	/*
	 * build adjacent mapping.
	 */
	std::vector<std::unordered_set<int32_t>> adjacent(numberofnodes);
	mwSize vertex_1, vertex_2, vertex_3;

	for (int32_t i = 0; i < numberofelems; i++) {
		vertex_1 = pelems[numberofnodesperelem * i    ] - 1;
		vertex_2 = pelems[numberofnodesperelem * i + 1] - 1;
		vertex_3 = pelems[numberofnodesperelem * i + 2] - 1;

		adjacent[vertex_1].insert(vertex_2);
		adjacent[vertex_1].insert(vertex_3);
		adjacent[vertex_2].insert(vertex_1);
		adjacent[vertex_2].insert(vertex_3);
		adjacent[vertex_3].insert(vertex_1);
		adjacent[vertex_3].insert(vertex_2);
	}

#ifdef DEBUG
	for (int32_t i = 0; i < numberofnodes; i++) {
		for (auto it = adjacent[i].begin(); it != adjacent[i].end(); it++) {
			std::cout << *it << " ";
		}
		std::cout << std::endl;
	}
#endif
	/*
	 * adjacent mapping of nodes built.
	 */
#endif

	for (int32_t i = 0; i < nAngle; i++) {
		std::cout << "the " << i << "th run" << std::endl;
		theta = initAngle + 2 * i * M_PIl / nAngle;
		visited.clear();
		for (int32_t j = 0; j < numberofelems; j++) {
			for (int32_t k = 0; k < numberofnodesperelem; k++){
				/*
				 * index of vertex in nodes.
				 */
				vertex = pelems[numberofnodesperelem * j + k] - 1;
				if (visited.find(vertex) == visited.end()){
					visited.insert(vertex);
					/*
					 * calculate the ray intersects with the boundary.
					 */
					a = pnodes[2 * vertex    ];
					b = pnodes[2 * vertex + 1];
					for (int32_t l = 0; l < numberofedges; l++) {
						e_vl = pedges[l * numberofnodesperedge    ] - 1;
						e_vr = pedges[l * numberofnodesperedge + 1] - 1;
						x1 = pnodes[2 * e_vl    ];
						x2 = pnodes[2 * e_vr    ];
						y1 = pnodes[2 * e_vl + 1];
						y2 = pnodes[2 * e_vr + 1];
						if (fabs( (y1 - y2)*(a - x2) - (x1 - x2) * (b - y2)) < MEX_EPS) {
							// colinear.
						}
						else {
							// non-colinear.
							if (fabs((x1 - x2) * sin(theta) - (y1- y2) * cos(theta)) < MEX_EPS) {
								// ray along edge.
								continue;
							}
							else {
								t =((y1 - y2)*(a - x2) - (x1 - x2) * (b - y2))/
										((x1 - x2) * sin(theta) - (y1- y2) * cos(theta));

								eta = (sin(theta)*(a - x2) - cos(theta) * (b - y2))/
										(sin(theta) * (x1 - x2) - cos(theta) * (y1 - y2));
								if (t >= 0 && eta >= 0 && eta <= 1) {
#ifdef DEBUG
									std::cout << "vertex " << vertex << " at [" << a << ", " << b  << "]"<< " at distance of " << t << " from boundary"
											<< "["  << x1 * eta + (1 - eta) * x2 << ", " << y1 * eta + (1 - eta)*y2
											<< "]" << std::endl;
#endif




								}// end if
							} // end else
						}// end else
					} //end for
				} // end if
			} // end for
		} // end for
	} // end for
} // end

} /* namespace Core */


using namespace Core;
template class mexplus::Session<DiscreteOrinates>;

namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<DiscreteOrinates>::create(new DiscreteOrinates(const_cast<mxArray*>(input.get(0)))));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<DiscreteOrinates>::destroy(input.get(0));
}

MEX_DEFINE(rayint) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 1);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	plhs[0] = mxCreateNumericMatrix(1, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);
	DOM->RayInt(pI, CAST(prhs[1]), CAST(prhs[2]),
			CAST(prhs[3]), CAST(prhs[4]),
			CAST(prhs[5]), CAST(prhs[6]));
}

}

MEX_DISPATCH

