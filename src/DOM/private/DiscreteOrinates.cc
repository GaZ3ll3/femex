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
	Real_t q_x1, q_y1, q_x2, q_y2, q_x3, q_y3;
	Real_t q_eta, q_t;
	int32_t index;

	std::vector<Real_t> tmp;

	bool flag, intersect;
	std::queue<int32_t> q_elem;
	std::unordered_set<int32_t> s_elem;

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

					// clean up the queue for elements.(??)
					// it is not needed to clean up, since the queue
					// will be empty before next round.

					for (int32_t l = 0; l < numberofedges; l++) {
						e_vl = pedges[l * numberofnodesperedge    ] - 1;
						e_vr = pedges[l * numberofnodesperedge + 1] - 1;
						x1 = pnodes[2 * e_vl    ];
						x2 = pnodes[2 * e_vr    ];
						y1 = pnodes[2 * e_vl + 1];
						y2 = pnodes[2 * e_vr + 1];
						//**********************************************
						//** rewrite this part.
						//**********************************************
						if (fabs( (y1 - y2)*(a - x2) - (x1 - x2) * (b - y2)) < SCALE * MEX_EPS) {
							// colinear.
							// colinear must sit in between, otherwise, the polygon is not convex.
						}
						else {
							// non-colinear.
							if (fabs((x1 - x2) * sin(theta) - (y1- y2) * cos(theta)) < SCALE * MEX_EPS) {
								// ray along edge.
								continue;
							}
							else {
								t =((y1 - y2)*(a - x2) - (x1 - x2) * (b - y2))/
										((x1 - x2) * sin(theta) - (y1- y2) * cos(theta));

								eta = (sin(theta)*(a - x2) - cos(theta) * (b - y2))/
										(sin(theta) * (x1 - x2) - cos(theta) * (y1 - y2));
								if (t >= 0 && eta >= -MEX_EPS && eta <= 1 + MEX_EPS) {
									flag = true;
									break;
								}// end if
								else {
									flag = false;
								}
							} // end else
						}// end else
						//************************************************
					} //end for

					if (flag){
#ifdef DEBUG
						std::cout
							<< "vertex " << vertex << " at [" << a << ", " << b  << "]"
							<< " at distance of " << t << " from boundary"
							<< "["  << x1 * eta + (1 - eta) * x2 << ", " << y1 * eta + (1 - eta)*y2
							<< "]" << std::endl;
#endif
						/*
						 * search along the ray, find the elements which passed through.
						 *
						 * use queue to do it.
						 *
						 * stopping criteria:
						 *
						 * queue is empty.
						 *
						 * Adding into queue if barycentric point is pushing forward.
						 *
						 */
						// clean up set.
						s_elem.clear();

						// push current element
						q_elem.push(j);
						s_elem.insert(j);

						/*
						 * check the jth.
						 */


						index = j;

						q_x1 = pnodes[2 * (pelems[numberofnodesperelem * index] - 1)];
						q_y1 = pnodes[2 * (pelems[numberofnodesperelem * index] - 1) + 1];

						q_x2 = pnodes[2 * (pelems[numberofnodesperelem * index + 1] - 1)];
						q_y2 = pnodes[2 * (pelems[numberofnodesperelem * index + 1] - 1) + 1];

						q_x3 = pnodes[2 * (pelems[numberofnodesperelem * index + 2] - 1)];
						q_y3 = pnodes[2 * (pelems[numberofnodesperelem * index + 2] - 1) + 1];

						RayTrace(tmp, intersect, q_t, q_eta,
								q_x1, q_y1, q_x2, q_y2, q_x3, q_y3,
								a, b, theta);

						RayTrim(tmp, a, b);

#ifdef DEBUG
						if (tmp.size()){
							std::cout << index << " "<<a<< " " << b << " "<< tmp << std::endl;
						}
#endif
						while (!q_elem.empty()){
							auto top = q_elem.front();
							q_elem.pop();
							// loop over neighbors.
							// fail safe strategy, if it is possible to intersect, must push.
							for (int32_t q_elem_i = 0; q_elem_i < 3; q_elem_i ++){

								index = pneighbors[3 * top + q_elem_i] - 1;

								// test if it is going to push them in.
								// first : valid index of element
								if (index > -1 && s_elem.find(index) == s_elem.end()) {

									/*
									 * calculate the open angle limit. And find the suitable one.
									 *
									 * Target function is max(min(abs(langle - angle), abs(rangle - angle)))
									 *
									 * if both gives 0, then push both into queue.
									 */

									q_x1 = pnodes[2 * (pelems[numberofnodesperelem * index] - 1)];
									q_y1 = pnodes[2 * (pelems[numberofnodesperelem * index] - 1) + 1];

									q_x2 = pnodes[2 * (pelems[numberofnodesperelem * index + 1] - 1)];
									q_y2 = pnodes[2 * (pelems[numberofnodesperelem * index + 1] - 1) + 1];

									q_x3 = pnodes[2 * (pelems[numberofnodesperelem * index + 2] - 1)];
									q_y3 = pnodes[2 * (pelems[numberofnodesperelem * index + 2] - 1) + 1];


									/*
									 * [a , b] 's angle for the element,
									 */

									RayTrace(tmp, intersect, q_t, q_eta,
											q_x1, q_y1, q_x2, q_y2, q_x3, q_y3,
											a, b, theta);

									if (intersect){
										q_elem.push(index);
										s_elem.insert(index);
										/*
										 * calculate integral over this element
										 *
										 * first step: calculate the intersection,
										 * second step: integral
										 */
										RayTrim(tmp, a, b);
#ifdef DEBUG
										if (tmp.size()){
											std::cout << index << " "<<a<< " " << b << " "<< tmp << std::endl;
										}
#endif
									}//end if
								}// end if
							}// end for
						}// end while
					}// end if(flag)
				} // end if
			} // end for
		} // end for
	} // end for
} // end

void DiscreteOrinates::RayTrace(std::vector<Real_t>& tmp, bool& intersect, Real_t& q_t, Real_t& q_eta,
			Real_t& q_x1, Real_t& q_y1, Real_t& q_x2,
			Real_t& q_y2, Real_t& q_x3, Real_t& q_y3,
			Real_t& a, Real_t& b, Real_t& theta){

	intersect = false;
	tmp.clear();
	// 1 - 2 intersects.
	if (fabs(INTERSECT_DET(q_x1, q_y1, q_x2, q_y2, theta)) < SCALE * MEX_EPS){
		// ray parallel for edge.
		if (fabs(INTERSECT_CROSS(q_x1, q_y1, q_x2, q_y2, a, b)) < SCALE * MEX_EPS){
			// it is lucky to be colinear.
			if (fabs(q_x1 - q_x2) + MEX_EPS < fabs(q_x1 - a) + fabs(q_x2 - a) ||
					fabs(q_y1 - q_y2) + MEX_EPS < fabs(q_y1 - b) + fabs(q_y2 - b) ){

				// outside
			}
			else {
				intersect = true;
				tmp.push_back(a);
				tmp.push_back(b);
			}
		}
		else {
			// intersect = false;
		}
	}
	else{
		// not parallel, then there is a intersect.
		q_t =INTERSECT_CROSS(q_x1, q_y1, q_x2, q_y2, a, b)/
				INTERSECT_DET(q_x1, q_y1, q_x2, q_y2, theta);

		q_eta = INTERSECT_DET(a, b, q_x2, q_y2, theta)/
				INTERSECT_DET(q_x1, q_y1, q_x2, q_y2, theta);
		/*
		 * q_t = 0 means colinear
		 *
		 */
		if (q_t >= 0 && q_eta >= -MEX_EPS && q_eta <= 1 + MEX_EPS) {
			intersect = true;
			tmp.push_back(q_eta * q_x1 + (1- q_eta) * q_x2);
			tmp.push_back(q_eta * q_y1 + (1 -q_eta) * q_y2);
		}
	}
	// 2 - 3
	if (fabs(INTERSECT_DET(q_x2, q_y2, q_x3, q_y3, theta)) < SCALE * MEX_EPS){
		// ray parallel for edge.
		if (fabs(INTERSECT_CROSS(q_x2, q_y2, q_x3, q_y3, a, b)) < SCALE * MEX_EPS){
			// it is lucky to be colinear.
			if (fabs(q_x2 - q_x3) + MEX_EPS < fabs(q_x2 - a) + fabs(q_x3 - a) ||
					fabs(q_y2 - q_y3) + MEX_EPS < fabs(q_y2 - b) + fabs(q_y3 - b) ){

				// outside
			}
			else {
				intersect = true;
				tmp.push_back(a);
				tmp.push_back(b);

			}
		}
		else {
			// intersect = false;
		}
	}
	else{
		// not parallel, then there is a intersect.
		q_t =INTERSECT_CROSS(q_x2, q_y2, q_x3, q_y3, a, b)/
				INTERSECT_DET(q_x2, q_y2, q_x3, q_y3, theta);

		q_eta = INTERSECT_DET(a, b, q_x3, q_y3, theta)/
				INTERSECT_DET(q_x2, q_y2, q_x3, q_y3, theta);

		/*
		 * q_t = 0 means colinear
		 *
		 */
		if (q_t >= 0 && q_eta >= -MEX_EPS && q_eta <= 1 + MEX_EPS) {
			intersect = true;
			tmp.push_back(q_eta * q_x2 + (1 - q_eta) * q_x3);
			tmp.push_back(q_eta * q_y2 + (1 - q_eta) * q_y3);
		}
	}
	// 3 - 1
	if (fabs(INTERSECT_DET(q_x3, q_y3, q_x1, q_y1, theta)) < SCALE * MEX_EPS){
		// ray parallel for edge.
		if (fabs(INTERSECT_CROSS(q_x3, q_y3, q_x1, q_y1, a, b)) < SCALE * MEX_EPS){
			// it is lucky to be colinear.
			if (fabs(q_x3 - q_x1) + MEX_EPS < fabs(q_x3 - a) + fabs(q_x1 - a) ||
					fabs(q_y3 - q_y1) + MEX_EPS < fabs(q_y3 - b) + fabs(q_y1 - b) ){

				// outside
			}
			else {
				intersect = true;
				tmp.push_back(a);
				tmp.push_back(b);
			}
		}
		else {
			// intersect = false;
		}
	}
	else{
		// not parallel, then there is a intersect.
		q_t =INTERSECT_CROSS(q_x3, q_y3, q_x1, q_y1, a, b)/
				INTERSECT_DET(q_x3, q_y3, q_x1, q_y1, theta);

		q_eta = INTERSECT_DET(a, b, q_x1, q_y1, theta)/
				INTERSECT_DET(q_x3, q_y3, q_x1, q_y1, theta);
		/*
		 * q_t = 0 means colinear
		 *
		 */
		if (q_t >= 0 && q_eta >= -MEX_EPS && q_eta <= 1 + MEX_EPS) {
			intersect = true;
			tmp.push_back(q_eta * q_x3 + (1 - q_eta) * q_x1);
			tmp.push_back(q_eta * q_y3 + (1 - q_eta) * q_y1);
		}
	}
}

void DiscreteOrinates::RayTrim(std::vector<Real_t>& tmp, Real_t &a, Real_t &b){
	for (int32_t tmp_i = 0; tmp_i < tmp.size()/2; tmp_i ++){

		if (fabs(a - tmp[2 * tmp_i]) + fabs(b - tmp[2 * tmp_i + 1]) < MEX_EPS){
			tmp.erase(tmp.begin() + 2 * tmp_i, tmp.begin() + 2 * tmp_i + 2);
			tmp_i --;
		}
	}

	if (tmp.size() == 6) {

	/*
	 * remove duplicate coordinates
	 */
	// 2 - 3 duplicates
		if (fabs(tmp[2] - tmp[4]) + fabs(tmp[3] - tmp[5]) < MEX_EPS){
			tmp.erase(tmp.begin() + 2, tmp.begin() + 4);
		}
	// 3 - 1 duplicates or 1 - 2 duplicates
		else {
			tmp.erase(tmp.begin(), tmp.begin() + 2);
		}
	}

	if (tmp.size() == 4) {
	/*
	 * remove duplicates
	 */
		if (fabs(tmp[0] - tmp[2]) + fabs(tmp[1] - tmp[3]) < MEX_EPS) {
			tmp.clear();
		}
		else {
			// reorder
			if (pow(tmp[0] - a, 2) + pow(tmp[1] - b, 2) > pow(tmp[2] - a, 2) + pow(tmp[3] - b, 2)){
				swap(tmp[0], tmp[2]);
				swap(tmp[1], tmp[3]);
			}
		}
	}
}

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

