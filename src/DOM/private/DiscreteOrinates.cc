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
	Ray.resize(nAngle);

}

DiscreteOrinates::DiscreteOrinates(MatlabPtr _nAngle, MatlabPtr _initAngle) noexcept {
	nAngle = static_cast<int32_t>(*mxGetPr(_nAngle));
	initAngle = static_cast<int32_t>(*mxGetPr(_initAngle));
	Ray.resize(nAngle);
}

DiscreteOrinates::~DiscreteOrinates() {

	Ray.clear();
	RHS.clear();
	Source.clear();
	Average.clear();
	Sigma_t.clear();
	Sigma_s.clear();

#ifdef DEBUG
	mexPrintf("Discrete Orinates Method(DOM) detached.\n");
#endif
}

void DiscreteOrinates::RayInt(MatlabPtr nodes, MatlabPtr elems,
		MatlabPtr neighbors){


	auto numberofnodes        = mxGetN(nodes);
	auto numberofelems        = mxGetN(elems);
	auto numberofnodesperelem = mxGetM(elems);

	auto pnodes               = mxGetPr(nodes);
	auto pelems               = (int32_t*)mxGetPr(elems);
	auto pneighbors           = (int32_t*)mxGetPr(neighbors);


	Real_t theta ;
	std::vector<bool> visited;
	/*
	 * resize Ray to hold raylets
	 */
	for (int32_t i = 0; i < nAngle; i++) {
		Ray[i].resize(numberofnodes);
	}
	/*
	 * main part
	 */
	for (int32_t i = 0; i < nAngle; i++) {
		theta = initAngle + 2 * i * M_PIl / nAngle;
		visited.clear();
		visited.resize(numberofnodes, false);
		RayIntHelper(numberofelems, numberofnodesperelem,
				numberofnodes,
				pelems, pnodes, pneighbors, i , theta, visited);
	}
}

void DiscreteOrinates::RayIntHelper(size_t& numberofelems, size_t& numberofnodesperelem,
		size_t& numberofnodes,
		int32_t* pelems, Real_t* pnodes, int32_t* pneighbors,int32_t& i, Real_t& theta,
		std::vector<bool>& visited){


	mwSize vertex, e_vl, e_vr;
	Real_t x1, x2, y1, y2,  a, b;
	Real_t ret, t, eta;
	Real_t q_x1, q_y1, q_x2, q_y2, q_x3, q_y3;
	Real_t q_eta, q_t;
	int32_t index, j, k;
	bool intersect;
	Raylet tmp_ray;

	std::vector<Real_t> tmp;
	std::unordered_set<int32_t> s_elem;
	std::queue<int32_t> q_elem;


	auto nproc = omp_get_num_procs();

	omp_lock_t lock;
	omp_init_lock(&lock);

	/*
	 * for loop, main part
	 */
#pragma omp parallel private(vertex, e_vl, e_vr,x1, x2, y1, y2,  a, b , ret, t,\
		eta,q_x1, q_y1, q_x2, q_y2, q_x3, q_y3, q_eta, q_t, index, j, k,intersect ,\
		tmp_ray, tmp, s_elem, q_elem) num_threads(nproc)
{

	for (j = 0; j < numberofelems; j++) {
		for (k = 0; k < numberofnodesperelem; k++){
			/*
			 * index of vertex in nodes.
			 */
			vertex = pelems[numberofnodesperelem * j + k] - 1;


			omp_set_lock(&lock);
			if (visited[vertex] == false){
				visited[vertex] = true;
				omp_unset_lock(&lock);
			}
			else {
				omp_unset_lock(&lock);
				continue;
			}

			if (visited[vertex] == true){
				/*
				 * calculate the ray intersects with the boundary.
				 */
				a = pnodes[2 * vertex    ];
				b = pnodes[2 * vertex + 1];

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

				if (tmp.size()) {
						tmp_ray.elem = index;
					if (tmp.size() == 2) {
						tmp_ray.first[0] = a;
						tmp_ray.first[1] = b;
						tmp_ray.second[0] = tmp[0];
						tmp_ray.second[1] = tmp[1];
					}
					else {
						tmp_ray.first[0] = tmp[0];
						tmp_ray.first[1] = tmp[1];
						tmp_ray.second[0] = tmp[2];
						tmp_ray.second[1] = tmp[3];
					}

					/*
					 * remove duplicates, since sorted, only test first two.
					 */

					if (!Ray[i][vertex].size() || fabs(tmp_ray.first[0] - Ray[i][vertex].back().first[0]) > MEX_EPS ||
						fabs(tmp_ray.first[1] - Ray[i][vertex].back().first[1]) > MEX_EPS){

						Ray[i][vertex].push_back(tmp_ray);
					}
				}
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

								if (tmp.size()) {
										tmp_ray.elem = index;
									if (tmp.size() == 2) {
										tmp_ray.first[0] = a;
										tmp_ray.first[1] = b;
										tmp_ray.second[0] = tmp[0];
										tmp_ray.second[1] = tmp[1];
									}
									else {
										tmp_ray.first[0] = tmp[0];
										tmp_ray.first[1] = tmp[1];
										tmp_ray.second[0] = tmp[2];
										tmp_ray.second[1] = tmp[3];
									}

									/*
									 * remove duplicates, since sorted, only test first two.
									 */
									if (!Ray[i][vertex].size() || fabs(tmp_ray.first[0] - Ray[i][vertex].back().first[0]) > MEX_EPS ||
										fabs(tmp_ray.first[1] - Ray[i][vertex].back().first[1]) > MEX_EPS){

										Ray[i][vertex].push_back(tmp_ray);

									}
								}
							}
						}
					}
				}
				Ray[i][vertex].shrink_to_fit();
			}
		}
	}// main part
	/*
	 * for loop
	 */
	omp_destroy_lock(&lock);
}
}

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
		//remove duplicates
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

/*
 * depreciated later.
 */
void DiscreteOrinates::RayShow(){

	int32_t tmp_i, tmp_j;
	size_t tmp_total = 0;
	if (nAngle != 0) {
		for (int32_t i = 0 ; i < nAngle; i++){
			tmp_i = Ray[i].size();
			for (int32_t j = 0; j < tmp_i; j++){
				tmp_j = Ray[i][j].size();
				tmp_total += tmp_j * 40;

//				for (int32_t k = 0; k < tmp_j; k++) {
//					std::cout << i << "th Angle, "
//							<< j << "th node, "
//							<< k << "th raylet: passes through "
//							<< Ray[i][j][k].elem << ", starting from "
//							<< Ray[i][j][k].first << " --> "
//							<< Ray[i][j][k].second << std::endl;
//				}
			}
		}
	}
	std::cout
	<< tmp_total / 1024.0/ 1024.0/ 1024.0
	<< " GBytes used in Ray storage."
	<< std::endl;
}

/*
 * Source Iteration.
 *
 */

void DiscreteOrinates::SourceIteration_init(){
	/*
	 * consider source as f(x).
	 */

	mxAssert(nAngle > 0,
			"DiscreteOrinates::SourceIteration_init::nAngle non positive.\n");

	auto numberofnodes = Ray[0].size();

	Source.resize(numberofnodes);
	Sigma_t.resize(numberofnodes);
	Sigma_s.resize(numberofnodes);

	RHS.resize(numberofnodes);
	Average.resize(numberofnodes);

	Output.resize(nAngle);

	for (int32_t s_i; s_i < nAngle; s_i++) {
		Output[s_i].resize(numberofnodes);
	}
}

void DiscreteOrinates::SourceIteration_port(MatlabPtr Fcn,
		MatlabPtr Sigma_t_Fcn, MatlabPtr Sigma_s_Fcn){


	auto interp        = mxGetPr(Fcn);
	auto interp_t      = mxGetPr(Sigma_t_Fcn);
	auto interp_s      = mxGetPr(Sigma_s_Fcn);

	mxAssert(nAngle > 0,
			"DiscreteOrinates::SourceIteration_port::nAngle non positive.\n");

	auto numberofnodes = Ray[0].size();

	mxAssert(mxGetElementSize(Fcn) == numberofnodes,
			"DiscreteOrinates::SourceIteration_port::Dimension does not match.\n");
	mxAssert(mxGetElementSize(Sigma_t_Fcn) == numberofnodes,
			"DiscreteOrinates::SourceIteration_port::Dimension does not match.\n");
	mxAssert(mxGetElementSize(Sigma_s_Fcn) == numberofnodes,
			"DiscreteOrinates::SourceIteration_port::Dimension does not match.\n");

	memcpy(&Source[0], interp, sizeof(Real_t) * numberofnodes);
	memcpy(&Sigma_t[0], interp_t, sizeof(Real_t) * numberofnodes);
	memcpy(&Sigma_s[0], interp_s, sizeof(Real_t) * numberofnodes);

}

void DiscreteOrinates::SourceIteration_iter(MatlabPtr nodes, MatlabPtr elems){

	auto pnodes        = mxGetPr(nodes);
	auto pelems        = (int32_t *)mxGetPr(elems);
	auto numberofnodesperelem = mxGetM(elems);

	auto numberofnodes = mxGetN(nodes);

	for (int32_t s_j = 0; s_j < numberofnodes; s_j++){

		RHS[s_j] = Sigma_s[s_j] * Average[s_j];
		RHS[s_j] += Source[s_j];
	}


	mwSize vertex_1, vertex_2, vertex_3;
	Real_t x1, y1, x2, y2, x3, y3, det, lambda, eta, length, accum_s, accum_v;

	Real_t lv, rv, ls, rs;

	for (int32_t s_i = 0; s_i < nAngle; s_i++) {
		for (int32_t s_j = 0; s_j < numberofnodes; s_j++){
			accum_s = 0.;
			accum_v = 0.;
			if (Ray[s_i][s_j].size()){
				for (auto it : Ray[s_i][s_j]){
					vertex_1 = pelems[it.elem * numberofnodesperelem ] - 1;
					vertex_2 = pelems[it.elem * numberofnodesperelem + 1] - 1;
					vertex_3 = pelems[it.elem * numberofnodesperelem + 2] - 1;

					x1 = pnodes[2 * vertex_1];
					y1 = pnodes[2 * vertex_1 + 1];
					x2 = pnodes[2 * vertex_2];
					y2 = pnodes[2 * vertex_2 + 1];
					x3 = pnodes[2 * vertex_3];
					y3 = pnodes[2 * vertex_3 + 1];

					/*
					 * first node
					 */
					det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);

					eta = ((y3 - y1) * (it.first[0] - x3) + (x1 - x3) * (it.first[1] - y3));
					eta /= det;

					lambda = (y2 - y3) * (it.first[0] - x3) + (x3 -  x2) * (it.first[1] - y3);
					lambda /= det;

					lv = lambda * RHS[vertex_1] + eta * RHS[vertex_2] +
							(1 - lambda - eta) * RHS[vertex_3];
					ls = lambda * Sigma_t[vertex_1] + eta * Sigma_t[vertex_2] +
							(1 - lambda - eta) * Sigma_t[vertex_3];

					/*
					 * second node
					 */
					eta = ((y3 - y1) * (it.second[0] - x3) + (x1 - x3) * (it.second[1] - y3));
					eta /= det;

					lambda = (y2 - y3) * (it.second[0] - x3) + (x3 -  x2) * (it.second[1] - y3);
					lambda /= det;

					rv = lambda * RHS[vertex_1] + eta * RHS[vertex_2] +
							(1 - lambda - eta) * RHS[vertex_3];
					rs = lambda * Sigma_t[vertex_1] + eta * Sigma_t[vertex_2] +
							(1 - lambda - eta) * Sigma_t[vertex_3];

					/*
					 * length
					 */
					length = sqrt(pow(it.first[0] - it.second[0], 2) + pow(it.first[1] - it.second[1], 2));

					/*
					 *
					 *
					 * Note: since exp(-x) function is decreasing too fast, the early stage of
					 * integral should be very careful. Now using Simpson formula.
					 *
					 * first part only accurate if sigma_a is linear.
					 * second part never accurate, can obtain O(h^2) accuracy locally.
					 *
					 */
					accum_v += exp(-accum_s) * lv * length/6.0;

					accum_s += (0.5 * rs + 1.5 * ls) * length/ 4.0;

					accum_v += exp(-accum_s) * (lv + rv) * length / 3.0;

					accum_s += (1.5 * rs + 0.5 * ls) * length/ 4.0;

					accum_v += exp(-accum_s) * rv * length/6.0;
				}
			}
			else{
				accum_v = 0.;
			}
			Output[s_i][s_j] = accum_v;
		}
	}

	/*
	 * update Average
	 */
	for (int32_t s_j = 0; s_j < numberofnodes; s_j++){
		Average[s_j] = 0.;
		for (int32_t s_i = 0; s_i < nAngle; s_i++) {
			Average[s_j] += Output[s_i][s_j];
		}
		Average[s_j] /= nAngle;
	}
}

void DiscreteOrinates::SourceIteration_accl(MatlabPtr delta){

	auto delta_ptr = mxGetPr(delta);
	auto numberofnodes = mxGetM(delta);

	mxAssert(numberofnodes == Average.size(), "DiscreteOrinates::SourceIteration_accl::Dimension does not match.\n");

	for (int32_t delta_i  = 0; delta_i < numberofnodes; delta_i ++) {
		Average[delta_i] += delta_ptr[delta_i];
	}
}

void DiscreteOrinates::SourceIteration_set(MatlabPtr ave){
	auto ave_ptr = mxGetPr(ave);
	auto numberofnodes = mxGetM(ave);

	mxAssert(numberofnodes == Average.size(), "DiscreteOrinates::SourceIteration_set::Dimension does not match.\n");

	for (int32_t delta_i  = 0; delta_i < numberofnodes; delta_i ++) {
		Average[delta_i] = ave_ptr[delta_i];
	}
}

void DiscreteOrinates::SourceIteration_beam(MatlabPtr incoming){

	/*
	 * incoming is the boundary condition.
	 *
	 * On each edge nodes, there are Discrete Orinates for them.
	 *
	 * Core algorithm is easy.
	 */

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
	InputArguments input(nrhs, prhs, 4);
	OutputArguments output(nlhs, plhs, 0);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	DOM->RayInt(CAST(prhs[1]), CAST(prhs[2]),
			CAST(prhs[3]));
}

MEX_DEFINE(rayshow) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	DOM->RayShow();
}

MEX_DEFINE(ray_build) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
// better than si_build. non-linear building
// depreciated
}

MEX_DEFINE(ray_build_omp) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
// better than si_build. non-linear building
// faster building

	InputArguments input(nrhs, prhs, 5);
	OutputArguments output(nlhs, plhs, 2);

	auto pnodes        = mxGetPr(prhs[1]);
	auto pelems        = (int32_t *)mxGetPr(prhs[2]);
	auto numberofnodesperelem = mxGetM(prhs[2]);

	auto numberofnodes = mxGetN(prhs[1]);

	auto Sigma_t       = mxGetPr(prhs[3]);
	auto Sigma_s       = mxGetPr(prhs[4]);

	plhs[0] = mxCreateNumericMatrix(numberofnodes, numberofnodes, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(numberofnodes, numberofnodes, mxDOUBLE_CLASS, mxREAL);

	auto ptr = mxGetPr(plhs[0]);
	auto mptr = mxGetPr(plhs[1]);


	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	auto nAngle = DOM->nAngle;


	mwSize vertex_1, vertex_2, vertex_3;
	Real_t x1, y1, x2, y2, x3, y3, det, lambda1,lambda2, eta1, eta2, length, accum_s;
	Real_t lv, rv, ls, rs;
	Real_t common1, common2, common3, common4;
	int32_t s_i, s_j;

	omp_set_num_threads(omp_get_num_procs());

#pragma omp parallel for private(s_i, s_j, vertex_1, vertex_2, vertex_3,x1, y1, x2, y2, x3, y3, det,\
		lambda1,lambda2, eta1, eta2, length, accum_s, lv, rv, ls, rs, common1, common2, common3, common4) schedule(dynamic,1) collapse(2)
	for (s_i = 0; s_i < nAngle; s_i++) {
		for (s_j = 0; s_j < numberofnodes; s_j++) {
			accum_s = 0.;
			if (DOM->Ray[s_i][s_j].size()){
				for (auto it : DOM->Ray[s_i][s_j]){

					vertex_1 = pelems[it.elem * numberofnodesperelem ] - 1;
					vertex_2 = pelems[it.elem * numberofnodesperelem + 1] - 1;
					vertex_3 = pelems[it.elem * numberofnodesperelem + 2] - 1;

					x1 = pnodes[2 * vertex_1];
					y1 = pnodes[2 * vertex_1 + 1];
					x2 = pnodes[2 * vertex_2];
					y2 = pnodes[2 * vertex_2 + 1];
					x3 = pnodes[2 * vertex_3];
					y3 = pnodes[2 * vertex_3 + 1];

					/*
					 * length
					 */
					length = sqrt(pow(it.first[0] - it.second[0], 2) + pow(it.first[1] - it.second[1], 2));
					/*
					 * first node
					 */
					det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);

					eta1 = ((y3 - y1) * (it.first[0] - x3) + (x1 - x3) * (it.first[1] - y3));
					eta1 /= det;

					lambda1 = (y2 - y3) * (it.first[0] - x3) + (x3 -  x2) * (it.first[1] - y3);
					lambda1 /= det;

					ls = lambda1 * Sigma_t[vertex_1] + eta1 * Sigma_t[vertex_2] +
							(1 - lambda1 - eta1) * Sigma_t[vertex_3];
					lv = lambda1 * Sigma_s[vertex_1] + eta1 * Sigma_s[vertex_2] +
							(1 - lambda1 - eta1) * Sigma_s[vertex_3];


					/*
					 * second node
					 */
					eta2 = ((y3 - y1) * (it.second[0] - x3) + (x1 - x3) * (it.second[1] - y3));
					eta2 /= det;

					lambda2 = (y2 - y3) * (it.second[0] - x3) + (x3 -  x2) * (it.second[1] - y3);
					lambda2 /= det;

					rs = lambda2 * Sigma_t[vertex_1] + eta2 * Sigma_t[vertex_2] +
							(1 - lambda2 - eta2) * Sigma_t[vertex_3];

					rv = lambda2 * Sigma_s[vertex_1] + eta2 * Sigma_s[vertex_2] +
							(1 - lambda2 - eta2) * Sigma_s[vertex_3];


					/*
					 * assume linear function, lose high order information.
					 * should lose some radiation from nearby.
					 */

					common1 = (ZeroOrder(ls, rs, length));
					common2 = (FirstOrder(ls, rs, length));
					common3 = (SecondOrder(ls, rs, length));
					common4 = length * exp(-accum_s);

					mptr[numberofnodes * s_j + vertex_1] +=
							common4 *
							(common1 * lambda1 + common2 * (lambda2 - lambda1));

					mptr[numberofnodes * s_j + vertex_2] +=
							common4 *
							(common1 * eta1 + common2 * (eta2 - eta1));

					mptr[numberofnodes * s_j + vertex_3] +=
							common4 *
							(common1 * (1 - eta1 - lambda1) +
							 common2 * (lambda1 - lambda2 + eta1 - eta2));

					ptr[numberofnodes * s_j + vertex_1] +=
							common4 *
							(lv *
							(common1 * lambda1 + common2 * (lambda2 - lambda1))
							+
							(rv - lv) *
							(common2 * lambda1 + common3 *(lambda2 - lambda1))
							);
					ptr[numberofnodes * s_j + vertex_2] +=
							common4 *
							(lv *
							(common1 * eta1 + common2 * (eta2 - eta1))
							 +
							 (rv - lv) *
							 (common2 * eta1 + common3 * (eta2 - eta1))
							);

					ptr[numberofnodes * s_j + vertex_3] +=
							common4 *
							(lv *
							(common1 * (1 - eta1 - lambda1) + common2 * (lambda1 - lambda2 + eta1 - eta2))
							+
							(rv - lv) *
							(common2 * (1 - eta1 - lambda1) + common3 * (lambda1 - lambda2 + eta1 - eta2))
							);

					accum_s += length * (ls + rs)/2;
				}
			}
		}
	}
}

MEX_DEFINE(ray_scatter_grad)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 6);
	OutputArguments output(nlhs, plhs, 1);

	auto pnodes        = mxGetPr(prhs[1]);
	auto pelems        = (int32_t *)mxGetPr(prhs[2]);
	auto numberofnodesperelem = mxGetM(prhs[2]);

	auto numberofnodes = mxGetN(prhs[1]);

	auto uptr          = mxGetPr(prhs[4]);
	auto vptr          = mxGetPr(prhs[3]);
	auto Sigma_t       = mxGetPr(prhs[5]);

	plhs[0] = mxCreateNumericMatrix(numberofnodes, 1, mxDOUBLE_CLASS, mxREAL);

	auto ptr = mxGetPr(plhs[0]);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	auto nAngle = DOM->nAngle;

	mwSize vertex_1, vertex_2, vertex_3;
	Real_t x1, y1, x2, y2, x3, y3, det, lambda1,lambda2, eta1, eta2, length, accum_s;
	Real_t lv, rv, ls, rs;
	Real_t common1, common2, common3, common4;
	Real_t coeff1, coeff2, coeff3;
	int32_t s_i, s_j;

	std::unordered_map<int32_t, Real_t> record;

	omp_set_num_threads(omp_get_num_procs());

#pragma omp parallel for private(s_i, s_j, vertex_1, vertex_2, vertex_3,x1, y1, x2, y2, x3, y3, det,\
		lambda1,lambda2, eta1, eta2, length, accum_s, lv, rv, ls, rs, common1, common2, common3, common4,\
		record, coeff1, coeff2, coeff3) schedule(dynamic,1) collapse(2)
	for (s_i = 0; s_i < nAngle; s_i++) {
		for (s_j = 0; s_j < numberofnodes; s_j++) {
			accum_s = 0.;
			record.clear();
			if (DOM->Ray[s_i][s_j].size()){
				// records all information along the ray
				for (auto it : DOM->Ray[s_i][s_j]){

					vertex_1 = pelems[it.elem * numberofnodesperelem ] - 1;
					vertex_2 = pelems[it.elem * numberofnodesperelem + 1] - 1;
					vertex_3 = pelems[it.elem * numberofnodesperelem + 2] - 1;

					x1 = pnodes[2 * vertex_1];
					y1 = pnodes[2 * vertex_1 + 1];
					x2 = pnodes[2 * vertex_2];
					y2 = pnodes[2 * vertex_2 + 1];
					x3 = pnodes[2 * vertex_3];
					y3 = pnodes[2 * vertex_3 + 1];

					/*
					 * length
					 */
					length = sqrt(pow(it.first[0] - it.second[0], 2) + pow(it.first[1] - it.second[1], 2));
					/*
					 * first node
					 */
					det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);

					eta1 = ((y3 - y1) * (it.first[0] - x3) + (x1 - x3) * (it.first[1] - y3));
					eta1 /= det;

					lambda1 = (y2 - y3) * (it.first[0] - x3) + (x3 -  x2) * (it.first[1] - y3);
					lambda1 /= det;

					ls = lambda1 * Sigma_t[vertex_1] + eta1 * Sigma_t[vertex_2] +
							(1 - lambda1 - eta1) * Sigma_t[vertex_3];



					/*
					 * second node
					 */
					eta2 = ((y3 - y1) * (it.second[0] - x3) + (x1 - x3) * (it.second[1] - y3));
					eta2 /= det;

					lambda2 = (y2 - y3) * (it.second[0] - x3) + (x3 -  x2) * (it.second[1] - y3);
					lambda2 /= det;

					rs = lambda2 * Sigma_t[vertex_1] + eta2 * Sigma_t[vertex_2] +
							(1 - lambda2 - eta2) * Sigma_t[vertex_3];


					coeff1 = exp(-accum_s) * lambda1 * length/6.0;
					coeff2 = exp(-accum_s) * eta1 * length/ 6.0;
					coeff3 = exp(-accum_s) * (1 - lambda1 - eta1) * length/6.0;


					for (auto record_it : record) {
						ptr[record_it.first] -= uptr[s_j] * record_it.second * vptr[vertex_1] * coeff1;
						ptr[record_it.first] -= uptr[s_j] * record_it.second * vptr[vertex_2] * coeff2;
						ptr[record_it.first] -= uptr[s_j] * record_it.second * vptr[vertex_3] * coeff3;
					}

					accum_s += (0.5 * rs + 1.5 * ls) * length/ 4.0;

					if (record.find(vertex_1) != record.end()) {
						record[vertex_1] += (0.5 * lambda2 + 1.5 * lambda1) * length / 4.0;
					}
					else {
						record[vertex_1] = (0.5 * lambda2  + 1.5 * lambda2) * length / 4.0;
					}

					if (record.find(vertex_2) != record.end()) {
						record[vertex_2] += (0.5 * eta2 + 1.5 * eta1)* length /4.0;
					}
					else {
						record[vertex_2] = (0.5 * eta2 + 1.5 * eta1)* length /4.0;
					}

					if (record.find(vertex_3) != record.end()) {
						record[vertex_3] += (0.5 * (1 - lambda2 - eta2) + 1.5 * (1 - lambda1 - eta1)) * length / 4.0;
					}
					else {
						record[vertex_3] = (0.5 * (1 - lambda2 - eta2) + 1.5 * (1 - lambda1 - eta1)) * length / 4.0;
					}


					coeff1 = exp(-accum_s) * (lambda1 + lambda2) * length/3.0;
					coeff2 = exp(-accum_s) * (eta1 + eta2) * length/ 3.0;
					coeff3 = exp(-accum_s) * (2 - lambda1 - eta1 - lambda2 - eta2) * length/3.0;

					for (auto record_it : record) {
						ptr[record_it.first] -= uptr[s_j] * record_it.second * vptr[vertex_1] * coeff1;
						ptr[record_it.first] -= uptr[s_j] * record_it.second * vptr[vertex_2] * coeff2;
						ptr[record_it.first] -= uptr[s_j] * record_it.second * vptr[vertex_3] * coeff3;
					}

//
//					ptr[numberofnodes * s_j + vertex_1] += exp(-accum_s) * (lambda1 + lambda2) * length/3.0;
//					ptr[numberofnodes * s_j + vertex_2] += exp(-accum_s) * (eta1 + eta2) * length/ 3.0;
//					ptr[numberofnodes * s_j + vertex_3] += exp(-accum_s) * (2 - lambda1 - eta1 - lambda2 - eta2) * length/3.0;

					accum_s += (1.5 * rs + 0.5 * ls) * length / 4.0;

					if (record.find(vertex_1) != record.end()) {
						record[vertex_1] += (1.5 * lambda2 + 0.5 * lambda1) * length / 4.0;
					}
					else {
						record[vertex_1] = (1.5 * lambda2  + 0.5 * lambda2) * length / 4.0;
					}

					if (record.find(vertex_2) != record.end()) {
						record[vertex_2] += (1.5 * eta2 + 0.5 * eta1)* length /4.0;
					}
					else {
						record[vertex_2] = (1.5 * eta2 + 0.5 * eta1)* length /4.0;
					}

					if (record.find(vertex_3) != record.end()) {
						record[vertex_3] += (1.5 * (1 - lambda2 - eta2) + 0.5 * (1 - lambda1 - eta1)) * length / 4.0;
					}
					else {
						record[vertex_3] = (1.5 * (1 - lambda2 - eta2) + 0.5 * (1 - lambda1 - eta1)) * length / 4.0;
					}

					coeff1 = exp(-accum_s) * lambda2 * length/6.0;
					coeff2 = exp(-accum_s) * eta2 * length/ 6.0;
					coeff3 = exp(-accum_s) * (1 - lambda2 - eta2) * length/6.0;

					for (auto record_it : record) {
						ptr[record_it.first] -= uptr[s_j] * record_it.second * vptr[vertex_1] * coeff1;
						ptr[record_it.first] -= uptr[s_j] * record_it.second * vptr[vertex_2] * coeff2;
						ptr[record_it.first] -= uptr[s_j] * record_it.second * vptr[vertex_3] * coeff3;
					}
//					ptr[numberofnodes * s_j + vertex_1] += exp(-accum_s) * lambda2 * length/6.0;
//					ptr[numberofnodes * s_j + vertex_2] += exp(-accum_s) * eta2 * length/ 6.0;
//					ptr[numberofnodes * s_j + vertex_3] += exp(-accum_s) * (1 - lambda2 - eta2) * length/6.0;
				}
			}
		}
	}

}


MEX_DEFINE(si_init) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	DOM->SourceIteration_init();
}

MEX_DEFINE(si_import) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 4);
	OutputArguments output(nlhs, plhs, 0);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	DOM->SourceIteration_port(CAST(prhs[1]), CAST(prhs[2]), CAST(prhs[3]));
}

MEX_DEFINE(si_iter) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 0);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));


	DOM->SourceIteration_iter(CAST(prhs[1]), CAST(prhs[2]));
}

MEX_DEFINE(si_output) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 1);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	plhs[0] = mxCreateNumericMatrix(DOM->Average.size(),1, mxDOUBLE_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[0]), &(DOM->Average[0]), DOM->Average.size()*sizeof(Real_t));
}

MEX_DEFINE(si_dsa) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 0);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	DOM->SourceIteration_accl(CAST(prhs[1]));
}

MEX_DEFINE(si_set) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 0);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	DOM->SourceIteration_set(CAST(prhs[1]));
}

MEX_DEFINE(si_build)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	/*
	 * input: id, nodes, elems, sigma_t
	 */
	InputArguments input(nrhs, prhs, 4);
	OutputArguments output(nlhs, plhs, 1);

	auto pnodes        = mxGetPr(prhs[1]);
	auto pelems        = (int32_t *)mxGetPr(prhs[2]);
	auto numberofnodesperelem = mxGetM(prhs[2]);

	auto numberofnodes = mxGetN(prhs[1]);

	auto Sigma_t       = mxGetPr(prhs[3]);

	plhs[0] = mxCreateNumericMatrix(numberofnodes, numberofnodes, mxDOUBLE_CLASS, mxREAL);
	auto ptr = mxGetPr(plhs[0]);

	mwSize vertex_1, vertex_2, vertex_3;
	Real_t x1, y1, x2, y2, x3, y3, det, lambda1,lambda2, eta1, eta2, length, accum_s;
	Real_t lv, rv, ls, rs;

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	auto nAngle = DOM->nAngle;

	/*
	 * building the matrix will take around 1x time of one iteration.
	 *
	 * Since the matrix is a small valued matrix, various methods can be
	 * used to make the converging process faster.
	 *
	 * Now use pcg or gmres.
	 */

	for (int32_t s_i = 0; s_i < nAngle; s_i++) {
		for (int32_t s_j = 0; s_j < numberofnodes; s_j++){
			accum_s = 0.;
			if (DOM->Ray[s_i][s_j].size()){
				for (auto it : DOM->Ray[s_i][s_j]){
					vertex_1 = pelems[it.elem * numberofnodesperelem ] - 1;
					vertex_2 = pelems[it.elem * numberofnodesperelem + 1] - 1;
					vertex_3 = pelems[it.elem * numberofnodesperelem + 2] - 1;

					x1 = pnodes[2 * vertex_1];
					y1 = pnodes[2 * vertex_1 + 1];
					x2 = pnodes[2 * vertex_2];
					y2 = pnodes[2 * vertex_2 + 1];
					x3 = pnodes[2 * vertex_3];
					y3 = pnodes[2 * vertex_3 + 1];

					/*
					 * length
					 */
					length = sqrt(pow(it.first[0] - it.second[0], 2) + pow(it.first[1] - it.second[1], 2));
					/*
					 * first node
					 */
					det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);

					eta1 = ((y3 - y1) * (it.first[0] - x3) + (x1 - x3) * (it.first[1] - y3));
					eta1 /= det;

					lambda1 = (y2 - y3) * (it.first[0] - x3) + (x3 -  x2) * (it.first[1] - y3);
					lambda1 /= det;

					ls = lambda1 * Sigma_t[vertex_1] + eta1 * Sigma_t[vertex_2] +
							(1 - lambda1 - eta1) * Sigma_t[vertex_3];


					/*
					 * second node
					 */
					eta2 = ((y3 - y1) * (it.second[0] - x3) + (x1 - x3) * (it.second[1] - y3));
					eta2 /= det;

					lambda2 = (y2 - y3) * (it.second[0] - x3) + (x3 -  x2) * (it.second[1] - y3);
					lambda2 /= det;

					rs = lambda2 * Sigma_t[vertex_1] + eta2 * Sigma_t[vertex_2] +
							(1 - lambda2 - eta2) * Sigma_t[vertex_3];

					/*
					 * inserting
					 */
					ptr[numberofnodes * s_j + vertex_1] += exp(-accum_s) * lambda1 * length/6.0;
					ptr[numberofnodes * s_j + vertex_2] += exp(-accum_s) * eta1 * length/ 6.0;
					ptr[numberofnodes * s_j + vertex_3] += exp(-accum_s) * (1 - lambda1 - eta1) * length/6.0;

					accum_s += (0.5 * rs + 1.5 * ls) * length/ 4.0;

					ptr[numberofnodes * s_j + vertex_1] += exp(-accum_s) * (lambda1 + lambda2) * length/3.0;
					ptr[numberofnodes * s_j + vertex_2] += exp(-accum_s) * (eta1 + eta2) * length/ 3.0;
					ptr[numberofnodes * s_j + vertex_3] += exp(-accum_s) * (2 - lambda1 - eta1 - lambda2 - eta2) * length/3.0;

					accum_s += (1.5 * rs + 0.5 * ls) * length / 4.0;

					ptr[numberofnodes * s_j + vertex_1] += exp(-accum_s) * lambda2 * length/6.0;
					ptr[numberofnodes * s_j + vertex_2] += exp(-accum_s) * eta2 * length/ 6.0;
					ptr[numberofnodes * s_j + vertex_3] += exp(-accum_s) * (1 - lambda2 - eta2) * length/6.0;
				}
			}
		}
	}
}

MEX_DEFINE(si_build_omp)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	/*
	 * input: id, nodes, elems, sigma_t
	 */
	InputArguments input(nrhs, prhs, 4);
	OutputArguments output(nlhs, plhs, 1);

	auto pnodes        = mxGetPr(prhs[1]);
	auto pelems        = (int32_t *)mxGetPr(prhs[2]);
	auto numberofnodesperelem = mxGetM(prhs[2]);

	auto numberofnodes = mxGetN(prhs[1]);

	auto Sigma_t       = mxGetPr(prhs[3]);

	plhs[0] = mxCreateNumericMatrix(numberofnodes, numberofnodes, mxDOUBLE_CLASS, mxREAL);
	auto ptr = mxGetPr(plhs[0]);

//	mwSize vertex_1, vertex_2, vertex_3;
//	Real_t x1, y1, x2, y2, x3, y3, det, lambda1,lambda2, eta1, eta2, length, accum_s;
//	Real_t lv, rv, ls, rs;

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	auto nAngle = DOM->nAngle;

	/*
	 * building the matrix will take around 1x time of one iteration.
	 *
	 * Since the matrix is a small valued matrix, various methods can be
	 * used to make the converging process faster.
	 *
	 * Now use pcg or gmres.
	 */


	mwSize vertex_1, vertex_2, vertex_3;
	Real_t x1, y1, x2, y2, x3, y3, det, lambda1,lambda2, eta1, eta2, length, accum_s;
	Real_t lv, rv, ls, rs;
	int32_t s_i, s_j;

	omp_set_num_threads(omp_get_num_procs());

#pragma omp parallel for private(s_i, s_j, vertex_1, vertex_2, vertex_3,x1, y1, x2, y2, x3, y3, det,\
		lambda1,lambda2, eta1, eta2, length, accum_s, lv, rv, ls, rs) schedule(dynamic,1) collapse(2)
	for (s_i = 0; s_i < nAngle; s_i++) {
		for (s_j = 0; s_j < numberofnodes; s_j++){
			accum_s = 0.;
			if (DOM->Ray[s_i][s_j].size()){
				for (auto it : DOM->Ray[s_i][s_j]){
					vertex_1 = pelems[it.elem * numberofnodesperelem ] - 1;
					vertex_2 = pelems[it.elem * numberofnodesperelem + 1] - 1;
					vertex_3 = pelems[it.elem * numberofnodesperelem + 2] - 1;

					x1 = pnodes[2 * vertex_1];
					y1 = pnodes[2 * vertex_1 + 1];
					x2 = pnodes[2 * vertex_2];
					y2 = pnodes[2 * vertex_2 + 1];
					x3 = pnodes[2 * vertex_3];
					y3 = pnodes[2 * vertex_3 + 1];

					/*
					 * length
					 */
					length = sqrt(pow(it.first[0] - it.second[0], 2) + pow(it.first[1] - it.second[1], 2));
					/*
					 * first node
					 */
					det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);

					eta1 = ((y3 - y1) * (it.first[0] - x3) + (x1 - x3) * (it.first[1] - y3));
					eta1 /= det;

					lambda1 = (y2 - y3) * (it.first[0] - x3) + (x3 -  x2) * (it.first[1] - y3);
					lambda1 /= det;

					ls = lambda1 * Sigma_t[vertex_1] + eta1 * Sigma_t[vertex_2] +
							(1 - lambda1 - eta1) * Sigma_t[vertex_3];


					/*
					 * second node
					 */
					eta2 = ((y3 - y1) * (it.second[0] - x3) + (x1 - x3) * (it.second[1] - y3));
					eta2 /= det;

					lambda2 = (y2 - y3) * (it.second[0] - x3) + (x3 -  x2) * (it.second[1] - y3);
					lambda2 /= det;

					rs = lambda2 * Sigma_t[vertex_1] + eta2 * Sigma_t[vertex_2] +
							(1 - lambda2 - eta2) * Sigma_t[vertex_3];

					/*
					 * inserting
					 */
					ptr[numberofnodes * s_j + vertex_1] += exp(-accum_s) * lambda1 * length/6.0;
					ptr[numberofnodes * s_j + vertex_2] += exp(-accum_s) * eta1 * length/ 6.0;
					ptr[numberofnodes * s_j + vertex_3] += exp(-accum_s) * (1 - lambda1 - eta1) * length/6.0;

					accum_s += (0.5 * rs + 1.5 * ls) * length/ 4.0;

					ptr[numberofnodes * s_j + vertex_1] += exp(-accum_s) * (lambda1 + lambda2) * length/3.0;
					ptr[numberofnodes * s_j + vertex_2] += exp(-accum_s) * (eta1 + eta2) * length/ 3.0;
					ptr[numberofnodes * s_j + vertex_3] += exp(-accum_s) * (2 - lambda1 - eta1 - lambda2 - eta2) * length/3.0;

					accum_s += (1.5 * rs + 0.5 * ls) * length / 4.0;

					ptr[numberofnodes * s_j + vertex_1] += exp(-accum_s) * lambda2 * length/6.0;
					ptr[numberofnodes * s_j + vertex_2] += exp(-accum_s) * eta2 * length/ 6.0;
					ptr[numberofnodes * s_j + vertex_3] += exp(-accum_s) * (1 - lambda2 - eta2) * length/6.0;
				}
			}
		}
	}
}

}

MEX_DISPATCH

