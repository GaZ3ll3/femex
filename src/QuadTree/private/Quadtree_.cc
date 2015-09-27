/*
 * quadtree.cc
 *
 *  Created on: Sep 18, 2015
 *      Author: lurker
 */

#include "QuadTree.h"

#include <mexplus.h>
#include "utils.h"

using namespace std;
using namespace mexplus;


template class mexplus::Session<QuadTree>;

namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 1);

	auto posptr = mxGetPr(prhs[0]);
	std::vector<double> pos {posptr[0], posptr[1]};
	auto szptr = mxGetPr(prhs[1]);

	output.set(0, Session<QuadTree>::create(new QuadTree(
			pos,
			*szptr)));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	auto root = Session<QuadTree>::get(input.get(0));

	root->clearPhotons();

	Session<QuadTree>::destroy(input.get(0));
}

MEX_DEFINE(import) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])  {
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 0);

	auto particles = mxGetPr(prhs[1]);
	auto numofparticle = mxGetN(prhs[1]);

	auto root = Session<QuadTree>::get(input.get(0));

	std::vector<double> pos{0., 0.};
	for (int i = 0; i < numofparticle; i++) {
		/*
		 * photon is (x, y, value) form
		 */
		pos[0] = particles[3 * i];
		pos[1] = particles[3 *  i + 1];
		auto ph = new Photon(pos, particles[3 * i + 2]);
		ph->id = i;
		root->addPhoton(ph);
	}

	root->updateAttribute();
}

MEX_DEFINE(split) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);

	auto root = Session<QuadTree>::get(input.get(0));
	root->setStatus(Status::ROOT);
	root->populate();
}

MEX_DEFINE(buildmatrix) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 1);

	auto root = Session<QuadTree>::get(input.get(0));
	auto num = root->getPhotonSize();
	auto sigma_t_ptr = mxGetPr(prhs[1]);
	auto theta_ptr = mxGetPr(prhs[2]);
	if (num == 0) {
		std::cout << "empty cell, stopped" << std::endl;
		return;
	}
	plhs[0] = mxCreateNumericMatrix(num, num,mxDOUBLE_CLASS, mxREAL);
	auto Radptr  = mxGetPr(plhs[0]);

	omp_set_num_threads(omp_get_num_procs());

	int i;

#pragma omp parallel for shared(sigma_t_ptr, theta_ptr, root, Radptr, num) private(i)
	for (i = 0; i < num; i++) {
		trasverse(*sigma_t_ptr, *theta_ptr, root->photons[i], root, Radptr, num);
	}
}
}

MEX_DISPATCH
