/*
 * MCRT.cpp
 *
 *  Created on: Jul 20, 2016
 *      Author: lurker
 */

#include "MCRT.h"

MCRT::MCRT(size_t n_photons_, size_t n,  double mu_s_, double mu_t_) {
	n_nodes = n;
	n_photons = n_photons_;
	albedo = mu_s_ / mu_t_;
	mu_t = mu_t_;
	assert(albedo <= 1);

	intensity.resize(n);
	charge.resize(n);

	for (size_t i = 0; i < n; i++) {
		intensity[i].resize(n);
		charge[i].resize(n);
	}
}

MCRT::~MCRT() {
	// TODO Auto-generated destructor stub
}

double MCRT::random() {
	return double(rand())/RAND_MAX;
}


void MCRT::set_charge() {
	double side = 1.0/n_nodes;
	double x;
	size_t i;
#pragma omp parallel for shared(side) private (i, x)
	for (i = 0; i <n_nodes; i++) {
		x = side/2.0 + i * side;
		for (size_t j = 0; j < n_nodes; j++) {
			double y = side/1.0 + j * side;
			double r2 =  ((x - 0.6) * (x - 0.6) + (y - 0.4) * (y - 0.4));
			if (r2 >= 0.15*0.15 && r2<=0.25*0.25) {
				charge[i][j] = 1.0;
			}
		}
	}
}

void MCRT::scatter(photon& p) {
	double tau, s;
	while(p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1) {
		tau = -1.0 * log(random());
		s = tau / mu_t;
		p.x += s * cos(p.theta);
		p.y += s * sin(p.theta);
		// particle moves to new location within mean-free path.
		if (p.x < 0 || p.x > 1 || p.y < 0 || p.y > 1) {
			// if moving to outside, then the lifetime of particle is done.
			break;
		}
		else {
			// location to be illuminated.
			int x_ind = int((p.x) * n_nodes);
			int y_ind = int((p.y) * n_nodes);
			intensity[x_ind][y_ind] +=1;

			// scattering.
			if (random() < albedo) {
				// scattering probability = mu_s / mu_t
				p.theta = 2.0 * M_PI * random();
			}
			else {
				// absorption, the life-cycle is done here.
				break;
			}
		}
	}
}

void MCRT::simulate() {
	double side = 1.0/n_nodes;
	size_t i;
	omp_set_num_threads(omp_get_num_procs());
#pragma omp parallel for private(i) schedule(static)
	for (i = 0; i < n_nodes; i++) {
		for (size_t j = 0; j < n_nodes; j++) {
			if (charge[i][j] != 0) {
				photon p(0.,0.,0.);
				for (size_t k = 0; k < n_photons; k++) {
					p.x = side/2.0 + i * side;
					p.y = side/2.0 + j * side;
					p.theta = random() * 2 * M_PI;
					scatter(p);
				}
			}
		}
	}

	double divider = double(n_photons) * mu_t;
#pragma omp parallel for private(i)
	for (i = 0; i < n_nodes; i++) {
		for (size_t j = 0; j < n_nodes; j++) {
			intensity[i][j] /= divider;
		}
	}
}

//void MCRT::write(char* filename) {
//	ofstream file;
//	file.open(filename);
//	for (size_t i = 0; i < n_nodes; i++) {
//		for (size_t j = 0; j < n_nodes; j++) {
//			file << intensity[i][j] << " ";
//		}
//		file <<"\n";
//	}
//	file.close();
//}


template class mexplus::Session<MCRT>;

namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 4);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<MCRT>::create(new MCRT(
			size_t(*mxGetPr(input.get(0))), size_t(*mxGetPr(input.get(1))),
			*mxGetPr(input.get(2)), *mxGetPr(input.get(3)))));

}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<MCRT>::destroy(input.get(0));
}

MEX_DEFINE(simulate)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 1);
	auto mcrt = Session<MCRT>::get(input.get(0));
	mcrt->set_charge();
	mcrt->simulate();
	plhs[0] = mxCreateNumericMatrix(mcrt->n_nodes * mcrt->n_nodes,1, mxDOUBLE_CLASS, mxREAL);
	double* ptr = mxGetPr(plhs[0]);
	for (size_t i = 0; i < mcrt->n_nodes; i++) {
		for (size_t j = 0; j < mcrt->n_nodes; j++) {
			*(ptr++) = mcrt->intensity[i][j];
		}
	}
}
}
MEX_DISPATCH

