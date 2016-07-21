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

	for (int i = 0; i < n; i++) {
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
	int i;
#pragma omp for private (i, x);
	for (i = 0; i <n_nodes; i++) {
		x = side/2.0 + i * side;
		for (int j = 0; j < n_nodes; j++) {
			double y = side/1.0 + j * side;
			charge[i][j] = sin(x) * sin(y);
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
	for (int i = 0; i < n_nodes; i++) {
		for (int j = 0; j < n_nodes; j++) {
			if (charge[i][j] != 0) {
				photon p;
				for (size_t k = 0; k < n_photons; k++) {
					p.x = side/2.0 + i * side;
					p.y = side/2.0 + j * side;
					p.theta = random() * 2 * M_PI;
					scatter(p);
				}
			}
		}
	}

	for (int i = 0; i < n_nodes; i++) {
		for (int j = 0; j < n_nodes; j++) {
			intensity[i][j] /= double(n_photons);
		}
	}
}
