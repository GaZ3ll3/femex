/*
 * MCRT.h
 *
 * Monte Carlo for 2D radiative transport equation in [0, 1]^2.
 *
 *  Created on: Jul 20, 2016
 *      Author: lurker
 */

#ifndef SRC_MCRT_PRIVATE_MCRT_H_
#define SRC_MCRT_PRIVATE_MCRT_H_

#include <mexplus.h>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include "utils.h"

using namespace mexplus;
using namespace std;

typedef struct photon {
	double x;
	double y;
	double theta;
	photon(double x_, double y_, double theta_) {
		x = x_;
		y = y_;
		theta = theta_;
	}
} photon;

class MCRT {
public:
	MCRT(size_t n_photons, size_t n, double mu_s, double mu_t);
	virtual ~MCRT();
	vector<vector<double> > intensity;
	vector<vector<double> > charge;
	size_t n_nodes;
	size_t n_photons;
	double albedo;
	double mu_t;

	void set_charge();
	double random();
	void simulate();
	void scatter(photon& p);
	void write(char* filename);
};

#endif /* SRC_MCRT_PRIVATE_MCRT_H_ */

