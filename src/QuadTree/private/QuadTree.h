/*
 * QuadTree.h
 *
 *  Created on: Sep 18, 2015
 *      Author: lurker
 */

#ifndef QUADTREE_H_
#define QUADTREE_H_

#include <cstdlib>
#include <vector>
#include <cmath>

#include <iostream>
#include <iterator>
#include <string.h>

#include <pprint.h>

using namespace std;

using scalar_t = double;
using int_t    = int;

class Photon {
public:
	Photon(const vector<scalar_t>& _position_, scalar_t _intensity_) :
		position (_position_),
		intensity(_intensity_),
		id(-1){};
	virtual ~Photon() {};
	vector<scalar_t> position;
	scalar_t intensity;
	int_t id;
};

enum class Status {
	ROOT,
	BRANCH,
	LEAF,
	EMPTY,
	UNSET
};

class QuadTree {
public:
	QuadTree(const vector<scalar_t>& _position_, scalar_t _size_);
	virtual ~QuadTree();

	/*
	 * members
	 */
	vector<Photon*>                      photons;
	vector<QuadTree*>                    children;
	vector<scalar_t>                     position;

	QuadTree*                            parent;
	Status                               status;
	scalar_t                             size;
	scalar_t                             intensity;


	/*
	 * member functions
	 */
	void addPhoton(Photon*)              noexcept;
	void clearPhotons()                  noexcept;
	void populate()                      noexcept;
	void updateAttribute()               noexcept;
	/*
	 * getters, setters
	 */
	QuadTree*            getParent()     noexcept;
	size_t               getPhotonSize() noexcept;
	Status               getStatus()     noexcept;

	void setParent(QuadTree*)            noexcept;
	void setStatus(Status)               noexcept;

};


static std::vector<double> x {
	-0.90617984593866396370032134655048,
	 -0.53846931010568310771446931539685,
	                                   0,
	  0.53846931010568310771446931539685,
	  0.90617984593866396370032134655048

};
static std::vector<double> w {
	 0.23692688505618908489935847683228,
	 0.47862867049936647090291330641776,
	 0.56888888888888888888888888888889,
	 0.47862867049936647090291330641776,
	 0.23692688505618908489935847683228
};


inline double distance(std::vector<double>& x, std::vector<double>& y) noexcept {
	return sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]));
}

inline double eval(double sigma_t, std::vector<double>& x, std::vector<double>& y) {
	return exp(-sigma_t * distance(x, y))/distance(x, y);
}

inline double visit(double sigma_t, Photon* photon, QuadTree* cell) {

	auto s = cell->size;

	double sum = 0;

	std::vector<double> center {
		cell->position[0] + s/2.,
		cell->position[1] + s/2.
	};

	for (size_t i = 0, l = x.size(); i < l; i++) {
		center[0] += x[i] * s/2;
		for (size_t j = 0, k = x.size(); j < k; j++) {
			center[1] += x[j] * s/2;
			sum += eval(sigma_t, photon->position, center) * w[i] * w[j]/4.0;
			center[1] -= x[j] * s/2;
		}
		center[0] -= x[i] * s/2;
	}
	return sum * s * s;
}

inline void trasverse(double sigma_t, double theta, Photon* photon, QuadTree* cell, double* ptr, size_t n) {
	double d = distance(photon->position, cell->position);
	if (cell != nullptr && (cell->getStatus()== Status::LEAF ||
			cell->size/d< theta) ) {
		auto result = visit(sigma_t , photon, cell);
		for (auto _photon : cell->photons) {
			if (photon->id != _photon->id) {
				ptr[photon->id * n + _photon->id] = result/cell->getPhotonSize();
			}
			else {
				ptr[photon->id * n + _photon->id] = 2* M_PI * (1 - exp(-sigma_t * cell->getPhotonSize()))/(sigma_t);
			}
		}
	}
	else{
		for(auto child : cell->children) {
			if (child->getStatus() != Status::EMPTY) {
				trasverse(sigma_t, theta, photon, child, ptr, n);
			}
		}
	}
}
#endif /* QUADTREE_H_ */
