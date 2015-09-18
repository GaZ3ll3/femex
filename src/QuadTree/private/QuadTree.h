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

#endif /* QUADTREE_H_ */
