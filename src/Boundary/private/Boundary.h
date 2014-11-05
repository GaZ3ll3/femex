/*
 * Boundary.h
 *
 *  Created on: Oct 30, 2014
 *      Author: lurker
 */

#ifndef BOUNDARY_PRIVATE_BOUNDARY_H_
#define BOUNDARY_PRIVATE_BOUNDARY_H_

#include <cstdlib>
#include <vector>
#include <cmath>

#include <iostream>
#include <iterator>
#include <string.h>

#include <unordered_set>

#include <mexplus.h>
#include <pprint.h>

#include "utils.h"

using namespace std;
using namespace mexplus;

namespace MEX {

enum Boundary_t {DIRICHLET = 0, NEUMANN, ROBIN};

class Boundary {
public:
	Boundary(MatlabPtr);
	virtual ~Boundary();

	vector<int32_t> b_edges;
	unordered_set<int32_t> b_edge_set;


	// apply boundary condition on LHS and RHS

};

class DirichletBC:public Boundary {
public:
	explicit DirichletBC(MatlabPtr _edges) : Boundary(_edges) {}
};


} /* namespace MEX */

#endif /* BOUNDARY_PRIVATE_BOUNDARY_H_ */
