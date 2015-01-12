/*
 * Solver.h
 *
 *  Created on: Nov 5, 2014
 *      Author: lurker
 */

#ifndef SRC_SOLVER_PRIVATE_SOLVER_H_
#define SRC_SOLVER_PRIVATE_SOLVER_H_

#include <cstdlib>
#include <cstdint>
#include <vector>
#include <cmath>

#include <iostream>
#include <iterator>
#include <string.h>

#include <mexplus.h>
#include <pprint.h>

#include "utils.h"

using namespace std;
using namespace mexplus;

namespace MEX {

class Solver {
public:
	Solver();
	virtual ~Solver();
	/*
	 * public methods
	 */
	// return reference
	void Reference(MatlabPtr &DX, MatlabPtr &DY,
			MatlabPtr Points);
	// return gradient on each element
	void Gradient(MatlabPtr& GradX, MatlabPtr& GradY,MatlabPtr Solution, MatlabPtr Nodes, MatlabPtr Elems,
			MatlabPtr DX, MatlabPtr DY);
	// return Neumann data which matches the boundaries
	void Neumann(MatlabPtr& Neumann, MatlabPtr GradX, MatlabPtr GradY, MatlabPtr Nodes, MatlabPtr Elems,
			MatlabPtr Boundaries, MatlabPtr BoundaryId);
};

} /* namespace MEX */

#endif /* SRC_SOLVER_PRIVATE_SOLVER_H_ */
