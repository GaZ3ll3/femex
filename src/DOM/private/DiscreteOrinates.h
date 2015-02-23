/*
 * DiscreteOrinates.h
 *
 *  Created on: Feb 21, 2015
 *      Author: lurker
 */

#ifndef SRC_DOM_PRIVATE_DISCRETEORINATES_H_
#define SRC_DOM_PRIVATE_DISCRETEORINATES_H_

#include <cstdlib>
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

#define INTERSECT_DET(X1, Y1, X2, Y2, THETA) (cos((THETA)) * ((Y1) - (Y2))) - (sin((THETA)) * ((X1) - (X2)))

#define INTERSECT_CROSS(X1, Y1, X2, Y2, A, B) (((X1) - (X2)) * ((B) - (Y2))) -(((A) - (X2)) * ((Y1)- (Y2)))

namespace Core {

class DiscreteOrinates {
private:
	int32_t nAngle;
	Real_t  initAngle;
public:
	DiscreteOrinates(MatlabPtr _nAngle) noexcept;
	DiscreteOrinates(MatlabPtr _nAngle, MatlabPtr _initAngle) noexcept;
	virtual ~DiscreteOrinates();


	void RayInt(Real_t*& output, MatlabPtr nodes, MatlabPtr elems,
			MatlabPtr neighbors, MatlabPtr edges, MatlabPtr weights, MatlabPtr Fcn);



};

} /* namespace Core */

#endif /* SRC_DOM_PRIVATE_DISCRETEORINATES_H_ */
