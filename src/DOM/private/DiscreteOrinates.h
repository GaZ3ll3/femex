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
			MatlabPtr neighbors, MatlabPtr weights, MatlabPtr Fcn);



};

} /* namespace Core */

#endif /* SRC_DOM_PRIVATE_DISCRETEORINATES_H_ */
