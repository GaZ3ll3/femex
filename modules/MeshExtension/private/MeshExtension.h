/*
 * MeshExtension.h
 *
 *  Created on: Feb 7, 2015
 *      Author: lurker
 */

#ifndef MODULES_MESHEXTENSION_PRIVATE_MESHEXTENSION_H_
#define MODULES_MESHEXTENSION_PRIVATE_MESHEXTENSION_H_

#include <cstdlib>
#include <vector>
#include <cmath>

#include <iostream>
#include <iterator>
#include <string.h>

#include <unordered_set>

#include <mexplus.h>
#include <pprint.h>
#include "exprtk.hpp"

#include "utils.h"

using namespace std;
using namespace mexplus;


namespace Extension {

class MeshExtension {
public:
	MeshExtension();
	virtual ~MeshExtension();
	/*
	 * extract indices for the thin layer
	 */
	std::vector<size_t> layer;

	void LayerElementsIndex(MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr Interior);
};

} /* namespace Extension */

#endif /* MODULES_MESHEXTENSION_PRIVATE_MESHEXTENSION_H_ */
