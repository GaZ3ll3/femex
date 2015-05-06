/*
 * DG.h
 *
 *  Created on: May 3, 2015
 *      Author: lurker
 */

#ifndef SRC_DG_PRIVATE_DG_H_
#define SRC_DG_PRIVATE_DG_H_

#include <cstdlib>
#include <vector>
#include <queue>
#include <cmath>

#include <iostream>
#include <iterator>
#include <string.h>

#include <mexplus.h>
#include <pprint.h>

#include "utils.h"

/*
 * todo:
 *
 * 1. numerical flux options.
 * 2. given flux, assemble interface.
 * 3. average operator, and jump operator on edges.
 * 4. limiter.
 * 5. neighbor structure, from mesh.cc, each triangle has 3 neighbor triangles.
 * 	  need to adjust order to match interface.
 */

class DG {
public:
	DG();
	virtual ~DG();
};

#endif /* SRC_DG_PRIVATE_DG_H_ */
