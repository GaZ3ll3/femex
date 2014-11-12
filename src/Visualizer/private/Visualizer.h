/*
 * Visualizer.h
 *
 *  Created on: Nov 4, 2014
 *      Author: lurker
 */

#ifndef VISUALIZER_PRIVATE_VISUALIZER_H_
#define VISUALIZER_PRIVATE_VISUALIZER_H_


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

namespace MEX {

/*
 * visualization of solution, figures can be suppressed
 *
 * and output solution's gradient
 */

class Visualizer {
public:
	Visualizer();
	Visualizer(MatlabPtr& i_data, MatlabPtr& i_domain);
	virtual ~Visualizer();


	/*
	 * assorted visualization methods
	 */
	/*
	 * surf figure
	 */
	void TriSurf();
	void TriMesh();
	/*
	 * Vector figure
	 */
	void Quiver();
	/*
	 * surf figure
	 */
	void QuiverX();
	void QuiverY();


private:


};

} /* namespace MEX */

#endif /* VISUALIZER_PRIVATE_VISUALIZER_H_ */
