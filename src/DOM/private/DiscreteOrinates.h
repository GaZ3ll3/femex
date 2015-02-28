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
#include <queue>
#include <cmath>

#include <iostream>
#include <iterator>
#include <string.h>

#include <mexplus.h>
#include <pprint.h>

#include "utils.h"

using namespace std;
using namespace mexplus;

#define SCALE 1.0L


typedef struct Raylet {
	int32_t elem;
	Real_t first[2];
	Real_t second[2];
} Raylet;

namespace Core {

class DiscreteOrinates {
private:
	int32_t nAngle;
	Real_t  initAngle;
public:
	DiscreteOrinates(MatlabPtr _nAngle) noexcept;
	DiscreteOrinates(MatlabPtr _nAngle, MatlabPtr _initAngle) noexcept;
	virtual ~DiscreteOrinates();


	/*
	 *  variables, caution, this might be slow.
	 */

	std::vector<std::vector<std::vector<Raylet>>> Ray;

	std::vector<std::vector<Real_t>> Output;
//	std::vector<std::vector<Real_t>> Source;

	std::vector<Real_t> RHS;
	std::vector<Real_t> Source;
	std::vector<Real_t> Average;

	std::vector<Real_t> Sigma_t;
	std::vector<Real_t> Sigma_s;

	/*
	 *  methods
	 */

	void RayInt(MatlabPtr nodes, MatlabPtr elems,
			MatlabPtr neighbors);

	void RayTrace(std::vector<Real_t>& tmp, bool& intersect, Real_t& q_t, Real_t& q_eta,
			Real_t& q_x1, Real_t& q_y1, Real_t& q_x2,
			Real_t& q_y2, Real_t& q_x3, Real_t& q_y3,
			Real_t& a, Real_t& b, Real_t& theta);

	void RayTrim(std::vector<Real_t>& tmp, Real_t &a, Real_t &b);

	void RayShow();


	/*
	 * solve the transport equation.
	 */
	void SourceIteration_init();
	void SourceIteration_port(MatlabPtr Fcn, MatlabPtr Sigma_t_Fcn, MatlabPtr Sigma_s_Fcn);
	void SourceIteration_iter(MatlabPtr nodes, MatlabPtr elems);
	void SourceIteration_accl(MatlabPtr delta);
	void SourceIteration_set(MatlabPtr delta);

};

} /* namespace Core */

#endif /* SRC_DOM_PRIVATE_DISCRETEORINATES_H_ */
