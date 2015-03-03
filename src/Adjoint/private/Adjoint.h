/*
 * Adjoint.h
 *
 *  Created on: Feb 28, 2015
 *      Author: lurker
 *
 *      Disclaim:
 *
 *      The adjoint class is originally designed to use template.
 *     	For simplicity, I did not use this.
 *
 */

#ifndef SRC_ADJOINT_PRIVATE_ADJOINT_H_
#define SRC_ADJOINT_PRIVATE_ADJOINT_H_

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

namespace Core {

class Adjoint {
public:
	Adjoint(std::size_t param_size, std::size_t solution_size) noexcept;
	Adjoint() noexcept;
	virtual ~Adjoint();

	std::vector<Real_t> parameter;
	std::vector<Real_t> solution;
	std::vector<Real_t> adjoint_solution;

	std::vector<Real_t> gradient;
	Real_t              J;

	/*
	 * The process is:
	 *
	 * 1. initial value: m_0.
	 * 2. Update_Solution and Update_J.
	 * 3. Update_Adj and Update_Gradient
	 * 4. Linear search algorithm
	 * 	- 4.1 Along the gradient, Update_Param
	 * 	- 4.2 Update_Solution and Update_J
	 * 	- 4.3 repeat until J does not decrease.
	 * 5. repeat 2.
	 *
	 * ----------------------------------------
	 *
	 * The adjoint model typically returns
	 *
	 * [J, Gradient] = Adjoint(param);
	 *
	 * Which will call forward solver multiple
	 * times and Adjoint solver once.
	 *
	 * ----------------------------------------
	 */
	void Update_Gradient();

	void Update_Solution(MatlabPtr _solution);
	void Update_Param(MatlabPtr _parameter);

	/*
	 * this might introduce some serious overhead.
	 *
	 * Update_J is "not necessary".
	 *
	 * Update_Adj will call another solver, which will calculate solution
	 * and store the solution in adjoint_solution.
	 */
	void Update_J();
	void Update_Adj();

};

} /* namespace Core */

#endif /* SRC_ADJOINT_PRIVATE_ADJOINT_H_ */
