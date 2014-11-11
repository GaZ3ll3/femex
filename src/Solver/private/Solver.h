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

enum class Solver_type : int32_t {iterative = 0, direct};
enum class Smoother : int32_t {jac = 0, gsf, gsb};
enum class AMG : int32_t {none = 0, ilu = 1, amli, amg};

typedef struct Solver_option {
	Solver_type _type;
	int32_t     _maxIter;
	int32_t     _nrestart;
	Real_t      _restol;
	bool        _multilevel;
	Smoother    _smoother;
	AMG         _amg;
} Solver_option;

class Solver {
public:
	Solver();
	virtual ~Solver();
	Solver_option option;

	/*
	 * public methods
	 */
	void set_type(int32_t);
	void set_maxIter(int32_t);
	void set_nrestart(int32_t);
	void set_restol(Real_t);
	void set_amg(int32_t);
	void set_smoother(int32_t);

	int32_t get_type();
	int32_t get_maxIter();
	int32_t get_nrestart();
	Real_t  get_restol();
	int32_t get_amg();
	int32_t get_smoother();
};

} /* namespace MEX */

#endif /* SRC_SOLVER_PRIVATE_SOLVER_H_ */
