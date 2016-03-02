/*
 * Mesh3.h
 *
 *  Created on: Mar 1, 2016
 *      Author: lurker
 */

#ifndef SRC_MESH3_PRIVATE_MESH3_H_
#define SRC_MESH3_PRIVATE_MESH3_H_



#include "utils.h"

#include <iostream>
#include <iterator>
#include <string.h>

#include <unordered_map>

#include <mexplus.h>
#include <pprint.h>


#include "tetgen/tetgen.h"


using namespace std;
using namespace mexplus;


class Mesh3 {
public:
	Mesh3(MatlabPtr _Vertices,\
			MatlabPtr _PML, \
			MatlabPtr _Vol)noexcept;
	virtual ~Mesh3() noexcept;

	void Info() noexcept;

	void Promote(int32_t _deg) noexcept;

	Real_t min_vol;
	tetgenio _meshdata;

};

#endif /* SRC_MESH3_PRIVATE_MESH3_H_ */
