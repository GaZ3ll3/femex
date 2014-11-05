/*
 * Mesh.h
 *
 *  Created on: Oct 17, 2014
 *      Author: lurker
 */

#ifndef MESH_PRIVATE_MESH_H_
#define MESH_PRIVATE_MESH_H_

#include "utils.h"

#include <iostream>
#include <iterator>
#include <string.h>

#include <unordered_map>

#include <mexplus.h>
#include <pprint.h>

extern "C"
{
#include "triangle/triangle.h"
}


using namespace std;
using namespace mexplus;

namespace MEX {

typedef	struct Topology {
		vector<Real_t> nodes;
		std::unordered_map<std::string, vector<int32_t>> edges;
		vector<int32_t>    elems;
		vector<int32_t>    boundary;

		void clear(){
			nodes.clear(); edges.clear();
			elems.clear(); boundary.clear();
		}

		void resize(int32_t _nsize,int32_t _esize,int32_t _bsize){
			nodes.resize(_nsize);
			elems.resize(_esize);
			boundary.resize(_bsize);
		}
} Topology;

typedef Topology* TopologyPtr;

class Mesh {
public:
	/*
	 *  Mesh(MatlabPtr, MatlabPtr)
	 *
	 *  Construct new Mesh(triangle) object from data.
	 */
	Mesh(MatlabPtr, MatlabPtr) noexcept;
	/*
	 *  virtual ~Mesh()
	 *
	 *  Destroy Mesh object.
	 */
	virtual ~Mesh();

	/*
	 *  Real_t min_area
	 *
	 *  Stores minimal triangle area of Mesh object.
	 */
	Real_t min_area;

	/*
	 * Pointer to topology
	 *
	 * stores all necessary information
	 */
	Topology topology;

	/*
	 * Refine Mesh object by shrinking minimal area
	 */
	void Refine() noexcept;

	/*
	 * Promote Mesh to higher ordered Lagrange Mesh,
	 * Base linear Mesh won't change.
	 */
	void Promote(int32_t) noexcept;

	/*
	 *  Report all information
	 */

	void Info() noexcept;
	int32_t Info(const std::string&) noexcept;

private:
	struct triangulateio _meshdata;

	inline Real_t _Point_X(int32_t _index){
		return _meshdata.pointlist[2*_index ];
	}

	inline Real_t _Point_Y(int32_t _index){
		return _meshdata.pointlist[2*_index + 1];
	}

	inline int32_t _Edge_L(int32_t _index){
		return _meshdata.edgelist[2*_index];
	}

	inline int32_t _Edge_R(int32_t _index){
		return _meshdata.edgelist[2* _index + 1];
	}

	inline int32_t _Tri_U(int32_t _index){
		return _meshdata.trianglelist[3 * _index];
	}

	inline int32_t _Tri_V(int32_t _index){
		return _meshdata.trianglelist[3 * _index + 1];
	}

	inline int32_t _Tri_W(int32_t _index){
		return _meshdata.trianglelist[3 * _index + 2];
	}

	inline int32_t _Seg_L(int32_t _index){
		return _meshdata.segmentlist[2 * _index ];
	}

	inline int32_t _Seg_R(int32_t _index){
		return _meshdata.segmentlist[2 * _index + 1];
	}



	void clear() noexcept;



};

} /* namespace MEX */

#endif /* MESH_PRIVATE_MESH_H_ */
