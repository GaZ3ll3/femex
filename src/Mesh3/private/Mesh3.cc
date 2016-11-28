/*
 * Mesh3.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: lurker
 */

#include "Mesh3.h"

Mesh3::Mesh3(MatlabPtr _Vertices, MatlabPtr _PML, MatlabPtr _Vol) noexcept {
	min_vol            = *Matlab_Cast<Real_t>(_Vol);
	auto _VerticesPtr  = Matlab_Cast<Real_t>(_Vertices);
	auto _VerticesSize = mxGetN(_Vertices);


	auto _PMLPtr      = Matlab_Cast<Real_t>(_PML);
	auto _PMLSize     = mxGetN(_PML);

	tetgenio input, mid;

	input.numberofpoints = _VerticesSize;
	input.pointlist = new Real_t[_VerticesSize * 3];

	/*
	 * fill pointlist, _Vertices is column primary.
	 */
	for (int i = 0; i < _VerticesSize; i++) {
		input.pointlist[3 * i    ] = _VerticesPtr[3 * i];
		input.pointlist[3 * i + 1] = _VerticesPtr[3 * i + 1];
		input.pointlist[3 * i + 2] = _VerticesPtr[3 * i + 2];
	}

	tetrahedralize("pzcQ", &input, &mid, (tetgenio*)nullptr, (tetgenio*)nullptr);
	tetrahedralize(
			const_cast<char*>(("prq1.20/18a" + std::to_string(min_vol) + "QzenBC").c_str()),\
			&mid, &_meshdata,(tetgenio*)nullptr, (tetgenio*)nullptr);

	return;
}

Mesh3::~Mesh3() noexcept {
#if DEBUG
	cout << "3D MESH detached." <<"\n";
#endif
}

void Mesh3::Info() noexcept{
  std::cout << "Delaunay information:" << std::endl;
  std::cout << "number of points: " << _meshdata.numberofpoints << std::endl;
  std::cout << "number of edges: " << _meshdata.numberofedges << std::endl;
  std::cout << "number of trifaces: " << _meshdata.numberoftrifaces << std::endl;
  std::cout << "number of tetrahedra: " << _meshdata.numberoftetrahedra << std::endl;
  std::cout << "number of corners: " << _meshdata.numberofcorners << std::endl;
}

void Mesh3::Promote(int32_t _deg) noexcept{
	/*
	 * add quadrature points
	 */
}


template class mexplus::Session<Mesh3>;

namespace {

MEX_DEFINE(new) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Mesh3>::create(new Mesh3(
			const_cast<mxArray*>(input.get(0)),
			const_cast<mxArray*>(input.get(1)),
			const_cast<mxArray*>(input.get(2))
	)));
}

// Delete the instance by id
MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<Mesh3>::destroy(input.get(0));
}

MEX_DEFINE(report)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	auto mesh3 = Session<Mesh3>::get(input.get(0));
	mesh3->Info();
}

MEX_DEFINE(meshdata)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 2);
	auto mesh3 = Session<Mesh3>::get(input.get(0));
	/*
	 * todo: export information
	 */
	plhs[0] = mxCreateNumericMatrix(3, mesh3->_meshdata.numberofpoints, mxDOUBLE_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[0]), mesh3->_meshdata.pointlist, 3 * mesh3->_meshdata.numberofpoints*sizeof(Real_t));
	plhs[1] = mxCreateNumericMatrix(4, mesh3->_meshdata.numberoftetrahedra, mxINT32_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[1]), mesh3->_meshdata.tetrahedronlist, 4 * mesh3->_meshdata.numberoftetrahedra*sizeof(int32_t));
	MatlabPtr temp[] = {plhs[1], mxCreateDoubleScalar(1.0)};
	mexCallMATLAB(1, &plhs[1], 2, temp, "plus");
}

}

MEX_DISPATCH
