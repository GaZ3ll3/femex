/*
 * Mesh.cpp
 *
 *  Created on: Oct 17, 2014
 *      Author: lurker
 */

#include "Mesh.h"

namespace MEX {

Mesh::Mesh(MatlabPtr _Boundary, MatlabPtr _Area) noexcept{
	min_area = *Matlab_Cast<Real_t>(_Area);
	auto _BoundaryPtr = Matlab_Cast<Real_t>(_Boundary);
	int32_t _Boundary_L = Msize(_Boundary)/2;


	/*
	 * Two step meshing
	 */
	struct triangulateio  input, mid;
	input.numberofpoints = _Boundary_L;
	input.numberofpointattributes = 0;
	input.pointlist = (Real_t *) malloc(_Boundary_L * 2 * sizeof(Real_t));

	for (int32_t i = 0; i < _Boundary_L; i++)
	{
		input.pointlist[2*i] = *(_BoundaryPtr++);
		input.pointlist[2*i + 1] = *(_BoundaryPtr++);
	}

	input.pointmarkerlist = (int *)nullptr;
	input.numberofsegments = 0;
	input.numberofholes = 0;
	input.numberofregions = 0;
	input.regionlist = (Real_t *)nullptr;

	mid.pointlist = (Real_t *) nullptr;
	mid.pointmarkerlist = (int *) nullptr;
	mid.trianglelist = (int *) nullptr;
	mid.triangleattributelist = (Real_t *) nullptr;
	mid.neighborlist = (int *) nullptr;
	mid.segmentlist = (int *) nullptr;
	mid.segmentmarkerlist = (int *) nullptr;
	mid.edgelist = (int *) nullptr;
	mid.edgemarkerlist = (int *) nullptr;

	/*
	 * Coarse meshing
	 */
	triangulate("pcz", &input, &mid, (struct triangulateio *) nullptr);
	/*
	 * Fine meshing
	 */
	mid.trianglearealist = nullptr;

	_meshdata.pointlist = (Real_t *) nullptr;
	_meshdata.pointattributelist = (Real_t *) nullptr;
	_meshdata.trianglelist = (int *) nullptr;
	_meshdata.triangleattributelist = (Real_t *) nullptr;
	_meshdata.edgelist = (int *) nullptr;
	_meshdata.edgemarkerlist = (int *) nullptr;
	_meshdata.segmentlist = (int *) nullptr;

	triangulate(const_cast<char*>(("prq30.0a" + std::to_string(min_area) + "ezBC").c_str()), &mid, &_meshdata, (struct triangulateio *) nullptr);

	free(input.pointlist);
	free(input.pointmarkerlist);
	free(input.regionlist);
	free(mid.pointlist);
	free(mid.pointmarkerlist);
	free(mid.trianglelist);
	free(mid.triangleattributelist);
	free(mid.trianglearealist);
	free(mid.neighborlist);
	free(mid.segmentlist);
	free(mid.segmentmarkerlist);
	free(mid.edgelist);
	free(mid.edgemarkerlist);
}

Mesh::~Mesh() noexcept{
	clear();
#ifdef DEBUG
	mexPrintf("Mesh detached\n");
#endif
}


void Mesh::Refine() noexcept{
/*
 * mid point refinement:
 *
 * each edges will produce a new node to insert
 */

}

void Mesh::Promote(int32_t _deg) noexcept{

	auto deg = _deg;
	if (deg <= 0) {mexWarnMsgTxt("Degree has to be positive. Please use first order instead.\n");}
	if (deg >= 20) {mexPrintf("Degree Too Large\n");}
	else{
		topology.clear();

		Real_t P1_Coord_X, P1_Coord_Y, P2_Coord_X,P2_Coord_Y;
		Real_t Delta_X, Delta_Y;

		int32_t id = (deg + 1) * (deg + 2)/2;

		topology.resize(
				2*( _meshdata.numberofpoints +
				(deg - 1)*_meshdata.numberofedges +
				(id - 3*deg)*_meshdata.numberoftriangles),
				((deg + 1) * (deg + 2)/2) * _meshdata.numberoftriangles,
				(deg + 1) * _meshdata.numberofsegments);

		/*
		 * Fill nodes
		 */
		for (int32_t i = 0; i < _meshdata.numberofpoints; i++){
			topology.nodes[2*i    ] = _Point_X(i);
			topology.nodes[2*i + 1] = _Point_Y(i);
		}

		for (int32_t i = 0; i < _meshdata.numberofedges; i++){
			P1_Coord_X = _Point_X(_Edge_L(i));
			P1_Coord_Y = _Point_Y(_Edge_L(i));
			P2_Coord_X = _Point_X(_Edge_R(i));
			P2_Coord_Y = _Point_Y(_Edge_R(i));

			Delta_X = (P2_Coord_X - P1_Coord_X)/(Real_t)deg;
			Delta_Y = (P2_Coord_Y - P1_Coord_Y)/(Real_t)deg;

			vector<int32_t> EdgePoints(deg - 1);

			for (int32_t j = 0; j < deg - 1; j++) {
				EdgePoints[j] = _meshdata.numberofpoints  + i*(deg - 1)  + j;

				topology.nodes[2 * EdgePoints[j]    ] = P1_Coord_X + Delta_X * (j + 1);
				topology.nodes[2 * EdgePoints[j] + 1] = P1_Coord_Y + Delta_Y * (j + 1);
			}

			topology.edges.insert(make_pair(
					std::to_string(_Edge_L(i)) + "-" +
					std::to_string(_Edge_R(i)), EdgePoints));

		}

		int32_t counter = 2*_meshdata.numberofpoints + 2*(deg - 1)*_meshdata.numberofedges;

		std::unordered_map<std::string, std::vector<int32_t>>::const_iterator EdgeIterator;


		for (int32_t index = 0; index < _meshdata.numberoftriangles; index++){
			topology.elems[(deg + 1)*(deg + 2)/ 2 * index ] = _Tri_U(index);
			topology.elems[(deg + 1)*(deg + 2)/ 2 * index + 1 ] = _Tri_V(index);
			topology.elems[(deg + 1)*(deg + 2)/ 2 * index + 2 ] = _Tri_W(index);

			auto tri_u = std::to_string(topology.elems[(deg + 1)*(deg + 2)/ 2 * index ]);
			auto tri_v = std::to_string(topology.elems[(deg + 1)*(deg + 2)/ 2 * index + 1]);
			auto tri_w = std::to_string(topology.elems[(deg + 1)*(deg + 2)/ 2 * index + 2]);


			/*
			 * edge 0 - 1
			 */
			EdgeIterator = topology.edges.find(tri_u + "-" + tri_v);

			if (EdgeIterator != topology.edges.end()){
				std::copy(EdgeIterator->second.begin(), EdgeIterator->second.end(), topology.elems.begin() +
						(deg + 1)*(deg + 2)/ 2 * index + 3);
			}

			if (EdgeIterator == topology.edges.end()){
				EdgeIterator = topology.edges.find(tri_v + "-" + tri_u);
				std::reverse_copy(EdgeIterator->second.begin(), EdgeIterator->second.end(), topology.elems.begin() +
						(deg + 1)*(deg + 2)/ 2 * index + 3);
			}

			/*
			 * edge 1 - 2
			 */
			EdgeIterator = topology.edges.find(tri_v + "-" + tri_w);

			if (EdgeIterator != topology.edges.end()){
				std::copy(EdgeIterator->second.begin(), EdgeIterator->second.end(), topology.elems.begin() +
						(deg + 1)*(deg + 2)/ 2 * index + 3 +  (deg - 1));
			}

			if (EdgeIterator == topology.edges.end()){
				EdgeIterator = topology.edges.find(tri_w + "-" + tri_v);
				std::reverse_copy(EdgeIterator->second.begin(), EdgeIterator->second.end(), topology.elems.begin() +
						(deg + 1)*(deg + 2)/ 2 * index + 3 +  (deg - 1));
			}

			/*
			 * edge 2 - 0
			 */
			EdgeIterator = topology.edges.find(tri_w + "-" + tri_u);

			if (EdgeIterator != topology.edges.end()){
				std::copy(EdgeIterator->second.begin(), EdgeIterator->second.end(), topology.elems.begin() +
						(deg + 1)*(deg + 2)/ 2 * index + 3 + 2*((deg - 1)));
			}

			if (EdgeIterator == topology.edges.end()){
				EdgeIterator = topology.edges.find(tri_u + "-" + tri_w);
				std::reverse_copy(EdgeIterator->second.begin(), EdgeIterator->second.end(), topology.elems.begin() +
						(deg + 1)*(deg + 2)/ 2 * index + 3 + 2*((deg - 1)));
			}



			int32_t internal_counter = 0;
			for (int32_t i  = 1; i < deg - 1; i++){
				for (int32_t j = 1; j < deg - i; j++){
					int32_t k = deg - i - j;

					topology.elems[(deg + 1)*(deg + 2)/ 2 * index + 3*deg + internal_counter] = counter/2;
					internal_counter += 1;

					topology.nodes[counter] =
							(double)k/(double)deg * _Point_X(_Tri_U(index)) +
							(double)i/(double)deg * _Point_X(_Tri_V(index)) +
							(double)j/(double)deg * _Point_X(_Tri_W(index));

					topology.nodes[counter + 1] =
							(double)k/(double)deg * _Point_Y(_Tri_U(index)) +
							(double)i/(double)deg * _Point_Y(_Tri_V(index)) +
							(double)j/(double)deg * _Point_Y(_Tri_W(index));
					counter += 2;
				}
			}
		} // End of loop over elems



		for (int32_t index = 0; index < _meshdata.numberofsegments; index++){
			/*
			 * Orientation reversed
			 */
			topology.boundary[(deg + 1)*index + 1] = _Seg_L(index);
			topology.boundary[(deg + 1)*index    ] = _Seg_R(index);

			auto bound_l = std::to_string(topology.boundary[(deg + 1)*index    ]);
			auto bound_r = std::to_string(topology.boundary[(deg + 1)*index + 1]);
			/*
			 * In order insert the interpolated nodes.
			 */
			EdgeIterator = topology.edges.find( bound_l + "-" + bound_r);

			if (EdgeIterator != topology.edges.end()){
				std::copy(EdgeIterator->second.begin(), EdgeIterator->second.end(), topology.boundary.begin() + (deg + 1)*index + 2);
			}

			if (EdgeIterator == topology.edges.end()){
				EdgeIterator = topology.edges.find(bound_r + "_" + bound_l);
				std::reverse_copy(EdgeIterator->second.begin(), EdgeIterator->second.end(), topology.boundary.begin() + (deg + 1)*index + 2);
			}
		}// End of for

	}
}

void Mesh::Info() noexcept{
	mexPrintf("Mesh Generated:\n");
	mexPrintf("\tPoints  : %15d\n", _meshdata.numberofpoints);
	mexPrintf("\tElements: %15d\n", _meshdata.numberoftriangles);
	mexPrintf("\tEdges   : %15d\n", _meshdata.numberofedges);
	mexPrintf("\tSegments: %15d\n", _meshdata.numberofsegments);
	mexPrintf("\tRegions : %15d\n", _meshdata.numberofregions);
	mexPrintf("\tHoles   : %15d\n", _meshdata.numberofholes);
	mexPrintf("\tCorners : %15d\n", _meshdata.numberofcorners);
}


int32_t Mesh::Info(const std::string& _string) noexcept{
	if (!_string.compare("nodes")){
		return _meshdata.numberofpoints;
	}
	else if (!_string.compare("segs")) {
		return _meshdata.numberofsegments;
	}

	else if (!_string.compare("elems")) {
		return _meshdata.numberoftriangles;
	}
	else {
		return 0;
	}
}

void Mesh::clear() noexcept{
	/*
	 *  Should not free same pointer twice.
	 */
	if (_meshdata.pointlist != nullptr) {free(_meshdata.pointlist); _meshdata.pointlist = nullptr;}
	if (_meshdata.pointattributelist != nullptr) {free(_meshdata.pointattributelist); _meshdata.pointattributelist = nullptr;}
	if (_meshdata.trianglelist != nullptr) {free(_meshdata.trianglelist); _meshdata.trianglelist = nullptr;}
	if (_meshdata.triangleattributelist != nullptr) {free(_meshdata.triangleattributelist); _meshdata.triangleattributelist = nullptr;}
	if (_meshdata.segmentlist != nullptr) {free(_meshdata.segmentlist); _meshdata.segmentlist = nullptr;}
	if (_meshdata.edgemarkerlist != nullptr) {free(_meshdata.edgemarkerlist); _meshdata.edgemarkerlist = nullptr;}
	if (_meshdata.edgelist != nullptr) {free(_meshdata.edgelist); _meshdata.edgelist = nullptr;}
} // clear


} /* namespace MEX */

using namespace MEX;

template class mexplus::Session<Mesh>;

namespace {

// Create a new instance of Mesh and return its session id.
MEX_DEFINE(new) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Mesh>::create(new Mesh(const_cast<mxArray*>(input.get(0)),
			const_cast<mxArray*>(input.get(1)))));
}

// Delete the instance by id
MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<Mesh>::destroy(input.get(0));
}

MEX_DEFINE(report)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	auto mesh = Session<Mesh>::get(input.get(0));

	bool verbose = false;

	mesh->Info();
}


MEX_DEFINE(promote)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 3);
	auto mesh = Session<Mesh>::get(input.get(0));
	auto deg  = input.get<int>(1);

//	mexPrintf("%d\n", deg);
	mesh->Promote(deg);

	plhs[0] = mxCreateNumericMatrix(2, mesh->topology.nodes.size()/2, mxDOUBLE_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[0]), &mesh->topology.nodes[0], mesh->topology.nodes.size()*sizeof(Real_t));
	plhs[1] = mxCreateNumericMatrix((deg + 1)*(deg + 2)/2, mesh->topology.elems.size()/((deg + 1)*(deg + 2)/2), mxINT32_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[1]), &mesh->topology.elems[0], mesh->topology.elems.size()*sizeof(int));
	MatlabPtr temp_1[] = {plhs[1],mxCreateDoubleScalar(1.0)};
	mexCallMATLAB(1, &plhs[1], 2, temp_1, "plus");
	plhs[2] = mxCreateNumericMatrix((deg + 1), mesh->topology.boundary.size()/((deg + 1)), mxINT32_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[2]), &mesh->topology.boundary[0], mesh->topology.boundary.size()*sizeof(int));
	MatlabPtr temp_2[] = {plhs[2],mxCreateDoubleScalar(1.0)};
	mexCallMATLAB(1, &plhs[2], 2, temp_2, "plus");
}


MEX_DEFINE(export)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 3);
	auto mesh = Session<Mesh>::get(input.get(0));
	plhs[0] = MxArray::from(mesh->Info("nodes"));
	plhs[1] = MxArray::from(mesh->Info("elems"));
	plhs[2] = MxArray::from(mesh->Info("segs"));
}
} // namespace


MEX_DISPATCH




