/*
 * Mesh.cpp
 *
 *  Created on: Oct 17, 2014
 *      Author: lurker
 */

#include "Mesh.h"

namespace MEX {

Mesh::Mesh(MatlabPtr _Boundary, MatlabPtr _PML, MatlabPtr _Area) noexcept{

	min_area            = *Matlab_Cast<Real_t>(_Area);
	auto _BoundaryPtr   = Matlab_Cast<Real_t>(_Boundary);
	auto _BoundarySize  = mxGetM(_Boundary)/2;

	auto _PMLPtr        = Matlab_Cast<Real_t>(_PML);
	auto _PMLSize       = mxGetM(_PML)/2;

	/*
	 * All segments
	 */
	int32_t _Boundary_L = (_BoundarySize  + _PMLSize);
	/*
	 * Two step meshing
	 */
	struct triangulateio  input, mid;
	input.numberofpoints = _Boundary_L;
	input.numberofpointattributes = 0;
	input.pointlist = (Real_t *) malloc(_Boundary_L * 2 * sizeof(Real_t));
	input.numberofsegments = _Boundary_L;
	input.segmentlist = (int32_t *) malloc(_Boundary_L * 2 * sizeof(int32_t));
	input.segmentmarkerlist = (int *) nullptr;


	for (int32_t i = 0; i < _BoundarySize; i++)
	{
		/*
		 * closed boundary input
		 */
		input.pointlist[2 * i] = _BoundaryPtr[2 * i];
		input.pointlist[2 * i + 1] = _BoundaryPtr[2 * i + 1];
		input.segmentlist[2 * i] = i ;
		input.segmentlist[2 * i + 1] = i + 1;
	}

	// getting a loop for closed boundary
	input.segmentlist[2 * _BoundarySize - 1] = 0;


	if (_PMLSize > 0){

		for (int32_t i = 0; i < _PMLSize; i++) {
			input.pointlist[2 * i + 2 * _BoundarySize] = _PMLPtr[2 * i];
			input.pointlist[2 * i + 2 * _BoundarySize +  1] = _PMLPtr[2 * i + 1];
			input.segmentlist[2 * i + 2 * _BoundarySize] = _BoundarySize + i;
			input.segmentlist[2 * i + 2 * _BoundarySize + 1] = _BoundarySize + i + 1;
		}

		input.segmentlist[2 * _Boundary_L - 1] = _BoundarySize;

	}
	// end of segments

	input.pointmarkerlist = (int *)nullptr;
	input.numberofholes = 0;
	input.numberofregions = 0;
	input.regionlist = (Real_t *)nullptr;

	mid.pointlist = (Real_t *) nullptr;
	mid.pointmarkerlist = (int *) nullptr;
	mid.trianglelist = (int *) nullptr;
	mid.triangleattributelist = (Real_t *) nullptr;
//	mid.neighborlist = (int *) nullptr;
	mid.segmentlist = (int *) nullptr;
	mid.segmentmarkerlist = (int *) nullptr;
	mid.edgelist = (int *) nullptr;
	mid.edgemarkerlist = (int *) nullptr;

	/*
	 * Coarse meshing
	 */
	triangulate("pczQ", &input, &mid, (struct triangulateio *) nullptr);
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
	_meshdata.neighborlist = (int *) nullptr;

	/*
	 * minimum angle chosen as 34.0 degrees now. should be adjusted.
	 */
	triangulate(const_cast<char*>(("prqQ34.0a" + std::to_string(min_area) + "eznBC").c_str()), &mid, &_meshdata, (struct triangulateio *) nullptr);

	free(input.pointlist);
	free(input.pointmarkerlist);
	free(input.regionlist);
	free(input.segmentlist);
	free(mid.pointlist);
	free(mid.pointmarkerlist);
	free(mid.trianglelist);
	free(mid.triangleattributelist);
	free(mid.trianglearealist);
//	free(mid.neighborlist);
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
				(deg + 1) * _meshdata.numberofsegments,
				_meshdata.numberofedges);


		/*
		 * Fill neighbors
		 */

		topology.neighbors.resize(3 * _meshdata.numberoftriangles);
		for (int32_t i = 0; i < _meshdata.numberoftriangles; i++) {
			topology.neighbors[3 * i] = _meshdata.neighborlist[3 * i] ;
			topology.neighbors[3 * i + 1] = _meshdata.neighborlist[3 * i + 1];
			topology.neighbors[3 * i + 2] = _meshdata.neighborlist[3 * i + 2];
		}

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

		// make a hash-set for boundary element
		std::unordered_set<std::string> _boundary_set;
		_boundary_set.reserve(_meshdata.numberofsegments);
		// temporary for storage
		std::unordered_map<std::string, int32_t>  boundary_elems;
		// setup _boundary_set
		for (int32_t index = 0; index < _meshdata.numberofsegments; index++){
			/*
			 * reverse the order
			 */
			auto bound_l = std::to_string(_Seg_R(index));
			auto bound_r = std::to_string(_Seg_L(index));
			_boundary_set.insert(bound_l + "-" + bound_r);
		}
		// end for

		int32_t counter = 2*_meshdata.numberofpoints + 2*(deg - 1)*_meshdata.numberofedges;

		std::unordered_map<std::string, std::vector<int32_t>>::const_iterator EdgeIterator;
		// iterator for boundary
		std::unordered_set<std::string>::const_iterator BoundaryIterator;
		std::unordered_map<std::string, int32_t>::const_iterator BoundaryElemIterator;


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


			/*
			 * check if boundary is in boundary set
			 */

			// 0 - 1
			BoundaryIterator = _boundary_set.find(tri_u + "-" + tri_v);
			if (BoundaryIterator != _boundary_set.end()) {
				boundary_elems.insert(make_pair(tri_u + "-" + tri_v, index));
			}
			// 1 - 2
			BoundaryIterator = _boundary_set.find(tri_v + "-" + tri_w);
			if (BoundaryIterator != _boundary_set.end()){
				boundary_elems.insert(make_pair(tri_v + "-" + tri_w, index));
			}
			// 2 - 0
			BoundaryIterator = _boundary_set.find(tri_w + "-" + tri_u);
			if (BoundaryIterator != _boundary_set.end()){
				boundary_elems.insert(make_pair(tri_w + "-" + tri_u, index));
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
			BoundaryElemIterator = boundary_elems.find(bound_l + "-" + bound_r);


			if (EdgeIterator != topology.edges.end()){
				std::copy(EdgeIterator->second.begin(), EdgeIterator->second.end(), topology.boundary.begin() + (deg + 1)*index + 2);
			}

			if (EdgeIterator == topology.edges.end()){
				EdgeIterator = topology.edges.find(bound_r + "-" + bound_l);
				std::reverse_copy(EdgeIterator->second.begin(), EdgeIterator->second.end(), topology.boundary.begin() + (deg + 1)*index + 2);
			}

			if (BoundaryElemIterator != boundary_elems.end()) {
				topology.boundary_index.push_back(BoundaryElemIterator->second);
			}
			else {
				BoundaryElemIterator = boundary_elems.find(bound_r + "-" + bound_l);
				if (BoundaryElemIterator != boundary_elems.end()) {
					topology.boundary_index.push_back(BoundaryElemIterator->second);
				}
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
	if (_meshdata.neighborlist != nullptr) {free(_meshdata.neighborlist); _meshdata.neighborlist = nullptr;}
} // clear


} /* namespace MEX */

using namespace MEX;

template class mexplus::Session<Mesh>;

namespace {

// Create a new instance of Mesh and return its session id.
MEX_DEFINE(new) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Mesh>::create(new Mesh(const_cast<mxArray*>(input.get(0)),
			const_cast<mxArray*>(input.get(1)), const_cast<mxArray*>(input.get(2)))));
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
	OutputArguments output(nlhs, plhs, 5);
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
	plhs[3] = mxCreateNumericMatrix(1, mesh->topology.boundary_index.size(), mxINT32_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[3]), &mesh->topology.boundary_index[0], mesh->topology.boundary_index.size()*sizeof(int));
	MatlabPtr temp_3[] = {plhs[3], mxCreateDoubleScalar(1.0)};
	mexCallMATLAB(1, &plhs[3], 2, temp_3, "plus");
	plhs[4] = mxCreateNumericMatrix(3, mesh->topology.neighbors.size()/3, mxINT32_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[4]), &mesh->topology.neighbors[0], mesh->topology.neighbors.size()*sizeof(int));
	MatlabPtr temp_4[] = {plhs[4], mxCreateDoubleScalar(1.0)};
	mexCallMATLAB(1, &plhs[4], 2, temp_4, "plus");
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




