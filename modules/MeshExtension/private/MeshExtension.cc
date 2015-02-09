/*
 * MeshExtension.cpp
 *
 *  Created on: Feb 7, 2015
 *      Author: lurker
 */

#include "MeshExtension.h"

namespace Extension {

MeshExtension::MeshExtension() {

}

MeshExtension::~MeshExtension() {
#ifdef DEBUG
	mexPrintf("MeshExtension detached...\n");
#endif
	layer.clear();
}

void MeshExtension::LayerElementsIndex(MatlabPtr Nodes,
		MatlabPtr Elems, MatlabPtr Interior){

	/*
	 * use naive linear search for each elements, only search based on its first 3 nodes' average
	 */

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  pedge_ptr            = mxGetPr(Interior);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	/*
	 * must be a polygon, input as [0 0 1 0 1 1 0 1]'
	 * or [0; 0; 1; 0; 1; 1; 0; 1] as column vector.
	 * length is twice the edges.
	 */
	auto numberofedges          = mxGetM(Interior)/2 ;

	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t center_x, center_y;

	for (size_t i =0; i < numberofelem; i++){

		vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


		center_x = (pnodes_ptr[2*vertex_1    ] + pnodes_ptr[2*vertex_2    ] + pnodes_ptr[2*vertex_3    ])/3.0;
		center_y = (pnodes_ptr[2*vertex_1 + 1] + pnodes_ptr[2*vertex_2 + 1] + pnodes_ptr[2*vertex_3 + 1])/3.0;


		/*
		 *
		 * Idea :
		 *
		 * If inside the interior, continue.
		 * Else pick out.
		 */

		/*
		 *
		 *	- Draw a horizontal line to the right of each point and extend it to infinity
		 *
		 *	- Count the number of times the line intersects with polygon edges.
		 *
		 *	- A point is inside the polygon if either count of intersections is odd or
		 *	  point lies on an edge of polygon.  If none of the conditions is true, then
		 *	  point lies outside.
		 */


		/*
		 * remark: the center of triangle won't be on a edge.
		 */

		Real_t f_x, f_y, s_x, s_y;
		size_t counter = 0;

		for (size_t j = 0; j < numberofedges; j++) {

			f_x = pedge_ptr[2 * j    ];
			f_y = pedge_ptr[2 * j + 1];

			if (j != numberofedges - 1){
				s_x = pedge_ptr[2 * (j + 1)];
				s_y = pedge_ptr[2 * (j + 1) + 1];
			}
			else{
				s_x = pedge_ptr[0];
				s_y = pedge_ptr[1];
			}


			if (fabs(f_y - s_y) > MEX_EPS) {

				if ((center_y - s_y) * (center_y - f_y) <= 0 &&
						(s_x*(f_y - center_y) + f_x*(center_y - s_y))/(f_y - s_y) >= center_x){
					counter++;
				}
			}
			else {
				if ((center_y - s_y) * (center_y - f_y) <= 0 &&
						(center_x - s_x) * (center_x - f_x) <= 0){
					counter++;
				}
			}

		}// end for
		if (!counter&1){
			layer.push_back(i);
		}
	}//end for
}//end all




} /* namespace Extension */

using namespace Extension;

template class mexplus::Session<MeshExtension>;

namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 0);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<MeshExtension>::create(new MeshExtension()));
}

MEX_DEFINE(delete)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<MeshExtension>::destroy(input.get(0));
}

MEX_DEFINE(layer)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 4);
	OutputArguments output(nlhs, plhs, 1);
	MeshExtension* mesh_ex = Session<MeshExtension>::get(input.get(0));

	mesh_ex->LayerElementsIndex(CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]));

	auto numberofindex = mesh_ex->layer.size();

	plhs[0] = mxCreateNumericMatrix(numberofindex, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* ptr = mxGetPr(plhs[0]);

	for (size_t i = 0; i < numberofindex; i++) {
		ptr[i] = mesh_ex->layer[i];
	}
}
}

MEX_DISPATCH
