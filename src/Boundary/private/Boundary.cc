/*
 * Boundary.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: lurker
 */

#include "Boundary.h"

namespace MEX {

Boundary::Boundary() {

}


Boundary::Boundary(MatlabPtr _edges) {

	b_edges.resize(mxGetM(_edges)*mxGetN(_edges));
	memcpy(&b_edges[0], mxGetPr(_edges), mxGetM(_edges)*mxGetN(_edges)*sizeof(int32_t));
	for (auto it = b_edges.begin(); it != b_edges.end(); it++) {
		b_edge_set.insert(*it);
	}

	// b_expr now is empty
}

Boundary::~Boundary() {

	b_edges.clear();
	b_edge_set.clear();
	b_expr.clear();
#ifdef DEBUG
	mexPrintf("Boundary detached\n");
#endif
}

void Boundary::setDirichlet(MatlabPtr _edges) {

	b_edges.resize(mxGetM(_edges)*mxGetN(_edges));
	memcpy(&b_edges[0], mxGetPr(_edges), mxGetM(_edges)*mxGetN(_edges)*sizeof(int32_t));
	for (auto it = b_edges.begin(); it != b_edges.end(); it++) {
		b_edge_set.insert(*it);
	}
}

} /* namespace MEX */



using namespace MEX;

template class mexplus::Session<DirichletBC>;


namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<DirichletBC>::create(new DirichletBC(
			CAST(input.get(0))
						)));
}

MEX_DEFINE(placeholder)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 0);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<DirichletBC>::create(new DirichletBC()));

}

MEX_DEFINE(set_dirichlet)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 0);
	auto boundary = Session<DirichletBC>::get(input.get(0));
	boundary->setDirichlet(CAST(input.get(1)));

}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<DirichletBC>::destroy(input.get(0));
}

MEX_DEFINE(report) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	auto boundary = Session<DirichletBC>::get(input.get(0));
	mexPrintf("Dirichlet Boundary Nodes Number: %d \n", boundary->b_edge_set.size());

}

// can be implemented by set_difference from <algorithm>
MEX_DEFINE(dofs) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 1);
	auto boundary = Session<DirichletBC>::get(input.get(0));
	auto _range   = input.get<int32_t>(1);

	plhs[0] = mxCreateNumericMatrix(_range - boundary->b_edge_set.size(), 1, mxINT32_CLASS, mxREAL);

	auto plhs_ptr = Matlab_Cast<int32_t>(plhs[0]);
	for (int32_t i = 1; i <= _range; ++i) {
		if (boundary->b_edge_set.find(i) != boundary->b_edge_set.end()) {
			continue;
		}
		else {
			*plhs_ptr++ = i;
		}
	}

}

// use exprTk for judging edges
// it is not very suitable to be with DirichletBC, since it could be also
// some other boundary conditions.
MEX_DEFINE(set_boundary)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 2);
	auto boundary = Session<DirichletBC>::get(input.get(0));
	std::string expr_string(mxArrayToString(prhs[1]));
	auto curr_size = boundary->b_expr.size();
	boundary->b_expr.push_back(expr_string);
}

MEX_DEFINE(get_boundary)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	/*
	 * Does not apply the input and output argument here,
	 * since the size of LHS is unclear.
	 */
	InputArguments input(nrhs, prhs, 3);
	auto boundary  = Session<DirichletBC>::get(input.get(0));

	auto edges_ptr = Matlab_Cast<int32_t>(CAST(prhs[1]));
	auto nodes_ptr = Matlab_Cast<Real_t>(CAST(prhs[2]));

	auto _num_edge = mxGetN(prhs[1]);
	auto _num_node = mxGetM(prhs[1]);


	if (_num_node < 2) {
		mexErrMsgTxt(
				"Error:Boundary:get_boundary::Each edge should contain at least 2 nodes.\n");
	}

	auto _num_expr = boundary->b_expr.size();

	if (_num_expr != nlhs) {
		mexErrMsgTxt("Error:Boundary:get_boundary:: RHS should contain as many as boundary conditions.\n");
	}

	/*
	 * stores edges for each boundary
	 */
	vector<vector<int32_t>> _global;
	_global.resize(_num_expr);

	exprtk::symbol_table<Real_t> symbol_table;

	Real_t x = 0., y = 0.;

	symbol_table.add_variable("x", x);
	symbol_table.add_variable("y", y);
	symbol_table.add_constants();

	vector<exprtk::expression<Real_t>> expressions(_num_expr);
	vector<exprtk::parser<Real_t>> parsers(_num_expr);

	for (size_t i = 0; i < _num_expr; i++) {
		expressions[i].register_symbol_table(symbol_table);
		parsers[i].compile(boundary->b_expr[i], expressions[i]);
	}

	Real_t _x_1, _y_1, _x_2, _y_2, result_1, result_2;

	for (size_t i = 0; i < _num_edge; i++){

		_x_1 = nodes_ptr[2*(edges_ptr[i*_num_node] - 1)];
		_y_1 = nodes_ptr[2*(edges_ptr[i*_num_node] - 1) + 1];
		_x_2 = nodes_ptr[2*(edges_ptr[i*_num_node + 1] - 1)];
		_y_2 = nodes_ptr[2*(edges_ptr[i*_num_node + 1] - 1) + 1];


		for (size_t j = 0 ; j < _num_expr; j++){
			x = _x_1; y= _y_1;
			result_1 = fabs(expressions[j].value());
			x = _x_2; y= _y_2;
			result_2 = fabs(expressions[j].value());

			if (result_1 < MEX_EPS && result_2 < MEX_EPS){
				for (size_t k = 0 ; k < _num_node; k++){
					_global[j].push_back(edges_ptr[i*_num_node + k]);
				}
			}
			else {
				continue;
			}
		}
	}
	for (size_t i = 0; i < _num_expr; i++) {
		plhs[i] = mxCreateNumericMatrix(_num_node, _global[i].size()/_num_node, mxINT32_CLASS,mxREAL);
		// implicit convert pointer
		memcpy(mxGetPr(plhs[i]),
				&_global[i][0], sizeof(int32_t)*_global[i].size());
	}
	_global.clear();
	expressions.clear();
	parsers.clear();
}

} // namespace

MEX_DISPATCH
