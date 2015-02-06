/*
 * AssemblerExtension.cpp
 *
 *  Created on: Feb 5, 2015
 *      Author: lurker
 */

#include "AssemblerExtension.h"

namespace Extension {

AssemblerExtension::AssemblerExtension(){

}

AssemblerExtension::~AssemblerExtension(){
#ifdef DEBUG
	mexPrintf("Assembler Extension detached\n");
#endif
}


void AssemblerExtension::AssembleGradXFunc(Real_t* &pI, Real_t* &pJ, Real_t* &pV,
		MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr Ref, MatlabPtr RefX, MatlabPtr RefY,
		MatlabPtr Weights, MatlabPtr Fcn){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area;
	Real_t Jacobian[2][2];


	// evaluate Fcn at all quadrature nodes before calculation
	if (mxGetNumberOfElements(Fcn) == numberofelem*numberofqnodes ){

		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + Interp[i*numberofqnodes + l] *
						           (Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] +
						            Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])
									* reference[k+ l*numberofnodesperelem]*
									weights[l];
					}
					*pV  = (*pV)/2.0;
					pI++; pJ++; pV++;
				}
			}
		}//end for
	}//end if
	else {
		mexErrMsgTxt("Error:AssemberExtension:GradXFunc:Dimension does not match.\n");
		mexErrMsgTxt("Not implemented yet...\n");
	}
}

void AssemblerExtension::AssembleGradYFunc(Real_t* &pI, Real_t* &pJ, Real_t* &pV,
		MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr Ref, MatlabPtr RefX,MatlabPtr RefY,
		MatlabPtr Weights, MatlabPtr Fcn){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area;
	Real_t Jacobian[2][2];


	// evaluate Fcn at all quadrature nodes before calculation
	if (mxGetNumberOfElements(Fcn) == numberofelem*numberofqnodes ){

		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + Interp[i*numberofqnodes + l] *
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] +
								 Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])
								* reference[k+ l*numberofnodesperelem]*
								weights[l];
					}
					*pV  = (*pV)/2.0;
					pI++; pJ++; pV++;
				}
			}
		}//end for
	}//end if
	else {
		mexErrMsgTxt("Error:AssemberExtension:GradXFunc:Dimension does not match.\n");
		mexErrMsgTxt("Not implemented yet...\n");
	}
}

// todo: try to do integration on int  A(\nabla phi) \psi to make the calls into one.

}



using namespace Extension;

template class mexplus::Session<AssemblerExtension>;

namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 0);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<AssemblerExtension>::create(new AssemblerExtension()));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<AssemblerExtension>::destroy(input.get(0));
}

MEX_DEFINE(assemex_gradfunc_x) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 8);
	OutputArguments output(nlhs, plhs, 3);
	AssemblerExtension* assembler_ex = Session<AssemblerExtension>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);


	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	Real_t* pV = mxGetPr(plhs[2]);

	assembler_ex->AssembleGradXFunc(pI, pJ, pV, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]), CAST(prhs[7]));
}

MEX_DEFINE(assemex_gradfunc_y)  (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 8);
	OutputArguments output(nlhs, plhs, 3);
	AssemblerExtension* assembler_ex = Session<AssemblerExtension>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);


	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	Real_t* pV = mxGetPr(plhs[2]);

	assembler_ex->AssembleGradYFunc(pI, pJ, pV, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]), CAST(prhs[7]));
}
}
MEX_DISPATCH
