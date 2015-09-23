/*
 * Assembler.h
 *
 *  Created on: Oct 10, 2014
 *      Author: lurker
 */

#ifndef ASSEMBLER_PRIVATE_ASSEMBLER_H_
#define ASSEMBLER_PRIVATE_ASSEMBLER_H_

#include <cstdlib>
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


#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

class Assembler {
public:
	Assembler();
	virtual ~Assembler();

	/*
	 * public methods
	 */
	/*
	 * Reference Matrix
	 */
	void Reference(MatlabPtr&, MatlabPtr&, MatlabPtr&,MatlabPtr,MatlabPtr);
	/*
	 * Mass Matrix
	 */
	void AssembleMass(Real_t*&, Real_t*&, Real_t*&, MatlabPtr, MatlabPtr,
			MatlabPtr, MatlabPtr, MatlabPtr);
	/*
	 * Stiffness Matrix
	 */
	// if kernel is a function
	void AssembleStiff(Real_t*&, Real_t*&, Real_t*&,MatlabPtr, MatlabPtr,
			MatlabPtr, MatlabPtr, MatlabPtr, MatlabPtr);

	// with a symmetric/diagonal matrix kernel
	void AssembleStiff(Real_t*&, Real_t*&, Real_t*&,MatlabPtr, MatlabPtr,
			MatlabPtr, MatlabPtr, MatlabPtr, MatlabPtr, MatlabPtr);
	/*
	 * Load Vector from Robin Boundary
	 */
	void AssembleLoad(Real_t*& pLoad, MatlabPtr Nodes,
			MatlabPtr QNodes, MatlabPtr Elems,MatlabPtr Ref,
			MatlabPtr Weights, MatlabPtr Fcn);

	/*
	 * Integral for point source < Fcn , delta(x)> = Fcn(x) in infinite dimensional space
	 * if x is not a node, it will be some weighted expression(localized) function.
	 */
	void AssembleLoad(Real_t*& pLoad, MatlabPtr _point, MatlabPtr Fcn);

	// 1d Integral
	void Reference(MatlabPtr&, MatlabPtr&, MatlabPtr, MatlabPtr);

	void AssembleBC(Real_t*& pNeumann,  MatlabPtr Nodes,
			MatlabPtr QNodes, MatlabPtr eNeumann,
			MatlabPtr Ref, MatlabPtr Weights,  MatlabPtr Fcn);
	void AssembleBC(Real_t* &pI, Real_t* &pJ, Real_t* &pV,
			MatlabPtr Nodes, MatlabPtr eRobin,
			MatlabPtr Ref, MatlabPtr Weights, MatlabPtr Fcn);

	// Auxiliary
	void Qnodes2D(Real_t*& Coords, MatlabPtr Nodes, MatlabPtr QNodes, MatlabPtr Elems);
	void Qnodes1D(Real_t*& Coords, MatlabPtr Nodes, MatlabPtr QNodes, MatlabPtr Edges);


	// modules
	/*
	 * integral on int \phi_i_x \phi_j
	 */
	void AssembleGradXFunc(Real_t* &pI, Real_t* &pJ, Real_t* &pV,
			MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr Ref, MatlabPtr RefX,MatlabPtr RefY,
			MatlabPtr Weights, MatlabPtr Fcn);
	/*
	 * integral on int \phi_i_y \phi_j
	 */
	void AssembleGradYFunc(Real_t* &pI, Real_t* &pJ, Real_t* &pV,
			MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr Ref, MatlabPtr RefX, MatlabPtr RefY,
			MatlabPtr Weights, MatlabPtr Fcn);
	/*
	 * combine them all
	 */
	void AssembleGradXYFunc(Real_t* &pI, Real_t* &pJ, Real_t* &pV,Real_t* &pW,
			MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr Ref, MatlabPtr RefX,MatlabPtr RefY,
			MatlabPtr Weights, MatlabPtr Fcn_X, MatlabPtr Fcn_Y);



	/*
	 * load matrix
	 */
	void AssembleLoadMatrix(Real_t*& pI, Real_t*& pJ, Real_t*& pV, MatlabPtr Nodes,
			 MatlabPtr Elems,MatlabPtr Ref,
			MatlabPtr Weights, MatlabPtr Fcn);


	void AssembleOverElement(Real_t*& w, MatlabPtr Nodes, MatlabPtr Elems,
			MatlabPtr Ref, MatlabPtr RefX,
			MatlabPtr RefY, MatlabPtr Weights, MatlabPtr Fcn_s, MatlabPtr Fcn_a,
			MatlabPtr u, MatlabPtr v);


};

#endif /* ASSEMBLER_PRIVATE_ASSEMBLER_C_ */
