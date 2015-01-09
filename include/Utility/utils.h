/*
 * utils.h
 *
 *  Created on: Oct 17, 2014
 *      Author: lurker
 */

#ifndef INCLUDE_UTILITY_UTILS_H_
#define INCLUDE_UTILITY_UTILS_H_

#include "mex.h"

#include <omp.h>

typedef mxArray* MatlabPtr;
typedef double   Real_t;

#undef DEBUG


using namespace std;

#define MEX_EXPECT(condition) if (!(condition)) \
    mexErrMsgTxt(#condition " not true.")

/*
 * make const MatlabPtr into nonconst MatlabPtr
 */
#define CAST(PTR) const_cast<MatlabPtr>(PTR)

#define MEX_EPS 1e-9

template <typename T> T* Matlab_Cast(MatlabPtr _ptr){
	return static_cast<T*> (mxGetData(_ptr));
}

inline size_t Msize(MatlabPtr _ptr){
	return mxGetM(_ptr);
}

inline size_t Nsize(MatlabPtr _ptr){
	return mxGetN(_ptr);
}


#endif /* INCLUDE_UTILITY_UTILS_H_ */
