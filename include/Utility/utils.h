/*
 * utils.h
 *
 *  Created on: Oct 17, 2014
 *      Author: lurker
 */

#ifndef INCLUDE_UTILITY_UTILS_H_
#define INCLUDE_UTILITY_UTILS_H_

#include "mex.h"
#include "math.h"
#include <omp.h>

typedef mxArray* MatlabPtr;
typedef double   Real_t;

//#define DEBUG


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

inline Real_t INTERSECT_DET(Real_t X1, Real_t Y1, Real_t X2, Real_t Y2, Real_t THETA){
		return ((cos((THETA)) * ((Y1) - (Y2))) - (sin((THETA)) * ((X1) - (X2))));
	}
inline Real_t INTERSECT_CROSS(Real_t X1,Real_t Y1,Real_t X2,Real_t Y2,Real_t A,Real_t B) {
		return ((((X1) - (X2)) * ((B) - (Y2))) -(((A) - (X2)) * ((Y1)- (Y2))));
	}


extern "C" bool mxUnshareArray(mxArray *array_ptr, bool noDeepCopy);

#endif /* INCLUDE_UTILITY_UTILS_H_ */
