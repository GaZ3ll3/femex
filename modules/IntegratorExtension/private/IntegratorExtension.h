/*
 * IntegratorExtension.h
 *
 *  Created on: Feb 18, 2015
 *      Author: lurker
 */


/*
 * 1. calculate ray points for integration from internal point to boundary.
 *
 * 2. calculate the integration with function handle.
 *
 */

#ifndef MODULES_INTEGRATOREXTENSION_PRIVATE_INTEGRATOREXTENSION_H_
#define MODULES_INTEGRATOREXTENSION_PRIVATE_INTEGRATOREXTENSION_H_

class IntegratorExtension {
public:
	IntegratorExtension();
	virtual ~IntegratorExtension();

	/*
	 * G is the function to be integrated on nodes.
	 * ray is calculated and interpolated.
	 */
	void Int_Along_Ray(MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr Edges,
			MatlabPtr Fcn, MatlabPtr Theta, MatlabPtr G);


};

#endif /* MODULES_INTEGRATOREXTENSION_PRIVATE_INTEGRATOREXTENSION_H_ */
