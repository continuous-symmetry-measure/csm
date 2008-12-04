/**
 * This file provides an interface to for the equation solver
 */
#ifndef _EQUATION_SET_H_
#define _EQUATION_SET_H_

#include "vm/vec_mat.h"	

template <typename T>
class EquationSet {

public: 

	typedef Vector<T> vec;

	/** 
	 * Update any relevant parameters prior to the next step
	 * 
	 * Default implementation is empty, assuming there are no varying parameters
	 * 
	 * @param state The current (pre-step) state
	 */
	virtual void updateParameters(const vec& state) { }

	/** 
 	 * Compute the error, given the previous state, the current state, and the used delta
  	 */
	virtual double computeError(const vec& state, const vec& prevState, double usedDelta) { 
		return scalarMax(abs(prevState - state) / state) / usedDelta;
	}

	/** 
	 * According to a current state, and the 
	 * relevant equations, compute the time derivative.
	 * 
	 * @param time The current time
	 * @param state The current state to compute derivative for
	 * @return The time derivative of the equations
	 */ 
	virtual vec compute_derivative(const double t, const vec& state) = 0;

	virtual ~EquationSet() {}
};

#endif