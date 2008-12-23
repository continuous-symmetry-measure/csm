/**
 * A simple solver for steady state of system of linear equations using
 * newton's method
 */
#ifndef _NEWTON_SOLVER_H_
#define _NEWTON_SOLVER_H_

#include "equation_set.h"
#include "rk_solver.h"
#include <iostream>

using namespace std;

template <typename T>
class NewtonSolver : RKResultProcessor<T> {
public:
	typedef typename EquationSet<T>::vec vec;
private:
	rk_params params;
	RKResultProcessor<T>& processor;
	EquationSet<T>& eqSet;
	vec state; 
		
public: 

	NewtonSolver(RKResultProcessor<T>& resProcessor, EquationSet<T>& equations) : processor(resProcessor), eqSet(equations) {}

	virtual void solutionComplete(const vec& finalState) { 
		state = finalState;
	}

	void solve(const rk_params& paramSet, vec& initialState, size_t initialSteps) {
		// Copy Params
		params = paramSet;
		
		std::cerr << "Running Runga-Kutta to get initial guess" << endl; 

		// Run 50 RK steps
		RKSolver<T> solver(*this, eqSet);
		solver.runSteps(paramSet, initialState, initialSteps); 

		std::cerr << "Running newton's method using guess" << endl;

		processor.solvingStarted(params, state);
		
		double max_err = params.max_error;
		double err = max_err + 1;
		size_t step = 0;
		while (err > max_err) {
			vec newState = state - (eqSet.compute_derivative(0.0,state) / eqSet.compute_derivative_derivative(state));
			processor.stepPerformed(step++, 0, newState, state, true);
			err = eqSet.computeError(newState, state, 1);
			state = newState;			
			cout << "Error: " << err << endl;
		}
		
		processor.solutionComplete(state);
	}
};
#endif
