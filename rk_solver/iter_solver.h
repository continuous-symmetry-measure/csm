/**
 * A simple solver for steady state of system of linear equations using
 * newton's method
 */
#ifndef _NEWTON_SOLVER_H_
#define _NEWTON_SOLVER_H_

#include "equation_set.h"
#include "rk_solver.h"
#include <iostream>

#define MAX_DER_LIMIT 1e-40
using namespace std;

template <typename T>
class IterSolver : RKResultProcessor<T> {
public:
	typedef typename EquationSet<T>::vec vec;
private:
	rk_params params;
	RKResultProcessor<T>& processor;
	EquationSet<T>& eqSet;
	vec state; 
		
public: 

	IterSolver(RKResultProcessor<T>& resProcessor, EquationSet<T>& equations) : processor(resProcessor), eqSet(equations) {}

	virtual void solutionComplete(const vec& finalState) { 
		state = finalState;
	}

	void solve(const rk_params& paramSet, vec& initialState, double a = 0.1, size_t initialSteps = 0) {
		// Copy Params
		params = paramSet;

		state = initialState.copy();
		
		std::cerr << "Running Runga-Kutta to get initial guess" << endl; 

		// Run 50 RK steps
		RKSolver<T> solver(*this, eqSet);
		solver.runSteps(paramSet, initialState, initialSteps); 

		std::cerr << "Running newton's method using guess" << endl;

		processor.solvingStarted(params, state);
		
		double max_err = params.limit;
		double err = max_err + 1;
		double max_der = 1;
		size_t step = 0;
		while (err > max_err || max_der < MAX_DER_LIMIT) {
			vec newState = a * state + (1 - a) * eqSet.comptueRnDivWn(state);
			processor.stepPerformed(step++, 0, eqSet.normalize(newState), eqSet.normalize(state), true);
			err = eqSet.computeError(newState, state, 1);			
			state = newState;
			max_der = fabs( scalarMax(eqSet.compute_derivative(0.0, state)));
			cout << "Error: " << err << ", Max prob: " << max_der << endl;
		}
		
		processor.solutionComplete(state);
	}
};
#endif
